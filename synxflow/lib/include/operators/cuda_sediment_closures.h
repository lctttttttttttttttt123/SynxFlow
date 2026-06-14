// ======================================================================================
// T5 step3: sediment closures — settling velocity + erosion/deposition rate.
// LICENCE: GPLv3
// ======================================================================================
//
// Design authority: T5_design.md §6 (enum-dispatch, __device__ inline) + §2 (replaceable
// interfaces). NO physical number is hard-coded — every coefficient comes from SedimentParams
// (read from input/sediment_setup.dat). Adding a NEW formula = add a switch case + recompile;
// swapping an already-compiled formula or changing parameters = config-only, no recompile.
//
// Included only by CUDA translation units (nvcc) — provides __device__ functions.
//
#ifndef CUDA_SEDIMENT_CLOSURES_H
#define CUDA_SEDIMENT_CLOSURES_H

#include "Scalar.h"
#include "cuda_sediment.h"   // SedimentParams (POD, device-usable)

namespace GC {

  // ---- settling velocity w_s [m/s] (enum-dispatch; const + Zhang implemented, rest fall back) ----
  // settling_id: 0=const(p.w_s)  1=Zhang(张瑞瑾)  2=Stokes  3=floc.
  // Implemented: 0 (const) and 1 (Zhang 张瑞瑾 1989, single grain in clear water). 2/3 fall back to
  // const(p.w_s) for now (reserved extension points — see §2 "先实现 const+Zhang, 余 stub+TODO").
  __device__ __forceinline__ Scalar sed_settling_velocity(const SedimentParams& p) {
    switch (p.settling_id) {
      case 1: {
        // Zhang Ruijin (张瑞瑾) formula: w_s = sqrt( (13.95 nu/D)^2 + 1.09 (rho_s/rho_w - 1) g D ) - 13.95 nu/D
        // nu = kinematic viscosity of water [m2/s]; constants are the published Zhang coefficients
        // (formula structure, not a tunable basin number). rho_s/D from config; rho_w, g, nu standard.
        const Scalar nu    = (Scalar)1.0e-6;   // water kinematic viscosity @ ~20C
        const Scalar rho_w = (Scalar)1000.0;
        const Scalar g     = (Scalar)9.81;
        Scalar D = p.D50;
        if (D <= (Scalar)0.0) return p.w_s;     // ill-defined grain size -> fall back to const
        Scalar a = (Scalar)13.95 * nu / D;
        Scalar b = (Scalar)1.09 * (p.rho_s / rho_w - (Scalar)1.0) * g * D;
        return sqrt(a*a + b) - a;
      }
      case 0:
      default:
        // 0=const; 2(Stokes)/3(floc) reserved -> const fallback (one-time host-side notice elsewhere)
        return p.w_s;
    }
  }

  // ---- erosion - deposition net rate [kg/m2/s in depth-integrated conc units] ----
  // closure_id=0: Partheniades two-threshold cohesive (default). [1=Shields reserved, unimplemented]
  //   deposition (tau_b < tau_cd): D = w_s * C * (1 - tau_b/tau_cd)        (>=0, settles out)
  //   erosion    (tau_b > tau_ce): E = M * (tau_b/tau_ce - 1)             (>=0, re-suspends)
  //   dead band  (tau_cd <= tau_b <= tau_ce): 0
  // Returns ED = E - D : >0 net re-suspension (C up), <0 net deposition (C down toward bed).
  // C is the local concentration (hC/h); w_s from sed_settling_velocity. All thresholds from config.
  __device__ __forceinline__ Scalar sed_ED_rate(const SedimentParams& p, Scalar tau_b, Scalar C, Scalar w_s) {
    // closure_id is reserved for future dispatch; only Partheniades (0) is implemented.
    Scalar D = (Scalar)0.0;
    Scalar E = (Scalar)0.0;
    if (p.tau_cd > (Scalar)0.0 && tau_b < p.tau_cd) {
      D = w_s * C * ((Scalar)1.0 - tau_b / p.tau_cd);
      if (D < (Scalar)0.0) D = (Scalar)0.0;
    }
    if (p.tau_ce > (Scalar)0.0 && tau_b > p.tau_ce) {
      E = p.M * (tau_b / p.tau_ce - (Scalar)1.0);
      if (E < (Scalar)0.0) E = (Scalar)0.0;
    }
    return E - D;
  }

}  // namespace GC

#endif
