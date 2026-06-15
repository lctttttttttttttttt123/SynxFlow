// ======================================================================================
// T6 §4: phosphorus sorption closures — Langmuir isotherm + sorption rate (口径 I).
// LICENCE: GPLv3
// ======================================================================================
//
// Design authority: T6_design.md §4 + T6_wiring_plan.md decision (B)=口径 I:
//   - SINGLE equilibrium reference = Langmuir q*(Pd); EPC0 NEVER enters the kernel rate
//     (EPC0 is a derived diagnostic only, handled in the Python layer).
//   - kinetics (mode 0, default): dq/dt = k_ads*(q*-q), q* = langmuir_qstar(Qmax,K,Pd).
//   - instant equilibrium (mode 1): handled at OPERATOR level by a coupled multi-group solve
//     (see cuPhosphorusSorption in cuda_phosphorus.cu); not driven from this per-group closure.
// Enum-dispatch, __device__ inline, ZERO hard-coded physical number (every coefficient comes from
// PhosphorusParams, read from input/phosphorus_setup.dat). Adding a NEW closure = add a case +
// recompile; swapping an implemented closure or changing parameters = config-only, no recompile.
//
// Mirrors cuda_sediment_closures.h style. Included only by CUDA translation units (nvcc).
//
#ifndef CUDA_PHOSPHORUS_CLOSURES_H
#define CUDA_PHOSPHORUS_CLOSURES_H

#include "Scalar.h"
#include "cuda_phosphorus.h"   // PhosphorusParams (POD, device-usable)

namespace GC {

  // ---- Langmuir equilibrium sorbed load q*(Pd) [mg-P/g] ----
  // q* = Qmax*K*Pd / (1 + K*Pd). Qmax=0 => inert group => q*=0 (short-circuit, T6_design §3).
  // Denominator 1+K*Pd >= 1 for K,Pd >= 0 => no singularity; this function NEVER divides by q*.
  __device__ __forceinline__ Scalar langmuir_qstar(const PhosphorusParams& p, Scalar Pd) {
    if (p.Qmax <= (Scalar)0.0) return (Scalar)0.0;          // inert group short-circuit (no phosphorus)
    Scalar kpd   = p.K * Pd;
    Scalar denom = (Scalar)1.0 + kpd;
    if (denom <= (Scalar)0.0) return (Scalar)0.0;           // defensive (K,Pd normally >= 0)
    return p.Qmax * kpd / denom;                            // Qmax*K*Pd/(1+K*Pd)
  }

  // ---- net sorbed-load change Δq over one step [mg-P/g] ----
  // Sign convention: Δq>0 net adsorption (Pd->particle, Pd down), Δq<0 net desorption (particle->Pd).
  // Direction comes from (q*-q) alone — NO external sign source (T6_design 修订①: forbids the voided
  // (q*-q)^2*sign form; the dead band q≈q* gives Δq->0 continuously, no jump).
  //   mode_id=0 (default, one-order kinetics): Δq = k_ads*(q*-q)*dt.
  //   reserved sub-switch (NOT implemented): pseudo-second-order, direction still from (q*-q) sign,
  //     e.g. dq/dt = k_ads*(q*-q)*|q*-q|. enum slot kept; falls through to one-order for now.
  // STIFF HARD GUARD (default, T6_design §4/§10): |Δq| <= |q*-q|. Capping to the gap magnitude (the
  // amount that pushes q exactly to q*) is the implicit-stable upper bound at ~zero cost: q reaches q*
  // at most, never overshoots -> ANY configured k_ads stays bounded even if equilibrium time << dt.
  // Qmax=0 inert group: forced zero transfer, never a /q* divide.
  // NOTE: this is the per-group kinetic rate (mode 0). mode_id=1 (instant equilibrium) is resolved by
  // the coupled multi-group solver in the operator to keep dissolved-phase positivity under several
  // groups sharing one Pd; callers must NOT use this function's output to drive a mode-1 update.
  __device__ __forceinline__ Scalar phos_sorption_rate(const PhosphorusParams& p, int mode_id,
                                                       Scalar Pd, Scalar q, Scalar dt) {
    if (p.Qmax <= (Scalar)0.0) return (Scalar)0.0;          // Qmax=0 inert: zero transfer
    Scalar qstar = langmuir_qstar(p, Pd);
    Scalar gap   = qstar - q;                               // drive; sign carries direction
    Scalar dq;
    switch (mode_id) {
      // mode 1 target (single-group view): straight to q*(Pd). The operator overrides this with a
      // coupled solve; kept here only so the enum dispatch is total.
      case 1:  dq = gap;                       break;
      // mode 0 (and the reserved pseudo-second-order sub-switch) one-order kinetics
      case 0:
      default: dq = p.k_ads * gap * dt;        break;
    }
    // stiff hard guard: |Δq| <= |gap|, sign preserved (for mode 0, k_ads,dt>=0 => sign(dq)==sign(gap))
    Scalar cap = fabs(gap);
    if (dq >  cap) dq =  cap;
    if (dq < -cap) dq = -cap;
    return dq;
  }

}  // namespace GC

#endif
