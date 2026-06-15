// ======================================================================================
// T6 §5/§7: phosphorus operators — sorption exchange (Pd<->hPp_k) + bed transfer (route ①).
// LICENCE: GPLv3
// ======================================================================================
//
// Design authority: T6_design.md §5/§7 + T6_wiring_plan.md decisions (A)=路线① / (B)=口径 I.
// Physics formulas live in cuda_phosphorus_closures.h (Langmuir) and cuda_sediment_closures.h
// (sed_ED_rate / sed_settling_velocity — REUSED by §7, not re-derived). NO phosphorus parameter is
// hard-coded here — all come from per-group PhosphorusParams / SedimentParams uploaded to device.
// These operators run ONLY on the phosphorus_on (n_phos>0) path.
//
#include "cuda_phosphorus.h"
#include "cuda_phosphorus_closures.h"
#include "cuda_sediment.h"
#include "cuda_sediment_closures.h"
#include "cuda_kernel_launch_parameters.h"
#include "Scalar.h"
#include "Vector.h"
#include "cuda_mapped_field.h"

namespace GC {
  namespace fv {

    // ===== §7 bed transfer (route ①) — recompute dmass identically, move hPp_k<->bedPp_k ============
    // INVARIANTS (T6_wiring_plan 决策块 A):
    //  ① ORDER: this runs BEFORE cuSedimentErosionDeposition, reading the SAME pre-E-D read-only
    //     fields (h, tau_b, hC_k, bed_k). Writes ONLY hPp_k / bedPp_k — NEVER hC_k / bed_k.
    //  ② SINGLE SOURCE: dmass recomputed via the SAME sed_settling_velocity + sed_ED_rate closures and
    //     the SAME R5 clamp as cuSedimentErosionDepositionKernel (cuda_sediment.cu:69–77). The physics
    //     formula is reused by CALL, not re-derived. (The 4 clamp/dmass lines below are byte-for-byte
    //     the E-D kernel's; the consistency cross-test asserts bit-exact equality.)
    //  ③ GUARD: dmass_dbg (optional) lets the cross-test compare this dmass to the E-D kernel's hC delta.
    __global__ void cuPhosphorusBedTransferKernel(Scalar* h, Scalar* tau_b, Scalar** hCs, Scalar** beds,
                                                  Scalar** hPps, Scalar** bedPps, SedimentParams* sparams,
                                                  int n_phos, Scalar dt, Scalar** dmass_dbg,
                                                  unsigned int phi_size) {
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      const Scalar h_small = 1e-10;
      while (index < phi_size) {
        Scalar h_this = h[index];
        Scalar tb     = tau_b[index];
        for (int k = 0; k < n_phos; ++k) {
          Scalar hc  = hCs[k][index];
          Scalar bed = beds[k][index];
          // --- recompute dmass EXACTLY as cuSedimentErosionDepositionKernel does (same closures + R5) ---
          Scalar C   = (h_this >= h_small) ? hc / h_this : (Scalar)0.0;
          Scalar ws  = sed_settling_velocity(sparams[k]);
          Scalar ed  = sed_ED_rate(sparams[k], tb, C, ws);
          Scalar ed_min = -hc / dt;          // R5: ED >= -hC/dt
          Scalar ed_max = bed / dt;          // R5: ED <=  bed/dt
          if (ed < ed_min) ed = ed_min;
          if (ed > ed_max) ed = ed_max;
          Scalar dmass = ed * dt;            // >0 erosion (bed->susp), <0 deposition (susp->bed)
          if (dmass_dbg) dmass_dbg[k][index] = dmass;
          // --- move sorbed P by the SAME sediment mass fraction (phosphorus rides the identical dmass) ---
          Scalar hpp = hPps[k][index];
          Scalar bpp = bedPps[k][index];
          if (dmass > (Scalar)0.0) {                       // erosion: bed sediment re-suspends w/ its P
            Scalar f = (bed > (Scalar)0.0) ? dmass / bed : (Scalar)0.0;
            if (f > (Scalar)1.0) f = (Scalar)1.0;          // defensive (clamp guarantees f<=1)
            Scalar move = f * bpp;
            hPps[k][index]   = hpp + move;
            bedPps[k][index] = bpp - move;
          } else if (dmass < (Scalar)0.0) {                // deposition: suspended sediment settles w/ its P
            Scalar f = (hc > (Scalar)0.0) ? (-dmass) / hc : (Scalar)0.0;
            if (f > (Scalar)1.0) f = (Scalar)1.0;
            Scalar move = f * hpp;
            hPps[k][index]   = hpp - move;
            bedPps[k][index] = bpp + move;
          }
          // dmass==0 dead band: no transfer (dmass_dbg already written 0 above)
        }
        index += blockDim.x * gridDim.x;
      }
    }

    void cuPhosphorusBedTransfer(cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& tau_b,
                                 Scalar** hCs_dev, Scalar** beds_dev, Scalar** hPps_dev, Scalar** bedPps_dev,
                                 SedimentParams* sparams_dev, int n_phos, Scalar dt, Scalar** dmass_dbg_dev) {
      cuPhosphorusBedTransferKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        h.data.dev_ptr(), tau_b.data.dev_ptr(), hCs_dev, beds_dev, hPps_dev, bedPps_dev,
        sparams_dev, n_phos, dt, dmass_dbg_dev, h.data.size());
    }

    // ===== §5 sorption — dissolved hPd <-> particulate hPp_k exchange (Langmuir, 口径 I) ===========
    // q_k = hPp_k / hC_k with sediment-tends-to-zero divide guard (R5; revision C: never /hC_k when
    // sediment ~0). Total depth-integrated P conserved: Δ(hPd) = -Σ_k Δ(hPp_k).
    //  mode 0 (kinetics, default): per group, sequential, Δq = phos_sorption_rate (stiff-clamped);
    //     dissolved-positivity cap so hPd stays >= 0 under several groups sharing one Pd.
    //  mode 1 (instant equilibrium): coupled monotone bisection for Pd_eq (h*Pd + Σ hC_k*q*(Pd) = T),
    //     then q_k = q*(p_k, Pd_eq). One step to the joint equilibrium; positivity by construction.
    __global__ void cuPhosphorusSorptionKernel(Scalar* h, Scalar* hPd, Scalar** hPps, Scalar** hCs,
                                               PhosphorusParams* p, int mode_id, int n_phos, Scalar dt,
                                               unsigned int phi_size) {
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      const Scalar h_small  = 1e-10;
      const Scalar hc_small = 1e-10;    // sediment-tends-to-zero guard for q_k = hPp_k/hC_k (R5)
      while (index < phi_size) {
        Scalar h_this = h[index];
        if (h_this < h_small) { index += blockDim.x * gridDim.x; continue; }   // dry cell: no water phase
        Scalar hpd = hPd[index];
        if (mode_id == 1) {
          // ---- instant equilibrium: conserve total P, solve increasing f(Pd_eq)=0 by bisection ----
          Scalar T = hpd;
          for (int k = 0; k < n_phos; ++k) {
            Scalar hc = hCs[k][index];
            Scalar q  = (hc > hc_small) ? hPps[k][index] / hc : (Scalar)0.0;
            T += hc * q;
          }
          Scalar lo = (Scalar)0.0, hi = T / h_this;       // Pd_eq in [0, all-dissolved]
          for (int it = 0; it < 60; ++it) {
            Scalar mid = (Scalar)0.5 * (lo + hi);
            Scalar f = h_this * mid - T;
            for (int k = 0; k < n_phos; ++k) f += hCs[k][index] * langmuir_qstar(p[k], mid);
            if (f > (Scalar)0.0) hi = mid; else lo = mid;
          }
          Scalar Pd_eq = (Scalar)0.5 * (lo + hi);
          hPd[index] = h_this * Pd_eq;
          for (int k = 0; k < n_phos; ++k)
            hPps[k][index] = hCs[k][index] * langmuir_qstar(p[k], Pd_eq);   // hPp_k = hC_k*q*(Pd_eq)
        } else {
          // ---- mode 0 kinetics: per group, sequential, dissolved-positivity capped ----
          for (int k = 0; k < n_phos; ++k) {
            Scalar hc = hCs[k][index];
            if (p[k].Qmax <= (Scalar)0.0 || hc <= hc_small) continue;       // inert / sediment-empty: no /hc
            Scalar q  = hPps[k][index] / hc;               // sorbed load [mg-P/g]
            Scalar Pd = hPd[index] / h_this;               // refresh (earlier groups changed hPd)
            Scalar dq = phos_sorption_rate(p[k], 0, Pd, q, dt);             // stiff-clamped, sign from (q*-q)
            Scalar dP = hc * dq;                           // depth-integrated P moved dissolved->particulate
            Scalar hpd_now = hPd[index];
            if (dP > hpd_now) dP = hpd_now;                // can't adsorb more dissolved P than present
            hPd[index]     = hpd_now - dP;
            hPps[k][index] = hPps[k][index] + dP;          // particulate gains exactly what dissolved loses
          }
        }
        index += blockDim.x * gridDim.x;
      }
    }

    void cuPhosphorusSorption(cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& hPd,
                              Scalar** hPps_dev, Scalar** hCs_dev, PhosphorusParams* pparams_dev,
                              int mode_id, int n_phos, Scalar dt) {
      cuPhosphorusSorptionKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        h.data.dev_ptr(), hPd.data.dev_ptr(), hPps_dev, hCs_dev, pparams_dev,
        mode_id, n_phos, dt, h.data.size());
    }

  }  // namespace fv
}  // namespace GC
