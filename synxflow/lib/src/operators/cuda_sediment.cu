// ======================================================================================
// T5 step3: sediment source/sink operators — bed shear tau_b + erosion/deposition (E-D).
// LICENCE: GPLv3
// ======================================================================================
//
// Design authority: T5_design.md §8 (tau_b coupling + E-D + R5 + bed conservation) and §4
// (morphology = flip-a-switch: bed_k store + Delta-z_accum computed-but-NOT-applied; z never touched).
// Physics formulas live in cuda_sediment_closures.h (enum-dispatch __device__ inline). NO sediment
// parameter is hard-coded here — all come from per-group SedimentParams uploaded to device.
// These operators run ONLY on the n_sed>0 (WithSediment) path; the n_sed=0 passive path never calls them.
//
#include "cuda_sediment.h"
#include "cuda_sediment_closures.h"
#include "cuda_kernel_launch_parameters.h"
#include "Scalar.h"
#include "Vector.h"
#include "cuda_mapped_field.h"

namespace GC {
  namespace fv {

    // ---- tau_b = rho_w * g * n^2 * |u|^2 / h^(1/3)  [Pa] ; |u|=|hU|/h, dry cell -> 0 ----
    __global__ void cuSedimentBedShearKernel(Scalar* manning, Scalar* gravity, Scalar* h, Vector* hU,
                                             Scalar* tau_b, unsigned int phi_size) {
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      const Scalar h_small = 1e-10;
      const Scalar rho_w = 1000.0;          // water density [kg/m3] (universal constant, not a basin param)
      while (index < phi_size) {
        Scalar h_this = h[index];
        if (h_this < h_small) {
          tau_b[index] = (Scalar)0.0;
        } else {
          Vector2 hu = hU[index];
          Scalar u2 = (hu.x*hu.x + hu.y*hu.y) / (h_this*h_this);
          Scalar n = manning[index];
          Scalar g = gravity[index];
          tau_b[index] = rho_w * g * n * n * u2 / cbrt(h_this);
        }
        index += blockDim.x * gridDim.x;
      }
    }

    void cuSedimentBedShear(cuFvMappedField<Scalar, on_cell>& manning, cuFvMappedField<Scalar, on_cell>& gravity,
                            cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU,
                            cuFvMappedField<Scalar, on_cell>& tau_b) {
      cuSedimentBedShearKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        manning.data.dev_ptr(), gravity.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(),
        tau_b.data.dev_ptr(), h.data.size());
    }

    // ---- E-D source/sink per group (advection already advanced hC; here add the reaction) ----
    // For each cell, each group k:
    //   C  = hC_k/h ; w_s = settling(params_k) ; ED = E - D = sed_ED_rate(params_k, tau_b, C, w_s)
    //   R5 positivity: ED clamped to [ -hC_k/dt , bed_k/dt ]  (can't deposit below 0, can't erode empty bed)
    //   hC_k += ED*dt ; bed_k -= ED*dt          (suspended + bed mass conserved: Sum(hC_k)+Sum(bed_k) const)
    //   Delta-z_accum += (-ED*dt)/(rho_s*(1-porosity))   (COMPUTED, never applied to z in T5 — §4)
    __global__ void cuSedimentErosionDepositionKernel(Scalar* h, Scalar* tau_b, Scalar** hCs, Scalar** beds,
                                                      Scalar* dz_accum, SedimentParams* params, int n_sed,
                                                      Scalar dt, unsigned int phi_size) {
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      const Scalar h_small = 1e-10;
      while (index < phi_size) {
        Scalar h_this = h[index];
        Scalar tb = tau_b[index];
        Scalar dz_cell = (Scalar)0.0;
        for (int k = 0; k < n_sed; ++k) {
          Scalar hc  = hCs[k][index];
          Scalar bed = beds[k][index];
          Scalar C   = (h_this >= h_small) ? hc / h_this : (Scalar)0.0;
          Scalar ws  = sed_settling_velocity(params[k]);
          Scalar ed  = sed_ED_rate(params[k], tb, C, ws);
          // R5 positivity clamp (dt>0): deposition limited by suspended mass, erosion by bed store
          Scalar ed_min = -hc / dt;          // ED >= -hC/dt  -> hC + ED*dt >= 0
          Scalar ed_max = bed / dt;          // ED <=  bed/dt -> bed - ED*dt >= 0
          if (ed < ed_min) ed = ed_min;
          if (ed > ed_max) ed = ed_max;
          Scalar dmass = ed * dt;            // mass moved bed->suspension (>0 erosion, <0 deposition)
          hCs[k][index]  = hc + dmass;
          beds[k][index] = bed - dmass;
          Scalar denom = params[k].rho_s * ((Scalar)1.0 - params[k].porosity);
          if (denom > (Scalar)0.0) dz_cell += (-dmass) / denom;   // bed-elevation increment (deposition -> +)
        }
        dz_accum[index] += dz_cell;          // running accumulator over time; NOT applied to z (extension pt)
        index += blockDim.x * gridDim.x;
      }
    }

    void cuSedimentErosionDeposition(cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& tau_b,
                                     Scalar** hCs_dev, Scalar** beds_dev, Scalar* dz_accum_dev,
                                     SedimentParams* params_dev, int n_sed, Scalar dt) {
      cuSedimentErosionDepositionKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        h.data.dev_ptr(), tau_b.data.dev_ptr(), hCs_dev, beds_dev, dz_accum_dev,
        params_dev, n_sed, dt, h.data.size());
    }

  }  // namespace fv
}  // namespace GC
