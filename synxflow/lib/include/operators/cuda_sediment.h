// ======================================================================================
// T5: pluggable cohesive-sediment transport framework — config + closures + operators.
// LICENCE: GPLv3
// ======================================================================================
//
// Design authority: T5_design.md (4 forks: enum-dispatch closures / sectioned sediment_setup.dat
// + Python writer / vector<field>+fused multi-group stencil / suspended-only with morphology
// flip-a-switch extension point). NO physical number is hard-coded here — all read from config.
//
#ifndef CUDA_SEDIMENT_H
#define CUDA_SEDIMENT_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "Scalar.h"
#include "cuda_mapped_field.h"
#include "Flag.h"
#include "Vector.h"

namespace GC {

  // Compile-time upper bound on sediment groups (sizes the WithSediment kernel's per-cell stack
  // arrays). Runtime n_sed must be <= this; the flood solver guards it. Raise + rebuild for more
  // grain-size classes (only the n_sed>0 WithSediment kernel pays the register/stack cost).
  static const int SED_MAX_GROUPS = 8;

  // -------- per-group parameters (POD; usable on host & device) --------
  // All values come from input/sediment_setup.dat. Defaults are literature PLACEHOLDERS
  // (so a missing key is sane), NOT modelling decisions — real values go in the config file.
  struct SedimentParams {
    int    settling_id = 0;      // 0=const(w_s)  1=Zhang(张瑞瑾)  2=Stokes  3=floc
    int    closure_id  = 0;      // 0=Partheniades two-threshold (cohesive)  [1=Shields reserved, unimplemented]
    int    cohesive    = 1;      // 1=cohesive (default path)  0=non-cohesive (reserved)
    Scalar rho_s       = 2650.0; // grain density [kg/m3]
    Scalar porosity    = 0.4;    // bed porosity [-]
    Scalar D50         = 1.0e-5; // median grain diameter [m]
    Scalar w_s         = 0.0;    // settling velocity [m/s] (const, or fallback for formula)
    Scalar tau_ce      = 0.0;    // critical shear for erosion [Pa]
    Scalar tau_cd      = 0.0;    // critical shear for deposition [Pa]
    Scalar M           = 0.0;    // Partheniades erosion coefficient [kg/m2/s]
    Scalar bed_init    = 0.0;    // initial erodible bed mass store [kg/m2] (uniform seed; 0=no bed)
  };

  struct SedimentConfig {
    int morphology_on = 0;       // flip-a-switch extension point (default OFF; T5 never applies Δz_accum to z)
    std::vector<SedimentParams> groups;
    int n_groups() const { return (int)groups.size(); }
  };

  // -------- host-side parser for input/sediment_setup.dat (sectioned key=value) --------
  // Returns an empty config (n_groups==0) if the file is absent -> case runs with no sediment.
  inline SedimentConfig read_sediment_setup(const char* path) {
    SedimentConfig cfg;
    std::ifstream fin(path);
    if (!fin.is_open()) return cfg;             // no sediment for this case
    int declared = 0, cur = -1;
    std::string line;
    while (std::getline(fin, line)) {
      // strip inline comments (#) and trailing whitespace
      std::size_t hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      std::istringstream ss(line);
      std::string key; if (!(ss >> key)) continue;     // blank line
      if (key == "$n_groups")       { ss >> declared; }
      else if (key == "$morphology_on") { ss >> cfg.morphology_on; }
      else if (key == "$group")     { cfg.groups.emplace_back(); cur = (int)cfg.groups.size() - 1; }
      else if (cur >= 0) {
        SedimentParams& p = cfg.groups[cur];
        if      (key == "settling_id") ss >> p.settling_id;
        else if (key == "closure_id")  ss >> p.closure_id;
        else if (key == "cohesive")    ss >> p.cohesive;
        else if (key == "rho_s")       ss >> p.rho_s;
        else if (key == "porosity")    ss >> p.porosity;
        else if (key == "D50")         ss >> p.D50;
        else if (key == "w_s")         ss >> p.w_s;
        else if (key == "tau_ce")      ss >> p.tau_ce;
        else if (key == "tau_cd")      ss >> p.tau_cd;
        else if (key == "M")           ss >> p.M;
        else if (key == "bed_init")    ss >> p.bed_init;
        // unknown keys ignored (forward-compatible)
      }
    }
    (void)declared;   // $n_groups is a hint; actual count = number of $group sections parsed
    return cfg;
  }

  namespace fv {
    // T5 step3 operators (defined in src/operators/cuda_sediment.cu). WithSediment path only.
    // tau_b = rho_w*g*n^2*|u|^2/h^(1/3) ; dry cell -> 0.
    void cuSedimentBedShear(cuFvMappedField<Scalar, on_cell>& manning, cuFvMappedField<Scalar, on_cell>& gravity,
                            cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU,
                            cuFvMappedField<Scalar, on_cell>& tau_b);
    // Per-group E-D reaction on top of advection: hC_k += ED*dt, bed_k -= ED*dt (R5-clamped),
    // dz_accum += (-ED*dt)/(rho_s*(1-porosity)) COMPUTED-but-unapplied (z never touched in T5).
    // params_dev = device array of n_sed SedimentParams; hCs_dev/beds_dev = device Scalar* arrays.
    void cuSedimentErosionDeposition(cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& tau_b,
                                     Scalar** hCs_dev, Scalar** beds_dev, Scalar* dz_accum_dev,
                                     SedimentParams* params_dev, int n_sed, Scalar dt);
  }  // namespace fv

}  // namespace GC

#endif
