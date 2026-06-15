// ======================================================================================
// T6: pluggable phosphorus (P) sorption-transport framework — config layer (step 1).
// LICENCE: GPLv3
// ======================================================================================
//
// Design authority: T6_design.md (3 forks: kinetics-default + equilibrium-switchable enum sorption /
// two-phase {reactive dissolved hPd + per-group particulate hPp_k mirroring T5 bed_k} / Qmax·K·EPC0·k_ads
// all config-driven). Mirrors the T5 sediment config discipline: NO physical number is hard-coded here —
// every phosphorus quantity is read from input/phosphorus_setup.dat. Missing file -> phosphorus_on=0 ->
// byte-level fallback to the T5 path (nothing in this header is compiled into the flood .so until the
// kernel-wiring steps; step 1 ships parser + POD + Python writer only).
//
// SCOPE (step 1, config layer ONLY): PhosphorusParams POD + PhosphorusConfig + host parser.
// NO kernel, NO closures, NO E-D touch, NO §4/§5/§7 implementation. Those are later commits (T6_design §11).
//
#ifndef CUDA_PHOSPHORUS_H
#define CUDA_PHOSPHORUS_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "Scalar.h"

namespace GC {

  // -------- per-group phosphorus sorption parameters (POD; usable on host & device) --------
  // Indices align 1:1 with the sediment groups (group k's phosphorus sorbs on sediment group k).
  // ALL defaults are 0 (an inert group): Qmax=0 => Langmuir q*=Qmax·K·Pd/(1+K·Pd)=0 => no sorption,
  // divide-safe (denominator 1+K·Pd >= 1, no singularity; any q/q* normalisation short-circuits on q*=0).
  // This makes "unspecified / coarse group" inert by construction (T6_design §3) AND keeps the kernel
  // free of any baked-in phosphorus number (T6 red line: zero hard-coded physics). Real Qmax/K/EPC0/
  // k_ads/q0 come from the config file (literature placeholders live in T6_literature.md, not here).
  struct PhosphorusParams {
    Scalar Qmax  = 0.0;   // Langmuir max sorption capacity [mg-P/g]  (0 => inert group, no phosphorus)
    Scalar K     = 0.0;   // Langmuir affinity constant [L/mg]
    Scalar EPC0  = 0.0;   // equilibrium P concentration of zero net sorption [mg/L]
    Scalar k_ads = 0.0;   // sorption rate constant [1/s] (kinetics mode_id=0)
    Scalar q0    = 0.0;   // initial sorbed load on this group's sediment [mg-P/g]
  };

  struct PhosphorusConfig {
    int    phosphorus_on    = 0;   // master switch (default OFF -> byte-level T5 fallback)
    int    sorption_mode_id = 0;   // 0=kinetics(default)  1=instant equilibrium  [global; §4]
    Scalar Pd_init          = 0.0; // dissolved-phase initial concentration [mg/L] (scalar; grid-file 'Pd' overrides)
    std::vector<PhosphorusParams> groups;   // index-aligned with sediment groups (n_phos == n_sed by design)
    int n_groups() const { return (int)groups.size(); }
  };

  // Read a header-directive scalar whose value may be on the SAME line ("$phosphorus_on 1")
  // OR on a following line ("$phosphorus_on\n1\n", the layout the Python writer emits, mirroring the
  // sediment file's look). Unlike T5's $morphology_on — a dormant flip-a-switch where the value being
  // ignored is harmless — phosphorus_on/sorption_mode are FUNCTIONAL, so the value must actually parse;
  // this helper reads it robustly either way and skips blank/comment-only lines.
  template <typename T>
  inline void read_directive_value(std::istringstream& ss, std::ifstream& fin, T& out) {
    if (ss >> out) return;                         // value on the same line
    std::string vline;
    while (std::getline(fin, vline)) {             // value on a following line
      std::size_t h = vline.find('#');
      if (h != std::string::npos) vline = vline.substr(0, h);
      std::istringstream vs(vline);
      if (vs >> out) return;
      // blank / comment-only line: keep scanning
    }
  }

  // -------- host-side parser for input/phosphorus_setup.dat (sectioned key=value) --------
  // Mirrors read_sediment_setup: unknown keys ignored (forward-compatible); absent file ->
  // phosphorus_on=0, no groups -> case runs exactly as T5 (byte-level fallback).
  inline PhosphorusConfig read_phosphorus_setup(const char* path) {
    PhosphorusConfig cfg;
    std::ifstream fin(path);
    if (!fin.is_open()) return cfg;                // no phosphorus for this case
    int cur = -1;                                  // current $pgroup index, -1 if none
    bool in_dissolved = false;                     // inside the $dissolved section
    std::string line;
    while (std::getline(fin, line)) {
      // strip inline comments (#) and trailing whitespace
      std::size_t hash = line.find('#');
      if (hash != std::string::npos) line = line.substr(0, hash);
      std::istringstream ss(line);
      std::string key; if (!(ss >> key)) continue; // blank line
      if      (key == "$phosphorus_on") { read_directive_value(ss, fin, cfg.phosphorus_on);    cur = -1; in_dissolved = false; }
      else if (key == "$sorption_mode") { read_directive_value(ss, fin, cfg.sorption_mode_id); cur = -1; in_dissolved = false; }
      else if (key == "$pgroup")        { cfg.groups.emplace_back(); cur = (int)cfg.groups.size() - 1; in_dissolved = false; }
      else if (key == "$dissolved")     { in_dissolved = true; cur = -1; }
      else if (in_dissolved) {
        if (key == "Pd_init") ss >> cfg.Pd_init;
        // unknown keys ignored (forward-compatible)
      }
      else if (cur >= 0) {
        PhosphorusParams& p = cfg.groups[cur];
        if      (key == "Qmax")  ss >> p.Qmax;
        else if (key == "K")     ss >> p.K;
        else if (key == "EPC0")  ss >> p.EPC0;
        else if (key == "k_ads") ss >> p.k_ads;
        else if (key == "q0")    ss >> p.q0;
        // unknown keys ignored (forward-compatible)
      }
    }
    return cfg;
  }

}  // namespace GC

#endif
