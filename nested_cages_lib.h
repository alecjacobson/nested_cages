#ifndef NESTED_CAGES_LIB_H
#define NESTED_CAGES_LIB_H

// Reusable, mesh-in / mesh-out implementation of the "Nested Cages"
// [Sacht et al. 2015] pipeline. This is the exact per-level pipeline that used
// to live inside main() in nested_cages.cpp, factored out so it can be reused
// by both the command-line tool and the Python bindings.

#include <Eigen/Core>
#include <functional>
#include <string>
#include <vector>

// Compute a sequence of nested cages for the input mesh (V0,F0).
//
// Inputs:
//   V0  #V0 by 3 list of input (fine) mesh vertex positions
//   F0  #F0 by 3 list of input mesh triangle indices into V0
//   quad_order  quadrature order for the shrinking flow (1, 2 or 3)
//   num_faces   #levels list, target number of faces of each cage
//                 (only used for a level whose V_decim entry is empty)
//   regular     #levels list, true => regular decimation, false => adaptive
//                 (only used for a level whose V_decim entry is empty)
//   V_decim     #levels list of optional input decimations (vertices). If
//                 V_decim[i] has rows, it is used directly as the i-th input
//                 decimation instead of decimating the previous cage.
//   F_decim     #levels list of optional input decimations (faces).
//   energy_expansion  energy minimized during re-inflation:
//                 "None","DispStep","DispInitial","Volume","SurfARAP","VolARAP"
//   energy_final      energy minimized after re-inflation (additional
//                 processing), same set of choices.
// Outputs:
//   Vout  #levels list of cage vertex positions (one entry per level)
//   Fout  #levels list of cage face indices (one entry per level)
//
// Optional:
//   per_level_callback  if set, called with (level_index, C, F) immediately
//     after each cage is successfully computed. This lets callers (e.g. the
//     command-line tool) stream results to disk as soon as each level finishes,
//     exactly as the original program did.
//
// Returns true on success. Returns false on any per-level failure (an input
// decimation could not be repaired, the flow failed to take the fine mesh
// inside the coarse mesh within 1000 iterations, or a sanity check found the
// output cage self-intersecting or intersecting the previous layer) -- these
// are exactly the conditions that used to make main() return early.
bool nested_cages(
  const Eigen::MatrixXd & V0,
  const Eigen::MatrixXi & F0,
  const int quad_order,
  const std::vector<int> & num_faces,
  const std::vector<bool> & regular,
  const std::vector<Eigen::MatrixXd> & V_decim,
  const std::vector<Eigen::MatrixXi> & F_decim,
  const std::string & energy_expansion,
  const std::string & energy_final,
  std::vector<Eigen::MatrixXd> & Vout,
  std::vector<Eigen::MatrixXi> & Fout,
  const std::function<void(int, const Eigen::MatrixXd &, const Eigen::MatrixXi &)>
    & per_level_callback = nullptr);

#endif
