#ifndef EDGE_COLLAPSE_IS_VALID_H
#define EDGE_COLLAPSE_IS_VALID_H
#include "collapse_edge.h"
#include <Eigen/Core>
// Assumes (V,F) is a closed manifold mesh (except for previouslly collapsed
// faces which should be set to: 
// [COLLAPSE_EDGE_NULL COLLAPSE_EDGE_NULL COLLAPSE_EDGE_NULL].
// Tests whether 
// Collapsing exactly two faces and exactly 3 edges from E (e and one side of
// each face gets collapsed to the other) will result in a mesh with the same
// topology.
//
// Inputs:
//   e  index into E of edge to try to collapse. E(e,:) = [s d] or [d s] so
//     that s<d, then d is collapsed to s.
//   F  #F by 3 list of face indices into V.
//   E  #E by 2 list of edge indices into V.
//   EMAP #F*3 list of indices into E, mapping each directed edge to unique
//     unique edge in E
//   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
//     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
//     e=(j->i)
//   EI  #E by 2 list of edge flap corners (see above).
// Returns true if edge collapse is valid
bool edge_collapse_is_valid(
  const int e,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI);
#ifndef IGL_STATIC_LIBRARY
#  include "edge_collapse_is_valid.cpp"
#endif
#endif
