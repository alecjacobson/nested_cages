#ifndef EDGE_FLAPS_H
#define EDGE_FLAPS_H
#include <Eigen/Core>
// Determine "edge flaps": two faces on either side of a unique edge (assumes
// edge-manifold mesh)
//
// Inputs:
//   F  #F by 3 list of face indices
//   E  #E by 2 list of edge indices into V.
//   EMAP #F*3 list of indices into E, mapping each directed edge to unique
//     unique edge in E
// Outputs:
//   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
//     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
//     e=(j->i)
//   EI  #E by 2 list of edge flap corners (see above).
void edge_flaps(
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI);
// Only faces as input
void edge_flaps(
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI);
#ifndef IGL_STATIC_LIBRARY
#  include "edge_flaps.cpp"
#endif

#endif
