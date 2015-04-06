#ifndef CIRCULATION_H
#define CIRCULATION_H
#include <Eigen/Core>
#include <vector>

// Return list of faces around the end point of an edge.
//
// Inputs:
//   e  index into E of edge to circulate
//   ccw  whether to _continue_ in ccw direction of edge (circulate around
//     E(e,1))
//   F  #F by 3 list of face indices
//   E  #E by 2 list of edge indices
//   EMAP #F*3 list of indices into E, mapping each directed edge to unique
//     unique edge in E
//   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
//     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
//     e=(j->i)
//   EI  #E by 2 list of edge flap corners (see above).
// Returns list of faces touched by circulation.
//   
std::vector<int> circulation(
  const int e,
  const bool ccw,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI);
// Wrapper with VectorXi output.
void circulation(
  const int e,
  const bool ccw,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & EMAP,
  const Eigen::MatrixXi & EF,
  const Eigen::MatrixXi & EI,
  Eigen::VectorXi & vN);

#ifndef IGL_STATIC_LIBRARY
#  include "circulation.cpp"
#endif
#endif
