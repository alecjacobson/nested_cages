#ifndef COLLAPSE_EDGE_H
#define COLLAPSE_EDGE_H
#include <Eigen/Core>
// Assumes (V,F) is a closed manifold mesh (except for previouslly collapsed
// faces which should be set to: 
// [COLLAPSE_EDGE_NULL COLLAPSE_EDGE_NULL COLLAPSE_EDGE_NULL].
// Collapses exactly two faces and exactly 3 edges from E (e and one side of
// each face gets collapsed to the other). This is implemented in a way that it
// can be repeatedly called until satisfaction and then the garbage in F can be
// collected by removing NULL faces.
//
// Inputs:
//   e  index into E of edge to try to collapse. E(e,:) = [s d] or [d s] so
//     that s<d, then d is collapsed to s.
///  p  dim list of vertex position where to place merged vertex
// Inputs/Outputs:
//   V  #V by dim list of vertex positions, lesser index of E(e,:) will be set
//     to midpoint of edge.
//   F  #F by 3 list of face indices into V.
//   E  #E by 2 list of edge indices into V.
//   EMAP #F*3 list of indices into E, mapping each directed edge to unique
//     unique edge in E
//   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
//     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
//     e=(j->i)
//   EI  #E by 2 list of edge flap corners (see above).
//   e1  index into E of edge collpased on left
//   e2  index into E of edge collpased on left
//   f1  index into E of edge collpased on left
//   f2  index into E of edge collpased on left
// Returns true if edge was collapsed
#define COLLAPSE_EDGE_NULL 0
bool collapse_edge(
  const int e,
  const Eigen::RowVectorXd & p,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI,
  int & e1,
  int & e2,
  int & f1,
  int & f2);
bool collapse_edge(
  const int e,
  const Eigen::RowVectorXd & p,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXi & E,
  Eigen::VectorXi & EMAP,
  Eigen::MatrixXi & EF,
  Eigen::MatrixXi & EI);
#ifndef IGL_STATIC_LIBRARY
#  include "collapse_edge.cpp"
#endif
#endif
