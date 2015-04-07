#ifndef LINPROG_H
#define LINPROG_H
#include <Eigen/Core>
// min  f'x
// s.t. A(    1:k,:) x <= b(1:k)
//      A(k+1:end,:) x = b(k+1:end)
//   ** x >= 0 **
//
// In contrast to other APIs the entries in b may be negative.
//
// maxit = -1 means inf
bool linprog_standard_form(
  const Eigen::VectorXd & c,
  const Eigen::MatrixXd & A,
  const Eigen::VectorXd & b,
  const int k,
  Eigen::VectorXd & f);

// min  f'x
// s.t. A x <= b
//      B x = c
//
bool linprog(
  const Eigen::VectorXd & f,
  const Eigen::MatrixXd & A,
  const Eigen::VectorXd & b,
  const Eigen::MatrixXd & B,
  const Eigen::VectorXd & c,
  Eigen::VectorXd & x);

#ifndef IGL_STATIC_LIBRARY
#  include "linprog.cpp"
#endif
#endif
