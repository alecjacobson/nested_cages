#ifndef LINPROG_H
#define LINPROG_H
#include <Eigen/Core>
// Solve a linear program given in "standard form"
//
// min  f'x
// s.t. A(    1:k,:) x <= b(1:k)
//      A(k+1:end,:) x = b(k+1:end)
//   ** x >= 0 **
//
// In contrast to other APIs the entries in b may be negative.
//
// Inputs:
//   c  #x list of linear coefficients
//   A  #A by #x matrix of linear constraint coefficients
//   b  #A list of linear constraint right-hand sides
//   k  number of inequality constraints as first rows of A,b
// Outputs:
//   x  #x solution vector
//
bool linprog(
  const Eigen::VectorXd & c,
  const Eigen::MatrixXd & A,
  const Eigen::VectorXd & b,
  const int k,
  Eigen::VectorXd & f);

// Wrapper in friendlier general form (no implicit bounds on x)
//
// min  f'x
// s.t. A x <= b
//      B x = c
//
// Inputs:
//   f  #x list of linear coefficients
//   A  #A by #x matrix of linear inequality constraint coefficients
//   b  #A list of linear constraint right-hand sides
//   B  #B by #x matrix of linear equality constraint coefficients
//   c  #B list of linear constraint right-hand sides
// Outputs:
//   x  #x solution vector
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
