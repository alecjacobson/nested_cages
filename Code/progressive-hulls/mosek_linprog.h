#ifndef MOSEK_LINPROG_H
#define MOSEK_LINPROG_H
#include <Eigen/Core>
#include <mosek.h>
// min c'x
// s.t. lc <= A x <= uc
//      lx <= x <= ux
//
bool mosek_linprog(
  const Eigen::VectorXd & c,
  const Eigen::SparseMatrix<double> & A,
  const Eigen::VectorXd & lc,
  const Eigen::VectorXd & uc,
  const Eigen::VectorXd & lx,
  const Eigen::VectorXd & ux,
  Eigen::VectorXd & x);
bool mosek_linprog(
  const Eigen::VectorXd & c,
  const Eigen::SparseMatrix<double> & A,
  const Eigen::VectorXd & lc,
  const Eigen::VectorXd & uc,
  const Eigen::VectorXd & lx,
  const Eigen::VectorXd & ux,
  const MSKenv_t & env,
  Eigen::VectorXd & x);

#ifndef IGL_STATIC_LIBRARY
#include "mosek_linprog.cpp"
#endif
#endif
