#ifndef GRADIENT_H
#define GRADIENT_H

#include <Eigen/Core>
#include <Eigen/SparseCore>


// add description here
void gradient_displacement(
  const Eigen::MatrixXi & C, 
  const Eigen::MatrixXd & C_hat,
  Eigen::MatrixXd & grad);

// add description here
void gradient(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const char* Energy,
  Eigen::MatrixXd & grad);

#endif 