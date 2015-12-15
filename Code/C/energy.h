#ifndef ENERGY_H
#define ENERGY_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

// add description
double energy_displacement(
  const Eigen::MatrixXi & C, 
  const Eigen::MatrixXd & C_hat);

// add description
double energy(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const char* Energy);

#endif 