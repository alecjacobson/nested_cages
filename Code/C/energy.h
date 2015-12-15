#ifndef ENERGY_H
#define ENERGY_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

// add description
double energy_displacement(
  const Eigen::MatrixXi & C, 
  const Eigen::MatrixXd & C_hat);

// add description
double energy_surface_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U);

// add description
double energy_volume(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F);

// add description
double energy(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const Eigen::MatrixXi & F,
  const char* Energy);

#endif 