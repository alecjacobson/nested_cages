#ifndef GRADIENT_H
#define GRADIENT_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <igl/arap.h>


// add description here
void gradient_displacement(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat,
  Eigen::MatrixXd & grad);

// add description here
void gradient_surface_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U, 
  const igl::ARAPData & data,
  Eigen::MatrixXd & grad);

void gradient_volume(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & grad);

// add description here
void gradient(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const Eigen::MatrixXi & F,
  const igl::ARAPData & data,
  const char* Energy,
  Eigen::MatrixXd & grad);

#endif 