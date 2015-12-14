#ifndef FILTER_H
#define FILTER_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

// Inputs:
//   Vf  fine mesh positions before stepping
//   T   fine mesh vertex indices into Vf 
//   Uf  Desired velocities of the fine mesh
//   C  coarse mesh positions before stepping
//   F_hat   coarse mesh vertex indices into C 
//   Uc  Desired velocities of the fine mesh
// Output:
//  Uc  Corrected coarse mesh velocities
void filter(
  const Eigen::MatrixXd & Vf, 
  const Eigen::MatrixXi & T, 
  const Eigen::MatrixXd & Uf, 
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXi & F_coarse, 
  Eigen::MatrixXd & Uc);

#endif 