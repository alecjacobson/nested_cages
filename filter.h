#ifndef FILTER_H
#define FILTER_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

// This file implements the filter subroutine of the algorithm in "Nested Cages"
// [Sacht et al. 2015]. It consists of subroutines:
//   - filter

void velocity_filter_ACM(
  const Eigen::MatrixXd & V0, 
  Eigen::MatrixXd & V1, 
  const Eigen::MatrixXi & F0, 
  double outer,
  double inner,
  int numinfinite);

double inflate_ACM(
  Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  double eps_distance,
  int numinfinite);

// Inputs:
//   Vf  fine mesh positions before stepping
//   T   fine mesh vertex indices into Vf 
//   Uf  Desired velocities of the fine mesh
//   C  coarse mesh positions before stepping
//   F_hat   coarse mesh vertex indices into C 
//   eps_prox   minimum separation distance
//   Uc  Desired velocities of the fine mesh
// Output:
//  Uc  Corrected coarse mesh velocities
void filter(
  const Eigen::MatrixXd & Vf, 
  const Eigen::MatrixXi & T, 
  const Eigen::MatrixXd & Uf, 
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXi & F_coarse, 
  double eps_prox,
  Eigen::MatrixXd & Uc);

#endif 