#ifndef REINFLATE_H
#define REINFLATE_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <stack>

// Given a history of fine flown meshes H, reinflates the fine
// mesh pushing the coarse mesh (C_hat,F_hat) away to new positions (C,F_hat).
// 
// Inputs:
//   H  history of flow mesh positions
//   T  connectivity of the fine meshes (same for all meshes in H)
//   C_hat  #C_hat by 3 list of coarse mesh initial positions
//   F_hat  #F_hat by 3 list of coarse mesh triangle indices C_hat
//   EnergyInflation  either "None", or "DispInitial", or "DispStep" or
//           "Volume" or "SurfARAP" or "VolARAP"
//   EnergyFinal  either "None", or "DispInitial", or "DispStep" or
//           "Volume" or "SurfARAP" or "VolARAP"
// Output:
//  C  #C by 3 list of final coarse mesh vertex positions (optimal cage)
void reinflate(
  std::stack<Eigen::MatrixXd> & H, 
  const Eigen::MatrixXi & T, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXi & F_hat, 
  const char* EnergyInflation,
  const char* EnergyFinal,
  Eigen::MatrixXd & C);

#endif 