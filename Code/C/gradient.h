#ifndef GRADIENT_H
#define GRADIENT_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <igl/arap.h>


// This file implements the gradient subroutine of the algorithm in "Nested Cages"
// [Sacht et al. 2015]. It consists of subroutines:
//   - gradient_displacement
//   - gradient_surface_arap
//   - gradient_volumetric_arap
//   - gradient_volume
//   - gradient

// computes the gradient of the squared norm of the difference between 2 meshes C and C_hat
//
// Inputs:
//   C  #C by 3 list of vertex positions of the 1st mesh
//   C_hat  #C by 3 list of vertex positions of the 2nd mesh
// Output:
//   grad   #C by 3 gradient of displacement energy
void gradient_displacement(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat,
  Eigen::MatrixXd & grad);

// Given rest-pose (surface) mesh (V,F) calculates the gradient of 
// the ARAP energy ('spokes-and-rims') of a deformed mesh (U,F)
//
// Inputs:
//   V  #V by 3 list of vertex positions of the rest pose mesh mesh
//   F  #F by 3 list of vertex indices into V
//   U  #V by 3 list of vertex positions of the deforming mesh
// Output:
//   grad   Gradient of surface ARAP energy
void gradient_surface_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U, 
  Eigen::MatrixXd & grad);

// Given rest-pose (tetrahedral) mesh (V,F) calculates the gradient of
// the ARAP energy ('elements') of a deformed mesh (U,F)
//
// Inputs:
//   V  #V by 3 list of vertex positions of the rest pose mesh mesh
//   F  #F by 4 list of vertex indices into V
//   U  #V by 3 list of vertex positions of the deforming mesh
// Output:
//   grad   Gradient of volumetric ARAP energy
void gradient_volumetric_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U,
  Eigen::MatrixXd & grad);

// Calculates the gradient of the volume of a mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of vertex positions of the rest pose mesh mesh
//   F  #F by 3 list of vertex indices into V
// Output:
//   grad   Gradient of volume at (V,F)
void gradient_volume(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & grad);

// Choose between gradients above depending on char* Energy and returns 
// energy value
//
// Inputs:
//   C  #C by 3 list of vertex positions of current coarse mesh
//   C_hat  #C by 3 list of vertex positions of initial coarse mesh
//   C_hat  #C by 3 list of vertex positions of previous coarse mesh
//   F  #F by 3 list of vertex indices into V
//   Energy  char specifying the energy to be calculated
// Output:
//   grad   #C by 3 gradient
//   bool   false if found error, true otherwise
bool gradient(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const Eigen::MatrixXi & F,
  const char* Energy,
  Eigen::MatrixXd & grad);

#endif 