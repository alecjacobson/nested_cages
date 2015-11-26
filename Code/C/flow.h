#ifndef FLOW_H
#define FLOW_H

// This file implements the flow subroutine of the algorithm in "Nested Cages"
// [Sacht et al. 2015]. It consists of subroutines:
//   - gradQ_to_gradV
//   - signed_distance_direction
//   - grad_energy
//   - flow_one_step
//   - flow_fine_inside_coarse
//

#include <Eigen/Core>
#include <Eigen/SparseCore>

// compute the matrix that converts gradient at quadrature 
// points to gradient at mesh vertices 
//
// Inputs:
//   V0  #V0 by 3 list of mesh vertex positions
//   F0  #F0 by 3 list of triangle indices into V0
//   area_0  #F0 list of triangle areas (precomputed)
//   quad_order  quadrature rule order: 1, 2, or 3
// Output:
//   A   #V0 by #A matrix taking quadrature points to gradients at vertices.
//     For quad_order=1,2,3 then #A=#F,3*#F,4*#F
//
void gradQ_to_gradV(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  const Eigen::MatrixXd & area_0, 
  const int quad_order,
  Eigen::SparseMatrix<double> & A);

// For a set of points P compute the direction of decrease according to the signed
// distance field of a mesh (V,F).
//
// Inputs:
//   P  #P by 3 list of points
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of triangle indices into V
// Outputs:
//   D  #P by 3 list of directions
//
void signed_distance_direction(
  const Eigen::MatrixXd & P, 
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  Eigen::MatrixXd & D);

// **Alec: I don't actually understand what this ones computing... Why is the
// output a vector and not a matrix?**
//
// Inputs:
//   V  #V by 3 list of fine mesh positions
//   F  #F by 3 list of fine mesh triangle indices into V
//   V_coarse  #V_coarse by 3 list of coarse mesh positions
//   F_coarse  #F_coarse by 3 list of coarse mesh triangle indices into
//     V_coarse
//   A   #V0 by #A matrix taking quadrature points to gradients at vertices.
//     For quad_order=1,2,3 then #A=#F,3*#F,4*#F
// Outputs:
//   grad  ?
void grad_energy(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv, 
  Eigen::MatrixXd & grad);

// Move the shrinking fine mesh one step along its flow inside the coarse mesh
// using a specified quadrature order.
//
// Inputs:
//   V  #V by 3 Previous positions of fine mesh vertex positions
//   F  #F by 3 fine mesh triangle indices into V
//   V_coarse  #V_coarse by 3 list of coarse mesh positions
//   F_coarse  #F_coarse by 3 list of coarse mesh triangle indices into
//   A   #V0 by #A matrix taking quadrature points to gradients at vertices.
//     For quad_order=1,2,3 then #A=#F,3*#F,4*#F
// Output:
//  V_new  #V by 3 list of new fine mesh vertex positions
//   
void flow_one_step(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv, 
  Eigen::MatrixXd & V_new);

// Flow the fine mesh inside the coarse mesh using a specified quadrature order
// 
// Inputs:
//   V0  #V0 by 3 Previous positions of fine mesh vertex positions
//   F0  #F0 by 3 fine mesh triangle indices into V0
//   V_coarse  #V_coarse by 3 list of coarse mesh positions
//   F_coarse  #F_coarse by 3 list of coarse mesh triangle indices into
//   A   #V0 by #A matrix taking quadrature points to gradients at vertices.
//     For quad_order=1,2,3 then #A=#F0,3*#F,4*#F0
// Output:
//  V  #V0 by 3 list of new fine mesh vertex positions
void flow_fine_inside_coarse(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv,
  Eigen::MatrixXd & V);
#endif 
