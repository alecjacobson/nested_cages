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
#include <stack>

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
  const Eigen::VectorXd & area_0, 
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

// Inputs:
//   V  #V by 3 list of fine mesh positions
//   F  #F by 3 list of fine mesh triangle indices into V
//   V_coarse  #V_coarse by 3 list of coarse mesh positions
//   F_coarse  #F_coarse by 3 list of coarse mesh triangle indices into
//     V_coarse
//   A_qv   #V0 by #A_qv matrix taking quadrature points to gradients at vertices.
//     For quad_order=1,2,3 then #A_qv=#F,3*#F,4*#F
//   M   mass matrix of the input fine matrix
// Outputs:
//   grad  #V by 3 list of gradient coordinates
void grad_energy(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv, 
  const Eigen::SparseMatrix<double> & M_inv, 
  Eigen::MatrixXd & grad);

// Calculates the diameter of the uninon of the vertex positions
// into V and V_coarse.
//
// Inputs:
//   V  #V by 3 list of fine mesh positions
//   V_coarse  #V_coarse by 3 list of coarse mesh positions
// Outputs:
//   diameter  
double diameter(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXd & V_coarse);

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
//   M_inv inverse of the mass matrix of the fine mesh
// Output:
//  V_new  #V by 3 list of new fine mesh vertex positions
//   
void flow_one_step(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv, 
  const double delta_t, 
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
//   H  history of flow mesh positions
//   a boolean, that says whether the flkow failed or succeeded
//      after 1000 iterations
bool flow_fine_inside_coarse(
  const Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  const Eigen::MatrixXd & V_coarse, 
  const Eigen::MatrixXi & F_coarse, 
  const Eigen::SparseMatrix<double> & A_qv,
  std::stack<Eigen::MatrixXd> & H);
#endif 
