#include "reinflate.h"
#include "filter.h"
#include <stdio.h>

// Need to include some IGL header to have igl namespace
#include <igl/signed_distance.h>

void reinflate(
  std::stack<Eigen::MatrixXd> & H, 
  const Eigen::MatrixXi & T, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXi & F_hat, 
  const char* EnergyInflation,
  const char* EnergyFinal,
  Eigen::MatrixXd & C)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;

  MatrixXd F = H.top();
  H.pop();
  C = C_hat;

  while (!H.empty()){

  	double beta = 1e-2;
  	// First, find a feasible state (before stepping and projecting)
  	MatrixXd Uc = MatrixXd::Zero(C_hat.rows(), 3);
  	MatrixXd Uf = H.top()-F;
  	H.pop();
  	filter(F,T,Uf,C,F_hat,Uc);
    // step size search

    // update
    C = C+Uc;

  }

}