#include "energy.h"
#include "gradient.h"
#include <stdio.h>
#include <igl/centroid.h>

double energy_displacement(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat)
{
  // gradient direction is the difference between current and reference meshes
	return ((C-C_hat).rowwise().norm()).norm();
}

double energy_surface_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U)
{

	using namespace Eigen;
  using namespace std;
  using namespace igl;

  // rest-pose mesh
  MatrixXd ref_V = V;
  MatrixXi ref_F = F;

  // precomputation of cotangent matrix of the rest-pose
  SparseMatrix<double> data_L;
  cotmatrix(ref_V,ref_F,data_L);

  // precomputation of stack of covariance matrices
  // note: for SurfARAP we use 'spokes_and_rims' energy
  SparseMatrix<double> data_CSM;
  covariance_scatter_matrix(ref_V,ref_F,ARAP_ENERGY_TYPE_SPOKES_AND_RIMS,data_CSM);

  // precomputation of RHS for ARAP solve
  SparseMatrix<double> data_K;
  arap_rhs(ref_V,ref_F,3,ARAP_ENERGY_TYPE_SPOKES_AND_RIMS,data_K);

  // repeat positions and mulitply by covraiance matrices
  MatrixXd U_rep;
  repmat(U,3,1,U_rep);
  MatrixXd S = data_CSM*U_rep;

  // rotation matrices (side-by-side)
  MatrixXd R;
  fit_rotations(S,false,R);

  // columnize R
  VectorXd Rcol;
  columnize(R,R.cols()/3,2,Rcol);

  // multiplication of R by RHS of ARAP
  MatrixXd dV;
  dV = data_K*Rcol;
  MatrixXd dV3(U.rows(),U.cols());
  dV3.col(0) = dV.block(0,0,U.rows(),1);
  dV3.col(1) = dV.block(U.rows(),0,U.rows(),1);
  dV3.col(2) = dV.block(2*U.rows(),0,U.rows(),1);

  // energy
  return (-U.transpose()*(0.5*data_L)*U - U.transpose()*dV3 - V.transpose()*(0.5*data_L)*V).trace();

}

double energy_volumetric_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U)
{

  using namespace Eigen;
  using namespace std;
  using namespace igl;

  // rest-pose mesh
  MatrixXd ref_V = V;
  MatrixXi ref_F = F;

  // precomputation of cotangent matrix of the rest-pose
  SparseMatrix<double> data_L;
  cotmatrix(ref_V,ref_F,data_L);

  // precomputation of stack of covariance matrices
  // note: for SurfARAP we use 'elements' energy
  SparseMatrix<double> data_CSM;
  covariance_scatter_matrix(ref_V,ref_F,ARAP_ENERGY_TYPE_ELEMENTS,data_CSM);

  // precomputation of RHS for ARAP solve
  SparseMatrix<double> data_K;
  arap_rhs(ref_V,ref_F,3,ARAP_ENERGY_TYPE_ELEMENTS,data_K);

  // repeat positions and mulitply by covraiance matrices
  MatrixXd U_rep;
  repmat(U,3,1,U_rep);
  MatrixXd S = data_CSM*U_rep;

  // rotation matrices (side-by-side)
  MatrixXd R;
  fit_rotations(S,false,R);

  // columnize R
  VectorXd Rcol;
  columnize(R,R.cols()/3,2,Rcol);

  // multiplication of R by RHS of ARAP
  MatrixXd dV;
  dV = data_K*Rcol;
  MatrixXd dV3(U.rows(),U.cols());
  dV3.col(0) = dV.block(0,0,U.rows(),1);
  dV3.col(1) = dV.block(U.rows(),0,U.rows(),1);
  dV3.col(2) = dV.block(2*U.rows(),0,U.rows(),1);

  // energy
  return (-U.transpose()*(0.5*data_L)*U - U.transpose()*dV3 - V.transpose()*(0.5*data_L)*V).trace();
}

double energy_volume(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F)
{

	using namespace Eigen;
  using namespace igl;

  VectorXd c(3,1);
  double volume;

  // call IGLS's function to calculate volume
  centroid(V,F,c,volume);

	return volume;

}

double energy(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev,
  const Eigen::MatrixXi & F, 
  const char* Energy)
{
  using namespace std;
  // call gradient depending on given energy
	if (strcmp(Energy,"DispStep")==0)
  {
		return energy_displacement(C,C_prev);
	}
	else if (strcmp(Energy,"DispInitial")==0)
	{
		return energy_displacement(C,C_hat);
	}
	else if (strcmp(Energy,"SurfARAP")==0)
	{
		return energy_surface_arap(C_hat,F,C);
	}
  else if (strcmp(Energy,"VolARAP")==0)
  {
    return energy_volumetric_arap(C_hat,F,C);
  }
	else if (strcmp(Energy,"Volume")==0)
	{
		return energy_volume(C,F);
	}
  else {
    cout << "ERROR: specify one of these energies:" << endl;
    cout << "DispStep, DispInitial, Volume, SurfARAP, VolARAP or None" << endl;
  }

	return 0.0;
}