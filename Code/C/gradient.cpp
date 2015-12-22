#include "gradient.h"
#include <igl/cotmatrix.h>
#include <igl/covariance_scatter_matrix.h>
#include <igl/arap_rhs.h>
#include <igl/repmat.h>
#include <igl/fit_rotations.h>
#include <igl/doublearea.h>
#include <igl/per_face_normals.h>
#include <stdio.h>

void gradient_displacement(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat,
  Eigen::MatrixXd & grad)
{
	grad = C-C_hat;
	return;
}

// add description here
void gradient_surface_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U,
  Eigen::MatrixXd & grad)
{

	using namespace Eigen;
  	using namespace std;
  	using namespace igl;

    MatrixXd ref_V = V; 
    MatrixXi ref_F = F;

    SparseMatrix<double> data_L;
    cotmatrix(ref_V,ref_F,data_L);

    SparseMatrix<double> data_CSM;
    covariance_scatter_matrix(ref_V,ref_F,ARAP_ENERGY_TYPE_SPOKES_AND_RIMS,data_CSM);

	SparseMatrix<double> data_K;
    arap_rhs(ref_V,ref_F,3,ARAP_ENERGY_TYPE_SPOKES_AND_RIMS,data_K);

    MatrixXd U_rep;
    repmat(U,3,1,U_rep);
    MatrixXd S = data_CSM*U_rep;

    MatrixXd R;
    fit_rotations(S,false,R);

    VectorXd Rcol;
	columnize(R,R.cols()/3,2,Rcol);

	MatrixXd dV;
	dV = data_K*Rcol;
	MatrixXd dV3(U.rows(),U.cols());
	dV3.col(0) = dV.block(0,0,U.rows(),1);
	dV3.col(1) = dV.block(U.rows(),0,U.rows(),1);
	dV3.col(2) = dV.block(2*U.rows(),0,U.rows(),1);

	grad = -(data_L*U + dV3);

}

// add description here
void gradient_volumetric_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U,
  Eigen::MatrixXd & grad)
{

	int n_boundary = grad.rows();

	using namespace Eigen;
  	using namespace std;
  	using namespace igl;

    MatrixXd ref_V = V; 
    MatrixXi ref_F = F;

    SparseMatrix<double> data_L;
    cotmatrix(ref_V,ref_F,data_L);

    SparseMatrix<double> data_CSM;
    covariance_scatter_matrix(ref_V,ref_F,ARAP_ENERGY_TYPE_ELEMENTS,data_CSM);

	SparseMatrix<double> data_K;
    arap_rhs(ref_V,ref_F,3,ARAP_ENERGY_TYPE_ELEMENTS,data_K);

    MatrixXd U_rep;
    repmat(U,3,1,U_rep);
    MatrixXd S = data_CSM*U_rep;

    MatrixXd R;
    fit_rotations(S,false,R);

    VectorXd Rcol;
	columnize(R,R.cols()/3,2,Rcol);

	MatrixXd dV;
	dV = data_K*Rcol;
	MatrixXd dV3(U.rows(),U.cols());
	dV3.col(0) = dV.block(0,0,U.rows(),1);
	dV3.col(1) = dV.block(U.rows(),0,U.rows(),1);
	dV3.col(2) = dV.block(2*U.rows(),0,U.rows(),1);

	grad = -(data_L*U + dV3);

    grad = grad.block(0,0,n_boundary,3);

}

void gradient_volume(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & grad)
{

	using namespace Eigen;
  	using namespace std;
  	using namespace igl;

	MatrixXd N;
	per_face_normals(V,F,N); // unit normals
	VectorXd area;
	doublearea(V,F,area);
	// weight normals by triangle areas
	for (int k=0; k<F.rows(); k++)
	{
		N(k,0) = (area(k)/2.0)*N(k,0);
		N(k,1) = (area(k)/2.0)*N(k,1);
		N(k,2) = (area(k)/2.0)*N(k,2);
	}

	grad.resize(V.rows(),3);
	for (int k=0; k<F.rows(); k++)
	{
		grad(F(k,0),0) = N(k,0);
		grad(F(k,0),1) = N(k,1);
		grad(F(k,0),2) = N(k,2);

		grad(F(k,1),0) = N(k,0);
		grad(F(k,1),1) = N(k,1);
		grad(F(k,1),2) = N(k,2);

		grad(F(k,2),0) = N(k,0);
		grad(F(k,2),1) = N(k,1);
		grad(F(k,2),2) = N(k,2);
	}

	return;
}

bool gradient(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const Eigen::MatrixXi & F,
  const char* Energy,
  Eigen::MatrixXd & grad)
{
	using namespace std;
	if (strcmp(Energy,"DispStep")==0)
	{
		gradient_displacement(C,C_prev,grad);
	}
	else if (strcmp(Energy,"DispInitial")==0)
	{
		gradient_displacement(C,C_hat,grad);
	}
	else if (strcmp(Energy,"SurfARAP")==0)
	{
		gradient_surface_arap(C_hat,F,C,grad);
	}
	else if (strcmp(Energy,"VolARAP")==0)
	{
		grad.resize(C_prev.rows(),C_prev.cols()); // C_prev is a surface mesh 
		gradient_volumetric_arap(C_hat,F,C,grad);
	}
	else if (strcmp(Energy,"Volume")==0)
	{
		gradient_volume(C,F,grad);
	}
	else {
    cout << "ERROR: specify one of these energies:" << endl;
    cout << "DispStep, DispInitial, Volume, SurfARAP, VolARAP or None" << endl;
    return false;
  }

	return true;
}