#include "gradient.h"
#include <igl/cotmatrix.h>
#include <igl/covariance_scatter_matrix.h>
#include <igl/arap_rhs.h>
#include <igl/repmat.h>
#include <igl/fit_rotations.h>
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
  const igl::ARAPData & data,
  Eigen::MatrixXd & grad)
{

	using namespace Eigen;
  	using namespace std;
  	using namespace igl;

    // MatrixXd ref_V = V;
    // MatrixXi ref_F = F;

    // ARAPData data_;
    // VectorXi b;
    // b.resize(0);
    // if (arap_precomputation(ref_V,ref_F,3,b,data_)){
    // 	cout << "ARAP precomputation successfully done " << endl;
    // }
    // else{
    // 	cout << "ARAP precomputation failed " << endl;
    // }

    SparseMatrix<double> data_L;
    cotmatrix(V,F,data_L);

 //    SparseMatrix<double> data_CSM;
 //    covariance_scatter_matrix(ref_V,ref_F,ARAP_ENERGY_TYPE_SPOKES_AND_RIMS,data_CSM);

	// SparseMatrix<double> data_K;
 //    arap_rhs(ref_V,ref_F,3,ARAP_ENERGY_TYPE_SPOKES_AND_RIMS,data_K);

    // MatrixXd S = MatrixXd::Zero(data_CSM.rows(), 3);
    // MatrixXd U_rep;
    // repmat(U,3,1,U_rep);
    // S = data_CSM*U_rep;

    MatrixXd S = data.CSM*U;

    MatrixXd R;
    fit_rotations(S,true,R);

    // Dirichlket energy with data_L (computed inside this function)
    // works fine
	grad = -(data_L*U);

	// the following makes Eltopo stuck (probably baqd direction)
    // grad = -(data.M*U);
    // the following makes Eltopo stuck (probably baqd direction)
    // grad = -(data.M*U + data.K*R);

    // grad = MatrixXd::Zero(V.rows(),3);

}

void gradient_volume(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & grad)
{

	// Is there an analogue to 'normals' in Matlab 

	//N = normals(CV,CF)/2;
    //grad_vol = full(sparse( ...
    //  repmat(CF(:),1,3),repmat(1:3,numel(CF),1),repmat(N,3,1),size(CV,1),3));

	return;
}

void gradient(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const Eigen::MatrixXi & F,
  const igl::ARAPData & data,
  const char* Energy,
  Eigen::MatrixXd & grad)
{
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
		gradient_surface_arap(C_hat,F,C,data,grad);
	}
	else if (strcmp(Energy,"Volume")==0)
	{
		gradient_volume(C,F,grad);
	}
	return;
}