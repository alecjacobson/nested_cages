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
  const igl::ARAPData & data,
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
    // cout << S.block(3950,0,100,3) << endl;
    // cout << S.rows() << " " << S.cols() << endl;

    MatrixXd R;
    fit_rotations(S,false,R);
    cout << R.block(0,0,9,3) << endl; // should be a stack of identities, but it isn't
    cout << R.rows() << " " << R.cols() << endl;

    // Dirichlket energy works fine
	grad = (data_L*U);

    // the following makes Eltopo stuck (probably bad direction)
    // grad = -(data.M*U + data_K*R);

    // grad = MatrixXd::Zero(V.rows(),3);

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