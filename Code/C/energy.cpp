#include "energy.h"
#include "gradient.h"
#include <stdio.h>
// Need to include some IGL header to have igl namespace
#include <igl/centroid.h>

double energy_displacement(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat)
{
	return ((C-C_hat).rowwise().norm()).norm();
}

// add description here
double energy_surface_arap(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & U)
{

	using namespace Eigen;
  	using namespace std;
  	using namespace igl;

    return 0.0;

}

// add description
double energy_volume(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F)
{

	using namespace Eigen;
  	using namespace igl;

  	VectorXd c(3,1);
  	double volume;

  	centroid(V,F,c,volume);

	return volume;

}

double energy(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev,
  const Eigen::MatrixXi &F, 
  const char* Energy)
{
	if (strcmp(Energy,"DispStep")==0){
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
	else if (strcmp(Energy,"Volume")==0)
	{
		return energy_volume(C,F);
	}

	return 0.0;
}