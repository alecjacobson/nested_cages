#include "energy.h"
#include <stdio.h>

double energy_displacement(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat)
{
	return ((C-C_hat).rowwise().norm()).norm();
}

double energy(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const char* Energy)
{
	if (strcmp(Energy,"DispStep")==0){
		return energy_displacement(C,C_prev);
	}
	else if (strcmp(Energy,"DispStep")==0)
	{
		return energy_displacement(C,C_hat);
	}
}