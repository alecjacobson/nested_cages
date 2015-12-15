#include "gradient.h"
#include <stdio.h>

void gradient_displacement(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat,
  Eigen::MatrixXd & grad)
{
	grad = C-C_hat;
	return;
}

void gradient(
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXd & C_hat, 
  const Eigen::MatrixXd & C_prev, 
  const char* Energy,
  Eigen::MatrixXd & grad)
{
	if (strcmp(Energy,"DispStep")==0)
	{
		gradient_displacement(C,C_prev,grad);
	}
	else if (strcmp(Energy,"DispStep")==0)
	{
		gradient_displacement(C,C_hat,grad);
	}
	return;
}