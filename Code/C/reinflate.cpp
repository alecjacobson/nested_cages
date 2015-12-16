#include "reinflate.h"
#include "filter.h"
#include "gradient.h"
#include "energy.h"
#include <igl/arap.h>
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
  MatrixXd C_prev = C; // needed for energies that depend on the previous step (e.g. DispStep)
  int step = 1;
  while (!H.empty()){

    // Intialize energy as infinity
    double current_energy = std::numeric_limits<double>::infinity();
    // step size search
    double beta = 1e-2;

    // Fine mesh velocity
    MatrixXd Uf = H.top()-F;
    H.pop();
    
    if (strcmp(EnergyInflation,"None")==0)
    {
    	// First, find a feasible state (before stepping and projecting)
    	MatrixXd Uc = MatrixXd::Zero(C_hat.rows(), 3);
      filter(F,T,Uf,C,F_hat,Uc);
      C = C+Uc;
      cout << "Reinflation step " << step << ": Energy = None (fesible state only)" << endl;
    }
    else if (strcmp(EnergyInflation,"None")!=0)
    {
      ARAPData data; // precompute ARAP data for ARAP energies
      if (strcmp(EnergyInflation,"SurfARAP")==0)
      {
        VectorXi b;
        b.resize(0);
        if (arap_precomputation(C_hat,F_hat,3,b,data))
        {
          cout << "ARAP precomputation successfully done " << endl;
        }
        else
        {
          cout << "ARAP precomputation failed, returning... " << endl;
          return;
        }
      }
      while (true)
      {
        cout << "Reinflation step " << step << ": Energy = " << EnergyInflation  << endl;
        MatrixXd grad;
        gradient(C,C_hat,C_prev,F_hat,data,EnergyInflation,grad);
        MatrixXd Uc = -beta*grad;
        filter(F,T,Uf,C,F_hat,Uc);
        double new_energy = energy(C+Uc,C_hat,C_prev,F_hat,data,EnergyInflation);
        if (new_energy>current_energy)
        {
          beta = 0.5*beta;
          cout << "energy increased, descreasing beta to " << beta << endl;
          if (beta<1e-3)
          {
            cout << "beta too small. Quitting line search loop " << endl;
            break;
          }
        }
        else{
          C = C+Uc;
          current_energy = new_energy;
          beta = 1.1*beta;
          cout << "energy decreased to " << current_energy << ", increasing beta to " << beta << endl;
        }
        // If tiny step, then brak
        if (((Uc).rowwise().norm()).norm()<1e-5) break;
      }
    }


    // update
    F = F+Uf;
    step = step+1;
    C_prev = C;

  }

  return;

}