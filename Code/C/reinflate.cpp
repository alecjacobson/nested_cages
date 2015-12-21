#include "reinflate.h"
#include "filter.h"
#include "gradient.h"
#include "energy.h"
#include <igl/arap.h>
#include <stdio.h>

// Need to include some IGL header to have igl namespace
#include <igl/signed_distance.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

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
  using namespace igl::copyleft::tetgen;

  MatrixXd F = H.top();
  H.pop();
  C = C_hat;
  MatrixXd C_prev = C; // needed for energies that depend on the previous step (e.g. DispStep)
  int step = 1;
  int total_steps = H.size();

  // reinflation
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
    	MatrixXd Uc = MatrixXd::Zero(C_hat.rows(), 3);
      // this eps=5e-3 is chosen so that there is space for the final optimization
      filter(F,T,Uf,C,F_hat,5e-3,Uc);
      F = F+Uf;
      C = C+Uc;
      cout << "Reinflation step " << step << "/" << total_steps << ": Energy = None (feasible state only)" << endl;
    }
    else if (strcmp(EnergyInflation,"None")!=0)
    {
      // // First, find a feasible state (before stepping and projecting)
      // MatrixXd Uc = MatrixXd::Zero(C_hat.rows(), 3);
      // filter(F,T,Uf,C,F_hat,1e-3,Uc);
      // F = F+Uf;
      // Uf = MatrixXd::Zero(F.rows(), 3);

      // double current_energy = energy(C+Uc,C_hat,C_prev,F_hat,EnergyInflation);
      // C = C+Uc;
      // cout << "Reinflation step " << step << "/" << total_steps << ". Feasible state only. With energy = " << current_energy << endl;
      
      MatrixXd TV;
      MatrixXi TT;
      MatrixXi TF;
      if (strcmp(EnergyInflation,"VolARAP")==0){
        string tetgen_flags ("-Cpg -q100");
        // tetrahedralize(C_hat,F_hat,tetgen_flags,TV,TT,TF);
      }

      while (true)
      {
        cout << "Reinflation step " << step << "/" << total_steps << ": Energy = " << EnergyInflation  << endl;
        MatrixXd grad;
        gradient(C,C_hat,C_prev,F_hat,EnergyInflation,grad);
        MatrixXd Uc = -beta*grad;
        filter(F,T,Uf,C,F_hat,1e-3,Uc);

        // update fine mesh (it has to be bone here, beacuse
        // otherwise the collision solver wilkl check intersections 
        // between new coarse mesh and old fine mesh - they rarely intersect,
        // but can be a problem)
        F = F+Uf;
        Uf = MatrixXd::Zero(F.rows(), 3);

        double new_energy = energy(C+Uc,C_hat,C_prev,F_hat,EnergyInflation);
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
        if (((Uc).rowwise().norm()).maxCoeff()<1e-5) 
        {
          cout << "Max change in postion = " << ((Uc).rowwise().norm()).maxCoeff() << " too small. Quitting line search loop " << endl;
          break;
        }
      }
    }

    step = step+1;
    C_prev = C;

  }

  // final (optional) minimization
  if (strcmp(EnergyFinal,"None")!=0)
  {

    // Intialize energy as infinity
    double current_energy = std::numeric_limits<double>::infinity();
    double beta = 1e-2;

    MatrixXd Uf = MatrixXd::Zero(F.rows(), 3);

    // to avoid tons of collisions
    bool final_minimization = true;
    // cout << "current volume = " << energy(C,C_hat,C_prev,F_hat,"Volume") << ", original volume = " <<  energy(F,C_hat,C_prev,T,"Volume") << endl;
    // if (strcmp(EnergyFinal,"Volume")==0)
    // {
    //   if (energy(C,C_hat,C_prev,F_hat,"Volume") < 1.05*energy(F,C_hat,C_prev,T,"Volume"))
    //   {
    //     cout << "current cage volume is too close to fine mesh volume. Skipping final minimization" << endl;
    //     final_minimization = false;
    //   }
    // }

    while (final_minimization)
      {
        cout << "Final optimization. Energy = " << EnergyFinal  << endl;
        MatrixXd grad;
        gradient(C,C_hat,C_prev,F_hat,EnergyFinal,grad);
        MatrixXd Uc = -beta*grad;
        // this eps=5e-4 is chosen so that tere's is space for final optimizarion
        filter(F,T,Uf,C,F_hat,5e-4,Uc); 
        double new_energy = energy(C+Uc,C_hat,C_prev,F_hat,EnergyFinal);
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

  return;

}