#include "reinflate.h"
#include "filter.h"
#include "gradient.h"
#include "energy.h"
#include "flow.h"
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

  // data needed for Volumetric ARAP
  MatrixXd TV;
  MatrixXd TV0;
  MatrixXi TT;
  MatrixXi TF;
  ARAPData data;
  VectorXi b(C.rows());
  for (int k=0;k<C.rows();k++) b(k) = k;

  // reinflation
  while (!H.empty()){

    // Intialize energy as infinity
    double current_energy = std::numeric_limits<double>::infinity();
    // Initial step size
    double beta = 1e-2;

    // Fine mesh velocity
    MatrixXd Uf = H.top()-F;
    H.pop();

    // calculate diameter to scale spearation between meshes
    double diam = diameter(F,C);
    
    // If no energy presecribed, find only a feasible state (by filtering Uc=0)
    if (strcmp(EnergyInflation,"None")==0)
    {
    	MatrixXd Uc = MatrixXd::Zero(C_hat.rows(), 3);
      // this eps=5e-3 is chosen so that there is space for the final optimization
      filter(F,T,Uf,C,F_hat,diam*5e-3,Uc);
      // Update fine mesh
      F = F+Uf;
      // Update coarse mesh with filtered velocities
      C = C+Uc;
      #ifdef VERBOSE_DEBUG
        cout << "Reinflation step " << step << "/" << total_steps << ": Energy = None (feasible state only)" << endl;
      #endif
    }
    // If energy for re-inflation is prescribed, minimize it
    else if (strcmp(EnergyInflation,"None")!=0)
    {
      // If Volumetric ARAP, tetrahedralize and precompute ARAP Data
      if (strcmp(EnergyInflation,"VolARAP")==0)
      {
        string tetgen_flags ("q2Y");
        tetrahedralize(C_hat,F_hat,tetgen_flags,TV0,TT,TF); 
        if (!arap_precomputation(TV0,TT,3,b,data))
        {
          cout << "ARAP Precomputation failed" << endl;
        }
        // Initialize deformed tet mesh as initial one
        TV = TV0;
      }
      // Stepping and projecting
      while (true)
      {
        #ifdef VERBOSE_DEBUG
          cout << "Reinflation step " << step << "/" << total_steps << ": Energy = " << EnergyInflation  << endl;
        #endif
        MatrixXd grad;
        // for Volumetric ARAP, calculate gradient after finding optimal 
        // deformed tet mesh (with local-globe solver)
        if (strcmp(EnergyInflation,"VolARAP")==0)
        {
          if (!arap_solve(C,data,TV))
          {
            cout << "ARAP Solve failed" << endl;
            return;
          }
          if (!gradient(TV,TV0,C_prev,TT,EnergyInflation,grad))
          {
            cout << "ERROR calculating gradient " << endl;
            return;
          }

        }
        // for other energies, just calculate the gradient
        else
        {
          if (!gradient(C,C_hat,C_prev,F_hat,EnergyInflation,grad))
          {
            cout << "ERROR calculating gradient " << endl;
            return;
          }
        }
        // multiply gradient by current step size
        MatrixXd Uc = -beta*grad;
        // filter coarse mesh velocities
        filter(F,T,Uf,C,F_hat,diam*1e-3,Uc);

        // update fine mesh (it has to be bone here, because
        // otherwise the collision solver will check intersections 
        // between new coarse mesh and old fine mesh - they rarely intersect,
        // but can be a problem)
        F = F+Uf;
        Uf = MatrixXd::Zero(F.rows(), 3);

        // update energy value
        double new_energy;
        if (strcmp(EnergyInflation,"VolARAP")==0)
        {
          new_energy = energy(TV,TV0,C_prev,TT,EnergyInflation);
        }
        else
        {
          new_energy = energy(C+Uc,C_hat,C_prev,F_hat,EnergyInflation);
        }

        // if energy increased, cut step size by half
        if (new_energy>current_energy)
        {
          beta = 0.5*beta;
          #ifdef VERBOSE_DEBUG
            cout << "energy increased, descreasing beta to " << beta << endl;
          #endif
          // If step size is too small, then it has converged. Break
          if (beta<1e-3)
          {
            #ifdef VERBOSE_DEBUG
              cout << "beta too small. Quitting line search loop " << endl;
            #endif
            break;
          }
        }
        // if energy decreased, update coarse mesh, update energy value and 
        // and increase beta 10%
        else{
          C = C+Uc;
          current_energy = new_energy;
          beta = 1.1*beta;
          #ifdef VERBOSE_DEBUG
            cout << "energy decreased to " << current_energy << ", increasing beta to " << beta << endl;
          #endif
        }
        // If tiny step, then it has converged. Break
        if (((Uc).rowwise().norm()).maxCoeff()<diam*1e-5) 
        {
          #ifdef VERBOSE_DEBUG
            cout << "Max change in postion = " << ((Uc).rowwise().norm()).maxCoeff() << " too small. Quitting line search loop " << endl;
          #endif
          break;
        }
      }
    }

    // Update step (only for printing purposes)
    step = step+1;
    // Update previous coarse mesh (currently used only for DispStep energy)
    C_prev = C;

  }

  // final (optional) minimization - Useful especially for volume minimization,
  // when it makes little sense to minimize volume when the fine mesh is not 
  // completely re-inflated
  if (strcmp(EnergyFinal,"None")!=0)
  {

    // Intialize energy as infinity
    double current_energy = std::numeric_limits<double>::infinity();
    // Intial step size
    double beta = 1e-2;

    // calculate diameter to scale spearation between meshes
    double diam = diameter(F,C);

    // At this moment the fine mesh is back to its original
    // (input) positions, so its velocity should be zero
    MatrixXd Uf = MatrixXd::Zero(F.rows(), 3);

    // If Volumetric ARAP, tetrahedralize and precompute ARAP Data
    if (strcmp(EnergyFinal,"VolARAP")==0){
      string tetgen_flags ("q2Y");
      tetrahedralize(C_hat,F_hat,tetgen_flags,TV0,TT,TF); 
      arap_precomputation(TV0,TT,3,b,data);
      TV = TV0;
    }
    // Stepping and projecting
    while (true)
      {
        #ifdef VERBOSE_DEBUG
          cout << "Final optimization. Energy = " << EnergyFinal  << endl;
        #endif
        MatrixXd grad;
        // for Volumetric ARAP, calculate gradient after finding optimal 
        // deformed tet mesh (with local-globe solver)
        if (strcmp(EnergyFinal,"VolARAP")==0)
        {
          if (!arap_solve(C,data,TV))
          {
            cout << "ARAP Solve failed " << endl;
            return;
          }
          if (!gradient(TV,TV0,C_prev,TT,EnergyFinal,grad))
          {
            cout << "ERROR calculating gradient " << endl;
            return;
          }
        }
        // for other energies, just calculate the gradient
        else
        {
          if (!gradient(C,C_hat,C_prev,F_hat,EnergyFinal,grad))
          {
            cout << "ERROR calculating gradient " << endl;
            return;
          }
        }
        // multiply gradient by current step size
        MatrixXd Uc = -beta*grad;
        // filter coarse mesh velocities
        // this eps=5e-4 is chosen so that tere's is space for final optimizarion
        filter(F,T,Uf,C,F_hat,diam*5e-4,Uc); 

        // update energy value
        double new_energy;
        if (strcmp(EnergyFinal,"VolARAP")==0)
        {
          new_energy = energy(TV,TV0,C_prev,TT,EnergyFinal);
        }
        else
        {
          new_energy = energy(C+Uc,C_hat,C_prev,F_hat,EnergyFinal);
        }
        // if energy increased, cut step size by half
        if (new_energy>current_energy)
        {
          beta = 0.5*beta;
          #ifdef VERBOSE_DEBUG
            cout << "energy increased, descreasing beta to " << beta << endl;
          #endif
          // If step size is too small, then it has converged. Break
          if (beta<1e-3)
          {
            #ifdef VERBOSE_DEBUG
              cout << "beta too small. Quitting line search loop " << endl;
            #endif
            break;
          }
        }
        // if energy decreased, update coarse mesh, update energy value and 
        // and increase beta 10%
        else{
          C = C+Uc;
          current_energy = new_energy;
          beta = 1.1*beta;
          #ifdef VERBOSE_DEBUG
            cout << "energy decreased to " << current_energy << ", increasing beta to " << beta << endl;
          #endif
        }
        // If tiny step, then it has converged. Break
        if (((Uc).rowwise().norm()).norm()<diam*1e-5) break;
      }

  }

  return;

}