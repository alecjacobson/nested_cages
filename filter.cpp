#include "filter.h"

// common
#include <array2.h>
#include <ccd_wrapper.h>
#include <collisionqueries.h>
#include <marching_tiles_hires.h>
#include <util.h>
#include <vec.h>
#include <wallclocktime.h>

// el topo
#include <collisionpipeline.h>
#include <eltopo.h>
#include <meshrenderer.h>
#include <runstats.h>
#include <surftrack.h>
#include <trianglequality.h>

// ACM
#include <VelocityFilter.h>
#include <Distance.h>

double inflate_ACM(
   Eigen::MatrixXd & V0, 
  const Eigen::MatrixXi & F0, 
  double eps_distance,
  int numinfinite)
{
  using namespace Eigen;
  using namespace std;

  Matrix3Xi F0_t(3,F0.rows());
  for (int k=0; k<F0.rows();k++){
      F0_t(0,k) = F0(k,0);
      F0_t(1,k) = F0(k,1);
      F0_t(2,k) = F0(k,2);
  }

  VectorXd q(3*V0.rows());
  for (int k=0; k<V0.rows();k++){
      q(3*k) = V0(k,0);
      q(3*k+1) = V0(k,1);
      q(3*k+2) = V0(k,2);
  }
  VectorXd invmasses(3*V0.rows());
  for (int k=0; k<3*numinfinite; k++){
      invmasses(k) = 0.0;
  }
  for (int k=3*numinfinite; k<3*V0.rows(); k++){
      invmasses(k) = 1.0;
  }

  set<int> fixedVerts;
  for(int i=0; i<numinfinite; i++)
      fixedVerts.insert(i);

  double distance = Distance::meshSelfDistance(q, F0_t, fixedVerts);
  double prev_distance = distance;
  while(distance < eps_distance)
  {
      cout << "Distance is " << distance << ". Inflating... " << endl;
      VectorXd qnew = q;
      VelocityFilter::velocityFilter(q, qnew, F0_t, invmasses, 2.0*distance, 0.5*distance);
      q = qnew;
      distance = Distance::meshSelfDistance(q, F0_t, fixedVerts);
      if (distance<prev_distance){
          cout << "Couldn't recah prescribed distance " << eps_distance << ". Returned " << prev_distance;
          break;
      }
      prev_distance = distance;
  }

  // overwrite V0
  for (int k=0; k<V0.rows(); k++){
      V0(k,0) = q(3*k);
      V0(k,1) = q(3*k+1);
      V0(k,2) = q(3*k+2);
  }

  if (distance>eps_distance) return eps_distance;
  else return distance;

}

void velocity_filter_ACM(
  const Eigen::MatrixXd & V0, 
  Eigen::MatrixXd & V1, 
  const Eigen::MatrixXi & F0, 
  double outer,
  double inner,
  int numinfinite)
{
  using namespace Eigen;
  using namespace std;

  Matrix3Xi F0_t(3,F0.rows());
  for (int k=0; k<F0.rows();k++){
      F0_t(0,k) = F0(k,0);
      F0_t(1,k) = F0(k,1);
      F0_t(2,k) = F0(k,2);
  }

  VectorXd qstart(3*V0.rows());
  VectorXd qend(3*V0.rows());
  for (int k=0; k<V0.rows();k++){
      qstart(3*k) = V0(k,0);
      qstart(3*k+1) = V0(k,1);
      qstart(3*k+2) = V0(k,2);
  }
  for (int k=0; k<V0.rows();k++){
      qend(3*k) = V1(k,0);
      qend(3*k+1) = V1(k,1);
      qend(3*k+2) = V1(k,2);
  }
  VectorXd invmasses(3*V0.rows());
  for (int k=0; k<3*numinfinite; k++){
      invmasses(k) = 0.0;
  }
  for (int k=3*numinfinite; k<3*V0.rows(); k++){
      invmasses(k) = 1.0;
  }

  int sim_status = VelocityFilter::velocityFilter(qstart, qend, F0_t, invmasses, outer, inner);
  cout << "simultaion_status = " << sim_status << endl;

  for (int k=0; k<V0.rows(); k++){
      V1(k,0) = qend(3*k);
      V1(k,1) = qend(3*k+1);
      V1(k,2) = qend(3*k+2);
  }

  return;
}

void filter(
  const Eigen::MatrixXd & Vf, 
  const Eigen::MatrixXi & T, 
  const Eigen::MatrixXd & Uf, 
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXi & F_hat, 
  double eps_prox,
  Eigen::MatrixXd & Uc)
{
  using namespace Eigen;
  using namespace std;
  // Define concataned meshes (treat fine and coarse meshes as a big single mesh)
  MatrixXd V0(Vf.rows()+C.rows(),3);
  MatrixXd V1(Vf.rows()+C.rows(),3);
  // Initial positions
  V0.block(0,0,Vf.rows(),3) = Vf;
  V0.block(Vf.rows(),0,C.rows(),3) = C;
  // Final positions
  V1.block(0,0,Vf.rows(),3) = Vf+Uf;
  V1.block(Vf.rows(),0,C.rows(),3) = C+Uc;
  
  // First part of F
  MatrixXi F_all(T.rows()+F_hat.rows(),3);
  // Second part of F
  F_all.block(0,0,T.rows(),3) = T;
  F_all.block(T.rows(),0,F_hat.rows(),3) = F_hat;
  for (int k=T.rows();k<F_all.rows();k++)
  {
    F_all(k,0) = F_all(k,0)+Vf.rows();
    F_all(k,1) = F_all(k,1)+Vf.rows();
    F_all(k,2) = F_all(k,2)+Vf.rows();
  }

  // set vertex masses (= infty for fine mesh vertices 
  // and =1.0 for the last ones)
  double *masses = new double[V0.rows()];
  for (int i=0; i<Vf.rows(); i++){
      masses[i] = std::numeric_limits<double>::infinity();
  }
  for (int i=Vf.rows(); i<V0.rows(); i++){
      masses[i] = 1.0;
  }

  // Convert V0 from Eigen matrix to array
  double *V0a = new double[3*V0.rows()];
  for (int k=0; k<V0.rows(); k++)
  {
    V0a[3*k] = V0(k,0);
    V0a[3*k+1] = V0(k,1);
    V0a[3*k+2] = V0(k,2);
  } 

  // Convert F_all from Eigen matrix to array
  int *F_alla = new int[3*F_all.rows()];
  for (int k=0; k<F_all.rows(); k++)
  {
    F_alla[3*k] = F_all(k,0);
    F_alla[3*k+1] = F_all(k,1);
    F_alla[3*k+2] = F_all(k,2);
  } 

  // encapsulate all data into an ElTopoMesh
  ElTopoMesh eltopo_time0;
  eltopo_time0.num_vertices = V0.rows();
  eltopo_time0.vertex_locations = V0a;
  eltopo_time0.num_triangles = F_all.rows();
  eltopo_time0.triangles = F_alla;
  eltopo_time0.vertex_masses = masses;

  // Convert V1 from Eigen matrix to array
  double *V1a = new double[3*V1.rows()];
  for (int k=0; k<V1.rows(); k++)
  {
    V1a[3*k] = V1(k,0);
    V1a[3*k+1] = V1(k,1);
    V1a[3*k+2] = V1(k,2);
  } 

  // Set general parameters
  ElTopoGeneralOptions sim_general_options;
  // do not print stuff to the console
  sim_general_options.m_verbose = 0;
  // do avoid self-intersections
  sim_general_options.m_collision_safety = 1;
  // separation between colliding meshes 
  sim_general_options.m_proximity_epsilon = eps_prox; 


  // Set Simulation parameters
  ElTopoIntegrationOptions sim_integration_options;
  sim_integration_options.m_friction_coefficient = 0.0;
  sim_integration_options.m_dt = 1.0;


  double* V_final;
  double out_dt = 0.0;
  // We start with 1.0 to step
  double rest_dt = 1.0;
  // While we haven't reached final positions
  while (rest_dt>1e-6)
  {
    // call Eltopo main function
    el_topo_integrate(&eltopo_time0, V1a, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
    // update the rest to go
    rest_dt = (1-out_dt)*rest_dt;
    // if we haven't reached final positions, print how much we have stepped
    // and update vertex positions 
    if (out_dt < 1.0)
    {
      cout << "out_dt = " << out_dt << endl;
      eltopo_time0.vertex_locations = V_final;
    }
    // if stepped too little, give up. To-do: call Asynchronous Contact Mechanincs
    if (out_dt<1e-3){
    	cout << "Eltopo couldn't reach final positions." << endl;
      cout << "It steped " << (1-rest_dt) << "< 1.0" << endl;
      cout << "Calling Asyncronous Contact Mechanics (slower)" << endl;
      double out_proximity = inflate_ACM(V0,F_all,eps_prox,Vf.rows());
      velocity_filter_ACM(V0,V1,F_all,out_proximity,0.01*out_proximity,Vf.rows());
      for (int k=0; k<V1.rows(); k++)
      {
        V_final[3*k] = V1(k,0);
        V_final[3*k+1] = V1(k,1);
        V_final[3*k+2] = V1(k,2);
      } 
    }
  }
  // output corrected velocities
  for (int k=0; k<C.rows(); k++)
  {
    Uc(k,0) = V_final[3*Vf.rows()+3*k]-C(k,0);
    Uc(k,1) = V_final[3*Vf.rows()+3*k+1]-C(k,1);
    Uc(k,2) = V_final[3*Vf.rows()+3*k+2]-C(k,2); 
  }

  delete[] masses;
  delete[] V0a;
  delete[] V1a;
  delete[] F_alla;
  
}
