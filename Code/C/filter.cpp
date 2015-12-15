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
#include <iomesh.h>
#include <meshrenderer.h>
#include <runstats.h>
#include <surftrack.h>
#include <trianglequality.h>

// talpa
#include <framestepper.h>
#include <meancurvature.h>
#include <meshdriver.h>
#include <scriptinit.h>
#include <simulation.h>

void filter(
  const Eigen::MatrixXd & Vf, 
  const Eigen::MatrixXi & T, 
  const Eigen::MatrixXd & Uf, 
  const Eigen::MatrixXd & C, 
  const Eigen::MatrixXi & F_hat, 
  Eigen::MatrixXd & Uc)
{
  using namespace Eigen;
  using namespace std;
  // Define concataned meshes (treat fine and coarse meshes as a big single mesh)
  MatrixXd V0(Vf.rows()+C.rows(),3);
  MatrixXd V1(Vf.rows()+C.rows(),3);
  V0.block(0,0,Vf.rows(),3) = Vf;
  V0.block(Vf.rows(),0,C.rows(),3) = C;
  V1.block(0,0,Vf.rows(),3) = Vf+Uf;
  V1.block(Vf.rows(),0,C.rows(),3) = C+Uc;
  MatrixXi F_all(T.rows()+F_hat.rows(),3);
  F_all.block(0,0,T.rows(),3) = T;

  F_all.block(T.rows(),0,F_hat.rows(),3) = F_hat;
  for (int k=T.rows();k<F_all.rows();k++)
  {
    F_all(k,0) = F_all(k,0)+Vf.rows();
    F_all(k,1) = F_all(k,1)+Vf.rows();
    F_all(k,2) = F_all(k,2)+Vf.rows();
  }

  // set vertex masses (=infty for fine mesh vertices 
  // and =1.0 for the last ones)
  double masses[V0.rows()];
  for (int i=0; i<Vf.rows(); i++){
      masses[i] = std::numeric_limits<double>::infinity();
  }
  for (int i=Vf.rows(); i<V0.rows(); i++){
      masses[i] = 1.0;
  }

  // Convert from Eigen matrix to array
  double V0a[3*V0.rows()];
  for (int k=0; k<V0.rows(); k++)
  {
    V0a[3*k] = V0(k,0);
    V0a[3*k+1] = V0(k,1);
    V0a[3*k+2] = V0(k,2);
  } 

  int F_alla[3*F_all.rows()];
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

  double V1a[3*V1.rows()];
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
  sim_general_options.m_proximity_epsilon = 1e-3; 


  // Set Simulation parameters
  ElTopoIntegrationOptions sim_integration_options;
  sim_integration_options.m_friction_coefficient = 0.0;
  sim_integration_options.m_dt = 1.0;


  double* V_final;
  double out_dt = 0.0;
  el_topo_integrate(&eltopo_time0, V1a, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
  if (out_dt<1.0){
  	cout << "Failed to reach final positions out_dt = " << out_dt << endl;
  	return;
  }
  // output corrected velocities
  for (int k=0; k<C.rows(); k++)
  {
    Uc(k,0) = V_final[3*Vf.rows()+3*k]-C(k,0);
    Uc(k,1) = V_final[3*Vf.rows()+3*k+1]-C(k,1);
    Uc(k,2) = V_final[3*Vf.rows()+3*k+2]-C(k,2); 
  }

  
}