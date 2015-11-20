#include <igl/doublearea.h>
#include <igl/cgal/signed_distance.h>
#include <igl/cgal/point_mesh_squared_distance.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace igl;
using namespace Eigen;
using namespace std;

SparseMatrix<double> gradQ_to_gradV(MatrixXd , MatrixXi , MatrixXd, int );
MatrixXd signed_distance_direction(MatrixXd , MatrixXd , MatrixXi );
void signed_distance(MatrixXd , MatrixXd , MatrixXi , MatrixXd*, MatrixXd* , VectorXi* , VectorXd* );
VectorXd grad_energy(MatrixXd , MatrixXi , MatrixXd , MatrixXi , int , SparseMatrix<double>);
MatrixXd flow_one_step(MatrixXd , MatrixXi , MatrixXd , MatrixXi , int , SparseMatrix<double>);
void flow_fine_inside_coarse(MatrixXd , MatrixXi , MatrixXd , MatrixXi , int , SparseMatrix<double>);