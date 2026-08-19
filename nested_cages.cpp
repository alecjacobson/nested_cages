// Command-line front-end for the "Nested Cages" [Sacht et al. 2015] pipeline.
// The actual algorithm lives in nested_cages_lib.{h,cpp}; this file only parses
// arguments and reads/writes meshes so that the exact same code path can be
// reused from the Python bindings.

#include "nested_cages_lib.h"
#include "io.h"

// libigl includes
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char * argv[])
{
  using namespace std;

  if (argc==1)
  {
    cout << R"(Usage:

    ./nested_cages input.off q L(1) L(2) ... L(k) EnergyExpansion EnergyFinal output

input: the program accepts files in the following formats: .off, .obj, .ply, .stl, .wrl .mesh

output: cages will be saved as output_1.obj, output_2.obj, ..., output_k.obj

q is the quadrature order for the shrinking flow

L(1) > L(2) > ... > L(k) is the number of faces for each cage.
If L(k) is followed by 'r' the initial decimation for this cage will be regular
(adaptive if no 'r').
Each L(k) can be replace by a file with an input decimation.

EnergyExpansion is the energy to be minimized for the re-inflation
Energies implemented: None, DispStep, DispInitial, Volume, SurfARAP, VolARAP

EnergyFinal is the energy to be minimized after the re-inflation (additional processing)
Energies implemented: None, DispStep, DispInitial, Volume, SurfARAP, VolARAP
)";
    return EXIT_FAILURE;
  }

  // number of layers
  const int k = argc-6;
  #ifdef VERBOSE_DEBUG
    cout << "number of layers = " << k << endl;
  #endif

  // quadrature order
  const int quad_order = atoi(argv[2]);

  // read input mesh
  Eigen::MatrixXd V0;
  Eigen::MatrixXi F0;
  if (!igl::read_triangle_mesh(argv[1],V0,F0)){
    cout << "unable to read input file"  << endl;
    return 0;
  }

  const string prefix = argv[argc-1];
  const string energy_expansion = argv[argc-3];
  const string energy_final     = argv[argc-2];

  // output input mesh as level output_0.obj
  igl::writeOBJ(prefix + "_0.obj", V0, F0);

  // Per-level inputs to nested_cages(): either an explicit input decimation
  // (read from a file) or a target face count + regular/adaptive flag.
  vector<int> num_faces(k,0);
  vector<bool> regular(k,false);
  vector<Eigen::MatrixXd> V_decim(k);
  vector<Eigen::MatrixXi> F_decim(k);

  for(int i = 0;i<k;i++){
    // check if argv[i+3] is a valid file. If so, it will be used as input decimation
    std::ifstream is_file(argv[i+3]);
    if (is_file)
    {
      if (!igl::read_triangle_mesh(argv[i+3],V_decim[i],F_decim[i])){
        cout << "error: input decimation is an existing file, but couldn't be read"  << endl;
        return 0;
      }
    }
    else
    {
      // first check if last character of argv[i+3] is r. If it is, drop the 'r'
      // and use regular (rather than adaptive) decimation.
      regular[i] = remove_all_chars_and_count(argv[i+3], 'r')!=0;
      // check if argv[i+3] is a valid integer (throw an error if it is not)
      if (!legal_int(argv[i+3])){
        cout << "you have to pass integer values or valid input deimatations"  << endl;
        cout << "the invalid argument you have passed is " << argv[i+3] << endl;
        return 0;
      }
      // specified number of faces for this cage
      num_faces[i] = atoi(argv[i+3]);
    }
  }

  // run the pipeline, writing each cage to output_{i+1}.obj as soon as it is
  // ready (matching the original program's behaviour on partial success).
  vector<Eigen::MatrixXd> Vout;
  vector<Eigen::MatrixXi> Fout;
  auto write_cage = [&](int i, const Eigen::MatrixXd & C, const Eigen::MatrixXi & Fc){
    igl::writeOBJ(prefix + "_" + std::to_string(i+1) + ".obj", C, Fc);
  };
  if(!nested_cages(
       V0,F0,quad_order,num_faces,regular,V_decim,F_decim,
       energy_expansion,energy_final,Vout,Fout,write_cage))
  {
    // nested_cages() already printed the reason it quit.
    return 0;
  }

  return 1;
}
