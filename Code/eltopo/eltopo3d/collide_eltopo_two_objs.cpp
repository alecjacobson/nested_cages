// ---------------------------------------------------------
//
//  collide_eltopo_two_objs.cpp
//  Leonardo Sacht, 2014.
//
//  Takes as input 2 OBJ files, and a number N that sets
//  the first N vertices to have infinite mass. Outputs
//  a new OBJ file.
//
//  To compile:
//     gcc collide_eltopo_two_objs.cpp robustleastsquares.cpp -I../common \
//       -I../eltopo3d -I../talpa -I../talpa/drivers -I../common/tunicate \
//       -llapack -lblas -lstdc++ \
//       /Users/Leo/PHD_Work/Cage_Generation_2013/code/starlab-mcfskel/core/external/cholmod-4.0.0/lib/osx64/libcholmod.a /Users/Leo/PHD_Work/Cage_Generation_2013/code/starlab-mcfskel/core/external/cholmod-4.0.0/lib/osx64/libamd.a /Users/Leo/PHD_Work/Cage_Generation_2013/code/starlab-mcfskel/core/external/cholmod-4.0.0/lib/osx64/libcolamd.a libeltopo_release.a ../talpa/obj/bfstream.o ../talpa/obj/iomesh.o -I/opt/local/include/eigen3 -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include/igl /opt/local/lib/gcc47/libgomp.a /Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/lib/libigl.a /opt/local/lib/libSuiteSparse.dylib -o collide_eltopo_two_objs
// 
//     g++ collide_eltopo_two_objs.cpp \
//       ../common/bfstream.cpp ../talpa/iomesh.cpp -I../common \
//       -I../eltopo3d -I../talpa -I../talpa/drivers -I../common/tunicate \
//       -llapack -lblas -lstdc++ -I/opt/local/include/eigen3 -lcholmod -L. \
//       -leltopo_release -L/opt/local/lib -lSuiteSparse -DNO_GUI -o \
//       collide_eltopo_two_objs
// 
//
// Example usage:
// ./collide_eltopo_two_objs V0.obj V1.obj 100 V_out.obj
//
//

// std
#include <cstdio>
#include <fenv.h>
#include <fstream>
#include <vector>
#include <queue>
#include <cfloat>

// common
#include <array2.h>
#include <ccd_wrapper.h>
#include <collisionqueries.h>
#include <expansion.h>
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

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include <chrono>
#include <thread>

// ---------------------------------------------------------
// Global interface declarations
// ---------------------------------------------------------

extern "C" {
    void exactinit();    // in predicates.c
}

//int main(int argc, char **argv);

int main(int argc, char **argv)
{
    if( argc != 5 )
    {
        std::cout << "Usage: <executable> <OBJ_file_1> <OBJ_file_2> <num_constrained_vertices> <OBJ_output_file>" << std::endl;
        std::cout << "e.g.: " << std::endl;
        std::cout << "$ ./collide_eltopo_two_objs meshes/V0_test1.obj meshes/V1_test1.obj 100 meshes/V_out_test1.obj" << std::endl;
        return 0;
    }

    // Mesh at time 0
    NonDestructiveTriMesh mesh_time0;
    std::vector<Vec3d> x0;
    // Read mesh_time0 from OBJ file
    read_objfile(mesh_time0, x0, argv[1]);
//    // Read mesh_time0 from BIN file
//    read_binary_pos_tris_file(mesh_time0, x0, argv[1]);
    // get number of vertices
    int num_vertices = 0;
    for (std::vector<Vec3d>::iterator it = x0.begin() ; it != x0.end(); ++it){
        num_vertices = num_vertices+1;
    }
    // get vertex positions
    double V0[3*num_vertices];
    int i = 0;
    for (std::vector<Vec3d>::iterator it = x0.begin() ; it != x0.end(); ++it){
        Vec3d vertex = *it;
        V0[3*i] = vertex[0];
        V0[3*i+1] = vertex[1];
        V0[3*i+2] = vertex[2];
        i = i + 1;
    }
    // get numver of triangles
    int num_triangles = mesh_time0.num_triangles();
    // get face indices
    int F0[3*num_triangles];
    std::vector<Vec3st> tri_pointer = mesh_time0.get_triangles();
    i = 0;
    for (std::vector<Vec3st>::iterator it = tri_pointer.begin(); it<tri_pointer.end(); ++it){
        Vec3st triangle = *it;
        F0[3*i] = triangle[0];
        F0[3*i+1] = triangle[1];
        F0[3*i+2] = triangle[2];
        i = i + 1;
    }
    // set vertex masses (=infty for the first vertices and =eps for the last ones)
    double masses[num_vertices];
    for (int i=0; i<atoi(argv[3]); i++){
        masses[i] = std::numeric_limits<double>::infinity();
    }
    for (int i=atoi(argv[3]); i<num_vertices; ++i){
        masses[i] = 1.0;
    }

    // encapsulate all data into an ElTopoMesh
    ElTopoMesh eltopo_time0;
    eltopo_time0.num_vertices = num_vertices;
    eltopo_time0.vertex_locations = V0;
    eltopo_time0.num_triangles = num_triangles;
    eltopo_time0.triangles = F0;
    eltopo_time0.vertex_masses = masses;

    // Mesh at time 1
    NonDestructiveTriMesh mesh_time1;
    std::vector<Vec3d> x1;
    // Read mesh_time0 from OBJ file
    read_objfile(mesh_time1, x1, argv[2]);
//    // Read mesh_time0 from BIN file
//    read_binary_pos_tris_file(mesh_time1, x1, argv[2]);
    // get vertex positions
    double V1[3*num_vertices];
    i = 0;
    for (std::vector<Vec3d>::iterator it = x1.begin() ; it != x1.end(); ++it){
        Vec3d vertex = *it;
        V1[3*i] = vertex[0];
        V1[3*i+1] = vertex[1];
        V1[3*i+2] = vertex[2];
        i = i + 1;
    }

    // Set general parameters
    ElTopoGeneralOptions sim_general_options;
    // do not print stuff to the console
    sim_general_options.m_verbose = false;
    // do avoid self-intersections
    sim_general_options.m_collision_safety = 1;
    sim_general_options.m_proximity_epsilon = 1e-4;

    // Set Simulation parameters
    ElTopoIntegrationOptions sim_integration_options;
    sim_integration_options.m_friction_coefficient = 0.0;
    sim_integration_options.m_dt = 1.0;

    // run simulation (cut time step if it does not satisfy constraints)
    double* V_final;
    double out_dt = 0.0;
//    el_topo_integrate(&eltopo_time0, V1, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
    double rest_dt = 1.0;
    int attempts = 0;
    while (rest_dt > 1e-1){
        el_topo_integrate(&eltopo_time0, V1, &sim_general_options, &sim_integration_options, &V_final, &out_dt);
        eltopo_time0.vertex_locations = V_final;
        std::cout << "out_dt = " << out_dt << std::endl;
        rest_dt = (1-out_dt)*rest_dt;
        std::cout << "rest_dt = " << rest_dt << std::endl;
        attempts = attempts+1;

        if (out_dt<1e-1){
            // it didn't find a time step, throw mex error
            return 0;
        }

    }

//    std::cout << "actual dt = " << out_dt << std::endl;

    // Output corrected vertices to x1
    i = 0;
    std::vector<Vec3d> x2(num_vertices);
    for (std::vector<Vec3d>::iterator it = x2.begin() ; it != x2.end(); ++it){
        Vec3d vertex;
        vertex[0] = V_final[3*i];
        vertex[1] = V_final[3*i+1];
        vertex[2] = V_final[3*i+2];
        *it = vertex;
//        std::cout << i << " " << "V_final = " << V_final[3*i] << " " << V_final[3*i+1] << " " << V_final[3*i+2] << " V1 = " << " " << V1[3*i] << " " << V1[3*i+1] << " " << V1[3*i+2]  << std::endl;
        i = i + 1;
    }

    // Save on argv[5] obj path
    write_objfile(mesh_time0,x2,argv[4]);

    return 0;

}
