Readme: Matlab interface for Etienne's Velocity Filter

1) Clone collisiondetection repo:
https://github.com/evouga/collisiondetection

2) Build libcollisions library by issuing 'make' in the 'collisiondetection' folder.

3) Mexify velocity_filter_mex (change libcollisions path and include below)
Leo''s Macbook
mex velocity_filter_mex.cpp -I/usr/local/include/eigen3 -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include -I/Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/include /Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/bin/libcollisions.a
mex self_distance_mex.cpp -I/usr/local/include/eigen3 -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include -I/Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/include /Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/bin/libcollisions.a
mex testNewSequence_mex.cpp -I/usr/local/include/eigen3 -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include -I/Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/include /Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/bin/libcollisions.a
mex inflate_mex.cpp -I/usr/local/include/eigen3 -I/Users/Leo/PHD_Work/Volume_Meshing_2012/libigl/include -I/Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/include /Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/bin/libcollisions.a

mex velocity_filter_mex.cpp -I/opt/local/include/eigen3 -I/usr/local/igl/libigl/include -I/Users/ajx/Documents/collisiondetection/include /Users/ajx/Documents/collisiondetection/bin/libcollisions.a
mex self_distance_mex.cpp   -I/opt/local/include/eigen3 -I/usr/local/igl/libigl/include -I/Users/ajx/Documents/collisiondetection/include /Users/ajx/Documents/collisiondetection/bin/libcollisions.a
mex testNewSequence_mex.cpp -I/opt/local/include/eigen3 -I/usr/local/igl/libigl/include -I/Users/ajx/Documents/collisiondetection/include /Users/ajx/Documents/collisiondetection/bin/libcollisions.a
mex inflate_mex.cpp         -I/opt/local/include/eigen3 -I/usr/local/igl/libigl/include -I/Users/ajx/Documents/collisiondetection/include /Users/ajx/Documents/collisiondetection/bin/libcollisions.a
