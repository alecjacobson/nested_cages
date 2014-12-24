Readme: Matlab interface for Etienne's Velocity Filter

1) Clone collisiondetection repo:
https://github.com/evouga/collisiondetection

2) Build libcollisions library by issuing 'make' in the 'collisiondetection' folder.

3) Mexify velocity_filter_mex (change libcollisions path and include below)
Leo''s Macbook
mex velocity_filter_mex.cpp -I/usr/local/include/eigen3 -I/Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/include /Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/bin/libcollisions.a
mex self_distance_mex.cpp -I/usr/local/include/eigen3 -I/Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/include /Users/Leo/PHD_Work/Cage_Generation_2013/code/collisiondetection/bin/libcollisions.a
