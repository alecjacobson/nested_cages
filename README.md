# Nested Cages 

This C++ project implements:

[Nested Cages] (http://www.cs.columbia.edu/cg/nested-cages/)  
ACM Transactions on Graphics, vol. 34, no. 6 (SIGGRAPH Asia 2015).  
Leonardo Sacht, Etienne Vouga and Alec Jacobson

> Get started with:
>
```bash
git clone --recursive https://github.com/alecjacobson/nested_cages.git
```

[![Bunny teaser from "Nested Cages"](http://www.cs.columbia.edu/cg/nested-cages/bunny-shelf-teaser.jpg)](http://www.cs.columbia.edu/cg/nested-cages/)

## Compilation

This code has been tested on Linux and Mac OS X. In theory this should also
work on Windows.

To compile, 

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
```

This will build all remaining dependencies and the `nested_cages` executable.

### Dependencies

The main dependencies libigl, eltopo and collisiondetection are included as  [git
submodules](https://git-scm.com/docs/git-submodule). If you clone this repo
using `git clone --recursive` then the dependency layout should be:

    nested_cages/
      collisiondetection/
      eltopo/
      eigen/
      libigl/

The cmake build of libigl will further download cgal and tetgen as external
dependencies.


## Example usages

Help information

    ./nested_cages

Obtain 2 volume minimizing nested cages for `../gargo.off`: one regular with
1000 faces and the other regular with 500 faces. Output resulting cages to
`../test_1.off` and `../test_2.off`

    ./nested_cages ../gargo.off 2 1000r 500r None Volume ../test

The same as above, but outputs adaptive decimations (instead of regular)

    ./nested_cages ../gargo.off 2 1000 500 None Volume ../test

Obtain 2 nested cages for `../gargo.off` that minimize surface ARAP energy,
using as input decimations `../gargo_1000.off` and `gargo_500.off`

    ./nested_cages ../gargo.off 2 ../gargo_1000.off ../gargo_500.off SurfARAP None ../test


## Contact

If you have any comments or questions, please contact Leonardo Sacht by e-mail:
leo@mtm.ufsc.br
