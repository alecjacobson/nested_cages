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

## Compilation

This code has been tested on Linux and Mac OS X. In theory this should also
work on Windows.

To compile, first _install CGAL_, then issue:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release  ..
make
```

This will build all remaining dependencies and the `nested_cages` executable.

### Dependencies

_Except CGAL_, all dependencies are included, either explicitly or as [git
submodules](https://git-scm.com/docs/git-submodule). If you clone this repo
using `git clone --recursive` then the dependency layout should be:

    nested_cages/
      collisiondetection/
      eltopo/
      eigen/
      libigl/
      meshfix/
        JMeshExt-1.0alpha_src/
          JMeshLib-1.2/
          OpenNL3.2.1/
            SuperLU/
      tetgen/


## Example usages

Help information

    ./nested_cages

Obtain 2 volume minimizing nested cages for `fertility.off`: one regular with
8000 faces and the other regular with 1000 faces. Output resulting cages to
`test_1.off` and `test_2.off`

    ./nested_cages fertility.off 2 8000r 1000r None Volume test

The same as above, but outputs adaptive decimations (instead of regular)

    ./nested_cages fertility.off 2 8000 1000 None Volume test

Obtain 2 nested cages for fertility.off that minimize surface ARAP energy,
using as input decimations `fertilitiy_8000.off` and `fertilitiy_1000.off`

    ./nested_cages fertility.off 2 fertilitiy_8000.off fertilitiy_1000.off SurfARAP None test

## Contact

If you have any comments or questions, please contact Leonardo Sacht by e-mail:
leo@mtm.ufsc.br
