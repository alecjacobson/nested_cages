cmake_minimum_required(VERSION 3.16)

# Some of the fetched third-party projects (e.g. MeshFix) still declare an old
# cmake_minimum_required(VERSION < 3.5). CMake >= 3.31/4.x refuses those by
# default; this lets them configure without editing their sources.
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING
      "Minimum policy version for old third-party CMake projects")
endif()

project(nested_cages LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# Global settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)      # nested_cages.cpp / eltopo rely on GNU extensions
# Build a PIC core library so it can be reused by the Python extension module.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected, defaulting to Release")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(NESTED_CAGES_BUILD_CLI    "Build the nested_cages command-line tool"      ON)
option(NESTED_CAGES_BUILD_PYTHON "Build the nanobind Python module"              OFF)
option(NESTED_CAGES_FETCH_BLAS   "Fetch & build OpenBLAS instead of find_package" ON)
option(WITH_MESHFIX  "Repair self-intersecting decimations with MeshFix" ON)
option(VERBOSE_DEBUG "Verbose debug output" OFF)

# ---------------------------------------------------------------------------
# Eigen (declared *before* libigl so this version wins the FetchContent race)
# ---------------------------------------------------------------------------
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(eigen)

# ---------------------------------------------------------------------------
# libigl v2.6.0 (transitively fetches CGAL and tetgen)
# ---------------------------------------------------------------------------
set(LIBIGL_INSTALL         OFF CACHE BOOL "" FORCE)
set(LIBIGL_COPYLEFT_CGAL   ON  CACHE BOOL "" FORCE)
set(LIBIGL_COPYLEFT_TETGEN ON  CACHE BOOL "" FORCE)
FetchContent_Declare(libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG v2.6.0)
FetchContent_MakeAvailable(libigl)

# ---------------------------------------------------------------------------
# BLAS / LAPACK for eltopo (USE_FORTRAN_BLAS).  Fetch reference BLAS/LAPACK by
# default so we do not depend on a system-wide install; fall back to
# find_package on request or on Apple (Accelerate).  Reference LAPACK builds the
# standard Fortran `blas` and `lapack` targets that eltopo links against; a
# system OpenBLAS can be substituted via NESTED_CAGES_FETCH_BLAS=OFF for speed.
# ---------------------------------------------------------------------------
if(APPLE)
  set(NESTED_CAGES_BLAS_LIBS "-framework Accelerate")
elseif(NESTED_CAGES_FETCH_BLAS)
  # Build reference BLAS/LAPACK in a fully isolated ExternalProject.  This is
  # deliberately *not* FetchContent/add_subdirectory: Eigen unconditionally
  # defines its own `blas`/`lapack` targets, which would collide.  The isolated
  # build also guarantees the full set of LAPACK routines eltopo needs
  # (dgelsd_, dgesdd_, dsyev_, dgetri_, ...).
  include(ExternalProject)
  set(_lapack_prefix "${CMAKE_BINARY_DIR}/lapack-install")
  set(_lapack_lib "${_lapack_prefix}/lib/liblapack.a")
  set(_blas_lib   "${_lapack_prefix}/lib/libblas.a")
  ExternalProject_Add(reference_lapack
    GIT_REPOSITORY https://github.com/Reference-LAPACK/lapack.git
    GIT_TAG v3.12.0
    GIT_SHALLOW TRUE
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=Release
      -DBUILD_TESTING=OFF
      -DCBLAS=OFF
      -DLAPACKE=OFF
      -DBUILD_SHARED_LIBS=OFF
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON
      -DCMAKE_INSTALL_PREFIX=${_lapack_prefix}
      -DCMAKE_INSTALL_LIBDIR=lib
    BUILD_BYPRODUCTS ${_lapack_lib} ${_blas_lib})
  file(MAKE_DIRECTORY "${_lapack_prefix}/lib")
  add_library(nc_lapack STATIC IMPORTED GLOBAL)
  set_target_properties(nc_lapack PROPERTIES IMPORTED_LOCATION "${_lapack_lib}")
  add_library(nc_blas STATIC IMPORTED GLOBAL)
  set_target_properties(nc_blas PROPERTIES IMPORTED_LOCATION "${_blas_lib}")
  # LAPACK/BLAS are Fortran; link the gfortran runtime.
  set(NESTED_CAGES_BLAS_LIBS nc_lapack nc_blas gfortran)
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  set(NESTED_CAGES_BLAS_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

# ---------------------------------------------------------------------------
# eltopo & collisiondetection (formerly git submodules, no CMake of their own)
# ---------------------------------------------------------------------------
FetchContent_Declare(eltopo
  GIT_REPOSITORY https://github.com/leokollersacht/eltopo.git
  GIT_TAG b2f1dc0667b23f9764cdb8df81d2007a1da7fff1)
FetchContent_Declare(collisiondetection
  GIT_REPOSITORY https://github.com/evouga/collisiondetection.git
  GIT_TAG d7122d61a137984ae12cfd24520c4feb65328449)
FetchContent_MakeAvailable(eltopo collisiondetection)

# --- build eltopo ---
set(ELTOPO_DIR "${eltopo_SOURCE_DIR}")
file(GLOB ELTOPO_SRCFILES "${ELTOPO_DIR}/eltopo3d/*.cpp"
                          "${ELTOPO_DIR}/common/*.cpp"
                          "${ELTOPO_DIR}/common/newsparse/*.cpp"
                          "${ELTOPO_DIR}/tunicate/*.cpp")
# Don't build the eltopo viewers that depend on OpenGL and GLUT
list(REMOVE_ITEM ELTOPO_SRCFILES "${ELTOPO_DIR}/common/gluvi.cpp")
list(REMOVE_ITEM ELTOPO_SRCFILES "${ELTOPO_DIR}/eltopo3d/meshrenderer.cpp")
add_library(eltopo STATIC ${ELTOPO_SRCFILES})
target_compile_definitions(eltopo PUBLIC
  -D__LITTLE_ENDIAN__ -DUSE_FORTRAN_BLAS -DNO_GUI -DGRID_100 -DEXTRA_PASSES -DREMOVE_RIZ)
# eltopo uses #include <myfile.h> instead of #include "myfile.h"
target_include_directories(eltopo PUBLIC
  "${ELTOPO_DIR}/common/"
  "${ELTOPO_DIR}/common/tunicate"
  "${ELTOPO_DIR}/common/newsparse"
  "${ELTOPO_DIR}/eltopo3d/")
target_link_libraries(eltopo PUBLIC ${NESTED_CAGES_BLAS_LIBS})
if(TARGET reference_lapack)
  # ensure the ExternalProject BLAS/LAPACK is built before anything links it
  add_dependencies(eltopo reference_lapack)
endif()

# --- build collisiondetection ---
set(COLLISIONDETECTION_DIR "${collisiondetection_SOURCE_DIR}")
file(GLOB COLLISIONDETECTION_SRCFILES "${COLLISIONDETECTION_DIR}/src/*.cpp")
add_library(collisiondetection STATIC ${COLLISIONDETECTION_SRCFILES})
target_include_directories(collisiondetection PUBLIC "${COLLISIONDETECTION_DIR}/include")
target_link_libraries(collisiondetection PUBLIC igl::core)

# ---------------------------------------------------------------------------
# MeshFix (has its own CMakeLists -> target `MeshFix`)
# ---------------------------------------------------------------------------
if(WITH_MESHFIX)
  set(MESHFIX_LIBRARY ON CACHE BOOL "Build MeshFix as library" FORCE)
  FetchContent_Declare(meshfix
    GIT_REPOSITORY https://github.com/alecjacobson/MeshFix-v2.1.git
    GIT_TAG 2538b17cd6070835f33d1c8bb834d85a38161596
    # guard closeLoadingSession()'s fclose(NULL) crash (see cmake/patch_meshfix.cmake)
    PATCH_COMMAND ${CMAKE_COMMAND} -DMESHFIX_IO_CPP=src/TMesh/io.cpp
                  -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/patch_meshfix.cmake")
  FetchContent_MakeAvailable(meshfix)
  if(TARGET MeshFix)
    # MeshFix is old C/C++ that relies on type punning; -fno-strict-aliasing
    # keeps modern GCC from miscompiling it (segfaults at -O3 otherwise).
    target_compile_options(MeshFix PRIVATE -fno-strict-aliasing)
  endif()
endif()

# ---------------------------------------------------------------------------
# nested_cages core library (algorithm) shared by the CLI and the Python module
# ---------------------------------------------------------------------------
set(NESTED_CAGES_CORE_SRC
  "${CMAKE_CURRENT_SOURCE_DIR}/cgal.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/energy.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/filter.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/flow.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/gradient.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/io.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/reinflate.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/nested_cages_lib.cpp")
add_library(nested_cages_core STATIC ${NESTED_CAGES_CORE_SRC})
target_include_directories(nested_cages_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(nested_cages_core PUBLIC
  igl::core
  igl_copyleft::cgal
  igl_copyleft::tetgen
  eltopo
  collisiondetection)

if(WITH_MESHFIX)
  target_link_libraries(nested_cages_core PUBLIC MeshFix)
  target_compile_definitions(nested_cages_core PUBLIC WITH_MESHFIX)
  # for meshfix_eigen.h
  target_include_directories(nested_cages_core PUBLIC "${meshfix_SOURCE_DIR}/examples")
endif()

if(VERBOSE_DEBUG)
  target_compile_definitions(nested_cages_core PUBLIC VERBOSE_DEBUG)
endif()

# ---------------------------------------------------------------------------
# Command-line tool
# ---------------------------------------------------------------------------
if(NESTED_CAGES_BUILD_CLI)
  add_executable(nested_cages "${CMAKE_CURRENT_SOURCE_DIR}/nested_cages.cpp")
  target_link_libraries(nested_cages PRIVATE nested_cages_core)
endif()

# ---------------------------------------------------------------------------
# Python bindings
# ---------------------------------------------------------------------------
if(NESTED_CAGES_BUILD_PYTHON)
  add_subdirectory(python)
endif()
