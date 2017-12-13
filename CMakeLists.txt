cmake_minimum_required(VERSION 2.8.6)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(APPLE)
  set(CMAKE_CXX_LINK_FLAGS "-framework Accelerate")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

# We need C++11. Put this directive after CGAL's include.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wno-deprecated-declarations -march=native")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif() 

SET(MESHFIX_WITH_LIBIGL OFF CACHE BOOL "Use some expat option")
add_subdirectory("${PROJECT_SOURCE_DIR}/meshfix")
#include_directories( "${PROJECT_SOURCE_DIR}/eigen/" )
include_directories( "${PROJECT_SOURCE_DIR}/libigl/include/" )
include_directories( "${PROJECT_SOURCE_DIR}/libigl/external/nanogui/ext/eigen/" )
include_directories( "${PROJECT_SOURCE_DIR}/cgal/../")
include_directories( "${PROJECT_SOURCE_DIR}/collisiondetection/include" )
include_directories( "${PROJECT_SOURCE_DIR}/eltopo/common/" )
include_directories( "${PROJECT_SOURCE_DIR}/eltopo/common/tunicate" )
include_directories( "${PROJECT_SOURCE_DIR}/eltopo/common/newsparse" )
include_directories( "${PROJECT_SOURCE_DIR}/eltopo/eltopo3d/" )
set(MESHFIX_DIR "${PROJECT_SOURCE_DIR}/meshfix/" )
set(JMESHEXT_DIR "${MESHFIX_DIR}/JMeshExt-1.0alpha_src/")
set(JMESH_DIR "${JMESHEXT_DIR}/JMeshLib-1.2/")
include_directories("${MESHFIX_DIR}")
include_directories("${JMESHEXT_DIR}/include")
include_directories("${JMESH_DIR}/include")
set(TETGEN_DIR "${PROJECT_SOURCE_DIR}/libigl/external/tetgen")
include_directories("${TETGEN_DIR}")
add_definitions(-DTETLIBRARY)

option(VERBOSE_DEBUG "VERBOSE_DEBUG" OFF)
if (VERBOSE_DEBUG)
  add_definitions(-DVERBOSE_DEBUG)
endif ()

file(GLOB SRCFILES                    "${PROJECT_SOURCE_DIR}/*.cpp")

# Tetgen, CollisionDetection, and Eltopo don't use CMAKE, so add their source
# files directly
file(GLOB TETGEN_SRCFILES             "${TETGEN_DIR}/*.cxx")
file(GLOB COLLISIONDETECTION_SRCFILES "${PROJECT_SOURCE_DIR}/collisiondetection/src/*.cpp")
file(GLOB ELTOPO_SRCFILES             "${PROJECT_SOURCE_DIR}/eltopo/eltopo3d/*.cpp"
                                      "${PROJECT_SOURCE_DIR}/eltopo/common/*.cpp"
                                      "${PROJECT_SOURCE_DIR}/eltopo/common/newsparse/*.cpp"
                                      "${PROJECT_SOURCE_DIR}/eltopo/tunicate/*.cpp")
# Don't build the eltopo viewers that depends on OpenGL and GLUT
list(REMOVE_ITEM ELTOPO_SRCFILES "${PROJECT_SOURCE_DIR}/eltopo/common/gluvi.cpp")
list(REMOVE_ITEM ELTOPO_SRCFILES "${PROJECT_SOURCE_DIR}/eltopo/eltopo3d/meshrenderer.cpp")
add_definitions(-D__LITTLE_ENDIAN__ -DUSE_FORTRAN_BLAS -DNO_GUI -DGRID_100 -DEXTRA_PASSES -DREMOVE_RIZ -fPIC)

add_executable(nested_cages 
  ${SRCFILES} 
  ${TETGEN_SRCFILES}
  ${COLLISIONDETECTION_SRCFILES}
  ${ELTOPO_SRCFILES})

target_link_libraries(
  nested_cages 
  meshfix 
  jmesh 
  jmeshext 
  nl 
  superlu
  ${CGAL_LIBRARIES} 
  ${BLAS_LIBRARIES} 
  ${LAPACK_LIBRARIES} 
  )
