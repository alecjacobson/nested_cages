cmake_minimum_required(VERSION 2.8.6)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


find_package(EIGEN REQUIRED)
#set(CGAL_BRANCH_BUILD ON)
#add_subdirectory("${PROJECT_SOURCE_DIR}/cgal/Installation")
find_library(BLAS REQUIRED)
if(APPLE)
  set(CMAKE_CXX_LINK_FLAGS "-framework Accelerate")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
# We need C++11. Put this directive after CGAL's include.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wextra -Wno-reorder -Wno-unknown-pragmas -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wno-deprecated-declarations -march=native")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif() 

add_subdirectory("${PROJECT_SOURCE_DIR}/meshfix")
include_directories( "${EIGEN_INCLUDE_DIR}" )
include_directories( "${PROJECT_SOURCE_DIR}/libigl/include/" )
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
include_directories("${PROJECT_SOURCE_DIR}/tetgen")
add_definitions(-DTETLIBRARY)

link_directories(
  /usr/local/lib
  /opt/local/lib
  {EIGEN_DIRS}
)

file(GLOB SRCFILES                    "${PROJECT_SOURCE_DIR}/*.cpp")

# Tetgen, CollisionDetection, and Eltopo don't use CMAKE, so add their source
# files directly
file(GLOB TETGEN_SRCFILES             "${PROJECT_SOURCE_DIR}/tetgen/*.cxx")
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
  ${BLAS_LIBRARIES} 
  )