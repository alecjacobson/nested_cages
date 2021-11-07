cmake_minimum_required(VERSION 3.1)
project(nested_cages)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl + cgal + tetgen
option(LIBIGL_USE_PREBUILT_LIBRARIES "Use prebuilt libraries"       OFF)
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ${LIBIGL_USE_PREBUILT_LIBRARIES})
option(LIBIGL_WITH_CGAL              "Use CGAL"                     ON)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   ON)
find_package(LIBIGL REQUIRED QUIET)

# MeshFix
option(WITH_MESHFIX "Use meshfix" ON)
if(WITH_MESHFIX)
  set(MESHFIX_LIBRARY ON CACHE BOOL "Build MeshFix as library" FORCE)
  set(MESHFIX_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MeshFix-v2.1")
  add_subdirectory("${MESHFIX_DIR}")
endif()

# build eltopo library
set(ELTOPO_DIR "${PROJECT_SOURCE_DIR}/eltopo")
file(GLOB ELTOPO_SRCFILES             "${ELTOPO_DIR}/eltopo3d/*.cpp"
                                      "${ELTOPO_DIR}/common/*.cpp"
                                      "${ELTOPO_DIR}/common/newsparse/*.cpp"
                                      "${ELTOPO_DIR}/tunicate/*.cpp")
# Don't build the eltopo viewers that depends on OpenGL and GLUT
list(REMOVE_ITEM ELTOPO_SRCFILES "${ELTOPO_DIR}/common/gluvi.cpp")
list(REMOVE_ITEM ELTOPO_SRCFILES "${ELTOPO_DIR}/eltopo3d/meshrenderer.cpp")
add_library(eltopo ${ELTOPO_SRCFILES})
target_compile_definitions(eltopo PUBLIC -D__LITTLE_ENDIAN__ -DUSE_FORTRAN_BLAS -DNO_GUI -DGRID_100 -DEXTRA_PASSES -DREMOVE_RIZ)
set_property(TARGET eltopo PROPERTY POSITION_INDEPENDENT_CODE ON)
# eltopo uses #include <myfile.h> instead of #include "myfile.h"
target_include_directories(eltopo PUBLIC "${ELTOPO_DIR}/common/" )
target_include_directories(eltopo PUBLIC "${ELTOPO_DIR}/common/tunicate" )
target_include_directories(eltopo PUBLIC "${ELTOPO_DIR}/common/newsparse" )
target_include_directories(eltopo PUBLIC "${ELTOPO_DIR}/eltopo3d/" )
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(APPLE)
  set(CMAKE_CXX_LINK_FLAGS "-framework Accelerate")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# Build collisiondetection library
set(COLLISIONDETECTION_DIR "${PROJECT_SOURCE_DIR}/collisiondetection")
file(GLOB COLLISIONDETECTION_SRCFILES "${COLLISIONDETECTION_DIR}/src/*.cpp")
add_library(collisiondetection ${COLLISIONDETECTION_SRCFILES})
target_include_directories(collisiondetection PUBLIC "${COLLISIONDETECTION_DIR}/include" )
target_link_libraries(collisiondetection igl::core)

file(GLOB SRCFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SRCFILES})

if(WITH_MESHFIX)
  target_link_libraries(${PROJECT_NAME} MeshFix)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "WITH_MESHFIX")
  # for meshfix_eigen.h
  target_include_directories(${PROJECT_NAME} PUBLIC "${MESHFIX_DIR}/examples")
endif()
target_link_libraries(${PROJECT_NAME} eltopo)
target_link_libraries(${PROJECT_NAME} collisiondetection)
# For each libigl module, link to target (and flip definition flag)
target_link_libraries(${PROJECT_NAME} igl::core)
list(APPEND modules "cgal" "comiso" "cork" "embree" "matlab" "mosek" "opengl" "opengl_glfw" "opengl_glfw_imgui" "png" "tetgen" "triangle" "predicates" "xml")
foreach(module_name IN LISTS modules)
  string(TOUPPER "LIBIGL_WITH_${module_name}" option_string)
  set(option_variable ${${option_string}})
  if(option_variable)
    target_compile_definitions(${PROJECT_NAME} PRIVATE -D${option_string})
    target_link_libraries(${PROJECT_NAME} "igl::${module_name}")
  endif()
endforeach()


#cmake_minimum_required(VERSION 2.8.6)
#
#SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#
#
#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)
#if(APPLE)
#  set(CMAKE_CXX_LINK_FLAGS "-framework Accelerate")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#endif()
#
#set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
#find_package(CGAL REQUIRED)
#include(${CGAL_USE_FILE})
#
## We need C++11. Put this directive after CGAL's include.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wno-deprecated-declarations -march=native")
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
#endif() 
#
#option(WITH_MESHFIX "Use meshfix" ON)
#
#include_directories( "${PROJECT_SOURCE_DIR}/libigl/include/" )
#include_directories( "${PROJECT_SOURCE_DIR}/libigl/external/nanogui/ext/eigen/" )
#include_directories( "${PROJECT_SOURCE_DIR}/cgal/../")
#include_directories( "${PROJECT_SOURCE_DIR}/collisiondetection/include" )
#include_directories( "${PROJECT_SOURCE_DIR}/eltopo/common/" )
#include_directories( "${PROJECT_SOURCE_DIR}/eltopo/common/tunicate" )
#include_directories( "${PROJECT_SOURCE_DIR}/eltopo/common/newsparse" )
#include_directories( "${PROJECT_SOURCE_DIR}/eltopo/eltopo3d/" )
#
#if(WITH_MESHFIX)
#  add_subdirectory("${PROJECT_SOURCE_DIR}/MeshFix-v2.1")
#endif()
#
#option(VERBOSE_DEBUG "VERBOSE_DEBUG" OFF)
#if (VERBOSE_DEBUG)
#  add_definitions(-DVERBOSE_DEBUG)
#endif ()
#
#file(GLOB SRCFILES                    "${PROJECT_SOURCE_DIR}/*.cpp")
#
## Tetgen, CollisionDetection, and Eltopo don't use CMAKE, so add their source
## files directly
#file(GLOB TETGEN_SRCFILES             "${TETGEN_DIR}/*.cxx")
#file(GLOB COLLISIONDETECTION_SRCFILES "${PROJECT_SOURCE_DIR}/collisiondetection/src/*.cpp")
#file(GLOB ELTOPO_SRCFILES             "${PROJECT_SOURCE_DIR}/eltopo/eltopo3d/*.cpp"
#                                      "${PROJECT_SOURCE_DIR}/eltopo/common/*.cpp"
#                                      "${PROJECT_SOURCE_DIR}/eltopo/common/newsparse/*.cpp"
#                                      "${PROJECT_SOURCE_DIR}/eltopo/tunicate/*.cpp")
## Don't build the eltopo viewers that depends on OpenGL and GLUT
#list(REMOVE_ITEM ELTOPO_SRCFILES "${PROJECT_SOURCE_DIR}/eltopo/common/gluvi.cpp")
#list(REMOVE_ITEM ELTOPO_SRCFILES "${PROJECT_SOURCE_DIR}/eltopo/eltopo3d/meshrenderer.cpp")
#add_definitions(-D__LITTLE_ENDIAN__ -DUSE_FORTRAN_BLAS -DNO_GUI -DGRID_100 -DEXTRA_PASSES -DREMOVE_RIZ -fPIC)
#
#add_executable(nested_cages 
#  ${SRCFILES} 
#  ${TETGEN_SRCFILES}
#  ${COLLISIONDETECTION_SRCFILES}
#  ${ELTOPO_SRCFILES})
#
#target_link_libraries(
#  nested_cages 
#  ${CGAL_LIBRARIES} 
#  )
#
#if(WITH_MESHFIX)
#  target_compile_definitions(nested_cages PUBLIC "WITH_MESHFIX")
#target_link_libraries(
#  nested_cages 
#  meshfix 
#  jmesh 
#  jmeshext 
#  nl 
#  superlu
#  ${BLAS_LIBRARIES} 
#  ${LAPACK_LIBRARIES} 
#  )
#endif()
