cmake_minimum_required(VERSION 2.8.6)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
# We need C++11. Put this directive after CGAL's include.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wextra -Wno-reorder -Wno-unknown-pragmas -march=native")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif() 

find_package(TETGEN REQUIRED) 
include_directories( "${TETGEN_INCLUDE_DIR}" )

find_package(EIGEN REQUIRED)
include_directories( "${EIGEN_INCLUDE_DIR}" )
# Must be found before trying to find libraries in libigl/external
find_package(LIBIGL REQUIRED)
include_directories( "${LIBIGL_INCLUDE_DIR}" )
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED) 

find_package(ELTOPO REQUIRED)
include_directories( "${LIBELTOPO_INCLUDE_DIR}" )
find_package(COMMON REQUIRED)
include_directories( "${LIBCOMMON_INCLUDE_DIR}" )
find_package(TALPA REQUIRED)
include_directories( "${LIBTALPA_INCLUDE_DIR}" "${LIBTALPA_INCLUDE_DIR}/drivers" )
 
link_directories(
  /usr/local/lib
  /opt/local/lib
  {EIGEN_DIRS}
)

set(LIBS 
  ${LIBIGL_LIBRARIES} 
  ${OPENGL_LIBRARIES} 
  ${GLUT_LIBRARIES}
  ${CGAL_LIBRARIES}
  ${SHARED_LIBRARIES}
  ${LIBELTOPO_LIBRARY} -llapack -lblas
  -ltet
)

FILE(GLOB SRCFILES *.cpp)
add_executable(nested_cages_demo ${SRCFILES})
target_link_libraries(nested_cages_demo ${LIBS})
