
.PHONY: all

all: obj nest_view

.PHONY: nest_view

#CXX=clang++ -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0
CXX=g++
CFLAGS+=-DIGL_STATIC_LIBRARY -fopenmp
CFLAGS+=-g -O3 -DNDEBUG
#CFLAGS+=-g -O0 

LIBIGL=/usr/local/igl/libigl/
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglmatlab -liglembree -liglcgal -liglboolean

CORK=$(LIBIGL)/external/cork/
CORK_INC=-I$(CORK)/include/
CORK_LIB=-L$(CORK)/lib/ -lcork

CGAL=/opt/local/
# Needs to come before matlab include because matlab has a bogus gmp
CGAL_LIB=-L$(CGAL)/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt
CGAL_INC=-I$(CGAL)/include -I/usr/include/
# This is absolutely necessary for Exact Construction
CGAL_FLAGS=-frounding-math -fsignaling-nans 

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit

MATLAB_INC=-I$(MATLAB)/extern/include/
MATLAB_LIB=-L$(MATLAB)/bin/maci64 -lmx -lmat -lmex -lstdc++

CARBON_LIB=-framework Carbon

# Use free glut for mouse scrolling
#FREE_GLUT=/opt/local/
#FREE_GLUT_INC=-I$(FREE_GLUT)/include
#FREE_GLUT_LIB=-L$(FREE_GLUT)/lib -lglut
OPENGL_LIB=-framework OpenGL
GLUT_LIB=-framework GLUT
GLUT_INC=

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(MATLAB_INC) $(GLUT_INC) $(CORK_INC) $(CGAL_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(CGAL_LIB) $(MATLAB_LIB) $(CARBON_LIB) $(CORK_LIB) 

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11

nest_view: obj $(OBJ_FILES)
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -o nest_view $(OBJ_FILES) $(LIB)

obj:
	mkdir -p obj

obj/%.o: %.cpp
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f nest_view
