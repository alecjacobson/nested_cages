
.PHONY: all

all: obj nested_cages

.PHONY: nested_cages

CXX=g++
CFLAGS+=-fopenmp
CFLAGS+=-g -O3 -DNDEBUG

LIBIGL=/home/leo/PHD_Work/Cage_Generation_2013/libigl/
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -liglcgal -ligl

EIGEN3_INC=-I/usr/include/Eigen -I/usr/include/Eigen/unsupported

GLFW_INC=-I$(LIBIGL)/external/glfw/include 
GLFW_LIB=-L$(LIBIGL)/external/glfw/lib -lglfw3 -framework Carbon -framework QuartzCore -framework IOKit

OPENGL_LIB=-framework OpenGL

INC=-I/usr/include/ $(LIBIGL_INC) $(EIGEN3_INC)
LIB=$(LIBIGL_LIB) -L/usr/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread -lboost_system

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11

nested_cages: obj obj/io.o obj/cgal.o obj/flow.o obj/nested_cages.o
	${CXX} $(CFLAGS) -o nested_cages obj/io.o obj/cgal.o obj/flow.o obj/nested_cages.o $(LIB)

obj:
	mkdir -p obj

obj/%.o: %.cpp
	${CXX} $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f nested_cages
