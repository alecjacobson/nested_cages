
.PHONY: all

all: obj progressive_hulls

.PHONY: progressive_hulls

#CXX=clang++ -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0
CXX=g++
#CFLAGS+=-DIGL_STATIC_LIBRARY -fopenmp
#CFLAGS+=-fopenmp
CFLAGS+=-g -O3 -DNDEBUG
#CFLAGS+=-g -O0 

LIBIGL=/usr/local/igl/libigl/
LIBIGL_INC=-I$(LIBIGL)/include 
#LIBIGL_INC=-I$(LIBIGL)/include -DIGL_STATIC_LIBRARY
#LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglmosek -liglviewer

EIGEN3_INC=-I/usr/local/include/eigen3 -I/usr/local/include/eigen3/unsupported

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include -I$(LIBIGL)/external/AntTweakBar/src
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar
GLFW_INC=-I$(LIBIGL)/external/nanogui/ext/glfw/include 
GLFW_LIB=-L$(LIBIGL)/external/glfw/src -lglfw3 -framework Carbon -framework QuartzCore -framework IOKit
NANOGUI_INC=-I$(LIBIGL)/external/nanogui/include  -I$(LIBIGL)/external/nanogui/ext/nanovg/src
NANOGUI_LIB=-L$(LIBIGL)/external/nanogui/build -lnanogui

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	MOSEKPLATFORM=osx64x86
	MOSEK=/usr/local/mosek
	MOSEK_INC=-I$(MOSEK)/7/tools/platform/$(MOSEKPLATFORM)/h
	MOSEK_LIB=-L$(MOSEK)/7/tools/platform/$(MOSEKPLATFORM)/bin -lmosek64 -liomp5 -lpthread
	ANTTWEAKBAR_LIB+=-framework AppKit
endif

OPENGL_LIB=-framework OpenGL

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(GLFW_INC) $(NANOGUI_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(GLFW_LIB) $(NANOGUI_LIB)

CPP_FILES=progressive_hulls.cpp
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11

progressive_hulls: obj $(OBJ_FILES)
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -o progressive_hulls $(OBJ_FILES) $(LIB)

obj:
	mkdir -p obj

obj/%.o: %.cpp
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f progressive_hulls
