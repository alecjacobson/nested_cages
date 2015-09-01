
.PHONY: all

all: obj nest_view

.PHONY: nest_view

#CXX=clang++ -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0
CXX=g++
CFLAGS+=-DIGL_STATIC_LIBRARY -Wno-deprecated-declarations
#CFLAGS+=-fopenmp
CFLAGS+=-g -O3 -DNDEBUG
#CFLAGS+=-g -O0 

LIBIGL=/usr/local/igl/libigl/
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglboolean -liglcgal -liglanttweakbar -liglopengl -liglopengl2 -ligltriangle -liglmatlab
TRIANGLE_INC=-I$(LIBIGL)/external/triangle
TRIANGLE_LIB=-ltriangle
MATLAB_INC=-I/Applications/MATLAB_R2015a.app/extern/include/
MATLAB_LIB=-L/Applications/MATLAB_R2015a.app/bin/maci64/ -lmx -lmat -lmex -leng

IGL_WITH_CSG=0
IGL_WITH_TETGEN=0
ifeq ($(IGL_WITH_CSG),1)
	CFLAGS+=-DWITH_CSG
ifeq ($(IGL_WITH_TETGEN),1)
	LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -ligltetgen
	CFLAGS+=-DWITH_TETGEN
else
  #CFLAGS+=-DIGL_STATIC_LIBRARY -fopenmp
  #CFLAGS+=-DIGL_STATIC_LIBRARY
	CORK=$(LIBIGL)/external/cork/
	CORK_INC=-I$(CORK)/include/
	CORK_LIB=-L$(CORK)/lib/ -lcork
	CGAL=/usr/local/
# Needs to come before matlab include because matlab has a bogus gmp
	CGAL_LIB=-L$(CGAL)/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt
	CGAL_INC=-I$(CGAL)/include -I/usr/include/
# This is absolutely necessary for Exact Construction
	CGAL_FLAGS=-frounding-math -fsignaling-nans 
endif
  TETGEN=$(LIBIGL)/external/tetgen
  TETGEN_LIB=-L$(TETGEN) -ltet 
  TETGEN_INC=-I$(TETGEN)
endif

EIGEN3_INC=-I/usr/local/include/eigen3 -I/usr/local/include/eigen3/unsupported

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	CARBON_LIB=-framework Carbon
	ANTTWEAKBAR_LIB+=-framework AppKit
endif

# Use free glut for mouse scrolling
#FREE_GLUT=/usr/local/
#FREE_GLUT_INC=-I$(FREE_GLUT)/include
#FREE_GLUT_LIB=-L$(FREE_GLUT)/lib -lglut
OPENGL_LIB=-framework OpenGL
GLUT_LIB=-framework GLUT
GLUT_INC=

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(GLUT_INC) $(CORK_INC) $(CGAL_INC) $(TETGEN_INC) $(TRIANGLE_INC) $(MATLAB_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(CGAL_LIB) $(CARBON_LIB) $(CORK_LIB) $(TETGEN_LIB) $(TRIANGLE_LIB) $(MATLAB_LIB)

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11

nest_view: obj $(OBJ_FILES)
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -o nest_view $(OBJ_FILES) $(LIB)

# If you get this error when launching the app:
#     The application cannot be opened because its executable is missing.
# Try running:
#     /System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister -f contact_2d.app
NestView.app: nest_view bunny-logo.icns
	build_app $^ $@
	dylibbundler -od -b -x ./$@/Contents/MacOS/$< -d ./$@/Contents/libs/
	mkdir -p ./$@/Contents/Frameworks/
	cp -r /System/Library/Frameworks/GLUT.framework ./$@/Contents/Frameworks/
	install_name_tool -change /System/Library/Frameworks/GLUT.framework/Versions/A/GLUT @executable_path/../Frameworks/GLUT.framework/Versions/A/GLUT $@/Contents/MacOS/$<
	cp bunny_*.obj ./$@/Contents/Resources/

obj:
	mkdir -p obj

obj/%.o: %.cpp
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	${CXX} $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f nest_view
