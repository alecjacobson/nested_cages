#-------------------------------------------------
#
# Project created by QtCreator 2014-09-03T16:42:53
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = collisiontest
TEMPLATE = app

INCLUDEPATH += $$_PRO_FILE_PWD_/eigen_3.2.2/ $$_PRO_FILE_PWD_/SOIL/src /usr/include/suitesparse
LIBS += -lGLU -L$$_PRO_FILE_PWD_/SOIL/lib -lSOIL

QMAKE_CXXFLAGS += -g

SOURCES += main.cpp\
        mainwindow.cpp \
    glpanel.cpp \
    controller.cpp \
    simulation.cpp \
    simparameters.cpp \
    mesh.cpp \
    rigidbodytemplate.cpp \
    rigidbodyinstance.cpp \
    vectormath.cpp \
    camera.cpp \
    rotator.cpp \
    collisions/CTCD.cpp \
    rigidbodyctcd.cpp \
    kdopbroadphase.cpp

HEADERS  += mainwindow.h \
    glpanel.h \
    controller.h \
    simulation.h \
    simparameters.h \
    mesh.h \
    rigidbodytemplate.h \
    rigidbodyinstance.h \
    vectormath.h \
    camera.h \
    rotator.h \
    collisions/Distance.h \
    collisions/CTCD.h \
    quadprog/eiquadprog.hpp \
    rigidbodyctcd.h \
    kdopbroadphase.h

FORMS    += mainwindow.ui

libSOIL.commands = cd $$_PRO_FILE_PWD_/SOIL/projects/makefile/ && make

QMAKE_EXTRA_TARGETS += libSOIL
PRE_TARGETDEPS += libSOIL
