#include "controller.h"
#include "mainwindow.h"
#include "simulation.h"
#include <QDebug>
#include <Eigen/Core>

using namespace Eigen;

Controller::Controller(int fps) : QThread(), mw_(NULL), fps_(fps)
{
}

Controller::~Controller()
{
    delete sim_;
}

void Controller::initialize(MainWindow *mw)
{
    mw_ = mw;
    sim_ = new Simulation(params_);
}

void Controller::initializeGL()
{
    sim_->initializeGL();
}

void Controller::run()
{
    reset();
    connect(&simtimer_, SIGNAL(timeout()), this, SLOT(simTick()));
    simtimer_.start(1000/fps_);
    exec();
}

void Controller::reset()
{
    params_ = SimParameters();
    QMetaObject::invokeMethod(mw_, "setUIFromParameters", Q_ARG(SimParameters, params_));
    clearScene();
}

void Controller::clearScene()
{
    sim_->clearScene();
}

void Controller::updateParameters(SimParameters params)
{
    params_ = params;
}

void Controller::renderPlanes(bool transparent)
{
    sim_->renderPlanes(transparent);
}

void Controller::renderObjects()
{
    sim_->renderObjects();
}

void Controller::mouseClicked(double x, double y, double z, double dx, double dy, double dz)
{
    Vector3d pos(x,y,z);
    Vector3d dir(dx,dy,dz);
    sim_->addRigidBody(pos, dir);
}

void Controller::simTick()
{
    if(params_.simRunning)
        sim_->takeSimulationStep();
}
