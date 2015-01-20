#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include "SOIL.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "vectormath.h"
#include <Eigen/Dense>
#include "mesh.h"
#include "quadprog/eiquadprog.hpp"
#include "rigidbodyctcd.h"

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0), wallTex_(0)
{
    loadRigidBodies();
}

Simulation::~Simulation()
{
    clearScene();
    for(vector<RigidBodyTemplate *>::iterator it = templates_.begin(); it != templates_.end(); ++it)
    {
        delete *it;
    }
}

void Simulation::initializeGL()
{
    loadFloorTexture();
    loadWallTexture();
}

void Simulation::loadRigidBodies()
{
    const int numobjs = 4;
    string objNames[numobjs] = {"resources/sphere.obj", "resources/2by4.obj", "resources/bunny.obj", "resources/octopus.obj"};
    string cageObjNames[numobjs] = {"resources/sphere.obj", "resources/2by4.obj", "resources/bunny.obj", "resources/octopus_cage.obj"};
    for(int i=0; i<numobjs; i++)
    {
        RigidBodyTemplate *rbt = new RigidBodyTemplate(objNames[i], cageObjNames[i]);
        templates_.push_back(rbt);
    }
}

void Simulation::loadFloorTexture()
{
    floorTex_ = SOIL_load_OGL_texture("resources/grid.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(floorTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}

void Simulation::loadWallTexture()
{
    wallTex_ = SOIL_load_OGL_texture("resources/wall.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(wallTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, wallTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}

void Simulation::renderPlanes(bool transparent)
{
    renderLock_.lock();

    glEnable(GL_CULL_FACE);
    if(transparent)
    {
        glCullFace(GL_FRONT);
        glColor4f(1.0, 1.0, 1.0, 0.1);
    }
    else
    {
        glCullFace(GL_BACK);
        glColor4f(1.0, 1.0, 1.0, 1.0);
    }

    for(vector<Plane>::iterator it = planes_.begin(); it != planes_.end(); ++it)
        renderPlane(*it, it == planes_.begin());

    glDisable(GL_CULL_FACE);

    renderLock_.unlock();
}

void Simulation::renderPlane(const Plane &p, bool isFloor)
{    
    if(isFloor && floorTex_)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else if(!isFloor && wallTex_)
    {
        glBindTexture(GL_TEXTURE_2D, wallTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else
        glColor3f(0.5, 0.5, 0.5);

    double texsize = 5.0;
    double gridsize = 1000.0;

    double texmax = gridsize/texsize;

    Vector3d tangent1 = VectorMath::perpToAxis(p.normal);
    Vector3d tangent2 = tangent1.cross(p.normal);

    Vector3d corner;

    glBegin(GL_QUADS);
    {
        glTexCoord2f(texmax, texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos + gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(texmax, -texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos + gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, -texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos - gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, texmax);
        glNormal3f(p.normal[0], p.normal[1], p.normal[2]);
        corner = p.pos - gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);
    }
    glDisable(GL_TEXTURE_2D);
    glEnd();
}

void Simulation::renderObjects()
{
    renderLock_.lock();
    {
        for(vector<RigidBodyInstance *>::iterator it = bodies_.begin(); it != bodies_.end(); ++it)
        {
            (*it)->render();
        }
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    time_ += params_.timeStep;

    vector<Matrix3d> Roldthetas;

    vector<Vector3d> newc(bodies_.size());
    vector<Vector3d> newcvel(bodies_.size());
    vector<Vector3d> newtheta(bodies_.size());
    vector<Vector3d> neww(bodies_.size());
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        newc[bodyidx] = bodies_[bodyidx]->c;
        newtheta[bodyidx] = bodies_[bodyidx]->theta;
        newcvel[bodyidx] = bodies_[bodyidx]->cvel;
        neww[bodyidx] = bodies_[bodyidx]->w;
        Matrix3d Rtheta = VectorMath::rotationMatrix(bodies_[bodyidx]->theta);
        Roldthetas.push_back(Rtheta);
    }
    vector<ContactInfo> contacts;
    while(true)
    {
        bool newcollisions = false;
        // step
        for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
        {
            RigidBodyInstance &body = *bodies_[bodyidx];
            newc[bodyidx] = body.c + params_.timeStep*newcvel[bodyidx];
            Matrix3d Rhw = VectorMath::rotationMatrix(params_.timeStep*neww[bodyidx]);
            Matrix3d Rtheta = VectorMath::rotationMatrix(newtheta[bodyidx]);
            newtheta[bodyidx] = VectorMath::axisAngle(Rhw*Rtheta);

            newcollisions |= RigidBodyCTCD::detectCollisions(bodies_, planes_, newc, newtheta, contacts, params_.useCage);
        }

        if(!newcollisions)
            break;


        // solve
        int numconstraints = contacts.size();
        VectorXd rhs(numconstraints);
        SparseMatrix<double> M(numconstraints, numconstraints);
        for(int i=0; i<numconstraints; i++)
        {
            ContactInfo &ci = contacts[i];
            Matrix3d Rtheta1 = VectorMath::rotationMatrix(bodies_[ci.body1]->theta);
            Matrix3d Rtheta2;
            if(ci.body2 != -1)
                Rtheta2 = VectorMath::rotationMatrix(bodies_[ci.body2]->theta);
            Vector3d relvel = bodies_[ci.body1]->cvel + (bodies_[ci.body1]->w).cross(Rtheta1*ci.pt1);
            if(ci.body2 != -1)
                relvel = relvel - bodies_[ci.body2]->cvel - (bodies_[ci.body2]->w).cross(Rtheta2*ci.pt2);
            rhs[i] = fabs(2.0*relvel.dot(ci.n));

            double m1 = bodies_[ci.body1]->density * bodies_[ci.body1]->getTemplate().getVolume();
            double m2;
            if(ci.body2 != -1)
                m2 = bodies_[ci.body2]->density * bodies_[ci.body2]->getTemplate().getVolume();

            Matrix3d MI1 = bodies_[ci.body1]->getTemplate().getInertiaTensor();
            Matrix3d MI2;
            if(ci.body2 != -1)
                MI2 = bodies_[ci.body2]->getTemplate().getInertiaTensor();

            Vector3d impulsew1 = Rtheta1*MI1.inverse()*ci.pt1.cross(Rtheta1.transpose()*ci.n)/bodies_[ci.body1]->density;
            Vector3d impulsew2;
            if(ci.body2 != -1)
                impulsew2 = Rtheta2*MI2.inverse()*ci.pt2.cross(Rtheta2.transpose()*-ci.n)/bodies_[ci.body2]->density;

            for(int j=0; j<numconstraints; j++)
            {
                if(ci.body1 == contacts[j].body1)
                {
                    double val = ci.n.dot(contacts[j].n)/m1;
                    val += impulsew1.cross(Rtheta1*contacts[j].pt1).dot(contacts[j].n);
                    M.coeffRef(j,i) += val;
                }
                if(ci.body1 == contacts[j].body2)
                {
                    double val = ci.n.dot(-contacts[j].n)/m1;
                    val += impulsew1.cross(Rtheta1*contacts[j].pt2).dot(-contacts[j].n);
                    M.coeffRef(j,i) += val;
                }
                if(ci.body2 != -1)
                {
                    if(ci.body2 == contacts[j].body1)
                    {
                        double val = -ci.n.dot(contacts[j].n)/m2;
                        val += impulsew2.cross(Rtheta2*contacts[j].pt1).dot(contacts[j].n);
                        M.coeffRef(j,i) += val;
                    }
                    if(ci.body2 == contacts[j].body2)
                    {
                        double val = -ci.n.dot(-contacts[j].n)/m2;
                        val += impulsew2.cross(Rtheta2*contacts[j].pt2).dot(-contacts[j].n);
                        M.coeffRef(j,i) += val;
                    }
                }
            }
        }
        M.makeCompressed();
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver(M);
        VectorXd lambdas = solver.solve(rhs);
        for(int i=0; i<(int)bodies_.size(); i++)
        {
            newcvel[i] = bodies_[i]->cvel;
            neww[i] = bodies_[i]->w;
        }
        for(int i=0; i<numconstraints; i++)
        {
            ContactInfo &ci = contacts[i];
            double m1 = bodies_[ci.body1]->density * bodies_[ci.body1]->getTemplate().getVolume();
            double m2;
            if(ci.body2 != -1)
                m2 = bodies_[ci.body2]->density * bodies_[ci.body2]->getTemplate().getVolume();

            Matrix3d Rtheta1 = VectorMath::rotationMatrix(bodies_[ci.body1]->theta);
            Matrix3d Rtheta2;
            if(ci.body2 != -1)
                Rtheta2 = VectorMath::rotationMatrix(bodies_[ci.body2]->theta);

            Matrix3d MI1 = bodies_[ci.body1]->getTemplate().getInertiaTensor();
            Matrix3d MI2;
            if(ci.body2 != -1)
                MI2 = bodies_[ci.body2]->getTemplate().getInertiaTensor();

            Vector3d impulsew1 = Rtheta1*MI1.inverse()*ci.pt1.cross(Rtheta1.transpose()*ci.n)/bodies_[ci.body1]->density;
            neww[ci.body1] += impulsew1*lambdas[i];
            newcvel[ci.body1] += ci.n/m1*lambdas[i];

            if(ci.body2 != -1)
            {
                Vector3d impulsew2 = Rtheta2*MI2.inverse()*ci.pt2.cross(Rtheta2.transpose()*-ci.n)/bodies_[ci.body2]->density;
                neww[ci.body2] -= impulsew2*lambdas[i];
                newcvel[ci.body2] -= ci.n/m2*lambdas[i];
            }
        }
    }

    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        body.c = newc[bodyidx];
        body.theta = newtheta[bodyidx];
        body.w = neww[bodyidx];
        body.cvel = newcvel[bodyidx];
    }


    VectorXd cForce;
    VectorXd thetaForce;

    computeForces(cForce, thetaForce);

    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        Matrix3d Mi = body.getTemplate().getInertiaTensor();

        body.cvel += params_.timeStep*cForce.segment<3>(3*bodyidx)/body.density/body.getTemplate().getVolume();

        Vector3d newwguess(body.w);
        Matrix3d &Roldtheta = Roldthetas[bodyidx];

        int iter = 0;
        for(iter=0; iter<params_.NewtonMaxIters; iter++)
        {
            Matrix3d Dw1 = -VectorMath::TMatrix(params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(-body.theta);
            Matrix3d Dw2 = VectorMath::TMatrix(-params_.timeStep*body.w).inverse() * VectorMath::TMatrix(-body.theta);
            Matrix3d DRw = VectorMath::DrotVector(-body.theta, newwguess);
            Matrix3d Rnewtheta = VectorMath::rotationMatrix(body.theta);
            Vector3d fval = body.density * Dw1.transpose()*Rnewtheta*Mi*Rnewtheta.transpose()*newwguess;
            fval += -params_.timeStep*body.density * DRw.transpose() * Mi * Rnewtheta.transpose() * newwguess;
            fval += body.density * Dw2.transpose() * Roldtheta * Mi * Roldtheta.transpose() * body.w;
            fval += params_.timeStep*thetaForce.segment<3>(3*bodyidx);

            if(fval.norm() < params_.NewtonTolerance)
                break;

            Matrix3d Df = body.density * Dw1.transpose()*Rnewtheta*Mi*Rnewtheta.transpose();// + -params_.timeStep*(*it)->density * DRw.transpose() * Mi * Rnewtheta.transpose();

            Vector3d deltaw = Df.inverse() * (-fval);
            newwguess += deltaw;
        }
        body.w = newwguess;
    }


    if(spawnTime_ < time_)
    {
        spawnTime_ += 0.5;
        Vector3d spawnpos(0,0,5.0);
        Vector3d orient;
        for(int i=0; i<3; i++)
            orient[i] = 2.0*VectorMath::randomUnitIntervalReal()-1.0;

        RigidBodyInstance *newbody = new RigidBodyInstance(*templates_[SimParameters::R_CUSTOM], spawnpos, orient, params_.bodyDensity);
        newbody->cvel.setZero();
        newbody->w.setZero();

        bodies_.push_back(newbody);
    }
}

void Simulation::clearScene()
{
    time_ = 0;
    spawnTime_ = 0;
    renderLock_.lock();
    {
        for(vector<RigidBodyInstance *>::iterator it = bodies_.begin(); it != bodies_.end(); ++it)
            delete *it;
        bodies_.clear();
        planes_.clear();
        Plane groundPlane;
        groundPlane.pos << 0,0,0;
        groundPlane.normal << 0,0,1;
        planes_.push_back(groundPlane);

        Plane side;
        side.pos << 3,0,0;
        side.normal << -1,0,0.1;
        planes_.push_back(side);

        side.pos << -3,0,0;
        side.normal << 1,0,0.1;
        planes_.push_back(side);

        side.pos << 0,3,0;
        side.normal << 0,-1,0.1;
        planes_.push_back(side);

        side.pos << 0,-3,0;
        side.normal << 0,1,0.1;
        planes_.push_back(side);
    }
    renderLock_.unlock();
}

void Simulation::addRigidBody(Vector3d pos, Vector3d lookdir)
{
    renderLock_.lock();
    {
        if(params_.launchBody == SimParameters::R_PLANE)
        {
            Plane newplane;
            newplane.normal = -lookdir;
            newplane.normal.normalize();
            newplane.pos = pos + 10.0*lookdir;
            planes_.push_back(newplane);
        }
        else
        {
            Vector3d orient(0,0,0);
            if(params_.randomLaunchOrientation)
            {
                for(int i=0; i<3; i++)
                    orient[i] = 2.0*VectorMath::randomUnitIntervalReal()-1.0;
            }
            RigidBodyInstance *newbody = new RigidBodyInstance(*templates_[params_.launchBody], pos, orient, params_.bodyDensity);
            newbody->cvel = lookdir;
            newbody->cvel.normalize();
            newbody->cvel *= params_.launchVel;

            Vector3d orientvel(0,0,0);
            if(params_.randomLaunchAngVel)
            {
                for(int i=0; i<3; i++)
                    orientvel[i] = (2.0*VectorMath::randomUnitIntervalReal()-1.0);
                orientvel.normalize();
                orientvel *= params_.randomLaunchVelMagnitude;
            }
            newbody->w = VectorMath::rotationMatrix(orient)*orientvel;

            bodies_.push_back(newbody);
        }
    }
    renderLock_.unlock();
}

void Simulation::computeForces(VectorXd &Fc, VectorXd &Ftheta)
{
    Fc.resize(3*bodies_.size());
    Ftheta.resize(3*bodies_.size());
    Fc.setZero();
    Ftheta.setZero();

    for(int bodyidx=0; bodyidx<(int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];

        if(params_.activeForces & SimParameters::F_GRAVITY)
            Fc[3*bodyidx+2] += params_.gravityG*body.density*body.getTemplate().getVolume();

    }
}
