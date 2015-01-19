#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <QMutex>
#include "simparameters.h"
#include <QGLWidget>

class RigidBodyTemplate;
class RigidBodyInstance;

typedef Eigen::Triplet<double> Tr;

class SimParameters;

struct ContactInfo
{
    int body1;
    int body2;
    Eigen::Vector3d pt1;
    Eigen::Vector3d pt2;
    Eigen::Vector3d n;
};

struct Plane
{
    Eigen::Vector3d pos;
    Eigen::Vector3d normal;
};

class Simulation
{
public:
    Simulation(const SimParameters &params);
    ~Simulation();

    void takeSimulationStep();
    void initializeGL();

    void renderPlanes(bool transparent);
    void renderObjects();
    void clearScene();
    void addRigidBody(Eigen::Vector3d pos, Eigen::Vector3d lookdir);

private:
    void loadFloorTexture();
    void loadWallTexture();
    void loadRigidBodies();

    void renderPlane(const Plane &p, bool isFloor);

    void computeForces(Eigen::VectorXd &Fc, Eigen::VectorXd &Ftheta);

    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    GLuint floorTex_;
    GLuint wallTex_;

    std::vector<RigidBodyTemplate *> templates_;
    std::vector<RigidBodyInstance *> bodies_;

    std::vector<Plane> planes_;
};

#endif // SIMULATION_H
