#ifndef RIGIDBODYINSTANCE_H
#define RIGIDBODYINSTANCE_H

#include <Eigen/Core>

class RigidBodyTemplate;

class RigidBodyInstance
{
public:
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta, double density);

    void render();

    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;

    double density;

    const RigidBodyTemplate &getTemplate() const {return rbtemplate_;}

private:

    const RigidBodyTemplate &rbtemplate_;
};

#endif // RIGIDBODYINSTANCE_H
