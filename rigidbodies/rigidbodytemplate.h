#ifndef RIGIDBODYTEMPLATE_H
#define RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>

class Mesh;
class SignedDistanceField;

class RigidBodyTemplate
{
public:
    RigidBodyTemplate(std::string &meshFilename);
    ~RigidBodyTemplate();

    const Mesh &getMesh() const {return *m_;}
    double getVolume() const {return volume_;}
    const Eigen::Vector3d getPrincipleAxis(int axis) const {return principAxes_[axis];}

    const Eigen::Matrix3d getInertiaTensor() const {return inertiaTensor_;}    
    void computeSDF(const char *filename);
    const SignedDistanceField *getSDF() const {return sdf_;}

private:
    RigidBodyTemplate(const RigidBodyTemplate &other);
    RigidBodyTemplate &operator=(const RigidBodyTemplate &other);

    void computeVolume();
    Eigen::Vector3d computeCenterOfMass();
    double boundingSphereRadius();
    void computeInertiaTensor();

    Mesh *m_;

    double volume_;
    Eigen::Matrix3d inertiaTensor_;
    Eigen::Vector3d principAxes_[3];

    SignedDistanceField *sdf_;
};

#endif // RIGIDBODYTEMPLATE_H
