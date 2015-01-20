#ifndef RIGIDBODYCTCD_H
#define RIGIDBODYCTCD_H

#include <vector>
#include <Eigen/Core>
#include "simulation.h"

class RigidBodyInstance;

class RigidBodyCTCD
{
public:
    RigidBodyCTCD();

    static bool detectCollisions(const std::vector<RigidBodyInstance *> &bodies, std::vector<Plane> &planes, std::vector<Eigen::Vector3d> &newc, std::vector<Eigen::Vector3d> &newtheta, std::vector<Eigen::Vector3d> &newcvel, std::vector<Eigen::Vector3d> &neww, std::vector<ContactInfo> &contacts, SimParameters::UseCage  useCage);
    static bool boundingCubesIntersect(const RigidBodyInstance &body1, const RigidBodyInstance &body2, const Eigen::Vector3d &newc1, const Eigen::Vector3d &newc2, const Eigen::Vector3d &newtheta1, const Eigen::Vector3d &newtheta2);
};

#endif // RIGIDBODYCTCD_H
