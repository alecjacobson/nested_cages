#ifndef MATH_H
#define MATH_H

#include <Eigen/Core>

class VectorMath
{
public:
    static const Eigen::Matrix3d crossProductMatrix(const Eigen::Vector3d &v);
    static const Eigen::Matrix3d rotationMatrix(const Eigen::Vector3d &axisAngle);
    static const Eigen::Vector3d axisAngle(const Eigen::Matrix3d &rotationMatrix);
    static const Eigen::Vector3d perpToAxis(const Eigen::Vector3d &v);
    static const Eigen::Matrix3d DrotVector(const Eigen::Vector3d &axisangle, const Eigen::Vector3d &rotatingVector);
    static const Eigen::Matrix3d TMatrix(const Eigen::Vector3d &v);
    static double randomUnitIntervalReal();
};

#endif // MATH_H
