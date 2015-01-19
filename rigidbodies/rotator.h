#ifndef ROTATOR_H
#define ROTATOR_H

#include <Eigen/Core>

class Camera;

class Rotator
{
public:
    Rotator(Camera &c);

    void startRotation(const Eigen::Vector2d &pos);
    void updateRotation(const Eigen::Vector2d &pos);
    void stopRotation();

private:
    void rotate(Eigen::Matrix3d &M, const Eigen::Vector3d &axis, double radians) const;

    Camera &c_;
    bool rotating_;
    Eigen::Vector2d startpos_;
};

#endif // ROTATOR_H
