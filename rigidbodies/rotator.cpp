#include "rotator.h"
#include "camera.h"

#include <Eigen/Geometry>

using namespace Eigen;

const double PI = 3.1415926535898;

Rotator::Rotator(Camera &c) : c_(c), rotating_(false)
{
}

void Rotator::rotate(Matrix3d &M, const Vector3d &axis, double radians) const
{
    Vector3d n = axis;
    n.normalize();
    double x = n[0];
    double y = n[1];
    double z = n[2];

    double c = cos(radians);
    double s = sin(radians);
    double t = 1-c;

    M(0,0) = t*x*x + c;
    M(0,1) = t*x*y - z*s;
    M(0,2) = t*x*z + y*s;
    M(1,0) = t*x*y + z*s;
    M(1,1) = t*y*y + c;
    M(1,2) = t*y*z - x*s;
    M(2,0) = t*x*z - y*s;
    M(2,1) = t*y*z + x*s;
    M(2,2) = t*z*z + c;
}

void Rotator::startRotation(const Vector2d &pos)
{
    rotating_ = true;
    startpos_ = pos;
 }

void Rotator::updateRotation(const Vector2d &pos)
{
    if(!rotating_)
        return;

    Matrix3d M;
    if(startpos_ != pos)
    {
        Vector3d right, up, center;
        c_.getSpanningSet(right, up, center);
        double coeff = PI/2.0;
        double lr = pos[0] - startpos_[0];
        double ud = pos[1] - startpos_[1];
        Vector3d z(0,0,1.0);
        rotate(M, z, coeff*lr);
        c_.orbitCenter(M);
        rotate(M, z.cross(c_.getCenter() - c_.getEye()), coeff*ud);
        c_.orbitCenter(M);

        startpos_ = pos;
    }
}

void Rotator::stopRotation()
{
    rotating_ = false;
}
