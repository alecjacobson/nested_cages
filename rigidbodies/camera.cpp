#include "camera.h"
#include <Eigen/Geometry>
#include <QGLWidget>
#include <GL/glu.h>
#include <iostream>

using namespace Eigen;

Camera::Camera()
{
    viewCenter_ << 1, 0, 0;
    eye_ << 0, 0, 0;
    up_ << 0, 0, 1;

    p_fovy_ = 60;
    p_aspect_ = 1;

    z_near_ = 0.1;
    z_far_ = 200;
}

void Camera::setPerpective(double fovy, double aspect)
{
    p_fovy_ = fovy;
    p_aspect_ = aspect;
}

void Camera::setViewport(int width, int height)
{
    v_width_ = width;
    v_height_ = height;
}

void Camera::getViewport(int &width, int &height) const
{
    width = v_width_;
    height = v_height_;
}

void Camera::setZClipping(double near, double far)
{
    z_near_ = near;
    z_far_ = far;
}

void Camera::getZClipping(double &near, double &far) const
{
    near = z_near_;
    far = z_far_;
}

void Camera::setDefault3D()
{
    viewCenter_ = Vector3d(1,0,5);
    eye_ = Vector3d(0,0,5);
}

void Camera::translateEye(const Eigen::Vector3d &v)
{
    eye_ += v;
}

void Camera::translateCenter(const Eigen::Vector3d &v)
{
    viewCenter_ += v;
}

void Camera::orbit(const Eigen::Matrix3d &M)
{
    const Vector3d e = eye_ - viewCenter_;
    eye_ = M.transpose() * e;
    eye_ += viewCenter_;

    up_ = M.transpose() * up_;
}

void Camera::orbitCenter(const Eigen::Matrix3d &M)
{
    const Vector3d e = viewCenter_ - eye_;
    viewCenter_ = M.transpose()*e;
    viewCenter_ += eye_;

    up_ = M.transpose() * up_;
}

void Camera::getSpanningSet(Vector3d &right, Vector3d &up, Vector3d &center)
{
    center = viewCenter_ - eye_;
    center.normalize();
    up = up_ - up_.dot(center) * center;
    up.normalize();
    right = center.cross(up);
}

void Camera::applyProjection() const
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(p_fovy_, p_aspect_*v_width_/v_height_, z_near_, z_far_);
    glMatrixMode(GL_MODELVIEW);
}

void Camera::getPixelAt(const Vector2d &pos, GLubyte *pixelbuf) const
{
    glReadPixels((1.0+pos[0])*v_width_/2.0, (1.0+pos[1])*v_height_/2.0, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, pixelbuf);
}

void Camera::applyViewport() const
{
    glViewport(0, 0, v_width_, v_height_);
}

void Camera::applyLookAt() const
{
    glLoadIdentity();
    gluLookAt(eye_[0], eye_[1], eye_[2], viewCenter_[0], viewCenter_[1], viewCenter_[2], up_[0], up_[1], up_[2]);
}

const Vector3d &Camera::getEye() const
{
    return eye_;
}

const Vector3d &Camera::getCenter() const
{
    return viewCenter_;
}

void Camera::setCenter(const Vector3d &center)
{
    viewCenter_ = center;
}

void Camera::project(const Vector3d &pos, double &x, double &y)
{
    double z;
    applyProjection();
    double model[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    double proj[16];
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    int view[4];
    view[0] = -1;
    view[1] = -1;
    view[2] = 2;
    view[3] = 2;
    gluProject(pos[0], pos[1], pos[2], model, proj, view, &x, &y, &z);
}
