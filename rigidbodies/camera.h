#ifndef CAMERA_H
#define CAMERA_H
#include <Eigen/Core>
#include <QtOpenGL>
class Camera
{
public:
    Camera();

    void setPerpective(double fovy, double aspect);
    void setViewport(int width, int height);
    void getViewport(int &width, int &height) const;
    void setZClipping(double near, double far);
    void getZClipping(double &near, double &far) const;

    void applyProjection() const;
    void applyViewport() const;
    void applyLookAt() const;

    void setDefault3D();
    void translateEye(const Eigen::Vector3d &v);
    void translateCenter(const Eigen::Vector3d &v);
    void orbit(const Eigen::Matrix3d &M);
    void orbitCenter(const Eigen::Matrix3d &M);

    void getSpanningSet(Eigen::Vector3d &right, Eigen::Vector3d &up, Eigen::Vector3d &center);

    const Eigen::Vector3d &getEye() const;
    const Eigen::Vector3d &getCenter() const;
    void setCenter(const Eigen::Vector3d &center);
    void project(const Eigen::Vector3d &pos, double &x, double &y);

    void getPixelAt(const Eigen::Vector2d &pos, GLubyte *pixelbuf) const;

private:
    Eigen::Vector3d viewCenter_;
    Eigen::Vector3d eye_;
    Eigen::Vector3d up_;

    double p_fovy_;
    double p_aspect_;
    double z_near_;
    double z_far_;
    double v_width_;
    double v_height_;
};
#endif // CAMERA_H
