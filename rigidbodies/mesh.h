#ifndef MESH_H
#define MESH_H

#include <string>
#include <Eigen/Core>
#include <vector>

class Mesh
{
public:
    Mesh();
    Mesh(const std::string &filename);

    int getNumVerts() const {return numverts_;}
    int getNumFaces() const {return numfaces_;}

    const Eigen::Vector3d getVert(int idx) const;
    const Eigen::Vector3d getVertNormal(int idx) const;
    const Eigen::Vector3d getFaceNormal(int idx) const;
    const Eigen::Vector3i getFace(int idx) const;
    double getFaceArea(int idx) const;

    const int *getFacePointer() const {return &faces_[0];}
    const double *getVertexPointer() const {return &verts_[0];}
    const double *getVertexNormalsPointer() const {return &vertNormals_[0];}

    void translate(const Eigen::Vector3d &vec);
    void scale(double s);

private:
    void computeNormals();

    int numverts_;
    int numfaces_;
    Eigen::VectorXd verts_;
    Eigen::VectorXi faces_;
    std::vector<Eigen::Vector3d> faceNormals_;
    Eigen::VectorXd vertNormals_;
    std::vector<double> faceAreas_;
};

#endif // MESH_H
