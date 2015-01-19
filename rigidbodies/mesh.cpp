#include "mesh.h"
#include "fstream"
#include <vector>
#include <iostream>
#include <cassert>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;

Mesh::Mesh() : numverts_(0), numfaces_(0) {};

Mesh::Mesh(const std::string &filename) : numverts_(0), numfaces_(0)
{
    vector<double> coords;
    vector<int> faces;

    ifstream ifs(filename.c_str());

    while(ifs)
    {
        string s;
        ifs >> s;
        if(s == "v")
        {
            numverts_++;
            for(int j=0; j<3; j++)
            {
                double dummy;
                ifs >> dummy;
                coords.push_back(dummy);
            }
        }
        else if(s == "f")
        {
            numfaces_++;
            for(int j=0; j<3; j++)
            {
                int face;
                ifs >> face;
                char dummy;
                ifs >> dummy;
                ifs >> dummy;
                int fnorm;
                ifs >> fnorm;
                faces.push_back(face-1);
            }
        }
    }

    assert((int)coords.size() == 3*numverts_);
    assert((int)faces.size() == 3*numfaces_);

    verts_.resize(coords.size());

    for(size_t i = 0; i<coords.size(); i++)
        verts_[i] = coords[i];

    faces_.resize(faces.size());

    for(size_t i =0; i < faces.size(); i++)
        faces_[i] = faces[i];

    computeNormals();
}

const Vector3d Mesh::getVert(int idx) const
{
    return verts_.segment<3>(3*idx);
}

const Vector3i Mesh::getFace(int idx) const
{
    return faces_.segment<3>(3*idx);
}

const Vector3d Mesh::getVertNormal(int idx) const
{
    return vertNormals_.segment<3>(3*idx);
}

const Vector3d Mesh::getFaceNormal(int idx) const
{
    return faceNormals_[idx];
}

double Mesh::getFaceArea(int idx) const
{
    return faceAreas_[idx];
}

void Mesh::computeNormals()
{
    faceNormals_.clear();    
    faceAreas_.clear();

    vector<double> verttotalarea;

    vertNormals_.resize(verts_.size());
    vertNormals_.setZero();

    for(int i=0; i < (int)verts_.size()/3; i++)
    {
        verttotalarea.push_back(0);
    }

    for(int i=0; i < (int)faces_.size()/3; i++)
    {
        Vector3i fverts = faces_.segment<3>(3*i);
        Vector3d pts[3];
        for(int j=0; j<3; j++)
            pts[j] = verts_.segment<3>(3*fverts[j]);

        Vector3d normal = (pts[1]-pts[0]).cross(pts[2]-pts[0]);
        double norm = normal.norm();
        faceAreas_.push_back(norm/2.0);
        faceNormals_.push_back(normal/norm);
        for(int j=0; j<3; j++)
        {
            vertNormals_.segment<3>(3*fverts[j]) += normal;
            verttotalarea[fverts[j]] += norm;
        }
    }

    for(int i=0; i<(int)verts_.size()/3; i++)
        vertNormals_.segment<3>(3*i) /= verttotalarea[i];
}

void Mesh::translate(const Vector3d &vec)
{
    for(int i=0; i<(int)verts_.size()/3; i++)
        verts_.segment<3>(3*i) += vec;
}

void Mesh::scale(double s)
{
    for(int i=0; i<(int)verts_.size(); i++)
        verts_[i] *= s;
    for(int i=0; i<(int)faces_.size()/3; i++)
        faceAreas_[i] *= s*s;
}
