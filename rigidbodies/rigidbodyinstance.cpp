#include "rigidbodyinstance.h"
#include "vectormath.h"
#include <QGLWidget>
#include "rigidbodytemplate.h"
#include "mesh.h"
#include <Eigen/Geometry>
#include <iostream>

using namespace Eigen;
using namespace std;

RigidBodyInstance::RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta, double density) : c(c), theta(theta), density(density), rbtemplate_(rbtemplate)
{
    cvel.setZero();
    w.setZero();
}

void RigidBodyInstance::render()
{
    Matrix3d rot = VectorMath::rotationMatrix(theta);

    glShadeModel(GL_SMOOTH);
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
    glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
    glEnable ( GL_COLOR_MATERIAL );
    Vector3d color(0.6, 0.9, 0.9);
    glColor4d(color[0], color[1], color[2], 1.0);

    glPushMatrix();
    {
        GLdouble xform[16];
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
                xform[4*j+i] = rot.coeff(i,j);
            xform[4*i+3] = 0;
            xform[12+i] = c[i];
        }
        xform[15] = 1.0;
        glMultMatrixd(xform);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);

        glVertexPointer(3, GL_DOUBLE, 0, rbtemplate_.getMesh().getVertexPointer());
        glNormalPointer(GL_DOUBLE, 0, rbtemplate_.getMesh().getVertexNormalsPointer());

        glDrawElements(GL_TRIANGLES, rbtemplate_.getMesh().getNumFaces()*3, GL_UNSIGNED_INT, rbtemplate_.getMesh().getFacePointer());

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    glPopMatrix();

    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    {
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(c[0], c[1], c[2]);
        Vector3d tip = c + 1.1 * rot*rbtemplate_.getPrincipleAxis(0);
        glVertex3f(tip[0], tip[1], tip[2]);

        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(c[0], c[1], c[2]);
        tip = c + 1.1 * rot*rbtemplate_.getPrincipleAxis(2);
        glVertex3f(tip[0], tip[1], tip[2]);
    }
    glEnd();
}
