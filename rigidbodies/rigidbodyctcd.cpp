#include "rigidbodyctcd.h"
#include "collisions/CTCD.h"
#include "rigidbodyinstance.h"
#include "rigidbodytemplate.h"
#include "mesh.h"
#include "vectormath.h"
#include "kdopbroadphase.h"
#include <iostream>

using namespace Eigen;
using namespace std;

RigidBodyCTCD::RigidBodyCTCD()
{
}

bool RigidBodyCTCD::detectCollisions(const std::vector<RigidBodyInstance *> &bodies, std::vector<Plane> &planes, std::vector<Eigen::Vector3d> &newc, std::vector<Eigen::Vector3d> &newtheta, std::vector<ContactInfo> &contacts, SimParameters::UseCage  useCage)
{
    int numbodies = bodies.size();
    vector<pair<int, ContactInfo> > candidates;

    for(int i=0; i<numbodies; i++)
    {
        RigidBodyInstance &body = *bodies[i];
        int nverts = body.getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getNumVerts();
        for(int j=0; j<nverts; j++)
        {
            Vector3d oldpos = body.c + VectorMath::rotationMatrix(body.theta)*body.getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(j);
            Vector3d newpos = newc[i] + VectorMath::rotationMatrix(newtheta[i])*body.getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(j);

            for(int p=0; p<planes.size(); p++)
            {

                double num = (planes[p].pos - oldpos).dot(planes[p].normal);
                double denom = (newpos-oldpos).dot(planes[p].normal);
                double t= num/denom;
                if(t >= 0 && t < 1)
                {
                    ContactInfo ci;
                    ci.body1 = i;
                    ci.n = planes[p].normal;
                    ci.body2 = -1;
                    ci.pt1 = body.getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(j);
                    candidates.push_back(pair<int, ContactInfo>(t, ci));
                }
            }
        }

        for(int k=i+1; k<numbodies; k++)
        {
            if(!boundingCubesIntersect(*bodies[i], *bodies[k], newc[i], newc[k], newtheta[i], newtheta[k]))
                continue;

            RigidBodyInstance &target = *bodies[k];
            set<VertexFaceStencil> vfs;
            KDOPBroadPhase broad;
            broad.findCollisionCandidates(body, target, newc[i], newc[k], newtheta[i], newtheta[k], vfs, useCage == SimParameters::C_ALWAYS);

            //std::cout << "found " << vfs.size() << " candidates" << std::endl;
            for(set<VertexFaceStencil>::iterator it = vfs.begin(); it != vfs.end(); ++it)
            {
                int body1 = (it->body == 0 ? i : k);
                int body2 = (it->body == 0 ? k : i);

                Vector3d starts[3];
                Vector3d ends[3];
                Vector3i verts = bodies[body2]->getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getFace(it->face);

                Vector3d vertp1 = bodies[body1]->c + VectorMath::rotationMatrix(bodies[body1]->theta)*bodies[body1]->getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(it->vert);
                Vector3d vertp2 = newc[body1] + VectorMath::rotationMatrix(newtheta[body1])*bodies[body1]->getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(it->vert);

                for(int n=0; n<3; n++)
                {
                    starts[n] = bodies[body2]->c + VectorMath::rotationMatrix(bodies[body2]->theta)*bodies[body2]->getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(verts[n]);
                    ends[n] = bodies[body2]->c + VectorMath::rotationMatrix(bodies[body2]->theta)*bodies[body2]->getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(verts[n]);
                }
                double t;
                if(CTCD::vertexFaceCTCD(vertp1, starts[0], starts[1], starts[2], vertp2, ends[0], ends[1], ends[2], t))
                {
                    Vector3d mids[3];
                    for(int n=0; n<3; n++)
                        mids[n] = starts[n]*(1-t) + ends[n]*t;
                    ContactInfo ci;
                    ci.body1 = body1;
                    ci.body2 = body2;
                    ci.n = (mids[1]-mids[0]).cross(mids[2]-mids[0]);
                    ci.n /= ci.n.norm();
                    ci.pt1 = bodies[body1]->getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(it->vert);
                    ci.pt2.setZero();
                    for(int n=0; n<3; n++)
                        ci.pt2 += bodies[body2]->getTemplate().selectMesh(useCage == SimParameters::C_ALWAYS).getVert(verts[n]);
                    ci.pt2 *= 1.0/3.0;
                    candidates.push_back(pair<int, ContactInfo>(t, ci));
                }
            }
        }
   }

    if(candidates.empty())
        return false;

    double earliest = std::numeric_limits<double>::infinity();
    ContactInfo ci = candidates[0].second;

    for(int i=0; i<candidates.size(); i++)
    {
        if(candidates[i].first < earliest)
        {
            earliest = candidates[i].first;
            ci = candidates[i].second;
        }
    }

    contacts.push_back(ci);

    return true;
}

bool RigidBodyCTCD::boundingCubesIntersect(const RigidBodyInstance &body1, const RigidBodyInstance &body2, const Vector3d &newc1, const Vector3d &newc2, const Vector3d &newtheta1, const Vector3d &newtheta2)
{
    double mins1[3];
    double maxs1[3];
    double mins2[3];
    double maxs2[3];
    for(int i=0; i<3; i++)
    {
        mins1[i] = std::numeric_limits<double>::infinity();
        mins2[i] = std::numeric_limits<double>::infinity();
        maxs1[i] = -std::numeric_limits<double>::infinity();
        maxs2[i] = -std::numeric_limits<double>::infinity();
    }

    for(int i=0; i<3; i++)
    {
        mins1[i] = std::min(body1.c[i] - 1, newc1[i] - 1);
        mins2[i] = std::min(body2.c[i] - 1, newc2[i] - 1);
        maxs1[i] = std::max(body1.c[i] + 1, newc1[i] + 1);
        maxs2[i] = std::max(body2.c[i] + 1, newc2[i] + 1);
    }

    for(int i=0; i<3; i++)
    {
        if(maxs1[i] < mins2[i] || maxs2[i] < mins1[i])
            return false;
    }
    return true;
}
