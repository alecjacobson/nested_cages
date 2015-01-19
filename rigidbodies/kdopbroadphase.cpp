#include "kdopbroadphase.h"
#include <set>
#include <iostream>

#include "mesh.h"
#include "vectormath.h"

using namespace std;
using namespace Eigen;

KDOPBroadPhase::KDOPBroadPhase()
{
  DOPaxis.push_back(Vector3d(1.0, 0, 0));
  DOPaxis.push_back(Vector3d(0, 1.0, 0));
  DOPaxis.push_back(Vector3d(0, 0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, 1.0, 0));
  DOPaxis.push_back(Vector3d(1.0, -1.0, 0));
  DOPaxis.push_back(Vector3d(0, 1.0, 1.0));
  DOPaxis.push_back(Vector3d(0, 1.0, -1.0));
  DOPaxis.push_back(Vector3d(1.0, 0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, 0, -1.0));
  DOPaxis.push_back(Vector3d(1.0, 1.0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, 1.0, -1.0));
  DOPaxis.push_back(Vector3d(1.0, -1.0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, -1.0, -1.0));
  if(K>DOPaxis.size())
    exit(0);

  for(int i=0; i<(int)DOPaxis.size(); i++)
    {
      DOPaxis[i] /= DOPaxis[i].norm();
    }
}

void KDOPBroadPhase::findCollisionCandidates(const RigidBodyInstance &body1, const RigidBodyInstance &body2, const Vector3d &newc1, const Vector3d &newc2, const Vector3d &newtheta1, const Vector3d &newtheta2, set<VertexFaceStencil> &vfs)
{
    vfs.clear();
    KDOPNode *tree1 = buildKDOPTree(0, body1, newc1, newtheta1);
    KDOPNode *tree2 = buildKDOPTree(1, body2, newc2, newtheta2);
    intersect(tree1, tree2, body1, body2, vfs);
    delete tree1;
    delete tree2;
}

KDOPNode *KDOPBroadPhase::buildKDOPTree(int bodyid, const RigidBodyInstance &body, const Eigen::Vector3d &newc, const Eigen::Vector3d &newtheta)
{
  vector<KDOPNode *> leaves;
  int nfaces = body.getTemplate().getMesh().getNumFaces();
  for(int i=0; i<nfaces; i++)
  {
      KDOPLeafNode *node = new KDOPLeafNode;
      node->body = bodyid;
      node->face = i;
      int verts[3];
      for(int j=0; j<3; j++)
      {
        verts[j] = body.getTemplate().getMesh().getFace(i)[j];
      }
      for(int j=0; j<K; j++)
      {
        node->mins[j] = std::numeric_limits<double>::infinity();
        node->maxs[j] = -std::numeric_limits<double>::infinity();
      }
      for(int j=0; j<3; j++)
      {
        Vector3d oldpos = body.c + VectorMath::rotationMatrix(body.theta)*body.getTemplate().getMesh().getVert(verts[j]);
        for(int k=0; k<K; k++)
            {
              node->mins[k] = min(oldpos.dot(DOPaxis[k]), node->mins[k]);
              node->maxs[k] = max(oldpos.dot(DOPaxis[k]), node->maxs[k]);
            }
        Vector3d newpos = newc + VectorMath::rotationMatrix(newtheta)*body.getTemplate().getMesh().getVert(verts[j]);
        for(int k=0; k<K; k++)
            {
              node->mins[k] = min(newpos.dot(DOPaxis[k]), node->mins[k]);
              node->maxs[k] = max(newpos.dot(DOPaxis[k]), node->maxs[k]);
            }
      }
      leaves.push_back(node);
  }
return buildKDOPInterior(leaves);
}

KDOPNode *KDOPBroadPhase::buildKDOPInterior(vector<KDOPNode *> &children)
{
    int nchildren = children.size();
    assert(nchildren > 0);
    if(nchildren == 1)
        return children[0];

    KDOPInteriorNode *node = new KDOPInteriorNode;

    for(int i=0; i<K; i++)
    {
        node->mins[i] = std::numeric_limits<double>::infinity();
        node->maxs[i] = -std::numeric_limits<double>::infinity();
    }

    for(vector<KDOPNode *>::iterator it = children.begin(); it != children.end(); ++it)
    {
      for(int j=0; j<K; j++)
        {
          node->mins[j] = min((*it)->mins[j], node->mins[j]);
          node->maxs[j] = max((*it)->maxs[j], node->maxs[j]);
        }
    }
    double lengths[K];
    for(int i=0; i<K; i++)
    {
      lengths[i] = node->maxs[i] - node->mins[i];
    }
    int greatest = -1;
    double greatestlen = 0;
    for(int i=0; i<K; i++)
      {
        if(lengths[i] > greatestlen)
          {
        greatestlen = lengths[i];
        greatest = i;
          }
      }
    node->splitaxis = greatest;

    sort(children.begin(), children.end(), NodeComparator(node->splitaxis));

    vector<KDOPNode *> left;
    int child=0;
    for(; child<nchildren/2; child++)
        left.push_back(children[child]);
    node->left = buildKDOPInterior(left);
    vector<KDOPNode *> right;
    for(; child<nchildren; child++)
        right.push_back(children[child]);
    node->right = buildKDOPInterior(right);
    return node;
}

void KDOPBroadPhase::intersect(KDOPNode *left, KDOPNode *right, const RigidBodyInstance &body1, const RigidBodyInstance &body2, std::set<VertexFaceStencil> &vfs)
{
    for(int axis=0; axis<K; axis++)
    {
      if(left->maxs[axis] < right->mins[axis] ||
         right->maxs[axis] < left->mins[axis])
        return;
    }
    if(!left->isLeaf())
    {
        KDOPInteriorNode *ileft = (KDOPInteriorNode *)left;
        intersect(ileft->left, right, body1, body2, vfs);
        intersect(ileft->right, right, body1, body2, vfs);
    }
    else if(!right->isLeaf())
    {
        KDOPInteriorNode *iright = (KDOPInteriorNode *)right;
        intersect(left, iright->left, body1, body2, vfs);
        intersect(left, iright->right, body1, body2, vfs);
    }
    else
    {
        KDOPLeafNode *lleft = (KDOPLeafNode *)left;
        KDOPLeafNode *lright = (KDOPLeafNode *)right;

        // 6 vertex-face and 9 edge-edge
        for(int i=0; i<3; i++)
        {
            VertexFaceStencil lefts;
            lefts.body = lleft->body;
            lefts.vert = body1.getTemplate().getMesh().getFace(lleft->face)[i];
            lefts.face = lright->face;
            vfs.insert(lefts);

            VertexFaceStencil rights;
            rights.body = lright->body;
            rights.vert = body2.getTemplate().getMesh().getFace(lright->face)[i];
            rights.face = lleft->face;
            vfs.insert(rights);
        }
    }
}
