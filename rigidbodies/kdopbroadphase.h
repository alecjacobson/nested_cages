#ifndef KDOPBROADPHASE_H
#define KDOPBROADPHASE_H
#include "rigidbodyinstance.h"
#include "rigidbodytemplate.h"
#include <set>
#include <vector>

#define K 13

struct KDOPNode
{
  virtual ~KDOPNode() {}
  virtual bool isLeaf() = 0;

  double mins[K];
  double maxs[K];
};

class NodeComparator
{
public:
    NodeComparator(int axis) : axis(axis) {}

    bool operator()(KDOPNode *left, KDOPNode *right)
    {
        return left->mins[axis] < right->mins[axis];
    }

private:
    int axis;
};

struct KDOPInteriorNode : public KDOPNode
{
    virtual ~KDOPInteriorNode()
    {
        delete left;
        delete right;
    }

    virtual bool isLeaf() {return false;}

    int splitaxis;
    KDOPNode *left, *right;
};

struct KDOPLeafNode : public KDOPNode
{
    virtual bool isLeaf() {return true;}

    int body;
    int face;
};

struct VertexFaceStencil
{
    int body;
    int vert;
    int face;

    bool operator<(const VertexFaceStencil &other) const
    {
        if(body == other.body)
        {
            if(vert == other.vert)
            {
                return face < other.face;
            }
            return vert < other.vert;
        }
        return body < other.body;
    }
};

class KDOPBroadPhase
{
public:
  KDOPBroadPhase();

  void findCollisionCandidates(const RigidBodyInstance &body1, const RigidBodyInstance &body2, const Eigen::Vector3d &newc1, const Eigen::Vector3d &newc2, const Eigen::Vector3d &newtheta1, const Eigen::Vector3d &newtheta2, std::set<VertexFaceStencil> &vfs, bool useCage);
 private:
  KDOPNode *buildKDOPTree(int bodyid, const RigidBodyInstance &body, const Eigen::Vector3d &newc, const Eigen::Vector3d &newtheta, bool useCage);
  KDOPNode *buildKDOPInterior(std::vector<KDOPNode *> &children);
  void intersect(KDOPNode *left, KDOPNode *right, const RigidBodyInstance &body1, const RigidBodyInstance &body2, std::set<VertexFaceStencil> &vfs, bool useCage);

  std::vector<Eigen::Vector3d> DOPaxis;
};

#endif
