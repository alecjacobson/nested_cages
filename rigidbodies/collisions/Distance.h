#ifndef DISTANCE_H
#define DISTANCE_H

#include <Eigen/Core>
#include <algorithm>

class Distance
{
 public:
  // Computes the vector between a point p and the closest point to p on the triangle (q0, q1, q2). Also returns the barycentric coordinates of this closest point on the triangle;
  // q0bary is the barycentric coordinate of q0, etc. (The distance from p to the triangle is the norm of this vector.)
  static Eigen::Vector3d vertexFaceDistance(const Eigen::Vector3d &p, 
					    const Eigen::Vector3d &q0, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, 
					    double &q0bary, double &q1bary, double &q2bary)
  {
  Eigen::Vector3d ab = q1-q0;
  Eigen::Vector3d ac = q2-q0;
  Eigen::Vector3d ap = p-q0;

  double d1 = ab.dot(ap);
  double d2 = ac.dot(ap);

  // corner and edge cases

  if(d1 <= 0 && d2 <= 0)
    {
      q0bary = 1.0;
      q1bary = 0.0;
      q2bary = 0.0;
      return q0-p;
    }

  Eigen::Vector3d bp = p-q1;
  double d3 = ab.dot(bp);
  double d4 = ac.dot(bp);
  if(d3 >= 0 && d4 <= d3)
    {
      q0bary = 0.0;
      q1bary = 1.0;
      q2bary = 0.0;
      return q1-p;
    }

  double vc = d1*d4 - d3*d2;
  if((vc <= 0) && (d1 >= 0) && (d3 <= 0))
    {
      double v = d1 / (d1-d3);
      q0bary = 1.0 - v;
      q1bary = v;
      q2bary = 0;
      return (q0 + v*ab)-p;
    }
  
  Eigen::Vector3d cp = p-q2;
  double d5 = ab.dot(cp);
  double d6 = ac.dot(cp);
  if(d6 >= 0 && d5 <= d6)
    {
      q0bary = 0;
      q1bary = 0;
      q2bary = 1.0;
      return q2-p;
    }

  double vb = d5*d2 - d1*d6;
  if((vb <= 0) && (d2 >= 0) && (d6 <= 0))
    {
      double w = d2/(d2-d6);
      q0bary = 1-w;
      q1bary = 0;
      q2bary = w;
      return (q0 + w*ac)-p;
    }

  double va = d3*d6 - d5*d4;
  if((va <= 0) && (d4-d3 >= 0) && (d5-d6 >= 0))
    {
      double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      q0bary = 0;
      q1bary = 1.0 - w;
      q2bary = w;
      
      return (q1 + w*(q2-q1))-p;
    }

  // face case
  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom;
  double w = vc * denom;
  double u = 1.0 - v - w;
  q0bary = u;
  q1bary = v;
  q2bary = w;
  return (u*q0 + v*q1 + w*q2)-p;
  }

  // Computes the shotest vector between a segment (p0, p1) and segment (q0, q1). Also returns the barycentric coordinates of the closest points on both segments; p0bary is the barycentric
  // coordinate of p0, etc. (The distance between the segments is the norm of this vector).
  static Eigen::Vector3d edgeEdgeDistance(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
					  const Eigen::Vector3d &q0, const Eigen::Vector3d &q1,
					  double &p0bary, double &p1bary,
					  double &q0bary, double &q1bary)
  {  
  Eigen::Vector3d d1 = p1-p0;
  Eigen::Vector3d d2 = q1-q0;
  Eigen::Vector3d r = p0-q0;
  double a = d1.squaredNorm();
  double e = d2.squaredNorm();
  double f = d2.dot(r);

  double s,t;

  double c = d1.dot(r);
  double b = d1.dot(d2);
  double denom = a*e-b*b;
  if(denom != 0.0) 
    {
      s = clamp( (b*f-c*e)/denom );
    }
  else 
    {
      //parallel edges and/or degenerate edges; values of s doesn't matter
      s = 0;
    }
  double tnom = b*s + f;
  if(tnom < 0 || e == 0)
    {
      t = 0;
      if(a == 0)
	s = 0;
      else
	s = clamp(-c/a);  
    }
  else if(tnom > e)
    {
      t = 1.0;
      if(a == 0)
	s = 0;
      else
	s = clamp( (b-c)/a );
    }
  else
    t = tnom/e;	    

  Eigen::Vector3d c1 = p0 + s*d1;
  Eigen::Vector3d c2 = q0 + t*d2;

  p0bary = 1.0-s;
  p1bary = s;
  q0bary = 1.0-t;
  q1bary = t;

  return c2-c1;
  }

 private:
  static double clamp(double u)
  {
    return std::min(1.0, std::max(u, 0.0));
  }
};

#endif
