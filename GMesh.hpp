
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/algorithm.hpp>
#include "Quaternion.hpp"

using namespace std;

typedef boost::array<double,3> Point3D;
typedef boost::array<int,3>    Array3I;
typedef boost::array<double,3> Vec3D;
typedef boost::array<float,4>  Color;
typedef vector<int>            ISequence;

template <class T>
inline double getDotProduct( const boost::array<T,3> &A, const boost::array<T,3> &B)
{
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

template <class T>
inline double getMagnitude( const boost::array<T,3> &A)
{
    return sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2] );
}

template <class T>
inline double getVecAngle( const boost::array<T,3> &A, const boost::array<T,3> &B)
{
    double AB = getDotProduct(A,B);
    double Am = getMagnitude(A);
    double Bm = getMagnitude(B);

    if( Am < 1.0E-15 || Bm < 1.0E-15) return 0.0;

    double x = AB/(Am*Bm);

    if( x > 1.0) x = 1.0;
    if( x < -1.0) x = -1.0;

    return acos(x);
}

template <class T>
inline boost::array<T,3> getVector( const boost::array<T,3> &head, const boost::array<T,3> &tail)
{
    boost::array<T,3> vec;
    vec[0] = head[0] - tail[0];
    vec[1] = head[1] - tail[1];
    vec[2] = head[2] - tail[2];
    return vec;
}

template <class T>
inline boost::array<T,3> getCrossProduct( const boost::array<T,3> &a, const boost::array<T,3> &b)
{
    boost::array<T,3> c;
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
}

template <class T>
inline void getUnitVector( boost::array<T,3> &vec)
{
    double dl  = getMagnitude(vec);
    if( dl > 1.0E-08) {
    vec[0] = vec[0]/dl;
    vec[1] = vec[1]/dl;
    vec[2] = vec[2]/dl;
    }
}

template <class T>
inline boost::array<T,3> getCentroid( const boost::array<T,3>  &pa, 
                                      const boost::array<T,3>  &pb, 
                                      const boost::array<T,3>  &pc)
{
    boost::array<T,3> pd;
    pd[0] = (pa[0] + pb[0] + pc[0] )/3.0;
    pd[1] = (pa[1] + pb[1] + pc[1] )/3.0;
    pd[2] = (pa[2] + pb[2] + pc[2] )/3.0;
    return pd;
}

template <class T>
inline boost::array<T,3> getVectorHead( const boost::array<T,3> &ptail, 
                                        const boost::array<T,3> &vec, double len)
{
    boost::array<T,3> phead;
    phead[0] = ptail[0] + len*vec[0];
    phead[1] = ptail[1] + len*vec[1];
    phead[2] = ptail[2] + len*vec[2];
    return phead;
}

struct GFaceMesh
{
   size_t getSize(int e) const {
          if( e == 0) return nodeCoords.size();
          if( e == 2) return triangles.size();
          return 0;
   }

   bool      displayNodes = 1;
   bool      displayFaces = 1;
   bool      filledFaces  = 1;
   int       groupID;
   Vec3D     normal;
   Point3D   center;
   Color     orgColor, faceColor;

   vector<Point3D>  nodeCoords;
   vector<Array3I>  triangles;
   vector<Point3D>  faceNormal;
   vector<Point3D>  faceCentroid;
   void updateGeometry();
   void setFacesNormal();
   void setGroupCenters();
};

struct GMesh
{
    vector<GFaceMesh*>   facemeshes;
    Point3D  getModelCenter();
    void readMesh( const string &f);
    void updateGeometry();
    void reOrient();
    void setCenterAt(int groupid);
    void saveAs(const string &s);
    void alignAlong( const Vec3D &src, const Vec3D &dst);
    void buildGeomTopology();
};
