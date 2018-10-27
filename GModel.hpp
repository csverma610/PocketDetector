#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <array>
#include <algorithm>
#include <boost/array.hpp>
#include <deque>
#include <set>

#include "GMesh.hpp"

using namespace std;
typedef boost::array<int,3> Array3I;
typedef boost::array<int,2> Array2I;
typedef boost::array<double,2> Point2D;
typedef boost::array<double,3> Point3D;

#include <boost/bind.hpp>

#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <Geom_Plane.hxx>
#include <TopoDS.hxx>
#include <TopOpeBRepTool_ShapeExplorer.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepMesh_FastDiscret.hxx>

#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>

#include <Standard_Version.hxx>
#if (OCC_VERSION_MAJOR < 7)
#  include <Handle_Standard_Transient.hxx>
#else
#  include <Standard_Transient.hxx>
#endif

#include <TColStd_SequenceOfTransient.hxx>
#include <TColStd_HSequenceOfTransient.hxx>
#include <TColgp_HArray1OfPnt.hxx>

#include <gp_Pnt.hxx>
#include <gp_Lin.hxx>
#include <gp_Vec.hxx>

#include <BRepTools.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_HCurve.hxx>
#include <BRepAdaptor_HCompCurve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepAlgo_Section.hxx>

#include <Geom_Plane.hxx>
#include <Geom_BoundedCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomLProp_SLProps.hxx>

class GNode;
class GEdge;
class GFace;

///////////////////////////////////////////////////////////////////////////////
class GNode
{
public:
    static GNode* newObject();

    void setID(int i) {
        id = i;
    }
    int  getID() const {
        return id;
    }
    void addOCCNode(TopoDS_Vertex v) {
        occNodes.push_back(v);
    }
    void addRelation(TopoDS_Edge  e);
    void addRelation(TopoDS_Face  e);
    void setCoords( const Point3D &p) {
        xyzCoords = p;
    }
    const Point3D &getCoords() const {
        return xyzCoords;
    }

private:
    int id;
    Point3D xyzCoords;
    vector<TopoDS_Vertex> occNodes;
    vector<TopoDS_Edge>   relations1;
    vector<TopoDS_Face>   relations2;
};

typedef vector<GNode*> GNodeSequence;
typedef set<GNode*>    GNodeSet;

class JDualNode
{
public:
    JDualNode() {}
    JDualNode( GFace *f) {
        primalFace = f;
    }

    GFace* getFace() {
        return primalFace;
    }
    int getID() const;
private:
    GFace *primalFace = nullptr;
};


///////////////////////////////////////////////////////////////////////////////

inline GNode* GNode::newObject() {
    GNode *newobj = new GNode;
    return newobj;
}

class GEdge
{
public:
    static GEdge* newObject(GNode *v0, GNode *v1);

    GEdge() {};

    GEdge(GNode *v0, GNode *v1) {
        nodes.resize(2);
        nodes[0] = v0;
        nodes[1] = v1;
    }

    void setID(int i) {
        id = i;
    }
    int  getID() const {
        return id;
    }

    void addOCCEdge(TopoDS_Edge &e) {
        occEdges.push_back(e);
    }

    void addRelation(GFace *f) {
        if( find(relations2.begin(), relations2.end(), f) == relations2.end() ) 
        relations2.push_back(f);
    }

    vector<GFace*> getRelations2() const {
        return relations2;
    }

    bool isSame( const GEdge* rhs) const {
        auto rhsnodes = rhs->getNodes();
        if( (nodes[0] == rhsnodes[0]) && (nodes[1] == rhsnodes[1]) ) return 1;
        if( (nodes[0] == rhsnodes[1]) && (nodes[1] == rhsnodes[0]) ) return 1;
        return 0;
    }

    bool isSame( const GNode* v0, const GNode *v1) const {
        if( (nodes[0] == v0) && (nodes[1] == v1) ) return 1;
        if( (nodes[0] == v1) && (nodes[1] == v0) ) return 1;
        return 0;
    }

    int  getOrientation( const GNode *v0, const GNode *v1) const
    {
        if( (nodes[0] == v0) && (nodes[1] == v1) ) return  1;
        if( (nodes[0] == v1) && (nodes[1] == v0) ) return -1;
        return 0;
    }

    const vector<GNode*> &getNodes() const {
        return nodes;
    }

    bool isClosed() const {
        if( nodes.empty() ) return 0;
        if( nodes.size() == 1) return 1;
        if( nodes[0] == nodes[1] ) return 1;
        return 0;
    }

    GNode* getNodeAt(int i) const {
      assert( i < 2);
      return nodes[i];
    }

    bool hasNode( const GNode *v) const 
    {
      if( nodes[0] == v) return 1;
      if( nodes[1] == v) return 1;
      return 0;
    }

    bool isAdjacent( const GEdge *nextedge) const
    {
      if( nextedge->hasNode( nodes[0] )) return 1;
      if( nextedge->hasNode( nodes[1] )) return 1;
      return 0;
    }

    bool hasOCCEdge( TopoDS_Edge e) const
    {
        for( auto oedge : occEdges)
            if( e == oedge ) return 1;
        return 0;
    }

private:
    int id;
    GNodeSequence   nodes;
    vector<GFace*>  relations2;
    vector<TopoDS_Edge> occEdges;
};

typedef vector<GEdge*> GEdgeSequence;
typedef set<GEdge*>    GEdgeSet;

inline GEdge* GEdge::newObject(GNode *v0, GNode *v1) {
    GEdge *newobj = new GEdge(v0,v1);
    return newobj;
}

class JDualEdge
{
public:
    JDualEdge() {}

    JDualEdge( JDualNode *v0, JDualNode *v1)
    {
        nodes[0] = v0;
        nodes[1] = v1;
    }
    JDualNode* getNodeAt(int i) {
        if( i < 2) return nodes[i];
        return nullptr;
    }


private:
    JDualNode* nodes[2];
};

///////////////////////////////////////////////////////////////////////////////

class GFace
{
public:
    void setID( int i ) {
        id = i;
    }

    int  getID() const {
        return id;
    }

    void setOCCFace(TopoDS_Face &f) {
        occFace = f;
    }

    void addEdge(GEdge *e) {
         if( find( boundEdges.begin(), boundEdges.end(), e) == boundEdges.end() )  
            boundEdges.push_back(e);
    }

    GEdgeSequence getEdges() const { 
        return boundEdges;
    }

    GNodeSequence  getNodes();

    void setDualNode( JDualNode *vtx) {
        dualNode = vtx;
    }

    JDualNode* getDualNode() const {
        return dualNode;
    }

    bool isPlanar() const;
    bool isCircular() const;
    GFaceMesh *getSurfMesh();

    void addNeighbor(GFace *f){
         if( find( adjFaces.begin(), adjFaces.end(), f) == adjFaces.end() )  
            adjFaces.push_back(f);
    };

    vector<GFace*> getAdjacentFaces() const {
           return adjFaces;
    }

    vector<GEdgeSequence> getWires() const;

private:
    int id;
    TopoDS_Face  occFace;
    vector<GFace*> adjFaces;
    GEdgeSequence boundEdges;
 
    JDualNode *dualNode = nullptr;
    vector<GEdgeSequence> wires;
    GEdgeSequence getNewLoop( list<GEdge*> &l) const;
};

typedef vector<GFace*> GFaceSequence;
typedef set<GFace*>    GFaceSet;

///////////////////////////////////////////////////////////////////////////////
class JGraph
{
public:
    void addNode( JDualNode *v) {
        nodes.push_back(v);
    }

    void addEdge( JDualEdge *e) {
        edges.push_back(e);
    }

    int getSize(int e ) const {
        if( e == 0) return nodes.size();
        if( e == 1) return edges.size();
        return 0;
    }

    void saveAs( const string &s);

private:
    vector<JDualNode*> nodes;
    vector<JDualEdge*> edges;
};

///////////////////////////////////////////////////////////////////////////////

inline int JDualNode :: getID() const
{
    return primalFace->getID();
}
///////////////////////////////////////////////////////////////////////////////
class GModel
{
public:
    int readSTEPModel (const string &s) {
        readOCCModel(s);
    }
    void setOCCShape(TopoDS_Shape shape);

    void getOCCInfo();

    JGraph* getGraph() {
        return graph;
    }

    size_t getSize(int e) const {
        if(e == 0) {
            return gNodes.size();
        }
        if(e == 1) {
            return gEdges.size();
        }
        if(e == 2) {
            return gFaces.size();
        }
    }
    GNode *getNodeAt(int i) const {return gNodes[i]; }
    GEdge *getEdgeAt(int i) const { return gEdges[i]; }
    GFace *getFaceAt(int i) const { return gFaces[i]; }

    GMesh* getSurfMesh();
private:
    struct ONode
    {
        TopoDS_Vertex occNode;
        Point3D       xyz;
        bool          active;
        double getDistance( const ONode &rhs) const
        {
            double dx = xyz[0] - rhs.xyz[0];
            double dy = xyz[1] - rhs.xyz[1];
            double dz = xyz[2] - rhs.xyz[2];
            return sqrt(dx*dx + dy*dy + dz*dz);
        }
    };

    TopoDS_Shape occShape;
    std::vector<TopoDS_Compound>  occCompounds;
    std::vector<TopoDS_CompSolid> occCompsolids;
    std::vector<TopoDS_Solid>     occSolids;
    std::vector<TopoDS_Shell>     occShells;
    std::vector<TopoDS_Wire>      occWires;
    std::vector<TopoDS_Face>      occFaces;
    std::vector<TopoDS_Edge>      occEdges;
    std::vector<TopoDS_Vertex>    occNodes;
    int readOCCModel( const string &s);

    GNodeSequence  gNodes;
    GEdgeSequence  gEdges;
    GFaceSequence  gFaces;

    JGraph *graph = nullptr;

    void extract_geometrical_shapes();
    void extract_compound_shapes();

    GNode*  getClosestNode( const Point3D &p) const;
    double  getDistance( const GNode *a, const GNode *b) const;
    double  getDistance( const Point3D  &p0, const Point3D &p1) const;
    GEdge*  getEdgeOf( const GNode *v0, const GNode *v1) const;
    GEdge*  getEdgeOf( TopoDS_Edge e) const;
    GFaceSequence  getNeighbors(GFace *f) const;

    void buildTopology();
    void buildNodes();
    void buildEdges();
    void buildFaces();
    gp_Pnt point(const Point3D &p);
    Point3D point(const gp_Pnt &p);

    GFaceMesh *mesh;
};
