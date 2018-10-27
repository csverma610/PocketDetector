#pragma once

#include <QGLViewer/qglviewer.h>
#include <QtWidgets>
#include <QGLViewer/qglviewer.h>
#include <qapplication.h>

#include "GMesh.hpp"
#include "CADPocketDetector.hpp"

//////////////////////////////////////////////////////////////////////

struct GMeshRender : public GMesh
{
    boost::dynamic_bitset<> nodeBit;
    boost::dynamic_bitset<> faceBit;
    vector<Color>  orgFaceColor, faceColor;
    void buildGeomTopology();
};

//////////////////////////////////////////////////////////////////////

class UfabViewer : public QGLViewer 
{
public:
    UfabViewer(QWidget *parent = nullptr);

    void readSTEPModel( const string &fname);
    void setMesh( GMeshRender *m) {
        mesh = m;
    }

protected:
    virtual void draw();
    virtual void drawWithNames();
    virtual void init();
    virtual void keyPressEvent( QKeyEvent *e);
    virtual void mousePressEvent( QMouseEvent *e);
    virtual void mouseReleaseEvent( QMouseEvent *e);

private:
    GMesh *mesh = nullptr;
    CADPocketDetector  cpktDetector;

    int faceSelected = -1;
    bool renderWires   = 1;
    bool renderMeshNormals = 0;
    bool renderMeshEdges    = 0;
    float faceNormalsLength = 0.1;

    void   initMesh();
    void   assignNodeColor();
    void   assignEdgeColor();
    void   assignFaceColor();

    void   drawTargetAxes( const Point3D &pc, const Vec3D &vec);
    void   drawAxes();
    void   drawCenters();
    void   drawNodes();
    void   drawEdges();
    void   drawFaceEdges();
    void   drawFaces();
    void   drawFacesNormal();

    void  drawNodesWithName();
    void  drawEdgesWithName();
    void  drawFacesWithName();


    void  setCenter( const Point3D &p);
    void  detectAll();
    void  selectFace( int id);
};
//////////////////////////////////////////////////////////////////////
