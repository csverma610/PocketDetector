#pragma once

#include "GModel.hpp"
#include "GMesh.hpp"

/////////////////////////////////////////////////////////////////////////////

struct CADPocket
{
    double dim[3];
    vector<int> faces;
};

/////////////////////////////////////////////////////////////////////////////

class CADPocketDetector
{
public:
    int  readSTEPModel( const string &s);
    void setOCCShape(TopoDS_Shape s);

    int detectAll();

    vector<CADPocket> getPockets() { return pockets; }

    bool isPocket(int faceid);
    int  reOrient(int faceid);
    void saveMesh( const string &s) {
        if(mesh) mesh->saveAs(s);
    }

    GMesh *getSurfMesh() const {
        return mesh;
    }

    GModel *getGeomModel() const {
        return gModel;
    }

    void genVoxels();
private:
    bool verbose = 1;
    GModel *gModel = nullptr;
    GMesh  *mesh   = nullptr;
    double  getMinYCoord(const GFaceMesh *m);
    double  getMaxYCoord(const GFaceMesh *m);
    bool rule1(int faceid);
    bool rule2(int faceid);
    bool rule3(int faceid);
    vector<CADPocket> pockets;
};
