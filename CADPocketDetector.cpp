#include "CADPocketDetector.hpp"
//#include "Visibility_HelperClasses.h"

//////////////////////////////////////////////////////////////////////

int CADPocketDetector :: readSTEPModel( const string &filename)
{
    gModel = new GModel;
    gModel->readSTEPModel(filename);
    mesh = gModel->getSurfMesh();
    return 0;
}
/////////////////////////////////////////////////////////////////////////////

void CADPocketDetector :: setOCCShape( TopoDS_Shape shape)
{
    gModel = new GModel;
    gModel->setOCCShape(shape);
    mesh = gModel->getSurfMesh();
}
/////////////////////////////////////////////////////////////////////////////
int CADPocketDetector :: detectAll()
{
    if( mesh == nullptr) return 0;

    pockets.clear();

    set<int> faceSet;
    for( int i = 0; i < gModel->getSize(2); i++)
        faceSet.insert(i);

    int numPockets = 0;
    int numGeomFaces = gModel->getSize(2);
    while(!faceSet.empty()) {
        int id = *faceSet.begin();
        faceSet.erase(id);
        reOrient(id);
        if( isPocket(id) ) {
            numPockets++;
            auto neighfaces = gModel->getFaceAt(id)->getAdjacentFaces();
            for( auto f: neighfaces)
                faceSet.erase(f->getID() );
        }
    }

    cout << "#Pockets detected " << numPockets << endl;

    return numPockets;
}

/////////////////////////////////////////////////////////////////////////////

int  CADPocketDetector :: reOrient( int faceID )
{
    mesh->setCenterAt(faceID);
    mesh->updateGeometry();

    GFaceMesh *submesh = mesh->facemeshes[faceID];
    if( submesh == nullptr) return 1;

    // The center must be located at the center of the planar face ....
    Point3D pc = submesh->center;
    if( fabs(pc[0]) > 1.0E-06) {
        cout << "Warning: model is not correctly centered at the origin" << endl;
        return 2;
    }

    if( fabs(pc[1]) > 1.0E-06) {
        cout << "Warning: model is not correctly centered at the origin" << endl;
        return 2;
    }

    if( fabs(pc[2]) > 1.0E-06) {
        cout << "Warning: model is not correctly centered at the origin" << endl;
        return 2;
    }

    // Rotate the object so that normal of the facegroup is in the "Y" direction ...

    Vec3D dstVec;
    Vec3D srcVec = submesh->normal;

    dstVec[0] = 0.0;
    dstVec[1] = 1.0;
    dstVec[2] = 0.0;
    mesh->alignAlong(srcVec, dstVec);
    mesh->updateGeometry();

    // Make sure again that after the transformations, the center is still at the origin ...
    pc = submesh->center;
    if( fabs(pc[0]) > 1.0E-06) {
        cout << "Warning: model is not correctly centered at the origin" << endl;
        return 2;
    }

    if( fabs(pc[1]) > 1.0E-06) {
        cout << "Warning: model is not correctly centered at the origin" << endl;
        return 2;
    }

    if( fabs(pc[2]) > 1.0E-06) {
        cout << "Warning: model is not correctly centered at the origin" << endl;
        return 2;
    }

    // Check the normal, it should be in the +y Direction..
    Vec3D vn = submesh->normal;
    if( vn[1] < 0.0) {
        for( int j = 0; j < mesh->facemeshes.size(); j++) {
            GFaceMesh *m = mesh->facemeshes[j];
            for( size_t i = 0; i < m->getSize(0); i++) {
                Point3D p3d = m->nodeCoords[i];
                p3d[1] *= -1;
                p3d[2] *= -1;
                m->nodeCoords[i] = p3d;
            }
            m->updateGeometry();
        }
    }

    vn = submesh->normal;
    if( vn[1] < 0.999 ) {
        cout << "Warning: Face not oriented in Y direction " << endl;
        cout << "Normal : " << vn[0] << " " << vn[1] << " " << vn[2] << endl;
        return 2;
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool CADPocketDetector :: isPocket( int faceID )
{
    bool rulePassed = 0;
    rulePassed = rule1(faceID);
    if( !rulePassed ) {
        if(verbose)  cout << "Verdict: rule1: probably not a pocket" << endl;
        return 0;
    }

    int stat = reOrient(faceID);
    if( stat ) return 0;

    rulePassed = rule2( faceID);
    if( !rulePassed ) {
        if(verbose) cout << "Verdict: rule2: probably not a pocket" << endl;
        return 0;
    }

    rulePassed = rule3( faceID );
    if( !rulePassed ) {
        if(verbose) cout << "Verdict: fule3: probably not a pocket" << endl;
        return 0;
    }

    GFace *f = gModel->getFaceAt(faceID);
    if( f->getWires().size() > 1) return 0;

    auto neighfaces = gModel->getFaceAt(faceID)->getAdjacentFaces();
    double depth = 0.0;
    for( auto gface : neighfaces) {
        int id = gface->getID();
        auto submesh = mesh->facemeshes[id];
        double ymax = getMaxYCoord(submesh);
        depth = max(depth, ymax);
    }
    CADPocket newPocket;
    newPocket.faces.push_back(faceID);
    for( auto gface : neighfaces) {
        int id = gface->getID();
        newPocket.faces.push_back(id);
    }
    pockets.push_back(newPocket);


    vector<Point3D> pckPoints;
    vector<Array3I> pckTriangles;
    for( auto xyz : mesh->facemeshes[faceID]->nodeCoords) {
        xyz[1] = depth;
        pckPoints.push_back(xyz);
    }

    for( auto tri: mesh->facemeshes[faceID]->triangles)
        pckTriangles.push_back(tri);

    Array3I newtri;
    for( int i = 1; i < newPocket.faces.size(); i++) {
        for( auto xyz : mesh->facemeshes[i]->nodeCoords) {
            pckPoints.push_back(xyz);
            for( auto tri: mesh->facemeshes[i]->triangles) {
                newtri[0] = tri[2];
                newtri[1] = tri[1];
                newtri[2] = tri[0];
                pckTriangles.push_back(newtri);
            }
        }
    }

    cout << "Detect Pocket at  " << faceID << " Depth " << depth << endl;

    return 1;
}

////////////////////////////////////////////////////////////////////////////////

bool CADPocketDetector :: rule1(int faceID)
{
    GFace *f = gModel->getFaceAt(faceID);
    assert(f != nullptr);
    if( !(f->isPlanar()) ) return 0;
    return 1;
}

////////////////////////////////////////////////////////////////////////////////
double CADPocketDetector :: getMinYCoord( const GFaceMesh *submesh)
{
    double ymin = std::numeric_limits<double>::max();
    for( auto xyz : submesh->nodeCoords)
        ymin = std::min(ymin, xyz[1]) ;
    return ymin;
}

////////////////////////////////////////////////////////////////////////////////
double CADPocketDetector :: getMaxYCoord( const GFaceMesh *submesh)
{
    double ymax = -0.99*std::numeric_limits<double>::max();
    for( auto xyz : submesh->nodeCoords)
        ymax = std::max(ymax, xyz[1]) ;
    return ymax;
}
////////////////////////////////////////////////////////////////////////////////

bool CADPocketDetector :: rule2(int faceID)
{
    // The center of all sides must be above zero in Y direction ...
    auto neighfaces = gModel->getFaceAt(faceID)->getAdjacentFaces();
    auto thisface = gModel->getFaceAt(faceID);
    if( neighfaces.empty() ) {
        cout << "Warning: the face has no neighbors " << endl;
        return 0;
    }

    assert( !neighfaces.empty() );
    for( auto gface: neighfaces) {
        if( gface != thisface ) {
            int id = gface->getID();
            auto submesh = mesh->facemeshes[id];
            if( getMinYCoord(submesh) < -1.0E-10) return 0;
        }
    }
    return 1;
}

////////////////////////////////////////////////////////////////////////////////

bool CADPocketDetector :: rule3(int faceID)
{
    auto neighfaces = gModel->getFaceAt(faceID)->getAdjacentFaces();
    auto thisface = gModel->getFaceAt(faceID);

    auto A = mesh->facemeshes[faceID]->normal;
    for( auto gface: neighfaces) {
        if( gface != thisface ) {
            int id = gface->getID();
            auto B = mesh->facemeshes[id]->normal;
            if( getMagnitude(B) > 1.0E-06) {
                double angle = 180.0*getVecAngle(A,B)/M_PI;
                if( fabs(angle-90) > 2.0) {
                    if(verbose) {
                        cout << "One of the side face is not normal" << endl;
                        cout << "Angle detected " << angle << endl;
                        cout << "Face " << faceID << "Vec A" << A[0] << " " << A[1] << " " << A[2] << endl;
                        cout << "Face " << id << "Vec B" << B[0] << " " << B[1] << " " << B[2] << endl;
                    }
                    return 0;
                }
            }
        }
    }

    return 1;
}

////////////////////////////////////////////////////////////////////////////////
