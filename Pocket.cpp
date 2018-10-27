#include "UfabViewer.hpp"

/////////////////////////////////////////////////////////////////////////////
void Pocket :: preProcess()
{
    for( auto keyVal : mesh->groupFaces) {
        int gid = keyVal.first;
        double maxVal = 0.0;
        Vec3D A = mesh->groupFaces[gid].normal;
        for( auto faceid: keyVal.second.faces) {
            Vec3D B = mesh->faceNormal[faceid];
            double angle = getVecAngle(A, B);
            maxVal = max(angle, maxVal);
        }
        mesh->groupFaces[gid].maxNormalAngleDev = 180.0*maxVal/M_PI;
    }
}
/////////////////////////////////////////////////////////////////////////////
bool Pocket :: isPocket( int face )
{
    baseFaceID = face;
    bool rulePassed = 0;

    sideFaces.clear();
    topFaces.clear();

    rulePassed = rule1();
    if( !rulePassed ) {
        cout << "Verdict: probably not a pocket" << endl;
        return 0;
    }
    rulePassed = rule2();
    if( !rulePassed ) {
        cout << "Verdict: probably not a pocket" << endl;
        return 0;
    }
    rulePassed = rule3();
    if( !rulePassed ) {
        cout << "Verdict: probably not a pocket" << endl;
        return 0;
    }
    cout << "Verdict: Strong candidate to be a Pocket" << endl;
    return 1;
}
////////////////////////////////////////////////////////////////////////////////
bool Pocket :: rule1()
{
    if( mesh->groupFaces[baseFaceID].maxNormalAngleDev > 1.0 ) {
        cout << "The seed face is not planar " << endl;
        return 0;
    }

    Vec3D vn = mesh->groupFaces[baseFaceID].normal;
    if( vn[1] < 0.999) {
        cout << "Normal" << vn[0] << " " << vn[1] << " " << vn[2] << endl;
        cout << "The base face is not oriented in the y direction" << endl;
        return 0;
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////
bool Pocket :: rule2()
{
    // The center of all sides must be above zero in Y direction ...
    set<int> neighfaces = mesh->groupFaces[baseFaceID].geomFaces;
    neighfaces.erase(baseFaceID);
    for( auto id: neighfaces) {
        double y = mesh->groupFaces[id].center[1];
        if( y < 0.0) {
            cout << "One of the side face is below zero " << endl;
            return 0;
        }
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////

bool Pocket :: rule3()
{
    set<int> neighfaces = mesh->groupFaces[baseFaceID].geomFaces;
    neighfaces.erase(baseFaceID);

    Vec3D A = mesh->groupFaces[baseFaceID].normal;
    for( auto id: neighfaces) {
        Vec3D B = mesh->groupFaces[id].normal;
        double angle = 180.0*getVecAngle(A,B)/M_PI;
        if( fabs(angle-90) > 2.0) {
            cout << "One of the side face is not normal" << endl;
            cout << "Angle: " << angle << endl;
            return 0;
        }
    }
    return 1;
}
