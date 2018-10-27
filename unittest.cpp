#include "GModel.hpp"

struct Point3D
{
    double x, y, z;
    bool operator < ( const Point3D &rhs) const
    {
        return x < rhs.x;
        return y < rhs.y;
        return z < rhs.z;
        return 0;
    }

};

bool pntCompare( Point3D &ap, Point3D &bp)
{
    double eps = 1.0E-15;

    double dz = fabs(ap.z - bp.z);
    double dy = fabs(ap.y - bp.y);
    double dx = fabs(ap.x - bp.x);
    return (dx + dy + dz > eps);
}
/////////////////////////////////////////////////////////////////////////

void unitTest()
{
    vector<Point3D> v(100);
    for( int i = 0; i < 100; i++) {
        v[i].x = i;
        v[i].y = i+1;
        v[i].z = i+2;
    }
    std::sort( v.begin(), v.end(), pntCompare);
    /*
        vector<Point3D> allpoints;

        size_t numNodes = 1000;

        allpoints.resize(numNodes);
        for( int i = 0; i < numNodes; i++) {
    //      const gp_Pnt vpos = BRep_Tool::Pnt( occNodes[i] );
            allpoints[i][0] = drand48();
            allpoints[i][1] = drand48();
            allpoints[i][2] = drand48();
        }

        cout << "HELLO Unit Test" << endl;
        std::sort( allpoints.begin(), allpoints.end(), isGreater);
        cout << "HELLO Unit Test" << endl;
    */
}

int main(int argc, char **argv)
{
    unitTest();
    return 0;


}
