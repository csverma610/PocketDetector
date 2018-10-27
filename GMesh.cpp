#include "GMesh.hpp"
#include "Quaternion.hpp"
#include <boost/math/quaternion.hpp>


/////////////////////////////////////////////////////////////////////////////////
Point3D GMesh :: getModelCenter()
{
    Point3D pcenter, pc;

    pcenter[0] = 0.0;
    pcenter[1] = 0.0;
    pcenter[2] = 0.0;

    size_t nSize = 0;
    for( auto mesh : facemeshes) {
        int numNodes = mesh->getSize(0);
        for( int i  = 0; i < numNodes; i++) {
            pc  = mesh->nodeCoords[i];
            pcenter[0] += pc[0];
            pcenter[1] += pc[1];
            pcenter[2] += pc[2];
        }
        nSize += numNodes;
    }
    pcenter[0] /= (double)nSize;
    pcenter[1] /= (double)nSize;
    pcenter[2] /= (double)nSize;
    return pcenter;
}

/////////////////////////////////////////////////////////////////////////////////

void GFaceMesh::updateGeometry()
{
    Point3D p0, p1, p2, pc;
    Vec3D   va,  vb, vn;

    size_t numfaces = triangles.size();
    faceCentroid.resize(numfaces);
    faceNormal.resize(numfaces);

    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;
    for( int i = 0; i < numfaces; i++) {
        int v0 = triangles[i][0];
        int v1 = triangles[i][1];
        int v2 = triangles[i][2];
        Point3D p0 = nodeCoords[v0];
        Point3D p1 = nodeCoords[v1];
        Point3D p2 = nodeCoords[v2];
        // Set the normal..
        va = getVector(p1,p0);
        vb = getVector(p2,p0);
        vn = getCrossProduct(va, vb);
        getUnitVector(vn);
        faceNormal[i]  = vn;
        normal[0] += vn[0];
        normal[1] += vn[1];
        normal[2] += vn[2];

        // Set the centroid
        pc  = getCentroid(p0,p1,p2);
        faceCentroid[i] = pc;
        center[0] += pc[0];
        center[1] += pc[1];
        center[2] += pc[2];
    }

    normal[0] /= (double)numfaces;
    normal[1] /= (double)numfaces;
    normal[2] /= (double)numfaces;
    getUnitVector(normal);

    center[0] /= (double)numfaces;
    center[1] /= (double)numfaces;
    center[2] /= (double)numfaces;

}
//////////////////////////////////////////////////////////////////////////////////////

void GMesh:: setCenterAt(int faceid)
{
    Point3D pc = facemeshes[faceid]->center;
    for( auto submesh: facemeshes) {
        int numnodes = submesh->nodeCoords.size();
        for( int i = 0; i < numnodes; i++) {
            submesh->nodeCoords[i][0] -= pc[0];
            submesh->nodeCoords[i][1] -= pc[1];
            submesh->nodeCoords[i][2] -= pc[2];
        }
        submesh->updateGeometry();
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void GMesh:: alignAlong( const Vec3D &srcVec, const Vec3D &dstVec)
{
    Vec3D  perpAxis = getCrossProduct( srcVec, dstVec);
    double dl = getMagnitude(perpAxis);
    if( dl < 1.0E-08) return;

    perpAxis[0] /= dl;
    perpAxis[1] /= dl;
    perpAxis[2] /= dl;
    double angle = getVecAngle(srcVec, dstVec);
    if( fabs(angle) < 1.0E-15) return;

    double qcos = cos(0.5*angle);
    double qsin = sin(0.5*angle);

    boost::math::quaternion<double> q(qcos, qsin*perpAxis[0], qsin*perpAxis[1], qsin*perpAxis[2]);
    boost::math::quaternion<double> q1 = boost::math::conj(q);

 size_t index = 0;
    for( auto submesh: facemeshes) {
        size_t numNodes = submesh->nodeCoords.size();
        for( size_t i = 0; i < numNodes; i++) {
            Point3D p3d = submesh->nodeCoords[i];

            boost::math::quaternion<double> v(0.0, p3d[0], p3d[1], p3d[2]);
            boost::math::quaternion<double> result = q*v*q1;
            p3d[0]  = result.R_component_2();
            p3d[1]  = result.R_component_3();
            p3d[2]  = result.R_component_4();
            submesh->nodeCoords[i] = p3d;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
void GMesh::updateGeometry()
{
    for( auto mesh : facemeshes)
        mesh->updateGeometry();
}
/////////////////////////////////////////////////////////////////////////////////////
void GMesh::saveAs( const string &filename)
{
    size_t numnodes = 0;
    size_t numfaces = 0;
    for( auto mesh : facemeshes) {
        numnodes += mesh->nodeCoords.size();
        numfaces += mesh->triangles.size();
    }
    ofstream ofile( filename.c_str(), ios::out);
    ofile << "OFF" << endl;
    ofile << numnodes << " " << numfaces << " 0 " << endl;

    for( auto mesh : facemeshes) {
        for( int i = 0; i < mesh->nodeCoords.size(); i++) {
            ofile << mesh->nodeCoords[i][0] << " "
                  << mesh->nodeCoords[i][1] << " "
                  << mesh->nodeCoords[i][2] << endl;
        }
    }
    size_t offset = 0;
    for( auto mesh : facemeshes) {
        for( int i = 0; i < mesh->triangles.size(); i++)  {
            Array3I ar = mesh->triangles[i];
            ofile << "3 " << ar[0] + offset << " "
                  << ar[1] + offset << " "
                  << ar[2] + offset << endl;
        }
        offset += mesh->nodeCoords.size();
    }
}
////////////////////////////////////////////////////////////////////////////////


void GMesh :: buildGeomTopology()
{
    vector<Point3D> uniquePoints;

    /*
        auto getIndex = []( Point3D &qPoint, vector<Point3D> &container)
        {
           double minDist = std::numeric_limits<double>::max();
           for( auto xyz = container) {
                double dx = xyz[0] - qPoint[0];
                double dy = xyz[1] - qPoint[1];
                double dx = xyz[2] - qPoint[2];
        };

        double eps = 1.0E-06;
        for( auto submesh : facemeshes) {
             for( int i = 0; i < submesh->getSize(0); i++) {
                  Point3D qPoint = submesh->nodeCoords[i];
                  auto gIndex = getClosestPoint( qPoint, uniqueCoords);
                  if( distance(uniqueCoords[gIndex], qPoint) > eps )  {
                       uniquePoints.push_back(qPoint);
                       submesh->globalIndex[i] = uniquePoints.size();
                  } else
                       submesh->globalIndex[i] = globalPoint.first;
              }
        }
    */



    /*
        size_t numnodes = nodeCoords.size();
        size_t numfaces = faces.size();

        vector<set<int>> nodeGeomFaces;
        nodeGeomFaces.resize(numnodes);

        for( int i = 0; i < numfaces; i++) {
            int v0 = faces[i][0];
            int v1 = faces[i][1];
            int v2 = faces[i][2];
            int gid = faceGeomID[i];
            nodeGeomFaces[v0].insert(gid);
            nodeGeomFaces[v1].insert(gid);
            nodeGeomFaces[v2].insert(gid);
        }

        size_t numEdges = geomEdges.size();
        vector<int> commfaces;
        for( int i = 0; i < numEdges; i++) {
            int v0 = geomEdges[i].nodes[0];
            int v1 = geomEdges[i].nodes[1];
            commfaces.clear();
            boost::set_intersection( nodeGeomFaces[v0], nodeGeomFaces[v1],
                                     back_inserter(commfaces));
            assert( commfaces.size() == 2);
            groupFaces[commfaces[0]].geomFaces.insert(commfaces[1]);
            groupFaces[commfaces[1]].geomFaces.insert(commfaces[0]);
        }
    */
}
/////////////////////////////////////////////////////////////////////////////////

void GMesh::readMesh(const string &filename)
{
    /*
        ifstream infile( filename.c_str(), ios::in);
        if( infile.fail() ) {
            cout << "Warning: Cann't open file " << filename << endl;
            return;
        }
        string str;
        size_t  numNodes, numElems;

        infile >> str;
        assert( str == "$MeshFormat");
        int  file_type, data_size, id;
        float version;
        infile >> version >> file_type >> data_size;
        infile >> str;
        assert( str == "$EndMeshFormat");

        infile >> str;
        assert( str == "$Nodes");

        infile >> numNodes;

        nodeCoords.resize(numNodes);
        Point3D p3d;
        double  x, y, z;

        for( size_t i = 0; i < numNodes; i++) {
            infile >> id >> x >> y >> z;
            p3d[0] = x;
            p3d[1] = y;
            p3d[2] = z;
            nodeCoords[i] = p3d;
        }
        infile >> str;
        assert( str == "$EndNodes");
        infile >> str;
        assert( str == "$Elements");

        int elm_number, elm_type, numTags, tagval1, tagval2;
        vector<int> facenodes;
        string line;

        Edge edge;
        infile >> numElems;
        size_t faceIndex = 0;
        size_t edgeIndex = 0;
        for( size_t i = 0; i < numElems; i++) {
            infile >> elm_number >> elm_type >> numTags;
            assert( numTags == 2);
            infile >> tagval1 >> tagval2;

            int valid_elem = 0;
            switch( elm_type ) {
            case 1:
                infile >> edge.nodes[0] >> edge.nodes[1];
                edge.nodes[0] -= 1;
                edge.nodes[1] -= 1;
                geomEdges.push_back(edge);
    //          groupEdges[tagval2].edges.push_back(edgeIndex);
                edgeIndex++;
                break;
            case 2:
                facenodes.resize(3);
                for( size_t j = 0; j < 3; j++)
                    infile >> facenodes[j];
                facenodes[0] -= 1;
                facenodes[1] -= 1;
                facenodes[2] -= 1;
                faces.push_back(facenodes);
                faceGeomID.push_back(tagval2);
                groupFaces[tagval2].faces.push_back(faceIndex);
                faceIndex++;
                break;
            case 3:
                facenodes.resize(4);
                for( size_t j = 0; j < 4; j++)
                    infile >> facenodes[j];
                facenodes[0] -= 1;
                facenodes[1] -= 1;
                facenodes[2] -= 1;
                facenodes[3] -= 1;
                faces.push_back(facenodes);
                faceGeomID.push_back(tagval2);
                groupFaces[tagval2].faces.push_back(faceIndex);
                faceIndex++;
                break;
            default:
                getline( infile, line );
            }
        }

        buildGeomTopology();
        updateGeometry();
    */
}

