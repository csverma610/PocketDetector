#include "UfabViewer.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

UfabViewer :: UfabViewer(QWidget *w) : QGLViewer(w)
{
}

///////////////////////////////////////////////////////////////////////////////

void UfabViewer :: readSTEPModel( const string &filename)
{
    if( mesh ) delete mesh;
    cpktDetector.readSTEPModel(filename);
    mesh =  cpktDetector.getSurfMesh();
    initMesh();
}

///////////////////////////////////////////////////////////////////////////////

void UfabViewer :: initMesh()
{
    if( mesh == nullptr) return;
    double xmin =  0.99*std::numeric_limits<double>::max();
    double xmax =  -0.99*std::numeric_limits<double>::max();
    double ymin = xmin;
    double zmin = xmin;
    double ymax = xmax;
    double zmax = zmax;

    for( auto submesh: mesh->facemeshes) {
        size_t numNodes = submesh->nodeCoords.size();
        size_t numFaces = submesh->triangles.size();
        for( size_t i = 0; i < numNodes; i++) {
            Point3D p = submesh->nodeCoords[i];
            xmin = std::min(p[0], xmin);
            ymin = std::min(p[1], ymin);
            zmin = std::min(p[2], zmin);

            xmax = std::max(p[0], xmax);
            ymax = std::max(p[1], ymax);
            zmax = std::max(p[2], zmax);
        }
    }

    double xcenter = 0.5*(xmax+xmin);
    double ycenter = 0.5*(ymax+ymin);
    double zcenter = 0.5*(zmax+zmin);

    for( auto submesh: mesh->facemeshes) {
        size_t numNodes = submesh->nodeCoords.size();
        for( size_t i = 0; i < numNodes; i++) {
            submesh->nodeCoords[i][0] -= xcenter;
            submesh->nodeCoords[i][1] -= ycenter;
            submesh->nodeCoords[i][2] -= zcenter;
        }
    }
    double maxlen = fabs(xmax-xmin);
    maxlen = max(maxlen, fabs(ymax-ymin));
    maxlen = max(maxlen, fabs(zmax-zmin));

    for( auto submesh: mesh->facemeshes) {
        size_t numNodes = submesh->nodeCoords.size();
        for( size_t i = 0; i < numNodes; i++) {
            submesh->nodeCoords[i][0] /= maxlen;
            submesh->nodeCoords[i][1] /= maxlen;
            submesh->nodeCoords[i][2] /= maxlen;
        }
    }
    mesh->updateGeometry();
    assignEdgeColor();
    assignFaceColor();
}

/////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::assignEdgeColor()
{
    /*
            Color color;
            color[0] = 0.0;
            color[1] = 0.0;
            color[2] = 0.0;
            color[3] = 1.0;

            size_t numEdges = mesh->geomEdges.size();
            for( size_t i = 0; i < numEdges; i++)
                mesh->geomEdges[i].color = color;
    */
}

//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::assignFaceColor()
{
    if( mesh == nullptr) return;
    Color color;
    for( auto submesh : mesh->facemeshes) {
        color[0] = max(0.2, drand48() );
        color[1] = max(0.2, drand48() );
        color[2] = max(0.2, drand48() );
        color[3] = 1.0;
        submesh->faceColor = color;
        submesh->orgColor = color;
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawCenters()
{
    glPointSize(5);

    for( auto submesh : mesh->facemeshes) {
        glBegin(GL_POINTS);
        glVertex3dv( &submesh->center[0] );
        glEnd();
    }
    glPointSize(1);
}

//////////////////////////////////////////////////////////////////////////////////////

void
UfabViewer::detectAll()
{
    int n = cpktDetector.detectAll();

    vector<CADPocket> pockets = cpktDetector.getPockets();

    Color color;
    color[0] = 0.2;
    color[1] = 0.2;
    color[2] = 0.2;
    color[3] = 1.0;

    for( auto submesh : mesh->facemeshes) {
        submesh->faceColor = color;
    }

    color[0] = 1.0;
    color[1] = 1.0;
    color[2] = 1.0;
    color[3] = 1.0;
    for( int i = 0; i < pockets.size(); i++) {
        int nface = pockets[i].faces.size();
        int id = pockets[i].faces[0];
        mesh->facemeshes[id]->faceColor = color;
    }

    update();
}
//////////////////////////////////////////////////////////////////////////////////////

void
UfabViewer::keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_W) {
        if( faceSelected >= 0) {
            mesh->facemeshes[faceSelected]->faceColor  =
                mesh->facemeshes[faceSelected]->orgColor;
            faceSelected = -1;
            setSelectedName(-1);
        }
        renderMeshEdges = !renderMeshEdges;
    }

    if( e->key() == Qt::Key_F) {
        if( faceSelected >= 0) {
            faceSelected = -1;
            setSelectedName(-1);
        }
        for( auto submesh : mesh->facemeshes) {
            submesh->faceColor = submesh ->orgColor;
        }
    }

    if( e->key() == Qt::Key_R) {
        for( auto submesh : mesh->facemeshes) {
            int numFaces = submesh->getSize(2);
            for( int i = 0; i < numFaces; i++) {
                Array3I ar = submesh->triangles[i];
                submesh->triangles[i][0] = ar[0];
                submesh->triangles[i][1] = ar[2];
                submesh->triangles[i][2] = ar[1];
            }
            submesh->updateGeometry();
        }
    }

    if( e->key() == Qt::Key_N)
        renderMeshNormals = !renderMeshNormals;

    if( e->key() == Qt::Key_P) {
        if( faceSelected >= 0) {
            cpktDetector.isPocket(faceSelected);
        }
    }

    if( e->key() == Qt::Key_O) {
        detectAll();
    }

    if( e->key() == Qt::Key_Q) {
        setSelectedName(-1);
        cout << "Provide the faceID " << endl;
        int id;
        cin >> id;
        selectFace(id);
    }
    update();
    QGLViewer::keyPressEvent(e);

}

//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawAxes()
{
    glDisable( GL_LIGHTING);
    Point3D p0, p1;
    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    glColor3f( 1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3dv( &p0[0] );
    p1[0] = p0[0] + 10000;
    p1[1] = p0[1];
    p1[2] = p0[2];
    glVertex3dv( &p1[0] );
    glEnd();

    glColor3f( 0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3dv( &p0[0] );
    p1[0] = p0[0];
    p1[1] = p0[1] + 10000;
    p1[2] = p0[2];
    glVertex3dv( &p1[0] );
    glEnd();

    glColor3f( 0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3dv( &p0[0] );
    p1[0] = p0[0];
    p1[1] = p0[1];
    p1[2] = p0[2] + 10000;
    glVertex3dv( &p1[0] );
    glEnd();
    glEnable( GL_LIGHTING);

}
//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawTargetAxes( const Point3D &center, const Vec3D &vec)
{
    glDisable( GL_LIGHTING);
    Point3D p1;

    glColor3f( 1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3dv( &center[0] );
    p1[0] = center[0]  + vec[0] + 10000;
    p1[1] = center[1]  + vec[1] ;
    p1[2] = center[2]  + vec[2];
    glVertex3dv( &p1[0] );
    glEnd();

    glColor3f( 0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3dv( &center[0] );
    p1[0] = center[0]  + vec[0];
    p1[1] = center[1]  + vec[1] + 10000;
    p1[2] = center[2]  + vec[2];
    glVertex3dv( &p1[0] );
    glEnd();

    glColor3f( 0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3dv( &center[0] );
    p1[0] = center[0]  + vec[0];
    p1[1] = center[1]  + vec[1];
    p1[2] = center[2]  + vec[2] + 10000;
    glVertex3dv( &p1[0] );
    glEnd();
    glEnable( GL_LIGHTING);
}

////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawNodes()
{
    glDisable( GL_LIGHTING);
    glPointSize(1);
    glColor3f( 1.0, 1.0, 1.0);
    for( auto submesh : mesh->facemeshes) {
        if( submesh->displayNodes) {
            size_t numNodes = submesh->nodeCoords.size();
            glBegin(GL_POINTS);
            for( size_t i = 0; i < numNodes; i++) {
                glVertex3dv(&submesh->nodeCoords[i][0]);
            }
            glEnd();
        }
    }
    glEnable( GL_LIGHTING);
}
////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawEdges()
{
    /*
        glDisable( GL_LIGHTING);
        glLineWidth(2.0);
        size_t numEdges = mesh->geomEdges.size();
        for( size_t i = 0; i < numEdges; i++) {
            if( mesh->geomEdges[i].renderBit ) {
                glColor3fv( &mesh->geomEdges[i].color[0] );
                glBegin(GL_LINES);
                int v0 = mesh->geomEdges[i].nodes[0];
                glVertex3dv(&mesh->nodeCoords[v0][0]);
                int v1 = mesh->geomEdges[i].nodes[1];
                glVertex3dv(&mesh->nodeCoords[v1][0]);
                glEnd();
            }
        }
        glEnable( GL_LIGHTING);
    */
}

//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawFaces()
{
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for( auto submesh : mesh->facemeshes) {
        if( submesh->displayFaces) {
            size_t numFaces = submesh->triangles.size();
            glColor4fv( &submesh->faceColor[0] );
            for( size_t i = 0; i < numFaces; i++) {
                glBegin(GL_TRIANGLES);
                int v0 = submesh->triangles[i][0];
                int v1 = submesh->triangles[i][1];
                int v2 = submesh->triangles[i][2];
                glNormal3dv(&submesh->faceNormal[i][0] );
                glVertex3dv(&submesh->nodeCoords[v0][0]);
                glVertex3dv(&submesh->nodeCoords[v1][0]);
                glVertex3dv(&submesh->nodeCoords[v2][0]);
                glEnd();
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawFaceEdges()
{
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glColor3f( 0.1, 0.1, 0.1);
    glLineWidth(2.0);

    for( auto submesh : mesh->facemeshes) {
        if( submesh->displayFaces ) {
            size_t numFaces = submesh->triangles.size();
            for( size_t i = 0; i < numFaces; i++) {
                glBegin(GL_TRIANGLES);
                int v0 = submesh->triangles[i][0];
                int v1 = submesh->triangles[i][1];
                int v2 = submesh->triangles[i][2];
                glNormal3dv(&submesh->faceNormal[i][0] );
                glVertex3dv(&submesh->nodeCoords[v0][0]);
                glVertex3dv(&submesh->nodeCoords[v1][0]);
                glVertex3dv(&submesh->nodeCoords[v2][0]);
                glEnd();
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////
void UfabViewer::drawFacesWithName()
{
    glDisable( GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    int index = 0;
    for( auto submesh : mesh->facemeshes) {
        if( submesh->displayFaces) {
            glPushName(index);
            for( int i = 0; i < submesh->getSize(2); i++) {
                glBegin(GL_TRIANGLES);
                int v0 = submesh->triangles[i][0];
                int v1 = submesh->triangles[i][1];
                int v2 = submesh->triangles[i][2];
                glVertex3dv(&submesh->nodeCoords[v0][0]);
                glVertex3dv(&submesh->nodeCoords[v1][0]);
                glVertex3dv(&submesh->nodeCoords[v2][0]);
                glEnd();
            }
            glPopName();
        }
        index++;
    }
    glEnable( GL_LIGHTING);
    glLineWidth(1.0);
}
//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawFacesNormal()
{
    glDisable( GL_LIGHTING);
    glLineWidth(1.0);
    glColor3f( 1.0, 1.0, 1.0);

    Point3D head;
    double len = faceNormalsLength;
    len = 0.05;
    for( auto submesh : mesh->facemeshes) {
        if( !submesh->displayFaces ) continue;
        size_t numFaces = submesh->triangles.size();
        for( size_t i = 0; i < numFaces; i++) {
            glBegin(GL_LINES);
            glVertex3dv(&submesh->faceCentroid[i][0]);
            head[0] = submesh->faceCentroid[i][0] + len*submesh->faceNormal[i][0];
            head[1] = submesh->faceCentroid[i][1] + len*submesh->faceNormal[i][1];
            head[2] = submesh->faceCentroid[i][2] + len*submesh->faceNormal[i][2];
            glVertex3dv(&head[0]);
            glEnd();
        }
    }
    glEnable( GL_LIGHTING);
}
//////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::draw()
{
    if( mesh == nullptr) return;
    glClearColor( 0.8, 0.8, 0.8, 1.0);
    glPushMatrix();
//  glEnable(GL_CULL_FACE);
    drawAxes();
//    drawNodes();
    drawEdges();
    if( renderMeshEdges ) drawFaceEdges();
    drawFaces();
    if(renderMeshNormals) drawFacesNormal();
    glPopMatrix();
}

/////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::drawWithNames()
{
    if( mesh == nullptr) return;
    drawFacesWithName();
}

/////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::init() {
    restoreStateFromFile();
}

/////////////////////////////////////////////////////////////////////////////////////

void UfabViewer:: mousePressEvent( QMouseEvent *e)
{
    QGLViewer::mousePressEvent(e);
}

/////////////////////////////////////////////////////////////////////////////////////
void UfabViewer:: selectFace( int faceID )
{
    Color white;
    white[0] = 1.0;
    white[1] = 1.0;
    white[2] = 1.0;
    white[3] = 1.0;

    Color gray;
    gray[0] = 0.2;
    gray[1] = 0.2;
    gray[2] = 0.2;
    gray[3] = 1.0;

    faceSelected = faceID;
    cpktDetector.reOrient(faceSelected);

    for( auto submesh : mesh->facemeshes)
        submesh->faceColor = gray;

    GModel *gModel  = cpktDetector.getGeomModel();
    GFace *thisface = gModel->getFaceAt(faceSelected);

    vector<GFace*> neighfaces = thisface->getAdjacentFaces();
    cout << neighfaces.size() << endl;
    for( auto gface : neighfaces) {
        int id =  gface->getID();
        mesh->facemeshes[id]->displayFaces = 1;
        mesh->facemeshes[id]->faceColor = mesh->facemeshes[id]->orgColor;
    }
    mesh->facemeshes[faceSelected]->faceColor = white;
    update();
}

/////////////////////////////////////////////////////////////////////////////////////

void UfabViewer::mouseReleaseEvent( QMouseEvent *e)
{
    int id = this->selectedName();
    if( id == faceSelected ) return;

    if( faceSelected >= 0) {
        mesh->facemeshes[faceSelected]->faceColor =
            mesh->facemeshes[faceSelected]->orgColor;
    }

    if( id >= 0) {
        selectFace(id);
    }
    QGLViewer::mouseReleaseEvent(e);
}
/////////////////////////////////////////////////////////////////////////////////////

