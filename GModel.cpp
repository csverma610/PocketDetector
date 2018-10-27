#include "GModel.hpp"
#include <BRepMesh_IncrementalMesh.hxx>

///////////////////////////////////////////////////////////////////////////////
bool GFace :: isPlanar() const
{
    Handle(Geom_Surface) surf  = BRep_Tool::Surface(occFace);

    GeomAdaptor_Surface ga(surf);

    if( ga.GetType() == GeomAbs_Plane) return 1;

    return 0;
}
///////////////////////////////////////////////////////////////////////////////
bool GFace :: isCircular() const
{
}
///////////////////////////////////////////////////////////////////////////////

GFaceMesh* GFace :: getSurfMesh() {
    GFaceMesh *mesh = new GFaceMesh;

    // https://www.opencascade.com/content/triangulation-objects-step-files
    BRepMesh_IncrementalMesh(occFace, 0.1f);
    TopLoc_Location L;
    Handle (Poly_Triangulation) facing;

    // Generate triangulations for the dace
    facing = BRep_Tool::Triangulation(occFace, L);
    const Poly_Array1OfTriangle &triangles = facing->Triangles();
    const TColgp_Array1OfPnt &nodes = facing->Nodes();
    // double scale = 1.;// .2;
    int triangleID;
    bool isTesselationNull = facing.IsNull();
    int p1, p2, p3;
    Array3I tri;

    if(occFace.Orientation()!=TopAbs_REVERSED) {
        for(int i = 1; i <= facing->NbTriangles(); i++) {
            Poly_Triangle triangle = (facing->Triangles())(i);
            triangle.Get(p1, p2, p3);
            tri[0] = p1 - 1;
            tri[1] = p2 - 1;
            tri[2] = p3 - 1;
            mesh->triangles.push_back(tri);
        }
    }
    else {
        for (int i = 1; i <= facing->NbTriangles(); i++) {
            Poly_Triangle triangle = (facing->Triangles())(i);
            triangle.Get(p1, p2, p3);
            tri[0] = p3 - 1;
            tri[1] = p2 - 1;
            tri[2] = p1 - 1;
            mesh->triangles.push_back(tri);
        }
    }
    Point3D xyz;
    for(int i = 1; i <= facing->NbNodes(); i++) {
        gp_Pnt pp = (facing->Nodes())(i);
        xyz[0] = pp.X();
        xyz[1] = pp.Y();
        xyz[2] = pp.Z();
        mesh->nodeCoords.push_back(xyz);
    }
    return mesh;
}

///////////////////////////////////////////////////////////////////////////////
vector<GEdge*> GFace :: getNewLoop( std::list<GEdge*> &edgelist) const
{
    vector<GEdge*> newloop;

    if( edgelist.empty() ) return newloop;

    if( edgelist.size() == 1) {
        newloop.push_back(edgelist.front());
        edgelist.pop_front();
        return newloop;
    }

    GEdge *currEdge = edgelist.front();
    edgelist.pop_front();
    newloop.push_back(currEdge);

    while(1) {
        int nsize = edgelist.size();
        if( nsize == 0) break;
        int found = 0;
        for( size_t i = 0; i < nsize; i++) {
            GEdge *nextEdge = edgelist.front();
            edgelist.pop_front();
            if( nextEdge->isAdjacent(currEdge) ) {
                newloop.push_back(nextEdge);
                currEdge = nextEdge;
                found    = 1;
                break;
            } else
                edgelist.push_back(nextEdge);
        }
        if( !found ) break;
    }
    return newloop;
}
///////////////////////////////////////////////////////////////////////////////
vector<GEdgeSequence> GFace :: getWires() const
{
    
    vector<vector<GEdge*>> wires;

    if( boundEdges.size() == 1) {
        wires.resize(1);
        wires[0].resize(1);
        wires[0][0] = boundEdges[0];
        return wires;
    }

    list<GEdge*> edgelist;
    for( auto edge : boundEdges)
        edgelist.push_back(edge);

    cout << "Detectiing .... " << endl;

    while(1) {
        vector<GEdge*> newloop = getNewLoop(edgelist);
        if( newloop.empty() ) break;
        wires.push_back(newloop);
    }
    cout << "#Wires " << wires.size() << endl;

    return wires;
}
///////////////////////////////////////////////////////////////////////////////
int GModel:: readOCCModel( const string &filename)
{
    STEPControl_Reader reader;
    IFSelect_ReturnStatus stat;
    stat = reader.ReadFile(filename.c_str());
    assert( stat == IFSelect_RetDone);

    Standard_Boolean failsonly = Standard_False;
    IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
    reader.PrintCheckLoad (failsonly, mode);

    Standard_Integer nRoots = reader.TransferRoots();
    //selects all IGES entities (including non visible ones) in the
    //file and puts them into a list called MyList,

    assert( nRoots >0 );

    // Handle STEP Scale here.
    gp_Pnt Origin;
    gp_Trsf scale;
    scale.SetScale (Origin, 1.0);
//
    TopoDS_Shape sh = reader.OneShape();
    BRepBuilderAPI_Transform trans(sh, scale);
    occShape = trans.Shape();
//  occShape = sh;
    buildTopology();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void GModel:: setOCCShape( TopoDS_Shape shape)
{
    occShape = shape;
    buildTopology();
}

///////////////////////////////////////////////////////////////////////////////

void GModel:: getOCCInfo()
{
    cout << "Step Entities  : " << endl;
    cout << "   #Compounts  : " << occCompounds.size() << endl;
    cout << "   #CompSolids : " << occCompsolids.size() << endl;
    cout << "   #Solids     : " << occSolids.size() << endl;
    cout << "   #Shells     : " << occShells.size() << endl;
    cout << "   #wires      : " << occWires.size()  << endl;
    cout << "Step Geometric Entities  : " << endl;
    cout << "   #Faces  : " << occFaces.size() << endl;
    cout << "   #Edges  : " << occEdges.size() << endl;
    cout << "   #Nodes  : " << occNodes.size() << endl;
}
/////////////////////////////////////////////////////////////////////////

double GModel:: getDistance( const Point3D &p0, const Point3D &p1) const
{
    double dx = p0[0] - p1[0];
    double dy = p0[1] - p1[1];
    double dz = p0[2] - p1[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
}
/////////////////////////////////////////////////////////////////////////

double GModel:: getDistance( const GNode *anode, const GNode *bnode) const
{
    const Point3D p0 = anode->getCoords();
    const Point3D p1 = bnode->getCoords();
    return getDistance(p0, p1);
}
/////////////////////////////////////////////////////////////////////////

GNode* GModel:: getClosestNode( const Point3D &queryPoint) const
{
    GNode *qnode = nullptr;
    double mindist = numeric_limits<double>::max();
    for( auto currnode : gNodes) {
        const Point3D currPoint = currnode->getCoords();
        double dist  = getDistance(currPoint, queryPoint);
        if( dist < mindist) {
            mindist = dist;
            qnode = currnode;
        }
    }
    return qnode;
}
/////////////////////////////////////////////////////////////////////////

void GModel:: buildNodes()
{
    gNodes.clear();

    if( occNodes.empty() ) return;

    deque<ONode> nodesQ;

    size_t numNodes =  occNodes.size();

    ONode newnode;
    for( int i = 0; i < numNodes; i++) {
        newnode.active = 1;
        newnode.occNode  = occNodes[i];
        const gp_Pnt vpos = BRep_Tool::Pnt( occNodes[i] );
        newnode.xyz[0] = vpos.X();
        newnode.xyz[1] = vpos.Y();
        newnode.xyz[2] = vpos.Z();
        nodesQ.push_back(newnode);
    }

    double eps = 1.0E-15;
    GNode *jnode;
    int index = 0;
    while(1) {
        if( nodesQ.empty() ) break;
        ONode currnode = nodesQ.front();
        nodesQ.pop_front();
        if( currnode.active ) {
            jnode = new GNode;
            jnode->setID(index++);
            jnode->setCoords(currnode.xyz);
            jnode->addOCCNode(currnode.occNode);
            size_t nSize = nodesQ.size();
            for( size_t i = 0; i < nSize; i++) {
                double dist = currnode.getDistance(nodesQ[i]);
                if( dist < eps) {
                    jnode->addOCCNode( nodesQ[i].occNode);
                    nodesQ[i].active = 0;
                }
            }
            gNodes.push_back(jnode);
        }
    }

    cout << "#GeomNodes " << gNodes.size() << endl;
}

/////////////////////////////////////////////////////////////////////////
GEdge* GModel:: getEdgeOf( const GNode *v0, const GNode *v1) const
{
    for(auto edge : gEdges)
        if( edge->isSame(v0,v1) ) return edge;
    return nullptr;
}
/////////////////////////////////////////////////////////////////////////
GEdge* GModel:: getEdgeOf( TopoDS_Edge e) const
{
    for(auto edge : gEdges)
        if( edge->hasOCCEdge(e)) return edge;
    return nullptr;
}
/////////////////////////////////////////////////////////////////////////

void GModel:: buildEdges()
{
    gEdges.clear();
    cout  << "Building unique edges " << endl;

    TopoDS_Vertex  v0, v1;
    gp_Pnt   gPoint;
    size_t numEdges = occEdges.size();
    Point3D xyz;

    for( int i = 0; i < numEdges; i++) {
        TopoDS_Edge edge = occEdges[i];
        TopExp::Vertices(edge, v0, v1);
        gPoint = BRep_Tool::Pnt( v0 );
        xyz[0] =  gPoint.X();
        xyz[1] =  gPoint.Y();
        xyz[2] =  gPoint.Z();
        GNode *node1 = getClosestNode(xyz);

        gPoint = BRep_Tool::Pnt( v1 );
        xyz[0] =  gPoint.X();
        xyz[1] =  gPoint.Y();
        xyz[2] =  gPoint.Z();
        GNode *node2 = getClosestNode(xyz);

        GEdge *jedge  = getEdgeOf(node1,node2);
        if( jedge ==  nullptr) {
            jedge = GEdge::newObject(node1,node2);
            gEdges.push_back(jedge);
        }
        jedge->addOCCEdge(edge);
    }
    cout << "#GeomEdges " << gEdges.size() << endl;
}
/////////////////////////////////////////////////////////////////////////
void GModel:: buildFaces()
{
    gFaces.clear();
    TopExp_Explorer expl;

    int numFaces = occFaces.size();

    int index = 0;
    for( int i = 0; i < numFaces; i++) {
        TopoDS_Face oface = occFaces[i];
        GFace *jface = new GFace;
        jface->setOCCFace(oface);
        jface->setID( index++);
        for(expl.Init(oface, TopAbs_EDGE); expl.More(); expl.Next()) {
            TopoDS_Edge oedge = TopoDS::Edge(expl.Current());
            GEdge *jedge = getEdgeOf(oedge);
            if(jedge) {
                jedge->addRelation(jface);
                jface->addEdge(jedge);
            }
        }
        gFaces.push_back(jface);
    }

    for( auto iface : gFaces) {
        for( auto edge : iface->getEdges() ) {
            auto adjfaces = edge->getRelations2();
            for( auto jface : adjfaces) {
                if( jface != iface)
                    iface->addNeighbor(jface);
            }
        }
    }


    /*
        for( auto edge : gEdges) {
            vector<GFace*> faces = edge->getRelations2();
            if(faces.size() == 2) {
                faces[0]->addNeighbor(faces[1]);
                faces[1]->addNeighbor(faces[0]);
            } else {
                cout << "Warning: Edge has more faces " << endl;
                for( auto f : faces)
                     cout << f->getID() << " ";
                cout << endl;

            }
        }
    */

    cout << "#GeomFaces " << gFaces.size() << endl;
}

/////////////////////////////////////////////////////////////////////////

GFaceSequence GModel :: getNeighbors(GFace *seedface) const
{
    GFaceSequence neighfaces;
    if( seedface == nullptr) return neighfaces;

    GEdgeSequence boundedges = seedface->getEdges();

    GFaceSequence edgefaces;
    GFaceSet faceSet;

    for( auto edge: boundedges) {
        edgefaces = edge->getRelations2();
        for( auto f: edgefaces) faceSet.insert(f);
    }
    faceSet.erase(seedface);

    edgefaces.clear();
    for( auto f: faceSet)  edgefaces.push_back(f);

    return edgefaces;
}

/////////////////////////////////////////////////////////////////////////

void JGraph::saveAs( const string &filename)
{
    ofstream ofile(filename.c_str(), ios::out);

    ofile << "graph STEP {" << endl;
    for( int i = 0; i < nodes.size(); i++) {
        JDualNode *node = nodes[i];
        if( node->getFace()->isPlanar() )
            ofile << node->getID()+1 << " [shape=box, style=filled, fillcolor = red, bb = 10]" << endl;
        else
            ofile << node->getID()+1 << " [shape=circle]" << endl;
    }

    for( int i = 0; i < edges.size(); i++) {
        JDualEdge *edge = edges[i];
        int v0   = edge->getNodeAt(0)->getID();
        int v1   = edge->getNodeAt(1)->getID();
        ofile << v0+1 << " -- " << v1+1 << ";" << endl;
    }
    ofile << "}" << endl;
}

/////////////////////////////////////////////////////////////////////////
void GModel:: buildTopology()
{
    extract_geometrical_shapes();

    buildNodes();
    buildEdges();
    buildFaces();
}

/////////////////////////////////////////////////////////////////////////
GMesh *GModel:: getSurfMesh()
{
    GMesh *modelmesh = new GMesh;
    size_t numFaces = gFaces.size();
    for( int i = 0; i < numFaces; i++) {
        GFaceMesh *fmesh = gFaces[i]->getSurfMesh();
        fmesh->groupID = i;
        modelmesh->facemeshes.push_back(fmesh);
    }
    return modelmesh;
}
/////////////////////////////////////////////////////////////////////////

void GModel::extract_geometrical_shapes()
{
    occFaces.clear();
    occEdges.clear();
    occNodes.clear();

    TopExp_Explorer exp;
    for (exp.Init(occShape, TopAbs_FACE); exp.More(); exp.Next())
    {
        occFaces.push_back(TopoDS::Face(exp.Current()));
    }
    for (exp.Init(occShape, TopAbs_EDGE); exp.More(); exp.Next())
    {
        occEdges.push_back(TopoDS::Edge(exp.Current()));
    }
    for (exp.Init(occShape, TopAbs_VERTEX); exp.More(); exp.Next())
    {
        occNodes.push_back(TopoDS::Vertex(exp.Current()));
    }
}

/////////////////////////////////////////////////////////////////////////////////

void GModel::extract_compound_shapes()
{
    occCompounds.resize(0);
    occCompsolids.resize(0);
    occSolids.resize(0);
    occShells.resize(0);
    occWires.resize(0);

    TopExp_Explorer exp;
    for (exp.Init(occShape, TopAbs_COMPOUND); exp.More(); exp.Next())
    {
        occCompounds.push_back(TopoDS::Compound(exp.Current()));
    }
    for (exp.Init(occShape, TopAbs_COMPSOLID); exp.More(); exp.Next())
    {
        occCompsolids.push_back(TopoDS::CompSolid(exp.Current()));
    }
    for (exp.Init(occShape, TopAbs_SOLID); exp.More(); exp.Next())
    {
        occSolids.push_back(TopoDS::Solid(exp.Current()));
    }
    for (exp.Init(occShape, TopAbs_SHELL); exp.More(); exp.Next())
    {
        occShells.push_back(TopoDS::Shell(exp.Current()));
    }
    for (exp.Init(occShape, TopAbs_WIRE); exp.More(); exp.Next())
    {
        occWires.push_back(TopoDS::Wire(exp.Current()));
    }
}
/////////////////////////////////////////////////////////////////////////////////

gp_Pnt GModel::point(const Point3D &p)
{
    return gp_Pnt(p[0], p[1], p[2]);
}

/////////////////////////////////////////////////////////////////////////////////
Point3D GModel::point(const gp_Pnt &p)
{
    Point3D xyz;
    xyz[0] = p.X();
    xyz[1] = p.Y();
    xyz[2] = p.Z();
    return xyz;
}
/////////////////////////////////////////////////////////////////////////////////
