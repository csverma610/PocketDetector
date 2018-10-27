#include "GModel.hpp"

//////////////////////////////////////////////////////////////////////////////////////
void ReorientModel( JFace *seedface, GModel *model)
{
}

//////////////////////////////////////////////////////////////////////////////////////
void BFSearch( JDualNode *node, JGraph *graph, JFaceSequence &sidefaces)
{

}

//////////////////////////////////////////////////////////////////////////////////////
JPocket* SearchPocket( JFace *seedface, GModel *model )
{
    if( seedface == nullptr) return nullptr;
    if( !seedface->isPlanar() ) return nullptr;
    if( !seedface->isCircular() ) return nullptr;

    ReorientModel(seedface, model);

    JGraph *graph = model->getGraph();
    JDualNode *dualnode = seedface->getDualNode();

    JFaceSequence sidefaces;
    BFSearch(dualnode, graph, sidefaces);

    if( sidefaces.empty() ) return nullptr;

    JPocket *newpocket = new JPocket;
    newpocket->baseface = seedface;
    newpocket->sidefaces = sidefaces;
    return newpocket;
}
//////////////////////////////////////////////////////////////////////////////////////
vector<JPocket*> SearchPockets(GModel *gmodel)
{
    vector<JPocket*> pockets;
    size_t numfaces = gmodel->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        JFace *face = gmodel->getFaceAt(i);
        JPocket *pocket = SearchPocket(face, gmodel);
        if( pocket ) pockets.push_back(pocket);
    }
}
//////////////////////////////////////////////////////////////////////////////////////
/*
int main(int argc, char **argv)
{
    if( argc != 2) {
        cout << "Usage: " << argv[0] << " Stepfile" << endl;
        return 1;
    }

    GModel *gmodel = new GModel();
    assert( gmodel ) ;
    gmodel->readSTEP( argv[1] );
    gmodel->getOCCInfo();
}
*/
