#include <QMainWindow>

#include "UfabMainWindow.hpp"
#include "CADPocketDetector.hpp"

int main(int argc, char** argv)
{
/*
    if( argc != 2 ) {
        cout << "Usage : " << argv[0] << " StepFile " << endl;
        return 1;
    }
    CADPocketDetector cpocket;
    cpocket.readSTEPModel(argv[1] );
    cpocket.saveMesh("tmp.off");
    int n = cpocket.detectAll();
    cout << "#Pockets " << n << endl;
    exit(0);
*/
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QApplication application(argc,argv);

    UfabMainWindow  win;
    win.setWindowTitle("UfabViewer");

    win.show();
    return application.exec();
}


