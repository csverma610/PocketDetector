#include "UfabMainWindow.hpp"

UfabMainWindow :: UfabMainWindow( QWidget *parent) : QMainWindow(parent)
{
    setupUi(this);
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

UfabMainWindow :: ~UfabMainWindow()
{
//    clearup();
}
///////////////////////////////////////////////////////////////////////////////
void UfabMainWindow :: resizeEvent( QResizeEvent *e)
{
    int w = centralwidget->width();
    int h = centralwidget->height();
    viewer->resize(w,h);
}
///////////////////////////////////////////////////////////////////////////////

void UfabMainWindow :: clearup()
{
}
///////////////////////////////////////////////////////////////////////////////

void UfabMainWindow :: Quit()
{
    clearup();
//    cout << " The Coredump is because of Log4Cxx library  " << endl;
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

void UfabMainWindow :: openFileDialog()
{
    static QString lastSelectedDirectory;

    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Mesh File "),
                   lastSelectedDirectory,
                   *new QString( "Model Format (*.step)"));

    string fileName = qstr.toUtf8().constData();
    viewer->readSTEPModel( fileName );
}

///////////////////////////////////////////////////////////////////////////////

void UfabMainWindow :: makeConnections()
{
    connect( actionOpen,  SIGNAL(triggered() ), this, SLOT( openFileDialog() ) );
    connect( actionQuit,  SIGNAL(triggered()), this, SLOT( Quit() ));

}
///////////////////////////////////////////////////////////////////////////////
