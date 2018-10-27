#include <stdio.h>
#include <stdlib.h>

#include <QMainWindow>
#include <QMessageBox>
#include <QFileDialog>
#include <QAction>

#include <iostream>

#include "Ui_UfabMainWindow.hpp"

class UfabMainWindow  : public QMainWindow, public Ui::UfabMainWindow {
    Q_OBJECT

public:
    explicit UfabMainWindow( QWidget *parent = NULL );
    ~UfabMainWindow();

private slots:
    void openFileDialog();
    void resizeEvent( QResizeEvent *e);
    void Quit();

private:
    void clearup();
    void makeConnections();
};

