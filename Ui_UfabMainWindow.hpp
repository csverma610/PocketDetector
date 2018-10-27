/********************************************************************************
** Form generated from reading UI file 'ufabmainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_UFABMAINWINDOW_H
#define UI_UFABMAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QWidget>
#include "QGLViewer/qglviewer.h"

#include "UfabViewer.hpp"

QT_BEGIN_NAMESPACE

class Ui_UfabMainWindow
{
public:
    QAction *actionOpen;
    QAction *actionQuit;
    QAction *actionRender_Properties;
    QAction *actionFeatures_Detection;
    QWidget *centralwidget;
    UfabViewer *viewer;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuTools;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *UfabMainWindow)
    {
        if (UfabMainWindow->objectName().isEmpty())
            UfabMainWindow->setObjectName(QStringLiteral("UfabMainWindow"));
        UfabMainWindow->resize(800, 600);
        actionOpen = new QAction(UfabMainWindow);
        actionOpen->setObjectName(QStringLiteral("actionOpen"));
        actionQuit = new QAction(UfabMainWindow);
        actionQuit->setObjectName(QStringLiteral("actionQuit"));
        actionRender_Properties = new QAction(UfabMainWindow);
        actionRender_Properties->setObjectName(QStringLiteral("actionRender_Properties"));
        actionFeatures_Detection = new QAction(UfabMainWindow);
        actionFeatures_Detection->setObjectName(QStringLiteral("actionFeatures_Detection"));
        centralwidget = new QWidget(UfabMainWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        viewer = new UfabViewer(centralwidget);
        viewer->setObjectName(QStringLiteral("viewer"));
        viewer->setGeometry(QRect(8, 7, 782, 590));
        UfabMainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(UfabMainWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 800, 19));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        menuTools = new QMenu(menubar);
        menuTools->setObjectName(QStringLiteral("menuTools"));
        UfabMainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(UfabMainWindow);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        UfabMainWindow->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuTools->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionQuit);
        menuTools->addAction(actionRender_Properties);
        menuTools->addAction(actionFeatures_Detection);

        retranslateUi(UfabMainWindow);

        QMetaObject::connectSlotsByName(UfabMainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *UfabMainWindow)
    {
        UfabMainWindow->setWindowTitle(QApplication::translate("UfabMainWindow", "MainWindow", Q_NULLPTR));
        actionOpen->setText(QApplication::translate("UfabMainWindow", "Open", Q_NULLPTR));
        actionQuit->setText(QApplication::translate("UfabMainWindow", "Quit", Q_NULLPTR));
        actionRender_Properties->setText(QApplication::translate("UfabMainWindow", "Render Properties", Q_NULLPTR));
        actionFeatures_Detection->setText(QApplication::translate("UfabMainWindow", "Features Detection", Q_NULLPTR));
        menuFile->setTitle(QApplication::translate("UfabMainWindow", "File", Q_NULLPTR));
        menuTools->setTitle(QApplication::translate("UfabMainWindow", "Tools", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
class UfabMainWindow: public Ui_UfabMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_UFABMAINWINDOW_H
