OBJS =  UfabMainWindow.o UfabViewer.o moc_UfabMainWindow.o main.o  CADPocketDetector.o GModel.o GMesh.o Quaternion.o vec.o
#OBJS =  main.o CADPocketDetector.o GModel.o GMesh.o

CPPFLAGS = -pg -fPIC -std=c++11 -I$(DEAL2_DIR)/include -I$(OCC_DIR)/include/opencascade -I$(TBB_DIR)/include
CPPFLAGS += -DQT_NO_DEBUG -DQT_OPENGL_LIB -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_XML_LIB -DQT_CORE_LIB
CPPFLAGS += -DOCC_VERSION_MAJOR=7

CPPFLAGS += -I$(QTDIR)/include -I$(QTDIR)/include/QtCore -I$(QTDIR)/include/QtWidgets -I$(QTDIR)/include/QtXml -I$(QTDIR)/include/QtOpenGL -I$(QTDIR)/include/QtGui 

CPPFLAGS += -I$(BOOST_DIR)/include 
CPPFLAGS += -I$(QGLVIEWER_DIR)/include

$LIBS     = -L$(DEAL2_DIR)/lib -ldeal_II -ldeal_II.g
LIBS    += -L$(OCC_DIR)/lib -lTKBinL -lTKBin -lTKBinTObj -lTKBinXCAF -lTKBool -lTKBO -lTKBRep -lTKCAF -lTKCDF -lTKDCAF -lTKDraw -lTKernel -lTKFeat -lTKFillet -lTKG2d -lTKG3d -lTKGeomAlgo -lTKGeomBase -lTKHLR -lTKIGES -lTKLCAF -lTKMath -lTKMesh -lTKMeshVS -lTKOffset -lTKOpenGl -lTKPrim -lTKQADraw -lTKService -lTKShHealing -lTKStdL -lTKStd -lTKSTEP -lTKSTEP209 -lTKSTEPAttr -lTKSTEPBase  -lTKSTL -lTKTObjDRAW -lTKTObj -lTKTopAlgo -lTKTopTest -lTKV3d -lTKVCAF -lTKViewerTest -lTKVRML -lTKXCAF -lTKXDEDRAW -lTKXDEIGES -lTKXDESTEP -lTKXMesh -lTKXmlL -lTKXml -lTKXmlTObj -lTKXmlXCAF -lTKXSBase -lTKXSDRAW

LIBS += -L$(QTDIR)/lib -lQt5Core -lQt5Xml -lQt5OpenGL -lQt5Widgets -lQt5Gui -lGL -lGLU
LIBS += -L$(QGLVIEWER_DIR)/lib -lQGLViewer

fire:$(OBJS)
	g++ -o fire $(OBJS) $(LIBS)

.o:.cpp
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o fire

