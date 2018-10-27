#include <QGLViewer/qglviewer.h>

class UfabViewer : public QGLViewer {
public:
    UfabViewer(QWidget *parent = nullptr);
protected:
    virtual void draw();
    virtual void init();
    virtual QString helpString() const;
};
