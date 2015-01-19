#ifndef GLPANEL_H
#define GLPANEL_H

#include <QGLWidget>
#include "camera.h"
#include "rotator.h"

class Controller;

class GLPanel : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLPanel(QWidget *parent = 0);
    void setController(Controller *cont);

    void tick();

signals:

public slots:
    virtual void initializeGL();
    virtual void resizeGL(int w, int h);
    virtual void paintGL();

    virtual void mousePressEvent(QMouseEvent *me);
    virtual void mouseMoveEvent(QMouseEvent *me);
    virtual void mouseReleaseEvent(QMouseEvent *me);

    virtual void keyPressEvent(QKeyEvent *ke);
    virtual void keyReleaseEvent(QKeyEvent *ke);

private:
    enum MouseAction { MA_NONE, MA_TRANSLATE, MA_ROTATE, MA_ZOOM, MA_LAUNCH};
    enum TranslateDirection { TD_FWD=1, TD_BACK=2, TD_LEFT=4, TD_RIGHT=8};

    void scaleMousePos(int x, int y, double &scaledx, double &scaledy) const;
    MouseAction deduceAction(QMouseEvent *event);
    void multShadowMatrix();

    Controller *cont_;

    Camera c_;
    int translateDir_;
    Rotator rotator_;
    //Zoomer zoomer_;
    Eigen::Vector4d lightPos_;
};

#endif // GLPANEL_H
