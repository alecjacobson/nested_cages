#include "glpanel.h"
#include "controller.h"
#include <QMouseEvent>

using namespace Eigen;

GLPanel::GLPanel(QWidget *parent) :
    QGLWidget(parent), c_(), translateDir_(0), rotator_(c_)
{
    cont_ = NULL;
    c_.setDefault3D();
    lightPos_.setZero();
    lightPos_[2] = 1.0;
}

void GLPanel::setController(Controller *cont)
{
    cont_ = cont;
}

void GLPanel::initializeGL()
{
    assert(cont_);
    glShadeModel(GL_SMOOTH);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClearDepth(1.0);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    cont_->initializeGL();
}

void GLPanel::resizeGL(int w, int h)
{
    c_.setPerpective(60.0, 1.0);
    c_.setViewport(w, h);
}

void GLPanel::paintGL()
{
    assert(cont_);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3f (0.0, 0.0, 0.0);

    c_.applyViewport();
    c_.applyProjection();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    c_.applyLookAt();

    GLfloat lightColor0[] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat lightPosition[4] = { lightPos_[0], lightPos_[1], lightPos_[2], lightPos_[3] };
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    glDisable(GL_LIGHTING);
    cont_->renderPlanes(false);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    cont_->renderObjects();

    glPushMatrix();
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(-1.0, -1.0);
    glDisable(GL_LIGHT0);
    multShadowMatrix();
    cont_->renderObjects();
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPopMatrix();

    cont_->renderPlanes(true);
}

void GLPanel::multShadowMatrix()
{
    float light[4];
    float ground[4];

    for(int i=0; i<4; i++)
    {
        light[i] = lightPos_[i];
        ground[i] = 0;
    }
    ground[2] = 1.0;

    float  dot;
    float  shadowMat[4][4];

    dot = ground[0] * light[0] +
          ground[1] * light[1] +
          ground[2] * light[2] +
          ground[3] * light[3];

    shadowMat[0][0] = dot - light[0] * ground[0];
    shadowMat[1][0] = 0.0 - light[0] * ground[1];
    shadowMat[2][0] = 0.0 - light[0] * ground[2];
    shadowMat[3][0] = 0.0 - light[0] * ground[3];

    shadowMat[0][1] = 0.0 - light[1] * ground[0];
    shadowMat[1][1] = dot - light[1] * ground[1];
    shadowMat[2][1] = 0.0 - light[1] * ground[2];
    shadowMat[3][1] = 0.0 - light[1] * ground[3];

    shadowMat[0][2] = 0.0 - light[2] * ground[0];
    shadowMat[1][2] = 0.0 - light[2] * ground[1];
    shadowMat[2][2] = dot - light[2] * ground[2];
    shadowMat[3][2] = 0.0 - light[2] * ground[3];

    shadowMat[0][3] = 0.0 - light[3] * ground[0];
    shadowMat[1][3] = 0.0 - light[3] * ground[1];
    shadowMat[2][3] = 0.0 - light[3] * ground[2];
    shadowMat[3][3] = dot - light[3] * ground[3];

    glMultMatrixf((const GLfloat*)shadowMat);
}

void GLPanel::scaleMousePos(int x, int y, double &scaledx, double &scaledy) const
{
    int w, h;
    c_.getViewport(w,h);
    scaledx = 2.0 * x / double(w-1) - 1.0;
    scaledy = 2.0 * (h - y -1) / double(h-1) - 1.0;
}

GLPanel::MouseAction GLPanel::deduceAction(QMouseEvent *event)
{
    if(event->buttons() & Qt::LeftButton)
        return MA_LAUNCH;
    else if(event->buttons() & Qt::RightButton)
        return MA_ROTATE;
    return MA_NONE;
}

void GLPanel::mousePressEvent(QMouseEvent *event)
{
    int x = event->pos().x();
    int y = event->pos().y();
    Vector2d pos;
    scaleMousePos(x,y,pos[0],pos[1]);

    MouseAction ma = deduceAction(event);
    switch(ma)
    {
    case MA_ROTATE:
    {
        rotator_.startRotation(pos);
        break;
    }
    case MA_LAUNCH:
    {
        Vector3d pos = c_.getEye();
        Vector3d right, up, center;
        c_.getSpanningSet(right, up, center);
        QMetaObject::invokeMethod(cont_, "mouseClicked", Q_ARG(double, pos[0]), Q_ARG(double, pos[1]), Q_ARG(double, pos[2]), Q_ARG(double, center[0]), Q_ARG(double, center[1]), Q_ARG(double, center[2]));
        break;
    }
    default:
        break;
    }
}

void GLPanel::mouseMoveEvent(QMouseEvent *event)
{
    int x = event->pos().x();
    int y = event->pos().y();
    Vector2d pos;
    scaleMousePos(x,y,pos[0],pos[1]);
    rotator_.updateRotation(pos);
}

void GLPanel::mouseReleaseEvent(QMouseEvent *event)
{
    MouseAction ma = deduceAction(event);
    if(ma != MA_ROTATE)
    {
        int x = event->pos().x();
        int y = event->pos().y();
        Vector2d pos;
        scaleMousePos(x,y,pos[0],pos[1]);
        rotator_.stopRotation();
    }
}

void GLPanel::keyPressEvent(QKeyEvent *ke)
{
    switch(ke->key())
    {
    case 'W':
        translateDir_ |= TD_FWD;
        break;
    case 'A':
        translateDir_ |= TD_LEFT;
        break;
    case 'S':
        translateDir_ |= TD_BACK;
        break;
    case 'D':
        translateDir_ |= TD_RIGHT;
        break;
    }
}

void GLPanel::tick()
{

    Vector3d right, up, center;
    c_.getSpanningSet(right, up, center);
    center[2] = 0;
    center.normalize();
    right[2] = 0;
    right.normalize();
    int fwdamt = 0;
    if(translateDir_ & TD_FWD)
        fwdamt++;
    if(translateDir_ & TD_BACK)
        fwdamt--;
    Vector3d dir = fwdamt*center;
    int rghtamt = 0;
    if(translateDir_ & TD_LEFT)
        rghtamt--;
    if(translateDir_ & TD_RIGHT)
        rghtamt++;
    dir += rghtamt*right;
    double dnorm = dir.norm();
    if(dnorm != 0)
        dir /= dnorm;
    c_.translateCenter(dir);
    c_.translateEye(dir);
}

void GLPanel::keyReleaseEvent(QKeyEvent *ke)
{
    switch(ke->key())
    {
    case 'W':
        translateDir_ &= ~TD_FWD;
        break;
    case 'A':
        translateDir_ &= ~TD_LEFT;
        break;
    case 'S':
        translateDir_ &= ~TD_BACK;
        break;
    case 'D':
        translateDir_ &= ~TD_RIGHT;
        break;
    }
}
