/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Mon Jan 19 00:36:38 2015
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      22,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      19,   12,   11,   11, 0x0a,
      54,   11,   11,   11, 0x08,
      65,   11,   11,   11, 0x08,
      91,   11,   11,   11, 0x08,
     129,   11,   11,   11, 0x08,
     156,   11,   11,   11, 0x08,
     191,   11,   11,   11, 0x08,
     225,   11,   11,   11, 0x08,
     260,   11,   11,   11, 0x08,
     300,   11,   11,   11, 0x08,
     329,   11,   11,   11, 0x08,
     363,   11,   11,   11, 0x08,
     389,   11,   11,   11, 0x08,
     418,   11,   11,   11, 0x08,
     443,   11,   11,   11, 0x08,
     469,   11,   11,   11, 0x08,
     504,   11,   11,   11, 0x08,
     544,   11,   11,   11, 0x08,
     582,   11,   11,   11, 0x08,
     620,   11,   11,   11, 0x08,
     653,   11,   11,   11, 0x08,
     695,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0params\0"
    "setUIFromParameters(SimParameters)\0"
    "updateGL()\0on_actionExit_triggered()\0"
    "on_actionReset_Everything_triggered()\0"
    "on_actionReset_triggered()\0"
    "on_startSimulationButton_clicked()\0"
    "on_timeStepEdit_editingFinished()\0"
    "on_newtonTolEdit_editingFinished()\0"
    "on_newtonMaxItersEdit_editingFinished()\0"
    "on_gravityCheckBox_clicked()\0"
    "on_gravityGEdit_editingFinished()\0"
    "on_sphereButton_clicked()\0"
    "on_twoByFourButton_clicked()\0"
    "on_bunnyButton_clicked()\0"
    "on_customButton_clicked()\0"
    "on_launchVelEdit_editingFinished()\0"
    "on_randomOrienatationCheckBox_clicked()\0"
    "on_randomAngularVelCheckBox_clicked()\0"
    "on_randomVelMagEdit_editingFinished()\0"
    "on_densityEdit_editingFinished()\0"
    "on_penaltyStiffnessEdit_editingFinished()\0"
    "on_planeButton_clicked()\0"
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWindow *_t = static_cast<MainWindow *>(_o);
        switch (_id) {
        case 0: _t->setUIFromParameters((*reinterpret_cast< const SimParameters(*)>(_a[1]))); break;
        case 1: _t->updateGL(); break;
        case 2: _t->on_actionExit_triggered(); break;
        case 3: _t->on_actionReset_Everything_triggered(); break;
        case 4: _t->on_actionReset_triggered(); break;
        case 5: _t->on_startSimulationButton_clicked(); break;
        case 6: _t->on_timeStepEdit_editingFinished(); break;
        case 7: _t->on_newtonTolEdit_editingFinished(); break;
        case 8: _t->on_newtonMaxItersEdit_editingFinished(); break;
        case 9: _t->on_gravityCheckBox_clicked(); break;
        case 10: _t->on_gravityGEdit_editingFinished(); break;
        case 11: _t->on_sphereButton_clicked(); break;
        case 12: _t->on_twoByFourButton_clicked(); break;
        case 13: _t->on_bunnyButton_clicked(); break;
        case 14: _t->on_customButton_clicked(); break;
        case 15: _t->on_launchVelEdit_editingFinished(); break;
        case 16: _t->on_randomOrienatationCheckBox_clicked(); break;
        case 17: _t->on_randomAngularVelCheckBox_clicked(); break;
        case 18: _t->on_randomVelMagEdit_editingFinished(); break;
        case 19: _t->on_densityEdit_editingFinished(); break;
        case 20: _t->on_penaltyStiffnessEdit_editingFinished(); break;
        case 21: _t->on_planeButton_clicked(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 22)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 22;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
