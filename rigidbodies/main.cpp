#include "mainwindow.h"
#include "controller.h"
#include "simparameters.h"
#include <QApplication>
#include <QMetaType>

int main(int argc, char *argv[])
{
    qRegisterMetaType<SimParameters>("SimParameters");

    srand(0);

    int fps = 20;
    int simfps = 200;

    QApplication a(argc, argv);       
    Controller controller(simfps);
    MainWindow w(controller, fps);
    controller.initialize(&w);
    controller.start();
    w.show();

    int ret = a.exec();
    controller.quit();
    return ret;
}
