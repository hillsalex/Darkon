#include <QtGui/QApplication>
#include <iostream>
#include "qgl.h"
#include <pty.h>
#include "mainwindow.h"
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
    //Q_INIT_RESOURCE(resources);

    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
