#include <QtGui/QApplication>
#include "gui.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    smile::CSmileGUI w;
    w.show();
    a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));
    return a.exec();
}
