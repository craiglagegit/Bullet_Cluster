#include <QCoreApplication>
#include "core.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    smile::CSmileCore core;
    a.connect(&core, SIGNAL(scriptDone()), SLOT(quit()), Qt::QueuedConnection);
    core.runScript(argv[1]);
    return a.exec();
}
