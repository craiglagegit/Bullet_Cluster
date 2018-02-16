
TEMPLATE = app
INCLUDEPATH += . /home/eugvas/libs/glpk/include /home/eugvas/libs/qwt-5.2.1/include /home/eugvas/libs/qwtplot3d /home/eugvas/libs/gsl/include 
QMAKE_LIBDIR += /home/eugvas/libs/glpk/lib /home/eugvas/libs/qwt-5.2.1/lib /home/eugvas/libs/qwtplot3d /home/eugvas/libs/gsl/lib
# Input
HEADERS += common.h \
           core.h \
           fileio.h \
           gui.h \
           massmodel.h \
           optimization.h \
           orbit.h \
           orbitlib.h \
           potential.h \
           schwarzschild.h \
           stringconv.h 
FORMS += smile.ui
SOURCES += common.cpp \
           core.cpp \
           fileio.cpp \
           gui.cpp \
           main.cpp \
           massmodel.cpp \
           optimization.cpp \
           orbit.cpp \
           orbitlib.cpp \
           potential.cpp \
           schwarzschild.cpp 
RESOURCES = smile.qrc
QT += opengl
LIBS += -static-libgcc -Wl,-Bstatic -lgsl -lgslcblas -lqwt -lqwtplot3d -lglpk -Wl,-Bdynamic
macx:LIBS += -lz
macx:ICON = smile.icns
