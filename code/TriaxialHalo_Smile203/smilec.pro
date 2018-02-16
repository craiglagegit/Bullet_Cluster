
TEMPLATE = app
INCLUDEPATH += . /home/eugvas/libs/glpk/include /home/eugvas/libs/gsl/include 
QMAKE_LIBDIR += /home/eugvas/libs/glpk/lib /home/eugvas/libs/gsl/lib

# Input
HEADERS += common.h \
           core.h \
           fileio.h \
           massmodel.h \
           optimization.h \
           orbit.h \
           orbitlib.h \
           potential.h \
           schwarzschild.h \
           stringconv.h 
SOURCES += common.cpp \
           core.cpp \
           fileio.cpp \
           mainc.cpp \
           massmodel.cpp \
           optimization.cpp \
           orbit.cpp \
           orbitlib.cpp \
           potential.cpp \
           schwarzschild.cpp  
RC_FILE  = smile.rc
QT -= gui
CONFIG += console
LIBS += -static-libgcc -Wl,-Bstatic -lgsl -lgslcblas -lglpk -Wl,-Bdynamic
macx:LIBS += -lz
macx:ICON = smile.icns
