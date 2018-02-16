#include "gui.h"
#include "core.h"
#include "potential.h"
#include "orbit.h"
#include "orbitlib.h"
#include "schwarzschild.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include <QKeyEvent>
#include <QFileDialog>
#include <QMessageBox>
#include <QPrinter>
#include <QImage>
#include <QMap>
#include <QInputDialog>
#include <QMetaObject>

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_text.h>
#include <qwt_symbol.h>
#include <qwt_plot_picker.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_zoomer.h>
#include <qwt_scale_engine.h>
#include <qwt_plot_marker.h>

#ifdef USE3DFACETS
// for triangulation thread
#include <QVector>
#include <QProcess>
#include <QTextStream>
#endif

namespace smile{

static CSmileGUI* GUI;

static void my_message_gui(const std::string &message)
{
    std::cerr << message <<"\n";
    QMetaObject::invokeMethod(GUI, "info", Q_ARG(QString, QString::fromStdString(message)));
}

static void my_error_gui(const std::string &message)
{
    std::cerr << message <<"\n";
    QMetaObject::invokeMethod(GUI, "error", Q_ARG(QString, QString::fromStdString(message)));
}

CSmileGUI::CSmileGUI(QWidget *parent, Qt::WFlags flags)
    : QMainWindow(parent, flags)
{
    GUI=this;
    ui.setupUi(this);

    core=new CSmileCore();  // initialize the instance of "business logic" object
    my_error_ptr=&my_error_gui;
    my_message_ptr=&my_message_gui;
    // connect signals to the core
    connect(core, SIGNAL(signalOrbitFinished()), SLOT(CalcFinished()));
    connect(core, SIGNAL(signalFreqMapFinished()), SLOT(FreqMapFinished()));
    connect(core, SIGNAL(signalSchwOrbitLibraryFinished()), SLOT(SchwOrbitLibraryFinished()));
    connect(core, SIGNAL(signalSchwOptimizationFinished(const QString&)), SLOT(SchwOptimizationFinished(const QString&)));
    connect(core, SIGNAL(signalNbodyExportFinished(const QString&)), SLOT(NbodyExportFinished(const QString&)));
    connect(core, SIGNAL(signalInfo(const QString&)), SLOT(info(const QString&)));
    connect(core, SIGNAL(signalTimer()), SLOT(timerEventCore()));
    
    Eprev=0;
    FreqMapRunning=false;
    SchwModelRunning=false;
    ExportNbodyRunning=false;
    // fill potential, density, symmetry and SM types dropdown list
    ui.selectPotentialType->clear();
    for(PotentialNameMapType::iterator iterp=PotentialNames.begin(); iterp!=PotentialNames.end(); iterp++)
        ui.selectPotentialType->insertItem(100, QString::fromStdString(iterp->second));
    ui.selectDensityType->clear();
    for(DensityNameMapType::iterator iterd=DensityNames.begin(); iterd!=DensityNames.end(); iterd++)
        ui.selectDensityType->insertItem(100, QString::fromStdString(iterd->second));
    ui.selectSymmetryType->clear();
    for(SymmetryNameMapType::iterator iters=SymmetryNames.begin(); iters!=SymmetryNames.end(); iters++)
        ui.selectSymmetryType->insertItem(100, QString::fromStdString(iters->second));
    ui.select_smType->clear();
    ui.select_smType->insertItem(100, QString::fromStdString(CSchwModelClassic::myName()));
    ui.select_smType->insertItem(100, QString::fromStdString(CSchwModelSHGrid::myName()));
    ui.select_smType->insertItem(100, QString::fromStdString(CSchwModelSHBSE::myName()));
    
    // orbit plot
    PlotOrbit = new QwtPlot(QwtText("Orbit"), ui.frame_orbitarea);
    CurveOrbit= new QwtPlotCurve("");
    CurveOrbit->setPen(QColor(Qt::blue));
    CurveOrbit->attach(PlotOrbit);
    CurveEquipotential= new QwtPlotCurve("");
    CurveEquipotential->setPen(QColor(Qt::red));
    CurveEquipotential->attach(PlotOrbit);
    PlotOrbit->show();
    QVBoxLayout* LPlotOrbit = new QVBoxLayout(ui.frame_orbitarea);
    LPlotOrbit->addWidget(PlotOrbit);
    ZoomerOrbit= new QwtPlotZoomer2( QwtPlot::xBottom, QwtPlot::yLeft, 
        QwtPicker::DragSelection | QwtPicker::CornerToCorner, QwtPicker::ActiveOnly, PlotOrbit->canvas(), true);
    ZoomerOrbit->setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton);
    QwtPlotPanner* PannerOrbit = new QwtPlotPanner(PlotOrbit->canvas());
    PannerOrbit->setMouseButton(Qt::MidButton);

    PlotOrbitParam = new QwtPlot(QwtText("Distribution of orbital parameters"), ui.frame_orbitarea);
    CurveOrbitParam= new QwtPlotCurve("");
    CurveOrbitParam->attach(PlotOrbitParam);
    PlotOrbitParam->hide();
    LPlotOrbit->addWidget(PlotOrbitParam);
    ZoomerOrbitParam= new QwtPlotZoomer2( QwtPlot::xBottom, QwtPlot::yLeft, 
        QwtPicker::DragSelection | QwtPicker::CornerToCorner, QwtPicker::ActiveOnly, PlotOrbitParam->canvas(), true);
    ZoomerOrbitParam->setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton);
    QwtPlotPanner* PannerOrbitParam = new QwtPlotPanner(PlotOrbitParam->canvas());
    PannerOrbitParam->setMouseButton(Qt::MidButton);

#ifdef USEQWT3D
    // 3d orbit
    Plot3d=new Qwt3D::SurfacePlot(ui.frame_orbitarea);
    Plot3d->setCoordinateStyle(Qwt3D::FRAME);
    Plot3d->setRotation(90,0,0);
    QColor bkg=ui.tabWidget->palette().color(QPalette::Background);
    Plot3d->setBackgroundColor(Qwt3D::RGBA(bkg.redF(), bkg.greenF(), bkg.blueF()));
    Plot3d->setMeshColor(Qwt3D::RGBA(0.2,0,0.4));
    for (unsigned i=0; i!=Plot3d->coordinates()->axes.size(); ++i)
    {
      Plot3d->coordinates()->axes[i].setMajors(7);
      Plot3d->coordinates()->axes[i].setMinors(4);
    }
    LPlotOrbit->addWidget(Plot3d);
    ui.progressBarOrbit->hide();
#ifdef USE3DFACETS
    TriangThr = new CTriangThread(core->AppDir+"qdelaunay");
    connect(TriangThr, SIGNAL(finished(bool)), SLOT(redrawOrbit(bool)));
    Plot3d->setLightComponent(GL_DIFFUSE, 0.3);
    Plot3d->setMaterialComponent(GL_SPECULAR, 0.0);
    Plot3d->enableLighting();
    Plot3d->illuminate();
    Plot3d->updateGL();
#endif
#endif

    // poincare plot
    PlotPoincare = new QwtPlot(QwtText("Poincare section"), ui.frame_poincarearea);
    PlotPoincare->show();
    QVBoxLayout * LPlotPoincare = new QVBoxLayout(ui.frame_poincarearea);
    LPlotPoincare->addWidget(PlotPoincare);
    QwtPlotPicker *PickerPoincare = new QwtPlotPicker2(QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::PointSelection, QwtPlotPicker::NoRubberBand, QwtPicker::AlwaysOn, PlotPoincare->canvas());
    PickerPoincare->setMousePattern(QwtEventPattern::MouseSelect1, Qt::RightButton);
    connect(PickerPoincare, SIGNAL(selected(const QwtDoublePoint &)), SLOT(PSclicked(const QwtDoublePoint &)));
    ZoomerPoincare = new QwtPlotZoomer2( QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::DragSelection | QwtPicker::CornerToCorner, QwtPicker::AlwaysOn, PlotPoincare->canvas(), true);
    ZoomerPoincare->setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton, Qt::ControlModifier);
    ZoomerPoincare->setMousePattern(QwtEventPattern::MouseSelect2, Qt::RightButton, Qt::ControlModifier);
    ZoomerPoincare->setKeyPattern(QwtEventPattern::KeySelect1, Qt::Key_Space);
    QwtPlotPanner* PannerPoincare = new QwtPlotPanner(PlotPoincare->canvas());
    PannerPoincare->setMouseButton(Qt::MidButton);

    // frequency plot
    PlotFreq = new QwtPlot(QwtText("Frequencies"), ui.frame_frequencies);
    CurveFreqX= new QwtPlotCurve("X");
    CurveFreqX->setPen(QColor(Qt::blue));
    CurveFreqX->attach(PlotFreq);
    CurveFreqY= new QwtPlotCurve("Y");
    CurveFreqY->setPen(QColor(Qt::darkGreen));
    CurveFreqY->attach(PlotFreq);
    CurveFreqZ= new QwtPlotCurve("Z");
    CurveFreqZ->setPen(QColor(Qt::darkRed));
    CurveFreqZ->attach(PlotFreq);
    CurveFreqPeaks= new QwtPlotCurve();
    CurveFreqPeaks->setStyle(QwtPlotCurve::Sticks);
    CurveFreqPeaks->attach(PlotFreq);
    CurveFreqPeaksP= new QwtPlotCurve();
    CurveFreqPeaksP->setStyle(QwtPlotCurve::Sticks);
    CurveFreqPeaksP->setPen(QColor(Qt::white));
    CurveFreqPeaksP->attach(PlotFreq);
    PlotFreq->show();
    QVBoxLayout* LPlotFreq = new QVBoxLayout(ui.frame_frequencies);
    LPlotFreq->addWidget(PlotFreq);
    ZoomerFreq = new QwtPlotZoomer2( QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::DragSelection | QwtPicker::CornerToCorner, QwtPicker::AlwaysOn, PlotFreq->canvas(), true);
    ZoomerFreq->setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton);
    QwtPlotPanner* PannerFreq = new QwtPlotPanner(PlotFreq->canvas());
    PannerFreq->setMouseButton(Qt::MidButton);

    // Lyapunov plot
    PlotLyapunov = new QwtPlot(QwtText("Lyapunov number"), ui.tab_l);
    PlotLyapunov->setAxisScaleEngine(QwtPlot::xBottom, new QwtLog10ScaleEngine);
    PlotLyapunov->setAxisScaleEngine(QwtPlot::yLeft, new QwtLog10ScaleEngine);
    PlotLyapunov->setAxisScaleEngine(QwtPlot::yRight, new QwtLog10ScaleEngine);
    CurveDeviation = new QwtPlotCurve("");
    CurveDeviation->setPen(QColor(Qt::red));
    CurveDeviation->attach(PlotLyapunov);
    CurveDeviation->setYAxis(QwtPlot::yRight);
    PlotLyapunov->enableAxis(QwtPlot::yRight);
    CurveLyapunov = new QwtPlotCurve("");
    CurveLyapunov->setPen(QColor(Qt::blue));
    CurveLyapunov->attach(PlotLyapunov);
    PlotLyapunov->show();
    /*QwtPlotPicker *PickerLyapunov = */new QwtPlotPicker2(QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::PointSelection, QwtPlotPicker::NoRubberBand, QwtPicker::AlwaysOn, PlotLyapunov->canvas());
    QVBoxLayout* LPlotLyapunov = new QVBoxLayout(ui.tab_l);
    LPlotLyapunov->addWidget(PlotLyapunov);

    // Frequency Map plot
    PlotFreqMap = new QwtPlot(QwtText("Frequency map"), ui.frame_fmarea);
    QwtSymbol sym;
    sym.setStyle(QwtSymbol::Rect);
    sym.setSize(2);
    sym.setPen(Qt::NoPen);
    sym.setBrush(QColor(Qt::red));
    CurveFreqMapC = new QwtPlotCurve("");
    CurveFreqMapC->setSymbol(sym);
    CurveFreqMapC->setStyle(QwtPlotCurve::Dots);
    CurveFreqMapC->attach(PlotFreqMap);
    sym.setBrush(QColor(Qt::blue));
    CurveFreqMapR = new QwtPlotCurve("");
    CurveFreqMapR->setSymbol(sym);
    CurveFreqMapR->setStyle(QwtPlotCurve::Dots);
    CurveFreqMapR->attach(PlotFreqMap);
    drawFreqMapResonances();
    PlotFreqMap->show();
    QVBoxLayout* LPlotFreqMap = new QVBoxLayout(ui.frame_fmarea);
    LPlotFreqMap->addWidget(PlotFreqMap);
    QwtPlotPicker *PickerFreqMap = new QwtPlotPicker2(QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::PointSelection, QwtPlotPicker::NoRubberBand, QwtPicker::AlwaysOn, PlotFreqMap->canvas());
    PickerFreqMap->setMousePattern(QwtEventPattern::MouseSelect1, Qt::RightButton);
    connect(PickerFreqMap, SIGNAL(selected(const QwtDoublePoint &)), SLOT(FMclicked(const QwtDoublePoint &)));
    ZoomerFreqMap = new QwtPlotZoomer2( QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::DragSelection | QwtPicker::CornerToCorner, QwtPicker::AlwaysOn, PlotFreqMap->canvas(), true);
    ZoomerFreqMap->setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton, Qt::ControlModifier);
    ZoomerFreqMap->setMousePattern(QwtEventPattern::MouseSelect2, Qt::RightButton, Qt::ControlModifier);
    ZoomerFreqMap->setKeyPattern(QwtEventPattern::KeySelect1, Qt::Key_Space);
    QwtPlotPanner* PannerFreqMap = new QwtPlotPanner(PlotFreqMap->canvas());
    PannerFreqMap->setMouseButton(Qt::MidButton);

    PlotFreqMapHisto = new QwtPlot(QwtText("Histograms"), ui.frame_fmarea);
    PlotFreqMapHisto->setAxisScaleEngine(QwtPlot::xBottom, new QwtLog10ScaleEngine);
    CurveFreqMapHisto = new QwtPlotCurve("");
    CurveFreqMapHisto->setStyle(QwtPlotCurve::Steps);
    CurveFreqMapHisto->attach(PlotFreqMapHisto);
    LPlotFreqMap->addWidget(PlotFreqMapHisto);
    /*QwtPlotPicker *PickerFreqMapHisto = */new QwtPlotPicker2(QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::PointSelection, QwtPlotPicker::NoRubberBand, QwtPicker::AlwaysOn, PlotFreqMapHisto->canvas());
    PlotFreqMapHisto->hide();

    PlotFreqMapSS = new QwtPlot(QwtText("Start-space"), ui.frame_fmarea);
    sym.setStyle(QwtSymbol::Rect);
    sym.setSize(4);
    sym.setBrush(QColor(Qt::blue));
    CurveFreqMapSSR = new QwtPlotCurve("");
    CurveFreqMapSSR->setSymbol(sym);
    CurveFreqMapSSR->setStyle(QwtPlotCurve::Dots);
    CurveFreqMapSSR->attach(PlotFreqMapSS);
    sym.setBrush(QColor(Qt::red));
    CurveFreqMapSSC = new QwtPlotCurve("");
    CurveFreqMapSSC->setSymbol(sym);
    CurveFreqMapSSC->setStyle(QwtPlotCurve::Dots);
    CurveFreqMapSSC->attach(PlotFreqMapSS);
    QwtPlotPicker *PickerFreqMapSS = new QwtPlotPicker2(QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::PointSelection, QwtPlotPicker::NoRubberBand, QwtPicker::ActiveOnly, PlotFreqMapSS->canvas());
    PickerFreqMapSS->setMousePattern(QwtEventPattern::MouseSelect1, Qt::RightButton);
    connect(PickerFreqMapSS, SIGNAL(selected(const QwtDoublePoint &)), SLOT(FMclicked(const QwtDoublePoint &)));
    LPlotFreqMap->addWidget(PlotFreqMapSS);
    PlotFreqMapSS->hide();
    ui.progressBarFM->hide();

    // Schwarzschild plot
    QVBoxLayout* LPlotSchw = new QVBoxLayout(ui.frame_smarea);
    PlotSchw = new QwtPlot(QwtText("Schwarzschild model grid"), ui.frame_smarea);
    sym.setStyle(QwtSymbol::Diamond);
    sym.setSize(8);
    sym.setBrush(QColor(Qt::blue));
    CurveSchwCellF = new QwtPlotCurve("");
    CurveSchwCellF->setSymbol(sym);
    CurveSchwCellF->setStyle(QwtPlotCurve::Dots);
    CurveSchwCellF->attach(PlotSchw);
    sym.setBrush(QColor(Qt::red));
    CurveSchwCellI= new QwtPlotCurve("");
    CurveSchwCellI->setSymbol(sym);
    CurveSchwCellI->setStyle(QwtPlotCurve::Dots);
    CurveSchwCellI->attach(PlotSchw);
    QwtPlotPicker *PickerFreqMapSM = new QwtPlotPicker2(QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::PointSelection, QwtPlotPicker::NoRubberBand, QwtPicker::ActiveOnly, PlotSchw->canvas());
    PickerFreqMapSM->setMousePattern(QwtEventPattern::MouseSelect1, Qt::RightButton);
    connect(PickerFreqMapSM, SIGNAL(selected(const QwtDoublePoint &)), SLOT(SMclicked(const QwtDoublePoint &)));
    PlotSchw->show();
    LPlotSchw->addWidget(PlotSchw);
    ui.progressBarSM->hide();
    PlotSchwHisto = new QwtPlot(QwtText("Schwarzschild model orbit weights"), ui.frame_smarea);
    sym.setSize(2);
    sym.setStyle(QwtSymbol::Rect);
    sym.setBrush(QColor(Qt::blue));
    CurveSchwHistoR = new QwtPlotCurve("");
    CurveSchwHistoR->setStyle(QwtPlotCurve::Dots);
    CurveSchwHistoR->setSymbol(sym);
    CurveSchwHistoR->attach(PlotSchwHisto);
    sym.setBrush(QColor(Qt::red));
    CurveSchwHistoC = new QwtPlotCurve("");
    CurveSchwHistoC->setStyle(QwtPlotCurve::Dots);
    CurveSchwHistoC->setSymbol(sym);
    CurveSchwHistoC->attach(PlotSchwHisto);
    LPlotSchw->addWidget(PlotSchwHisto);
    QwtPlotPicker *PickerSchwHisto = new QwtPlotPicker2(QwtPlot::xBottom, QwtPlot::yLeft,
        QwtPicker::PointSelection, QwtPlotPicker::NoRubberBand, QwtPicker::AlwaysOn, PlotSchwHisto->canvas());
    PickerSchwHisto->setMousePattern(QwtEventPattern::MouseSelect1, Qt::RightButton);
    connect(PickerSchwHisto, SIGNAL(selected(const QwtDoublePoint &)), SLOT(SMHistoClicked(const QwtDoublePoint &)));
    PlotSchwHisto->hide();

    // connect various events to change corresponding parameters in core
    connect(ui.lineEdit_x, SIGNAL(textChanged(const QString&)), SLOT(ICchanged()));
    connect(ui.lineEdit_y, SIGNAL(textChanged(const QString&)), SLOT(ICchanged()));
    connect(ui.lineEdit_z, SIGNAL(textChanged(const QString&)), SLOT(ICchanged()));
    connect(ui.lineEdit_vx, SIGNAL(textChanged(const QString&)), SLOT(ICchanged()));
    connect(ui.lineEdit_vy, SIGNAL(textChanged(const QString&)), SLOT(ICchanged()));
    connect(ui.lineEdit_vz, SIGNAL(textChanged(const QString&)), SLOT(ICchanged()));
    connect(ui.lineEdit_E, SIGNAL(textChanged(const QString&)), SLOT(ICchanged()));
    connect(ui.lineEdit_q, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_p, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_Mbh, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_Rc, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_Gamma, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_Omega, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_Ncoefs_radial, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_Ncoefs_angular, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_Alpha, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_eps, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.lineEdit_tol, SIGNAL(textChanged(const QString&)), SLOT(potentialChanged()));
    connect(ui.selectPotentialType, SIGNAL(currentIndexChanged(int)), SLOT(potentialChanged()));
    connect(ui.selectDensityType, SIGNAL(currentIndexChanged(int)), SLOT(potentialChanged()));
    connect(ui.selectSymmetryType, SIGNAL(currentIndexChanged(int)), SLOT(potentialChanged()));
    connect(ui.buttonChooseNbodyFile, SIGNAL(clicked()), SLOT(chooseNbodyFile()));
    connect(ui.lineEdit_NbodyFile, SIGNAL(returnPressed()), SLOT(typeinNbodyFile()));
    connect(ui.radioButton_2dorbit, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_3dline, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_3dmesh, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_opPericenter, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_opLPeri, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_opL, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_otP, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_otL, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_otE, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_2d_xy, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_2d_xz, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_2d_yz, SIGNAL(toggled(bool)), SLOT(redrawOrbit(bool)));
    connect(ui.radioButton_fmmapdiff, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmmaplyapunov, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmhistodiff, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmhistolyapunov, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmssdiff, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmsslyapunov, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmspdiff, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmsplyapunov, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmsydiff, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_fmsylyapunov, SIGNAL(toggled(bool)), SLOT(redrawFreqMap(bool)));
    connect(ui.radioButton_smviewgrid, SIGNAL(toggled(bool)), SLOT(redrawSchwModel(bool)));
    connect(ui.radioButton_smviewweights, SIGNAL(toggled(bool)), SLOT(redrawSchwModel(bool)));
    connect(ui.radioButton_smviewhisto, SIGNAL(toggled(bool)), SLOT(redrawSchwModel(bool)));
    connect(ui.lineEdit_intTime, SIGNAL(textChanged(const QString&)), SLOT(intTime_changed()));
    connect(ui.lineEdit_timeStep,SIGNAL(textChanged(const QString&)), SLOT(intTime_changed()));
    connect(ui.lineEdit_fms, SIGNAL(textChanged(const QString&)), SLOT(fmparams_changed()));
    connect(ui.lineEdit_fmp, SIGNAL(textChanged(const QString&)), SLOT(fmparams_changed()));
    connect(ui.lineEdit_fmy, SIGNAL(textChanged(const QString&)), SLOT(fmparams_changed()));
    connect(ui.lineEdit_fmr, SIGNAL(textChanged(const QString&)), SLOT(fmparams_changed()));
    connect(ui.lineEdit_smNumOrbits, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.select_smType, SIGNAL(currentIndexChanged(int)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smNumLines, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smNumAngCoefs, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smNumRadCoefs, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smAlpha, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smNumShells, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smInnerShellMass, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smOuterShellMass, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smNumSamplingPoints, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smChaoticWeight, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.checkBox_smBeta, SIGNAL(toggled(bool)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smBetaIn, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smBetaOut, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_smMaxOrbitWeight, SIGNAL(textChanged(const QString&)), SLOT(smparams_changed()));
    connect(ui.lineEdit_MinFreqDiff, SIGNAL(textChanged(const QString&)), SLOT(ccparams_changed()));
    connect(ui.lineEdit_MinLambda, SIGNAL(textChanged(const QString&)), SLOT(ccparams_changed()));
    ui.lineEdit_x->installEventFilter(this);
    ui.lineEdit_y->installEventFilter(this);
    ui.lineEdit_z->installEventFilter(this);
    ui.lineEdit_vx->installEventFilter(this);
    ui.lineEdit_vy->installEventFilter(this);
    ui.lineEdit_vz->installEventFilter(this);
    ui.lineEdit_q->installEventFilter(this);
    ui.lineEdit_p->installEventFilter(this);
    ui.lineEdit_Mbh->installEventFilter(this);
    ui.lineEdit_Rc->installEventFilter(this);
    ui.lineEdit_intTime->installEventFilter(this);
    ui.lineEdit_timeStep->installEventFilter(this);
    ui.checkBox_Lyapunov->installEventFilter(this);
    ui.radioButton_2d->installEventFilter(this);
    ui.radioButton_3d->installEventFilter(this);
    ui.radioButton_2dorbit->installEventFilter(this);
    ui.radioButton_3dline->installEventFilter(this);
    ui.radioButton_3dmesh->installEventFilter(this);
    ui.ButtonClearPS->installEventFilter(this);
    ui.ButtonRandom->installEventFilter(this);
    ui.tab_f->installEventFilter(this);
    ui.tab_p->installEventFilter(this);
    ui.tab_o->installEventFilter(this);
    ui.tabWidget->installEventFilter(this);
    PlotOrbit->installEventFilter(this);
    PlotPoincare->installEventFilter(this);
    PlotFreq->installEventFilter(this);
    PlotFreqMap->installEventFilter(this);

    // load and apply settings
    loadApplySettings();
}

CSmileGUI::~CSmileGUI()
{
    if(ui.checkBox_saveSettings->isChecked())
        core->saveSettings();
    delete core;
}

bool CSmileGUI::loadApplySettings(QString fileName)
{
    ///nonpublic
    ui.label_Omega->setVisible(false);  // !!! no full support for rotating potential yet...
    ui.lineEdit_Omega->setVisible(false);

    bool result=core->loadSettings(fileName);
    InternalRecalc=true;
    // find potential name by its value
    ui.selectPotentialType->setCurrentIndex(ui.selectPotentialType->findText(QString::fromStdString(PotentialNames[configPotential.PotentialType])));
    ui.selectDensityType  ->setCurrentIndex(ui.selectDensityType  ->findText(QString::fromStdString(DensityNames  [configPotential.DensityType])));
    ui.selectSymmetryType ->setCurrentIndex(ui.selectSymmetryType ->findText(QString::fromStdString(SymmetryNames [configPotential.SymmetryType])));
    ui.select_smType->setCurrentIndex(ui.select_smType->findText(QString::fromStdString(getSchwModelNameByType(configSchw.ModelType))));
    ui.lineEdit_q->setText(QString::number(configPotential.q));
    ui.lineEdit_p->setText(QString::number(configPotential.p));
    ui.lineEdit_Mbh->setText(QString::number(configPotential.Mbh));
    ui.lineEdit_Rc->setText(QString::number(configPotential.Rc));
    ui.lineEdit_Gamma->setText(QString::number(configPotential.Gamma));
    ui.lineEdit_Omega->setText(QString::number(configPotential.Omega));
    ui.lineEdit_Ncoefs_radial->setText(QString::number(configPotential.Ncoefs_radial));
    ui.lineEdit_Ncoefs_angular->setText(QString::number(configPotential.Ncoefs_angular));
    ui.lineEdit_Alpha->setText(QString::number(configPotential.Alpha));
    ui.lineEdit_NbodyFile->setText(core->NbodyFile);
    ui.lineEdit_eps->setText(QString::number(configPotential.nbpotEps));
    ui.lineEdit_tol->setText(QString::number(configPotential.nbpotTol));
    ui.lineEdit_x->setText(QString::number(core->InitData.initCond.Pos[0]));
    ui.lineEdit_y->setText(QString::number(core->InitData.initCond.Pos[1]));
    ui.lineEdit_z->setText(QString::number(core->InitData.initCond.Pos[2]));
    ui.lineEdit_vx->setText(QString::number(core->InitData.initCond.Vel[0]));
    ui.lineEdit_vy->setText(QString::number(core->InitData.initCond.Vel[1]));
    ui.lineEdit_vz->setText(QString::number(core->InitData.initCond.Vel[2]));
    ui.lineEdit_E->setText(QString::number(core->initCondE));
    if(core->useICe)
    {
        ui.lineEdit_E->setEnabled(true);
        ui.groupBox_ICx->setEnabled(false);
    }
    else
    {
        ui.lineEdit_E->setEnabled(false);
        ui.groupBox_ICx->setEnabled(true);
    }
    ui.lineEdit_intTime->setText(QString::number(configCore.intTimeInPeriods));
    ui.lineEdit_timeStep->setText(QString::number(configCore.intTimeStepsPerPeriod));
    ui.checkBox_Lyapunov->setChecked(core->InitData.calcLyapunov);
    ui.checkBox_usePS->setChecked(configCore.usePS);
    ui.lineEdit_fms->setText(QString::number(configCore.fm_numOrbitsStationary));
    ui.lineEdit_fmp->setText(QString::number(configCore.fm_numOrbitsPrincipalPlane));
    ui.lineEdit_fmy->setText(QString::number(configCore.fm_numOrbitsYalpha));
    ui.lineEdit_fmr->setText(QString::number(configCore.fm_numOrbitsRandom));
    ui.lineEdit_smNumOrbits->setText(QString::number(configCore.sm_numOrbitsRandom));
    ui.lineEdit_smNumLines->setText(QString::number(configCore.sm_linesPerSegment));
    ui.lineEdit_smNumAngCoefs->setText(QString::number(configCore.sm_numAngularCoefs));
    ui.lineEdit_smNumRadCoefs->setText(QString::number(configCore.sm_numRadialCoefs));
    ui.lineEdit_smAlpha->setText(QString::number(configCore.sm_Alpha));
    ui.lineEdit_smNumShells->setText(QString::number(configCore.sm_numShells));
    ui.lineEdit_smInnerShellMass->setText(QString::number(configCore.sm_innerShellMass));
    ui.lineEdit_smOuterShellMass->setText(QString::number(configCore.sm_outerShellMass));
    ui.lineEdit_smNumSamplingPoints->setText(QString::number(configCore.sm_numSamplingPoints));
    ui.selectEnergyLevel->setMaximum(configCore.sm_numShells);
    ui.lineEdit_MinFreqDiff->setText(QString::number(configCore.chaoticMinFreqDiff));
    ui.lineEdit_MinLambda->setText(QString::number(configCore.chaoticMinLambda));
    ui.lineEdit_smChaoticWeight->setText(QString::number(configCore.chaoticWeightFactor));
    ui.lineEdit_smMaxOrbitWeight->setText(QString::number(configSchw.sm_maxWeight));
    ui.lineEdit_smBetaIn->setText(QString::number(configSchw.sm_betaIn));
    ui.lineEdit_smBetaOut->setText(QString::number(configSchw.sm_betaOut));
    ui.checkBox_smBeta->setChecked(configSchw.sm_constrainBeta);
    ui.lineEdit_smBetaIn->setEnabled(configSchw.sm_constrainBeta);
    ui.lineEdit_smBetaOut->setEnabled(configSchw.sm_constrainBeta);
    ui.ButtonSchwQP->setEnabled(configCore.sm_useBPMPD);
    if(core->useICe) ui.radioButton_ICe->setChecked(true); else ui.radioButton_ICe->setChecked(false);
    if(core->potential->N_dim==2) ui.radioButton_2d->setChecked(true); else ui.radioButton_3d->setChecked(true);
    InternalRecalc=false;
    potentialChanged();
    ui.buttonInitPotential->hide();   // this button will be shown only when parameters changed
    smparams_changed();
    ICchanged();
    return result;
}

bool CSmileGUI::eventFilter(QObject* target, QEvent* Event)
{
    if(Event->type() == QEvent::KeyPress) 
    {  
        QKeyEvent* keyEvent = static_cast<QKeyEvent*>(Event);
        if((keyEvent->key() == Qt::Key_Return) || (keyEvent->key() == Qt::Key_Enter))
        {
            on_ButtonStart_clicked();
            return true;
        }
    }
    return QMainWindow::eventFilter(target, Event);
}

void CSmileGUI::timerEventCore()
{
    if(FreqMapRunning)
    {
        ui.progressBarFM->setValue(core->orbitlib->numComplete());
        redrawFreqMap(true);
    }
    if(SchwModelRunning || ExportNbodyRunning)
        ui.progressBarSM->setValue(core->orbitlib->numComplete());
}

void CSmileGUI::timerEvent(QTimerEvent * Event)
{
    if(core->CalcThr!=NULL && Event->timerId()==myTimerOrbit)
    {
        ui.progressBarOrbit->setValue(core->orbit->getIntTime()*100 / (configCore.intTimeInPeriods*core->InitData.timeUnit));
        ui.progressBarOrbit->show();
        redrawOrbit(true);
    }
}

void CSmileGUI::info(const QString &message)
{
    ui.label_info->setPlainText(message);
    ui.label_info->repaint();  // to ensure immediate displaying
}

void CSmileGUI::error(const QString &message)
{
    QMessageBox::warning(this, "Error", message);
}

/// ----------------------------  single orbit integration --------- ///
void CSmileGUI::on_ButtonStart_clicked()
{
    if(core->CalcThr!=NULL) { core->orbit->halt(); return; }  // stop current computation
    if(core->InitData.timeUnit<=0 || configCore.intTimeInPeriods<=0 || configCore.intTimeStepsPerPeriod<=0) return;  // error
    ui.groupBox_N_dim->setEnabled(false);
    ui.groupBox_potential->setEnabled(false);
    ui.ButtonStart->setText("Stop!");
    ui.ButtonRandom->setEnabled(false);
    core->startOrbit();
    timeOI.start();
    myTimerOrbit = startTimer(2000);  // to refresh orbit plot
}

void CSmileGUI::CalcFinished()    // called from core->coreOrbitFinished()
{
    ui.groupBox_N_dim->setEnabled(true);
    ui.groupBox_potential->setEnabled(true);
    ui.ButtonStart->setText("Start");
    ui.ButtonRandom->setEnabled(true);
    timeOIelapsed=timeOI.elapsed();  // store to re-use in subsequent "refresh spectrum" events
    killTimer(myTimerOrbit);
    ui.progressBarOrbit->hide();
    const COrbitRuntimeTrajectory* trj = static_cast<const COrbitRuntimeTrajectory*>(core->orbit->getRuntimeFnc(0));   // 0th runtime function is the orbit trajectory analyzer
    assert(trj->FncType()==COrbitRuntimeTrajectory::FT_TRAJ_ANALYSIS);
    const COrbitInitData<double> InitData=core->orbit->getInitData();

    // Poincare section display (if it was used)
    const CPoincareInformation<float>* tps=NULL;
    for(size_t rf=0; rf<core->orbit->getInfoNum(); rf++)
        if(core->orbit->getRuntimeFnc(rf)->FncType() == CBasicOrbitRuntimeFnc::FT_POINCARE)
            tps = static_cast<const CPoincareInformation<float>* >(core->orbit->getNewInfo(rf));
    if(tps!=NULL && tps->size()>0)
    {
        QwtPlotCurve* CurvePS = new QwtPlotCurve();
        QwtSymbol sym;
        sym.setStyle(QwtSymbol::Rect);
        sym.setPen(Qt::NoPen);
        sym.setBrush(QColor(rand()%224, rand()%192, rand()%224));
        sym.setSize(2);
        CurvePS->setSymbol(sym);
        CurvePS->setStyle(QwtPlotCurve::Dots);
        CurvePS->attach(PlotPoincare);
        unsigned int nPS=tps->size();
        vectord psx(nPS), psvx(nPS);
        for(unsigned int i=0; i<nPS; i++)
        {
            psx[i]=tps->at(i).first;
            psvx[i]=tps->at(i).second;
        }
        // remove closest
        bool actiontaken;
        if(psx.size()>100)
        {
            double maxdistx=0, maxdistvx=0;
            for(unsigned int k1=1; k1<psx.size(); k1++)
            {
                double x1=psx[k1], y1=psvx[k1];
                for(unsigned int k2=0; k2<k1; k2++)
                {
                    double d=fabs(psx[k2]-x1);
                    if(d>maxdistx)
                        maxdistx=d;
                    d=fabs(psvx[k2]-y1);
                    if(d>maxdistvx) maxdistvx=d;
                }
            }
          do{
            double mindist=1000;
            unsigned int minp1=0, minp2=0;
            for(unsigned int k1=1; k1<psx.size(); k1++)
            {
                double x1=psx[k1], y1=psvx[k1];
                for(unsigned int k2=0; k2<k1; k2++)
                {
                    double d=fabs(psx[k2]-x1)/maxdistx + fabs(psvx[k2]-y1)/maxdistvx;
                    if(d<mindist)
                    {
                        mindist=d;
                        minp1=k1;
                        minp2=k2;
                    }
                }
            }
            if(mindist<1e-3)
            {
                psx.erase(psx.begin()+minp1);
                psvx.erase(psvx.begin()+minp1);
                actiontaken=true;
            }
            else
                actiontaken=false;
          }while(actiontaken);
        }
        CurvePS->setData(&(psx.front()), &(psvx.front()), static_cast<int>(psx.size()));
        PlotPoincare->replot();
        CurvesPoincare<<(CurvePS);
    }
    if(tps!=NULL) delete tps;

#ifdef USEQWT3D
#ifdef USE3DFACETS
    if(core->potential->N_dim==3)
    {
        // calc max distance between adjacent points (for triangulation purposes)
        maxdist=trj->getMaxDist() * 1.01;
        ui.radioButton_3dmesh_full->setChecked(true);
        ui.lineEdit_3dmesh_maxdist->setText(QString::number(maxdist, 'g', 3));
        on_pushButton_refresh3dmesh_clicked();
    }
#endif
#endif

    // perform frequency graphs and orbit classification (also set label text and redraw orbit)
    if(ui.spinBox_intervalCount->value()!=1)
        ui.spinBox_intervalCount->setValue(1);  // refreshes frequency graph and calls for orbit classification
    else
        on_spinBox_intervalCount_valueChanged(1);
    // lyapunov
    if(InitData.calcLyapunov)
    {
      size_t nActual=trj->getTrajSize();
      if(nActual<100)
      {
        CurveLyapunov->setData(NULL, NULL, 0);
        CurveDeviation->setData(NULL, NULL, 0);
        PlotLyapunov->replot();
      }
      else
      {
        vectord timesteps(nActual);
        vectord lambdatmp(nActual);
        double devvecmin=1e4;
        for(size_t s=100; s<nActual; s++) 
        {
            timesteps[s] = s*InitData.timeStep/InitData.timeUnit;
            double lnw=trj->getLnDevVec(s);
            lambdatmp[s] = std::max<double>(lnw/timesteps[s], 1e-4);
            devvecmin = std::min<double>(lnw - log(timesteps[s]), devvecmin);
        }
        double devvecadd=(devvecmin<0.1?0.1-devvecmin:0);
        CurveLyapunov->setData(&(timesteps.at(100)), &(lambdatmp.at(100)), static_cast<int>(nActual)-100);
        for(size_t s=100; s<nActual; s++) 
            lambdatmp[s] = trj->getLnDevVec(s) - log(timesteps[s]) + devvecadd;
        CurveDeviation->setData(&(timesteps.at(100)), &(lambdatmp.at(100)), static_cast<int>(nActual)-100);
        //for(int i=0; i<PlotLyapunov->axisCnt; i++) PlotLyapunov->setAxisAutoScale(i); 
        PlotLyapunov->setAxisScale(QwtPlot::yLeft, 1e-4, 10);
        PlotLyapunov->setAxisScale(QwtPlot::yRight, 0.1, 1e4);
        PlotLyapunov->setAxisScale(QwtPlot::xBottom, 1, timesteps.back()*1.05);
        PlotLyapunov->replot();
      }
    }
}

/// ----------------------------  orbit rendering --------- ///
void CSmileGUI::redrawOrbit(bool show)
{
    if(!show || !core->orbit || !core->orbit->getIntTime()) return;

    const COrbitRuntimeTrajectory* trj = static_cast<const COrbitRuntimeTrajectory*>(core->orbit->getRuntimeFnc(0));   // 0th runtime function is the orbit trajectory analyzer
    assert(trj->FncType()==COrbitRuntimeTrajectory::FT_TRAJ_ANALYSIS);
    size_t nActual=trj->getTrajSize();
    const COrbitInitData<double> InitData=core->orbit->getInitData();

    size_t intNum = ui.spinBox_IntervalNum->value()-1;
    size_t intCount = ui.spinBox_intervalCount->value();
    size_t start=nActual * intNum / intCount;
    size_t count=nActual / intCount; 
    if(start+count>nActual) count=nActual-start;
    ui.groupBox_2dplane->setEnabled(true);
    if(ui.radioButton_2dorbit->isChecked())
    {   // draw 2d orbit projection
        ui.groupBox_3dmesh->setEnabled(false);
        if(ui.radioButton_otP->isChecked())
        {
            vectord hor(count), ver(count);
            size_t indh=0, indv=0;
            if(InitData.potential->N_dim==2)
            {
                indh=0; indv=1;
            }
            else
            {
                if(ui.radioButton_2d_yz->isChecked()) { indh=1; indv=2; }
                if(ui.radioButton_2d_xz->isChecked()) { indh=0; indv=2; }
                if(ui.radioButton_2d_xy->isChecked()) { indh=0; indv=1; }
            }
            for(size_t i=0; i<count; i++)
            {
                const CPosVelPoint<double>& pt = trj->getTraj(i+start);
                hor[i]=pt.Pos[indh];
                ver[i]=pt.Pos[indv];
            }
            if(count>0)
                CurveOrbit->setData(&(hor.front()), &(ver.front()), static_cast<int>(count));
            else
                CurveOrbit->setData(NULL, NULL, 0);
            // equipotential surface
            double eqh[NUM_POINTS_PLOT_EQUIPOTENTIAL*4], eqv[NUM_POINTS_PLOT_EQUIPOTENTIAL*4];
            double E = core->initCondE;
            int num_points_plot_equipotential = (InitData.potential->PotentialType()==CDensity::PT_NB)?NUM_POINTS_PLOT_EQUIPOTENTIAL/5:NUM_POINTS_PLOT_EQUIPOTENTIAL;  
            CDensity::SYMMETRYTYPE sym=InitData.potential->symmetry();
            for(int i=0; i<num_points_plot_equipotential; i++)
            {
                double angle=(i+0.5)*M_PI/2/num_points_plot_equipotential;
                double scrcoords[4][3] = { // four directions(quadrants) and three coords: horizontal, vertical and unused
                    { cos(angle), sin(angle), 0}, 
                    {-cos(angle), sin(angle), 0}, 
                    {-cos(angle),-sin(angle), 0}, 
                    { cos(angle),-sin(angle), 0} };
                int coordindex[3]={0,0,0};  // X_k = scrcoords[quadrant][coordindex[k]];
                if(ui.radioButton_2d_xy->isChecked()) { coordindex[0]=0; coordindex[1]=1; coordindex[2]=2; }
                if(ui.radioButton_2d_yz->isChecked()) { coordindex[0]=2; coordindex[1]=0; coordindex[2]=1; }
                if(ui.radioButton_2d_xz->isChecked()) { coordindex[0]=0; coordindex[1]=2; coordindex[2]=1; }
                double r0=InitData.potential->findintersection(E, scrcoords[0][coordindex[0]], scrcoords[0][coordindex[1]], scrcoords[0][coordindex[2]]);
                double r1=r0, r2=r0, r3=r0;
                if((sym & CDensity::ST_TRIAXIAL) == 0)  // no symmetry w.r.t. changing sign of any coordinate
                {
                    r1=InitData.potential->findintersection(E, scrcoords[1][coordindex[0]], scrcoords[1][coordindex[1]], scrcoords[1][coordindex[2]]);
                    if((sym & CDensity::ST_REFLECTION) == 0)  // no symmetry w.r.t. reflection of all three coords
                    {
                        r2=InitData.potential->findintersection(E, scrcoords[2][coordindex[0]], scrcoords[2][coordindex[1]], scrcoords[2][coordindex[2]]);
                        r3=InitData.potential->findintersection(E, scrcoords[3][coordindex[0]], scrcoords[3][coordindex[1]], scrcoords[3][coordindex[2]]);
                    }
                    else {r2=r0; r3=r1; }
                }
                eqh[i]=scrcoords[0][0]*r0; 
                eqv[i]=scrcoords[0][1]*r0;
                eqh[num_points_plot_equipotential*2-1-i]=scrcoords[1][0]*r1;
                eqv[num_points_plot_equipotential*2-1-i]=scrcoords[1][1]*r1;
                eqh[num_points_plot_equipotential*2+i]  =scrcoords[2][0]*r2;
                eqv[num_points_plot_equipotential*2+i]  =scrcoords[2][1]*r2;
                eqh[num_points_plot_equipotential*4-1-i]=scrcoords[3][0]*r3;
                eqv[num_points_plot_equipotential*4-1-i]=scrcoords[3][1]*r3;
            }
            CurveEquipotential->setData(eqh, eqv, 4*num_points_plot_equipotential);
            for(int i=0; i<PlotOrbit->axisCnt; i++) PlotOrbit->setAxisAutoScale(i); 
            PlotOrbit->replot();
            ZoomerOrbit->setZoomBase(false);
        }
        else
        {
            vectord hor(count), ver(count);
            for(size_t i=0; i<count; i++)
            {
                const CPosVelPoint<double>& pt = trj->getTraj(i+start);
                double Lx=pt.Pos[1]*pt.Vel[2] - pt.Pos[2]*pt.Vel[1];
                double Ly=pt.Pos[2]*pt.Vel[0] - pt.Pos[0]*pt.Vel[2];
                double Lz=pt.Pos[0]*pt.Vel[1] - pt.Pos[1]*pt.Vel[0];
                if(ui.radioButton_otL->isChecked())
                {
                    if(ui.radioButton_2d_yz->isChecked()) { hor[i]=Ly; ver[i]=Lz; }
                    if(ui.radioButton_2d_xz->isChecked()) { hor[i]=Lx; ver[i]=Lz; }
                    if(ui.radioButton_2d_xy->isChecked()) { hor[i]=Lx; ver[i]=Ly; }
                }
                else
                {
                    double r=sqrt(pow_2(pt.Pos[0])+pow_2(pt.Pos[1])+pow_2(pt.Pos[2]));
                    double ex=pt.Vel[1]*Lz-pt.Vel[2]*Ly - pt.Pos[0]/r;
                    double ey=pt.Vel[2]*Lx-pt.Vel[0]*Lz - pt.Pos[1]/r;
                    double ez=pt.Vel[0]*Ly-pt.Vel[1]*Lx - pt.Pos[2]/r;
                    if(ui.radioButton_2d_yz->isChecked()) { hor[i]=ey; ver[i]=ez; }
                    if(ui.radioButton_2d_xz->isChecked()) { hor[i]=ex; ver[i]=ez; }
                    if(ui.radioButton_2d_xy->isChecked()) { hor[i]=ex; ver[i]=ey; }
                }
            }
            if(count>0)
                CurveOrbit->setData(&(hor.front()), &(ver.front()), static_cast<int>(count));
            else
                CurveOrbit->setData(NULL, NULL, 0);
            CurveEquipotential->setData (NULL, NULL, 0);
            for(int i=0; i<PlotOrbit->axisCnt; i++) PlotOrbit->setAxisAutoScale(i); 
            PlotOrbit->replot();
            ZoomerOrbit->setZoomBase(false);
        }
#ifdef USEQWT3D
        Plot3d->hide();
#endif
        PlotOrbitParam->hide();
        PlotOrbit->show();
    }
    else if(ui.radioButton_3dline->isChecked() || ui.radioButton_3dmesh->isChecked())
    {   // draw 3d orbit
        ui.groupBox_2dplane->setEnabled(false);
        ui.groupBox_3dmesh->setEnabled(false);
#ifdef USEQWT3D 
        if(ui.radioButton_3dline->isChecked())
        {
            Qwt3D::Cell trjcell;
            trjpos.clear();
            trjcells.clear();
            trjcell.clear();
            for(size_t i=0; i<count; i++)
            {
                const CPosVelPoint<double>& pt = trj->getTraj(i+start);
                trjpos.push_back(Qwt3D::Triple(pt.Pos[0], pt.Pos[1], pt.Pos[2]));
                trjcell.push_back(i);
            }
            trjcells.push_back(trjcell);
            // draw trajectory
            Plot3d->loadFromData(trjpos, trjcells);
            Plot3d->setPlotStyle(Qwt3D::WIREFRAME);
        }
#ifdef USE3DFACETS
        if(ui.radioButton_3dmesh->isChecked())
        {
            // draw mesh
            if(TriangThr->restart)   // triangulation is in process, so do not display anything
                Plot3d->loadFromData(Qwt3D::TripleField(), Qwt3D::CellField());
            else
                Plot3d->loadFromData(TriangThr->trjpos, TriangThr->facets);
            Plot3d->setPlotStyle(Qwt3D::FILLEDMESH);
            ui.groupBox_3dmesh->setEnabled(true);
        }
#endif
        Plot3d->updateData();
        Plot3d->updateGL();
        Plot3d->show();
        PlotOrbitParam->hide();
        PlotOrbit->hide();
#endif
    }
    else
    {   // orbit params distribution
        ui.groupBox_2dplane->setEnabled(false);
        ui.groupBox_3dmesh->setEnabled(false);
        vectord data;
        const COrbitRuntimePericenter* tpc=NULL;
        for(size_t rf=0; rf<core->orbit->getInfoNum(); rf++)
            if(core->orbit->getRuntimeFnc(rf)->FncType() == CBasicOrbitRuntimeFnc::FT_PERICENTER)
                tpc = static_cast<const COrbitRuntimePericenter* >(core->orbit->getRuntimeFnc(rf));
        if(ui.radioButton_opLPeri->isChecked())
        {
            if(tpc!=NULL)
                for(size_t i=0; i<tpc->size(); i++)
                    data.push_back(tpc->getl2(i));
        }
        else if(ui.radioButton_opPericenter->isChecked())
        {
            if(tpc!=NULL)
                for(size_t i=0; i<tpc->size(); i++)
                    data.push_back(tpc->getr(i));
        }
        else if(ui.radioButton_opL->isChecked())
        {
            data.resize(nActual);
            for(size_t i=0; i<nActual; i++)
            {
                const CPosVelPoint<double>& pt = trj->getTraj(i);
                data[i] = pow_2(pt.Pos[1]*pt.Vel[2] - pt.Pos[2]*pt.Vel[1]) +
                        + pow_2(pt.Pos[2]*pt.Vel[0] - pt.Pos[0]*pt.Vel[2]) +
                        + pow_2(pt.Pos[0]*pt.Vel[1] - pt.Pos[1]*pt.Vel[0]);
            }
        }
        size_t nump=data.size();
        if(nump>0)
        {
            std::sort(data.begin(), data.end());
            vectord frac(nump);
            for(size_t i=0; i<nump; i++)
                frac[i]=i*1.0/nump;
            CurveOrbitParam->setData(&(frac.front()), &(data.front()), static_cast<int>(nump));
        }
        else CurveOrbitParam->setData(NULL, NULL, 0);
        for(int i=0; i<PlotOrbitParam->axisCnt; i++) PlotOrbitParam->setAxisAutoScale(i); 
        //PlotOrbitParam->replot();
        ZoomerOrbitParam->setZoomBase(true);
#ifdef USEQWT3D
        Plot3d->hide();
#endif
        PlotOrbit->hide();
        PlotOrbitParam->show();
    }
}

void CSmileGUI::on_pushButton_refresh3dmesh_clicked()
{
#ifdef USE3DFACETS
    bool ok;
    double md = ui.lineEdit_3dmesh_maxdist->text().toDouble(&ok);   
    if(ok) maxdist=md;
    int meshstyle=0;
    if(ui.radioButton_3dmesh_x->isChecked()) meshstyle=1; else
    if(ui.radioButton_3dmesh_y->isChecked()) meshstyle=2; else
    if(ui.radioButton_3dmesh_z->isChecked()) meshstyle=3; 
    std::vector<CPosPoint<double> > trajp;
    assert(core->orbit);
    const COrbitRuntimeTrajectory* trj = static_cast<const COrbitRuntimeTrajectory*>(core->orbit->getRuntimeFnc(0));   // 0th runtime function is the orbit trajectory analyzer
    assert(trj->FncType()==COrbitRuntimeTrajectory::FT_TRAJ_ANALYSIS);
    size_t nActual=trj->getTrajSize();
    trajp.reserve(nActual);
    for(size_t i=0; i<nActual; i++)
        trajp.push_back(trj->getTraj(i));
    TriangThr->startTriangulation(trajp, maxdist, meshstyle);
#endif
}
/// ----------------------------  spectrum rendering --------- ///
void CSmileGUI::on_spinBox_intervalCount_valueChanged(int count)
{
    if(InternalRecalc) return;
    bool doNotRefresh = (ui.spinBox_IntervalNum->value() > count);
    ui.spinBox_IntervalNum->setMaximum(count);
    if(doNotRefresh) return;   // setting maximum more than current value led to change of value of IntervalNum and hence the refresh function was already called
    if(ui.spinBox_IntervalNum->value()!=1)
        ui.spinBox_IntervalNum->setValue(1);  // refreshes frequency graph and calls for orbit classification
    else
        refreshOrbitAndSpectrum();
}

void CSmileGUI::on_spinBox_IntervalNum_valueChanged(int)
{
    refreshOrbitAndSpectrum();
}

void CSmileGUI::refreshOrbitAndSpectrum() 
{
    const COrbitRuntimeTrajectory* trj = static_cast<const COrbitRuntimeTrajectory*>(core->orbit->getRuntimeFnc(0));   // 0th runtime function is the orbit trajectory analyzer
    assert(trj->FncType()==COrbitRuntimeTrajectory::FT_TRAJ_ANALYSIS);
    size_t nActual=trj->getTrajSize();
    const COrbitInitData<double> InitData=core->orbit->getInitData();

    size_t intNum = ui.spinBox_IntervalNum->value() - 1;
    size_t intCount = ui.spinBox_intervalCount->value();
    if (intNum>=intCount) return;
    size_t start=nActual * intNum / intCount;
    size_t count=nActual / intCount; 
    if(start+count>nActual) count=nActual-start;

    timeOI.start();
    // obtain spectrum for given interval
    CSpectrum sl[N_DIM];
    std::vector<double> fa[N_DIM];
    COrbitInformation<double>* infoOrbit = trj->analyzeOrbit<double>(start, count, sl, fa);
    int nPointsSpectrum = static_cast<int>(fa[0].size()); 
    if(nPointsSpectrum<=0) return;  // ??? shouldn't happen but who knows..
    vectord freq(nPointsSpectrum);
    for(int i=0; i<nPointsSpectrum; i++)
        freq[i]=i/core->orbit->getIntTime()*InitData.timeUnit * intCount;
    CurveFreqX->setData(&freq.front(), &fa[0].front(), nPointsSpectrum*std::min<double>(2*PLOT_FREQ_MAX_F/ui.lineEdit_timeStep->text().toDouble(),1));
    CurveFreqY->setData(&freq.front(), &fa[1].front(), nPointsSpectrum*std::min<double>(2*PLOT_FREQ_MAX_F/ui.lineEdit_timeStep->text().toDouble(),1));
    if(core->potential->N_dim==3)
        CurveFreqZ->setData(&freq.front(), &fa[2].front(), nPointsSpectrum*std::min<double>(2*PLOT_FREQ_MAX_F/ui.lineEdit_timeStep->text().toDouble(),1));
    else
        CurveFreqZ->setData(NULL, NULL, 0);
    double sf[MAX_SPECTRAL_LINES*N_DIM], sa[MAX_SPECTRAL_LINES*N_DIM];
    double sfp[MAX_SPECTRAL_LINES*N_DIM], sap[MAX_SPECTRAL_LINES*N_DIM];
    int si=0, sip=0;
    for(unsigned int d=0; d<core->potential->N_dim; d++)
      for(size_t s=0; s<sl[d].size() && s<MAX_SPECTRAL_LINES*N_DIM; s++) 
          if(sl[d][s].freq<5)
          {
            if(sl[d][s].precise)
            {
                sfp[sip]=sl[d][s].freq;
                sap[sip]=sl[d][s].ampl;
                sip++;
            }
            else
            {
                sf[si]=sl[d][s].freq;
                sa[si]=sl[d][s].ampl;
                si++;
            }
          }
    CurveFreqPeaks->setData(sf, sa, si);
    CurveFreqPeaksP->setData(sfp, sap, sip);
    for(int i=0; i<PlotFreq->axisCnt; i++) PlotFreq->setAxisAutoScale(i); 
    PlotFreq->replot();
    ZoomerFreq->setZoomBase(true);
    ///!!!
    /*{
        std::ofstream strm("d:\\document\\spec.txt", std::ios::out);
        if(strm) for(int i=0; i<nPointsSpectrum; i++)
        {
            strm << freq[i] << "\t" << fa[0][i] << "\t" << fa[1][i] << "\t" << fa[2][i] << "\n";
        }
        strm << "\n";
        for(size_t s=0; s<10 && s<MAX_SPECTRAL_LINES*N_DIM; s++) 
        {
            if(s<sl[0].size()) strm << sl[0][s].freq << "\t" << sl[0][s].ampl << "\t" << sl[0][s].precise << "\n";
            if(s<sl[1].size()) strm << sl[1][s].freq << "\t" << sl[1][s].ampl << "\t" << sl[1][s].precise << "\n";
            if(s<sl[2].size()) strm << sl[2][s].freq << "\t" << sl[2][s].ampl << "\t" << sl[2][s].precise << "\n";
        }
    }*/

    // set label text
    std::string strInfo = infoOrbit->toString();   // obtain information from orbit analysis - not for entire trajectory but for a segment of it.
    delete infoOrbit;
    for(size_t rf=1; rf<core->orbit->getInfoNum(); rf++)    // obtain the rest (starting from 1 since 0th runtime function is the orbit analysis which was called manually above)
    {
        const CBasicInformation* info=core->orbit->getNewInfo(rf);
        strInfo += info->toString();
        delete info;
    }
    strInfo+=core->orbit->toString(false);
    info(QString::fromStdString(strInfo)
        + "Time elapsed="+QString::number((timeOIelapsed+timeOI.elapsed())*0.001)+" s");
    redrawOrbit(true);  // redraw trajectory for given interval start--count
}

void CSmileGUI::on_checkBox_freqLines_stateChanged(int)
{
    if(ui.checkBox_freqLines->isChecked())
    {
        CurveFreqPeaksP->setStyle(QwtPlotCurve::Sticks);
        CurveFreqPeaks->setStyle(QwtPlotCurve::Sticks);
    }
    else
    {
        CurveFreqPeaksP->setStyle(QwtPlotCurve::NoCurve);
        CurveFreqPeaks->setStyle(QwtPlotCurve::NoCurve);
    }
    PlotFreq->replot();   
}

/// ----------------------------  update initial conditions --------- ///
void CSmileGUI::potentialChanged()
{
    if(InternalRecalc) return;
    InternalRecalc=true;    // temporary set flag so that changing options below does not call this function again
    CDensity::POTENTIALTYPE PotentialType = getPotentialTypeByName(ui.selectPotentialType->currentText().toStdString());
    CDensity::POTENTIALTYPE DensityType   = getDensityTypeByName  (ui.selectDensityType  ->currentText().toStdString());
    CDensity::SYMMETRYTYPE SymmetryType   = getSymmetryTypeByName (ui.selectSymmetryType ->currentText().toStdString());

    if(PotentialType==CDensity::PT_BSE || PotentialType==CDensity::PT_SPLINE)
    {
        ui.selectDensityType->show();
        ui.label_densityType->show();
        ui.lineEdit_Ncoefs_radial->show();
        ui.label_Ncoefs_radial->show();
    }
    else
    {
        ui.selectDensityType->hide();
        ui.label_densityType->hide();
        ui.lineEdit_Ncoefs_radial->hide();
        ui.label_Ncoefs_radial->hide();
    }
    if(PotentialType==CDensity::PT_BSE)
    {
        ui.lineEdit_Alpha->show();
        ui.label_Alpha->show();
    }
    else
    {
        ui.lineEdit_Alpha->hide();
        ui.label_Alpha->hide();
    }
    if(PotentialType==CDensity::PT_NB || ((PotentialType==CDensity::PT_BSE || PotentialType==CDensity::PT_SPLINE) && DensityType==CDensity::PT_NB))
    {
        ui.lineEdit_NbodyFile->show();
        ui.buttonChooseNbodyFile->show();
        ui.lineEdit_p->hide();
        ui.lineEdit_q->hide();
        ui.label_p->hide();
        ui.label_q->hide();
        if(PotentialType!=CPotential::PT_NB) {
            ui.selectSymmetryType->show(); 
            ui.label_symmetryType->show();
        }else{
            ui.selectSymmetryType->hide();
            ui.label_symmetryType->hide();
        }
    }
    else
    {
        ui.lineEdit_NbodyFile->hide();
        ui.buttonChooseNbodyFile->hide();
        ui.lineEdit_p->show();
        ui.lineEdit_q->show();
        ui.label_p->show();
        ui.label_q->show();
        ui.selectSymmetryType->hide();
        ui.label_symmetryType->hide();
    }
    if(PotentialType==CDensity::PT_BSE || PotentialType==CDensity::PT_SPLINE || PotentialType==CDensity::PT_SCALEFREESH)
    {
        ui.lineEdit_Ncoefs_angular->show();
        ui.label_Ncoefs_angular->show();
    }
    else
    {
        ui.lineEdit_Ncoefs_angular->hide();
        ui.label_Ncoefs_angular->hide();
    }
    if(PotentialType==CDensity::PT_LOG || 
        ((PotentialType==CDensity::PT_BSE || PotentialType==CDensity::PT_SPLINE) && DensityType==CDensity::PT_NFW) )
    {
        ui.lineEdit_Rc->show();
        ui.label_core->show();
    }
    else
    {
        ui.lineEdit_Rc->hide();
        ui.label_core->hide();
    }
    if(PotentialType==CDensity::PT_DEHNEN || PotentialType==CDensity::PT_SCALEFREE || PotentialType==CDensity::PT_SCALEFREESH || 
      ((PotentialType==CDensity::PT_BSE || PotentialType==CDensity::PT_SPLINE) && DensityType==CDensity::PT_DEHNEN))
    {
        ui.lineEdit_Gamma->show();
        ui.label_cusp->show();
    }
    else
    {
        ui.lineEdit_Gamma->hide();
        ui.label_cusp->hide();
    }
    if(PotentialType==CDensity::PT_NB)
    {
        ui.label_eps->show();
        ui.label_tol->show();
        ui.lineEdit_eps->show();
        ui.lineEdit_tol->show();
        // disable Lyapunov exponent
        core->InitData.calcLyapunov=false;
        ui.checkBox_Lyapunov->setChecked(false);
        ui.checkBox_Lyapunov->setEnabled(false);
    }
    else
    {
        ui.label_eps->hide();
        ui.label_tol->hide();
        ui.lineEdit_eps->hide();
        ui.lineEdit_tol->hide();
        ui.groupBox_2d3d->setEnabled(true);
        ui.checkBox_Lyapunov->setEnabled(true);
    }
    if(PotentialType==CDensity::PT_BSE || PotentialType==CDensity::PT_SPLINE || PotentialType==CDensity::PT_NB)
    {
        // set N_dim=3
        configPotential.N_dim=3;
        on_radioButton_2d_toggled(false);
        ui.radioButton_3d->setChecked(true);
        ui.groupBox_2d3d->setEnabled(false);
    }
    InternalRecalc=false;
    bool ok=true;
    double newq=(ui.lineEdit_q->text().toDouble(&ok));  if(!ok || (newq==0)) return;
    double newp=(ui.lineEdit_p->text().toDouble(&ok));  if(!ok || (newp==0)) return;
    double newRc=ui.lineEdit_Rc->text().toDouble(&ok);  if (!ok) return;
    double newGamma=ui.lineEdit_Gamma->text().toDouble(&ok);  if (!ok) return;
    double newOmega=ui.lineEdit_Omega->text().toDouble(&ok);  if (!ok) return;
    double newAlpha=ui.lineEdit_Alpha->text().toDouble(&ok);  if (!ok) return;
    double newMbh=ui.lineEdit_Mbh->text().toDouble(&ok);  if(!ok) return;
    unsigned int newNcoefs_radial=ui.lineEdit_Ncoefs_radial->text().toInt(&ok);  if(!ok) return;
    unsigned int newNcoefs_angular=ui.lineEdit_Ncoefs_angular->text().toInt(&ok);  if(!ok) return;
    double newnbpotEps=ui.lineEdit_eps->text().toDouble(&ok); if(!ok) return;
    double newnbpotTol=ui.lineEdit_tol->text().toDouble(&ok); if(!ok) return;
    if(PotentialType!=configPotential.PotentialType || 
        DensityType!=configPotential.DensityType || 
        SymmetryType!=configPotential.SymmetryType || 
        newq!=configPotential.q || newp!=configPotential.p || 
        newRc!=configPotential.Rc || newMbh!=configPotential.Mbh ||
        newGamma!=configPotential.Gamma || newAlpha!=configPotential.Alpha ||  
        newNcoefs_radial!=configPotential.Ncoefs_radial || 
        newNcoefs_angular!=configPotential.Ncoefs_angular || 
        newOmega!=configPotential.Omega || 
        newnbpotEps!=configPotential.nbpotEps || 
        newnbpotTol!=configPotential.nbpotTol)
    {
        configPotential.q=newq; configPotential.p=newp; 
        configPotential.Rc=newRc; 
        configPotential.Gamma=newGamma; 
        configPotential.Mbh=newMbh; 
        configPotential.Omega=newOmega;
        configPotential.Ncoefs_radial=newNcoefs_radial; 
        configPotential.Ncoefs_angular=newNcoefs_angular; 
        configPotential.Alpha=newAlpha;
        configPotential.nbpotEps=newnbpotEps; 
        configPotential.nbpotTol=newnbpotTol;
        configPotential.PotentialType=PotentialType;
        configPotential.DensityType=DensityType;
        configPotential.SymmetryType=SymmetryType;
        if(PotentialType!=CPotential::PT_BSE && PotentialType!=CPotential::PT_SPLINE && PotentialType!=CPotential::PT_NB)
            on_buttonInitPotential_clicked();   // do initialization for 'lightweight' potentials only; should manually press this button for 'heavy-weight' initialization of general-purpose potentials listed above
        else
            ui.buttonInitPotential->show();   // let the user press this button manually
    }
    ///nonpublic
/*    if(newOmega!=0)
    {   // calculate corotation radius
        info("Corotation radius="+QString::number(core->getPotential->corotationradius()));
    }*/
    ui.frame->repaint();
}

void CSmileGUI::on_buttonInitPotential_clicked()
{
    core->initPotential();
    ICchanged();
    if(configPotential.PotentialType==CPotential::PT_BSE || configPotential.PotentialType==CDensity::PT_SPLINE)
    {   // the symmetry type and parameters may have changed by potential constructor, especially if it loaded coefficients from file rather than initialized from Nbody file
        InternalRecalc=true;
        ui.selectSymmetryType ->setCurrentIndex(ui.selectSymmetryType ->findText(QString::fromStdString(SymmetryNames [configPotential.SymmetryType])));
        // ensure that N_radial, N_angular and Alpha input lines are consistent with corresponding variables (which may have been updated after loading a coef file)
        if(configPotential.Ncoefs_radial != static_cast<unsigned int>(ui.lineEdit_Ncoefs_radial->text().toInt()))
            ui.lineEdit_Ncoefs_radial->setText(QString::number(configPotential.Ncoefs_radial));
        if(configPotential.Ncoefs_angular != static_cast<unsigned int>(ui.lineEdit_Ncoefs_angular->text().toInt())) 
            ui.lineEdit_Ncoefs_angular->setText(QString::number(configPotential.Ncoefs_angular));
        if(configPotential.Alpha != ui.lineEdit_Alpha->text().toDouble())
            ui.lineEdit_Alpha->setText(QString::number(configPotential.Alpha));
        if(configPotential.Mbh != ui.lineEdit_Mbh->text().toDouble())
            ui.lineEdit_Mbh->setText(QString::number(configPotential.Mbh));
        InternalRecalc=false;
    }
    ui.buttonInitPotential->hide();   // this button will be shown only when parameters changed
}

void CSmileGUI::chooseNbodyFile()
{
    assignNbodyFile(QFileDialog::getOpenFileName(this, "Open Nbody file", core->WorkDir));
}

void CSmileGUI::typeinNbodyFile()
{
    QString fileName = ui.lineEdit_NbodyFile->text();
    if(!QFile::exists(fileName))
        fileName = QFileDialog::getOpenFileName(this, "Open Nbody file", core->WorkDir);
    assignNbodyFile(fileName);
}

void CSmileGUI::assignNbodyFile(QString fileName)
{
    if(fileName.isEmpty() || !QFile::exists(fileName)) return;
    core->WorkDir = getDirName(fileName);
    core->NbodyFile=fileName;
    ui.lineEdit_NbodyFile->setText(fileName);
    configPotential.PotentialType = getPotentialTypeByName(ui.selectPotentialType->currentText().toStdString());   // set explicitly because it may have been reset to default if previous load was unsuccessful
    on_buttonInitPotential_clicked();  // do initialization
}

void CSmileGUI::ICchanged()
{
    if(InternalRecalc) return;
    bool ok=true;
    double newE=ui.lineEdit_E->text().toDouble(&ok);  if(!ok) return;
    double X=ui.lineEdit_x->text().toDouble(&ok);  if(!ok) return;
    double Y=ui.lineEdit_y->text().toDouble(&ok);  if(!ok) return;
    double Z=(core->potential->N_dim==3)?ui.lineEdit_z->text().toDouble(&ok):0;  if(!ok) return;
    double Vx=ui.lineEdit_vx->text().toDouble(&ok);  if(!ok) return;
    double Vy=ui.lineEdit_vy->text().toDouble(&ok);  if(!ok) return;
    double Vz=(core->potential->N_dim==3)?ui.lineEdit_vz->text().toDouble(&ok):0;  if(!ok) return;
    if(core->useICe)
    {
        core->initCondE=newE;
    }
    else
    {
        core->InitData.initCond=CPosVelPoint<double>(X, Y, Z, Vx, Vy, Vz);
    }
    InternalRecalc=true;
    core->initIC();
    if(core->useICe)
    {
        ui.lineEdit_x->setText(QString::number(core->InitData.initCond.Pos[0]));
        ui.lineEdit_y->setText(QString::number(core->InitData.initCond.Pos[1]));
        ui.lineEdit_z->setText(QString::number(core->InitData.initCond.Pos[2]));
        ui.lineEdit_vx->setText(QString::number(core->InitData.initCond.Vel[0]));
        ui.lineEdit_vy->setText(QString::number(core->InitData.initCond.Vel[1]));
        ui.lineEdit_vz->setText(QString::number(core->InitData.initCond.Vel[2]));
    }
    else
    {
        ui.lineEdit_E->setText(QString::number(core->initCondE));
    }
    ui.lineEdit_Torb->setText(QString::number(core->InitData.timeUnit));
    InternalRecalc=false;
    if(fabs(core->initCondE-Eprev)> 1e-4)
        on_ButtonClearPS_clicked();
}

/// ----------------------------  modification of various parameters by GUI lineEdits ----- ///
void CSmileGUI::intTime_changed()
{
    if(InternalRecalc) return;
    bool ok;
    double t = ui.lineEdit_intTime->text().toDouble(&ok);
    if(ok)  configCore.intTimeInPeriods=t;
    t = ui.lineEdit_timeStep->text().toDouble(&ok);
    if(ok)  configCore.intTimeStepsPerPeriod=t;
}

void CSmileGUI::fmparams_changed()
{
    if(InternalRecalc) return;
    bool ok; int t;
    t = ui.lineEdit_fms->text().toInt(&ok);
    if(ok)  configCore.fm_numOrbitsStationary=t;
    t = ui.lineEdit_fmp->text().toInt(&ok);
    if(ok)  configCore.fm_numOrbitsPrincipalPlane=t;
    t = ui.lineEdit_fmy->text().toInt(&ok);
    if(ok)  configCore.fm_numOrbitsYalpha=t;
    t = ui.lineEdit_fmr->text().toInt(&ok);
    if(ok)  configCore.fm_numOrbitsRandom=t;
}

void CSmileGUI::smparams_changed()
{
    if(InternalRecalc) return;
    bool ok; int t; double d;
    // orbit library tab
    t = ui.lineEdit_smNumOrbits->text().toInt(&ok);
    if(ok)  configCore.sm_numOrbitsRandom=t;
    t = ui.lineEdit_smNumSamplingPoints->text().toInt(&ok);
    if(ok) configCore.sm_numSamplingPoints=t;
    // model params tab
    CBasicSchwModel::MODELTYPE ModelType = getSchwModelTypeByName(ui.select_smType->currentText().toStdString());
    if(ModelType==CBasicSchwModel::MT_CLASSIC)
    {
        ui.lineEdit_smNumLines->show();
        ui.label_smNumLines->show();
        ui.lineEdit_smNumAngCoefs->hide();
        ui.label_smNumAngCoefs->hide();
    }
    else
    {
        ui.lineEdit_smNumLines->hide();
        ui.label_smNumLines->hide();
        ui.lineEdit_smNumAngCoefs->show();
        ui.label_smNumAngCoefs->show();
    }
    if(ModelType==CBasicSchwModel::MT_SHBSE)
    {
        ui.lineEdit_smNumRadCoefs->show();
        ui.lineEdit_smAlpha->show();
        ui.label_smNumRadCoefs->show();
        ui.label_smAlpha->show();
    }
    else
    {
        ui.lineEdit_smNumRadCoefs->hide();
        ui.lineEdit_smAlpha->hide();
        ui.label_smNumRadCoefs->hide();
        ui.label_smAlpha->hide();
    }
    if(ModelType!=configSchw.ModelType)
    {
        ui.groupBox_smOptim->setEnabled(false);
        configSchw.ModelType=ModelType;
    }
    t = ui.lineEdit_smNumAngCoefs->text().toInt(&ok);
    if(ok&&(configCore.sm_numAngularCoefs!=t))
    {
        ui.groupBox_smOptim->setEnabled(false);
        configCore.sm_numAngularCoefs=t;
    }
    t = ui.lineEdit_smNumRadCoefs->text().toInt(&ok);
    if(ok&&(configCore.sm_numRadialCoefs!=t))
    {
        ui.groupBox_smOptim->setEnabled(false);
        configCore.sm_numRadialCoefs=t;
    }
    d = ui.lineEdit_smAlpha->text().toDouble(&ok);
    if(ok&&(configCore.sm_Alpha!=d))
    {
        ui.groupBox_smOptim->setEnabled(false);
        configCore.sm_Alpha=t;
    }
    t = ui.lineEdit_smNumLines->text().toInt(&ok);
    if(ok&&(configCore.sm_linesPerSegment!=t))
    {
        ui.groupBox_smOptim->setEnabled(false);
        configCore.sm_linesPerSegment=t;
    }
    t = ui.lineEdit_smNumShells->text().toInt(&ok);
    if(ok)
    {
        if(configCore.sm_numShells!=t) ui.groupBox_smOptim->setEnabled(false);
        ui.selectEnergyLevel->setMaximum(t);
        ui.selectEnergyLevel_fm->setMaximum(t);
        configCore.sm_numShells=t;
    }
    d = ui.lineEdit_smInnerShellMass->text().toDouble(&ok);
    if(ok&&(configCore.sm_innerShellMass!=d))
    {
        ui.groupBox_smOptim->setEnabled(false);
        configCore.sm_innerShellMass=d;
    }
    d = ui.lineEdit_smOuterShellMass->text().toDouble(&ok);
    if(ok&&(configCore.sm_outerShellMass!=d))
    {
        ui.groupBox_smOptim->setEnabled(false);
        configCore.sm_outerShellMass=d;
    }
    // optimization tab
    d = ui.lineEdit_smChaoticWeight->text().toDouble(&ok);
    if(ok) configCore.chaoticWeightFactor = d;
    d = ui.lineEdit_smMaxOrbitWeight->text().toDouble(&ok);
    if(ok) configSchw.sm_maxWeight = d;
    configSchw.sm_constrainBeta=ui.checkBox_smBeta->isChecked();
    ui.lineEdit_smBetaIn->setEnabled(configSchw.sm_constrainBeta);
    ui.lineEdit_smBetaOut->setEnabled(configSchw.sm_constrainBeta);
    d = ui.lineEdit_smBetaIn->text().toDouble(&ok);
    if(ok)  configSchw.sm_betaIn = d;
    d = ui.lineEdit_smBetaOut->text().toDouble(&ok);
    if(ok)  configSchw.sm_betaOut = d;
    // compute # of coefs
    int numCoefs=0;
    switch(ModelType)
    {
    case CBasicSchwModel::MT_CLASSIC: 
        numCoefs=3*configCore.sm_linesPerSegment*configCore.sm_linesPerSegment*configCore.sm_numShells;
        break;
    case CBasicSchwModel::MT_SHGRID: 
        numCoefs=1+configCore.sm_numShells*(configCore.sm_numAngularCoefs/2+1)*(configCore.sm_numAngularCoefs/2+2)/2;
        break;
    case CBasicSchwModel::MT_SHBSE: 
        numCoefs=(1+configCore.sm_numRadialCoefs)*(configCore.sm_numAngularCoefs/2+1)*(configCore.sm_numAngularCoefs/2+2)/2;
        break;
    default: info("Invalid type of Schwarzschild model");
    }
    ui.lineEdit_smNumCoefs->setText(QString::number(numCoefs));
}

void CSmileGUI::ccparams_changed()
{
    if(InternalRecalc) return;
    bool ok;
    double d = ui.lineEdit_MinFreqDiff->text().toDouble(&ok);
    if(ok)  configCore.chaoticMinFreqDiff=d;
    d = ui.lineEdit_MinLambda->text().toDouble(&ok);
    if(ok)  configCore.chaoticMinLambda=d;
    redrawFreqMap(true);
}

void CSmileGUI::on_checkBox_Lyapunov_stateChanged(int Checked)
{
    if(InternalRecalc) return;
    core->InitData.calcLyapunov = (Checked>0);
}

void CSmileGUI::on_checkBox_usePS_stateChanged(int Checked)
{
    if(InternalRecalc) return;
    configCore.usePS = (Checked>0);
}

void CSmileGUI::on_radioButton_2d_toggled(bool n2d)
{
    if(n2d)
    {
        ui.radioButton_2d_xy->setChecked(true);
        ui.radioButton_2d_xz->setEnabled(false);
        ui.radioButton_2d_yz->setEnabled(false);
        ui.radioButton_3dline->setEnabled(false);
        ui.radioButton_3dmesh->setEnabled(false);
    }
    else
    {
        ui.radioButton_2d_xz->setEnabled(true);
        ui.radioButton_2d_yz->setEnabled(true);
#ifdef USEQWT3D
        ui.radioButton_3dline->setEnabled(true);
#ifdef USE3DFACETS
        ui.radioButton_3dmesh->setEnabled(true);
#endif
#endif
    }
    if(InternalRecalc) return;
    configPotential.N_dim=n2d?2:3;
    core->initPotential();
    ICchanged();
}

void CSmileGUI::on_radioButton_ICe_toggled(bool useICe)
{
    if(InternalRecalc) return;
    core->useICe=useICe;
    if(useICe)
    {
        ui.lineEdit_E->setEnabled(true);
        ui.groupBox_ICx->setEnabled(false);
    }
    else
    {
        ui.lineEdit_E->setEnabled(false);
        ui.groupBox_ICx->setEnabled(true);
    }
}

/// ----------------------------  Random initial conditions --------- ///
void CSmileGUI::on_ButtonRandom_clicked()
{
    double Xmax=core->potential->longaxisradius(core->initCondE);
    if(Xmax<0) return; // error
    const int npsampling=1000;
    std::vector<double> Xs(npsampling), Ys(npsampling), Zs(npsampling), Ws(npsampling);
    double X, Y, Z, Et, totalWeight=0;
    for(int s=0; s<npsampling; s++)
    {   // sample accessible configuration space with appropriate weight
      do{
        X = (rand()*2.0/RAND_MAX-1) * Xmax;
        Y = (rand()*2.0/RAND_MAX-1) * Xmax;
        Z = core->potential->N_dim==3 ? (rand()*2.0/RAND_MAX-1) * Xmax : 0;
        Et = core->potential->Phi(X, Y, Z);
      }
      while(Et>core->initCondE);
      double W=sqrt(2*(core->initCondE-Et));
      Xs[s]=X; Ys[s]=Y; Zs[s]=Z; Ws[s]=W; totalWeight+=W;
    }
    // now select point from sampled
    double indw=rand() * totalWeight / RAND_MAX;
    int i=0;
    while(i<npsampling-1 && indw>=Ws[i]) 
    {
        indw-=Ws[i];
        i++;
    }
    X=Xs[i]; Y=Ys[i]; Z=Zs[i]; 
    double V=Ws[i];
    double costh=(core->potential->N_dim==3)? (rand()*2.0/RAND_MAX-1) : 0;
    double sinth=sqrt(1-pow_2(costh));
    double phi=rand()*2*M_PI/RAND_MAX;
    InternalRecalc=true;

    core->InitData.initCond=CPosVelPoint<double>(X, Y, Z, V*sinth*cos(phi), V*sinth*sin(phi), V*costh);
    //core->initIC();
    ui.lineEdit_x->setText(QString::number(core->InitData.initCond.Pos[0]));
    ui.lineEdit_y->setText(QString::number(core->InitData.initCond.Pos[1]));
    ui.lineEdit_z->setText(QString::number(core->InitData.initCond.Pos[2]));
    ui.lineEdit_vx->setText(QString::number(core->InitData.initCond.Vel[0]));
    ui.lineEdit_vy->setText(QString::number(core->InitData.initCond.Vel[1]));
    ui.lineEdit_vz->setText(QString::number(core->InitData.initCond.Vel[2]));
    InternalRecalc=false;
    on_ButtonStart_clicked();
}

/// ----------------------------  Poincare section operations --------- ///
void CSmileGUI::PSclicked(const QPointF &pos)
{
    double newx=pos.x();
    double newvx=pos.y();
    double newvy2=2*(core->initCondE - core->potential->Phi(newx, 0, 0)) - pow_2(newvx);
    if(newvy2>=0)
    {
        InternalRecalc=true;
        ui.lineEdit_x->setText(QString::number(newx, 'f', 7));
        ui.lineEdit_vx->setText(QString::number(newvx, 'f', 7));
        ui.lineEdit_vy->setText(QString::number(sqrt(newvy2), 'f', 7));
        ui.lineEdit_y->setText("0");
        ui.lineEdit_z->setText("0");
        ui.lineEdit_vz->setText("0");
        InternalRecalc=false;
        ICchanged();
    }
}

void CSmileGUI::on_ButtonClearPS_clicked()
{
    while(CurvesPoincare.size()>0)
    {
        delete CurvesPoincare[CurvesPoincare.size()-1];
        CurvesPoincare.remove(CurvesPoincare.size()-1);
    }
    double xm[NUM_POINTS_PLOT_POINCARE*4], vxm[NUM_POINTS_PLOT_POINCARE*4];
    double Xmax=core->potential->longaxisradius(core->initCondE);
    if(Xmax>0) 
    {
        for(unsigned int i=0; i<NUM_POINTS_PLOT_POINCARE; i++)
        {
            xm[i]=Xmax*(i+0.5)/NUM_POINTS_PLOT_POINCARE;
            double vx2=2*(core->initCondE - core->potential->Phi(xm[i], 0, 0));
            if(vx2>0) vxm[i]=sqrt(vx2); else vxm[i]=0;
            xm[NUM_POINTS_PLOT_POINCARE*2-1-i]=xm[i];
            vxm[NUM_POINTS_PLOT_POINCARE*2-1-i]=-vxm[i];
            xm[NUM_POINTS_PLOT_POINCARE*2+i]=-xm[i];
            vxm[NUM_POINTS_PLOT_POINCARE*2+i]=-vxm[i];
            xm[NUM_POINTS_PLOT_POINCARE*4-1-i]=-xm[i];
            vxm[NUM_POINTS_PLOT_POINCARE*4-1-i]=vxm[i];
        }
        QwtPlotCurve* CurveBoundary = new QwtPlotCurve();
        CurveBoundary->setPen(QColor(Qt::red));
        CurveBoundary->attach(PlotPoincare);
        CurveBoundary->setData(xm, vxm, 4*NUM_POINTS_PLOT_POINCARE);
        CurvesPoincare<<(CurveBoundary);
    }
    for(int i=0; i<PlotPoincare->axisCnt; i++) PlotPoincare->setAxisAutoScale(i); 
    PlotPoincare->replot();
    ZoomerPoincare->setZoomBase(true);
    Eprev=core->initCondE;
}

/// ----------------------------  Frequency map operations --------- ///
void CSmileGUI::on_ButtonFreqMap_clicked()
{
    if(FreqMapRunning)
    {
        emit KillFreqMapThreads();
        return;
    }
    if(SchwModelRunning || ExportNbodyRunning || core->SchwThr!=NULL) return;
    CurveFreqMapR->setData(NULL, NULL, 0);
    CurveFreqMapC->setData(NULL, NULL, 0);
    // start building freq map
    core->BuildFreqMapIC(core->initCondE, ui.checkBox_fmExistingSS->isChecked());
    core->StartOrbitLibrary(SLOT(coreFreqMapFinished()));

    ui.progressBarFM->setRange(0, core->orbitlib->size());
    ui.progressBarFM->setValue(0);
    ui.progressBarFM->show();
    timeFM.start();
    ui.ButtonFreqMap->setText("Stop");
    ui.groupBox_N_dim->setEnabled(false);
    ui.groupBox_potential->setEnabled(false);
    for(size_t t=0; t<core->CalcManyThr.size(); t++)
        connect(this, SIGNAL(KillFreqMapThreads()), core->CalcManyThr[t], SLOT(stopThread()));
    FreqMapRunning=true;
}

void CSmileGUI::FreqMapFinished()
{
    FreqMapRunning=false;
    ui.ButtonFreqMap->setText("Start");
    ui.progressBarFM->hide();
    redrawFreqMap(true);
    info(/*ui.label_info->toPlainText() + "\n"+*/QString::number(core->orbitlib->numComplete())+" orbits, " +
        QString::number(core->orbitlib->numComplete()*1000.0/timeFM.elapsed(),'g',3)+" orbits/second\n" );
    ui.groupBox_N_dim->setEnabled(true);
    ui.groupBox_potential->setEnabled(true);
}

bool comparePaird(paird elem1, paird elem2)
{
    return elem1.second < elem2.second;
}

void CSmileGUI::redrawFreqMap(bool show)
{
    if(!show || core->orbitlib==NULL) return;
    unsigned int nFMpoints = core->orbitlib->size();
    if(nFMpoints==0) return;
    int numShell = ui.selectEnergyLevel_fm->value();
    const CBasicShellSchwModel* model_sh = core->model!=NULL && (core->model->ModelType() & CBasicSchwModel::MT_SHELL) ? static_cast<const CBasicShellSchwModel*>(core->model) : NULL;
    if(numShell>0 && model_sh!=NULL)
        ui.lineEdit_E->setText(QString::number(model_sh->getShellEnergy(numShell-1)));
    bool useNonzeroOnly = ui.checkBox_fmConfigSpace->isChecked();
    const CShellOrbitFilteringFnc filter(model_sh, numShell);
    //const CChaosOrbitFilteringFnc chaos(1, configCore.chaoticMinLambda, configCore.chaoticMinFreqDiff);
    if(ui.radioButton_fmmapdiff->isChecked() || ui.radioButton_fmmaplyapunov->isChecked())
    {   // draw frequency map itself
        vectord w1w3r(nFMpoints),  w2w3r(nFMpoints),  w1w3c(nFMpoints),  w2w3c(nFMpoints);
        int nr=0, nc=0;
        for(unsigned int i=0; i<nFMpoints; i++)
        {
            bool useOrbit=true;
            useOrbit = filter.eval(core->orbitlib->getOrbitDesc(i))>=0.5 && (!useNonzeroOnly || core->orbitlib->getOrbitDesc(i)->getWeight()>0);
            if(core->orbitlib->getOrbitDesc(i)->getState()==COrbitDesc::OS_DONE && useOrbit)
            {
                const COrbitInformation<float>* info=static_cast<const COrbitInformation<float>*>(core->orbitlib->getOrbitDesc(i)->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
                if(info==NULL) continue;
                bool criterion;
                if(ui.radioButton_fmmapdiff->isChecked()) 
                    criterion = info->getlfccdiff() <= configCore.chaoticMinFreqDiff;
                else 
                    criterion = info->getlambda() <= configCore.chaoticMinLambda;
                double denom = (core->potential->N_dim==3)?info->getlf(2):1;
                if(denom==0) continue;
                if(criterion)
                {
                    w1w3r[nr] = info->getlf(0)/denom;
                    w2w3r[nr] = info->getlf(1)/denom;
                    nr++;
                }
                else
                {
                    w1w3c[nc] = info->getlf(0)/denom;
                    w2w3c[nc] = info->getlf(1)/denom;
                    nc++;
                }
            }
        }
        CurveFreqMapR->setData(&(w1w3r.front()), &(w2w3r.front()), nr);
        CurveFreqMapC->setData(&(w1w3c.front()), &(w2w3c.front()), nc);
        if(core->potential->N_dim==3)
        {
            PlotFreqMap->setAxisScale(QwtPlot::yLeft, 0.5, 1.1);
            PlotFreqMap->setAxisScale(QwtPlot::xBottom, 0.3, 1.1);
        }
        PlotFreqMap->replot();
        PlotFreqMapHisto->hide();
        PlotFreqMapSS->hide();
        PlotFreqMap->show();

        if(!FreqMapRunning)
        {
            info(getOrbitPopulation(numShell, useNonzeroOnly));
        }
    }
    else if(ui.radioButton_fmhistodiff->isChecked() || ui.radioButton_fmhistolyapunov->isChecked())
    {   // draw histogram of frequency diffusion rates or Lyapunov exponents
        std::vector<paird> values;
        double totalWeight=0;
        for(unsigned int i=0; i<nFMpoints; i++)
            if(core->orbitlib->getOrbitDesc(i)->getState()==COrbitDesc::OS_DONE && filter.eval(core->orbitlib->getOrbitDesc(i))>=0.5)
            {
                const COrbitInformation<float>* info=static_cast<const COrbitInformation<float>*>(core->orbitlib->getOrbitDesc(i)->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
                if(info==NULL) continue;
                double weight=(useNonzeroOnly ? core->orbitlib->getOrbitDesc(i)->getWeight() : 0);
                totalWeight+=weight;
                values.push_back(paird(weight, ui.radioButton_fmhistodiff->isChecked() ? info->getlfccdiff() : info->getlambda()));
            }
        std::sort(values.begin(), values.end(), comparePaird);
        int npoints=static_cast<int>(values.size());
        if(npoints>0)
        {
            vectord val(npoints), np(npoints);
            double minval=1, tw=0;
            for(int i=0; i<npoints; i++)
            {
                val[i]=values[i].second;
                if(totalWeight)
                {
                    np[i]=tw/totalWeight;
                    tw+=values[i].first;
                }
                else
                    np[i]=i*1.0/npoints;
                if((val[i]>0) && (val[i]<minval)) minval=val[i];
            }
            CurveFreqMapHisto->setData(&(val.front()), &(np.front()), npoints);
            PlotFreqMapHisto->setAxisScale(QwtPlot::xBottom, minval/1.5, val[npoints-1]*1.1);
            PlotFreqMapHisto->replot();
        }
        PlotFreqMap->hide();
        PlotFreqMapSS->hide();
        PlotFreqMapHisto->show();
    }
    else if(ui.radioButton_fmssdiff->isChecked() || ui.radioButton_fmsslyapunov->isChecked() ||
        ui.radioButton_fmspdiff->isChecked() || ui.radioButton_fmsplyapunov->isChecked() ||
        ui.radioButton_fmsydiff->isChecked() || ui.radioButton_fmsylyapunov->isChecked())
    {   // draw start spaces
        vectord Xr(nFMpoints), Yr(nFMpoints), Xc(nFMpoints), Yc(nFMpoints);
        int nr=0, nc=0;
        double rad=0;
        COrbitInitData<float> InitData;
        for(unsigned int i=0; i<nFMpoints; i++)
        {
            if(core->orbitlib->getOrbitDesc(i)->getState()==COrbitDesc::OS_DONE && filter.eval(core->orbitlib->getOrbitDesc(i))>=0.5)
            {
                bool criterion;
                const COrbitInformation<float>* info=static_cast<const COrbitInformation<float>*>(core->orbitlib->getOrbitDesc(i)->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
                if(info==NULL) continue;
                if(ui.radioButton_fmssdiff->isChecked() || ui.radioButton_fmspdiff->isChecked() || ui.radioButton_fmsydiff->isChecked()) 
                    criterion = info->getlfccdiff() <= configCore.chaoticMinFreqDiff;
                else
                    criterion = info->getlambda() <= configCore.chaoticMinLambda;
                const COrbitInitData<float> InitData=core->orbitlib->getOrbitDesc(i)->getInitData();
                double X=InitData.initCond.Pos[0];
                double Y=InitData.initCond.Pos[1];
                double Z=InitData.initCond.Pos[2];
                double Xp, Yp;
                if(ui.radioButton_fmsydiff->isChecked() || ui.radioButton_fmsylyapunov->isChecked())
                {
                    Xp=Y;
                    Yp=atan2(InitData.initCond.Vel[0], InitData.initCond.Vel[2]);
                }
                else
                    XYZtoXY(X,Y,Z,Xp,Yp);
                if(pow_2(Xp)+pow_2(Yp)>pow_2(rad)) rad=sqrt(pow_2(Xp)+pow_2(Yp));
                if(ui.radioButton_fmssdiff->isChecked() || ui.radioButton_fmsslyapunov->isChecked())
                {
                    if(pow_2(InitData.initCond.Vel[0])+pow_2(InitData.initCond.Vel[1])+pow_2(InitData.initCond.Vel[2])>0) continue;   // this is not a stationary orbit
                }
                else if(ui.radioButton_fmspdiff->isChecked() || ui.radioButton_fmsplyapunov->isChecked())
                {
                    if(X*Y*Z!=0 || (X==0 && Z==0)) continue;   // this is not a principal-plane orbit
                }
                else
                {
                    if(X!=0 || Z!=0) continue;  // this point is not on Y axis
                }
                if(criterion)
                {
                    Xr[nr] = Xp;
                    Yr[nr] = Yp;
                    nr++;
                }
                else
                {
                    Xc[nc] = Xp;
                    Yc[nc] = Yp;
                    nc++;
                }
            }
        }
        CurveFreqMapSSR->setData(&(Xr.front()), &(Yr.front()), nr);
        CurveFreqMapSSC->setData(&(Xc.front()), &(Yc.front()), nc);
        PlotFreqMapSS->setAxisScale(QwtPlot::yLeft, -rad, rad);
        PlotFreqMapSS->setAxisScale(QwtPlot::xBottom, -rad, rad);
        PlotFreqMapSS->replot();
        PlotFreqMap->hide();
        PlotFreqMapHisto->hide();
        PlotFreqMapSS->show();
    }
}

void CSmileGUI::on_selectEnergyLevel_fm_valueChanged(int /*numShell*/)
{
    redrawFreqMap(true);
    //if(numShell>0 && core->model!=NULL && (core->model->ModelType() & CBasicSchwModel::MT_SHELL))  // redundant
    //    ui.lineEdit_E->setText(QString::number(static_cast<CBasicShellSchwModel*>(core->model)->getShellEnergy(numShell-1)));
}

void CSmileGUI::on_checkBox_fmConfigSpace_stateChanged(int)
{
    redrawFreqMap(true);
}

QString CSmileGUI::getOrbitPopulation(int numShell, bool useWeights)
{
    if(core->orbitlib==NULL) return "Error!";
    // if model is initialized and numShell>0, then filter orbits according to energy shells they belong to
    const CBasicShellSchwModel* model_sh = core->model!=NULL && (core->model->ModelType() & CBasicSchwModel::MT_SHELL) ? static_cast<const CBasicShellSchwModel*>(core->model) : NULL;
    const CShellOrbitFilteringFnc filter(model_sh, numShell);
    const CChaosOrbitFilteringFnc chaos(1, configCore.chaoticMinLambda, configCore.chaoticMinFreqDiff);
    return QString::fromStdString(core->orbitlib->getOrbitPopulation(&filter, &chaos, useWeights));
}

void CSmileGUI::drawFreqMapResonances()
{
    const int NR=18;   // number of resonances depicted on the map
    double res[NR][5] = {
        {0,1,-1 ,0.80,1},
        {1,-1,0 ,0.55,0.55},
        {2,-1,0 ,0.45,0.90},
        {3,-2,0 ,0.34,0.51},
        {4,-3,0 ,0.39,0.52},
        {2,0,-1 ,0.50,1.08},
        {3,0,-2 ,0.66,0.52},
        {1,-3,2 ,0.34,0.78},
        {1,-2,1 ,0.80,0.90},
        {2,-3,1 ,0.33,0.55},
        {4,-2,-1,0.78,1.08},
        {3,-1,-1,0.69,1.08},
        {3,1,-3 ,0.63,1.08},
        {2,1,-2 ,0.74,0.52},
        {2,2,-3 ,0.80,0.70},
        {1,2,-2 ,0.74,0.63},
        {1,1,-1 ,0.47,0.53} };
    double xmin=0.3, xmax=1., ymin=0.4, ymax=1.1;
    double boundrect[4][3] = { {1, 0, -xmin}, {1, 0, -xmax}, {0, 1, -ymin}, {0, 1, -ymax} };
    //QString S;
    for(int i=0; i<NR; i++)
    {
        double x[2], y[2];
        if(res[i][1]==0)
        {
            x[0]=x[1]=-res[i][2]/res[i][0];
            y[0]=ymin;
            y[1]=ymax;
        }
        else if(res[i][0]==0)
        {
            x[0]=xmin;
            x[1]=xmax;
            y[0]=y[1]=-res[i][2]/res[i][1];
        }
        else
        {   // intersect with bounding rectangle
            std::vector<double> ax(0), ay(0);
            for(int c=0; c<4; c++)
            {
                double D = res[i][0]*boundrect[c][1] - res[i][1]*boundrect[c][0];
                ax.push_back( (-res[i][2]*boundrect[c][1] + res[i][1]*boundrect[c][2]) / D);
                ay.push_back(-(-res[i][2]*boundrect[c][0] + res[i][0]*boundrect[c][2]) / D);
            }
            std::sort(ax.begin(), ax.end());
            std::sort(ay.begin(), ay.end());
            int sgn = (res[i][0]*res[i][1]>0)?1:0;
            x[0]=ax[1]; x[1]=ax[2];
            y[0]=ay[1+sgn]; y[1]=ay[2-sgn];
        }
        QwtPlotCurve* c = new QwtPlotCurve();
        c->setPen(QPen(Qt::darkGray, 1, Qt::DotLine));
        c->setData(x, y, 2);
        c->attach(PlotFreqMap);
        QwtPlotMarker* m = new QwtPlotMarker();
        m->setValue(res[i][3], res[i][4]);
        m->setLabel(QString::number(res[i][0])+","+QString::number(res[i][1])+","+QString::number(res[i][2]));
        m->attach(PlotFreqMap);
    }
}

void CSmileGUI::FMclicked(const QPointF &pos)
{
    if(core->orbitlib==NULL || core->orbitlib->size()==0)
        return;
    int mode=-1;
    if(ui.radioButton_fmmapdiff->isChecked() || ui.radioButton_fmmaplyapunov->isChecked()) mode=0;
    else if(ui.radioButton_fmssdiff->isChecked() || ui.radioButton_fmsslyapunov->isChecked()) mode=1;
    else if(ui.radioButton_fmspdiff->isChecked() || ui.radioButton_fmsplyapunov->isChecked()) mode=2;
    if(mode<0) return;
    double px=pos.x();
    double py=pos.y();
    double mindist=1; 
    int np=0;
    int numShell = ui.selectEnergyLevel_fm->value();
    const CBasicShellSchwModel* model_sh = core->model!=NULL && (core->model->ModelType() & CBasicSchwModel::MT_SHELL) ? static_cast<const CBasicShellSchwModel*>(core->model) : NULL;
    const CShellOrbitFilteringFnc filter(model_sh, numShell);
    for(unsigned int i=0; i<core->orbitlib->size(); i++)
    {
        const COrbitInformation<float>* info=static_cast<const COrbitInformation<float>*>(core->orbitlib->getOrbitDesc(i)->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
        if(info==NULL) continue;
        double dist=1;
        if(mode==0)
        {
            double lf2 = (core->potential->N_dim==3)?info->getlf(2):1;
            double w1w3o = info->getlf(0)/lf2;
            double w2w3o = info->getlf(1)/lf2;
            dist = pow_2(w1w3o-px) + pow_2(w2w3o-py);
        }
        else
        {
            const COrbitInitData<float> InitData=core->orbitlib->getOrbitDesc(i)->getInitData();
            double X=InitData.initCond.Pos[0];
            double Y=InitData.initCond.Pos[1];
            double Z=InitData.initCond.Pos[2];
            double Xp, Yp;
            XYZtoXY(X,Y,Z,Xp,Yp);
            if((mode==1 && X*Y*Z==0) || (mode==2 && X*Y*Z!=0)) continue; 
            dist = pow_2(Xp-px) + pow_2(Yp-py);
        }
        if(dist<mindist && filter.eval(core->orbitlib->getOrbitDesc(i))>=0.5)
        {
            np=i;
            mindist=dist;
        }
    }
    if(mindist>=1) return;  // too far
    // load data from np'th orbit
    bool tmpLyapunov=core->InitData.calcLyapunov;
    core->InitData=core->orbitlib->getOrbitDesc(np)->getInitData();
    core->InitData.calcLyapunov=tmpLyapunov;
    // set exact values
    bool tmpUseE=core->useICe;
    core->useICe=false;  // need to put in values for IC, not for E
    core->initIC();
    core->useICe=tmpUseE;
    InternalRecalc=true;
    ui.lineEdit_x->setText(QString::number(core->InitData.initCond.Pos[0]));
    ui.lineEdit_y->setText(QString::number(core->InitData.initCond.Pos[1]));
    ui.lineEdit_z->setText(QString::number(core->InitData.initCond.Pos[2]));
    ui.lineEdit_vx->setText(QString::number(core->InitData.initCond.Vel[0]));
    ui.lineEdit_vy->setText(QString::number(core->InitData.initCond.Vel[1]));
    ui.lineEdit_vz->setText(QString::number(core->InitData.initCond.Vel[2]));
    ui.lineEdit_E->setText(QString::number(core->initCondE));
    ui.lineEdit_Torb->setText(QString::number(core->InitData.timeUnit));
    InternalRecalc=false;
    info(QString::fromStdString(core->orbitlib->getOrbitDesc(np)->toString()));
}

void CSmileGUI::XYZtoXY(const double X, const double Y, const double Z, double& x, double& y)
{
    x = (Y-X)/sqrt(2.0);
    y = (2*Z-X-Y)/sqrt(6.0);
}

/// ----------------------------  Schwarzschild modelling --------- ///
void CSmileGUI::on_ButtonSchwModel_clicked()
{
    if(FreqMapRunning || ExportNbodyRunning) return;
    if(core->SchwThr!=NULL)  // stop
    {
        error("Cannot stop optimization routine");
        return;
    }
    if(SchwModelRunning)
    {
        emit KillFreqMapThreads();
        return;
    }
    if(core->potential->N_dim!=3)
        return;

    // start building library
    core->BuildSchwIC(ui.checkBox_smExistingSS->isChecked());
    core->StartOrbitLibrary(SLOT(coreSchwOrbitLibraryFinished()));
    ui.groupBox_N_dim->setEnabled(false);
    ui.groupBox_potential->setEnabled(false);
    ui.groupBox_smOptim->setEnabled(false);
    ui.groupBox_smGrid->setEnabled(false);
    ui.ButtonSchwModel->setText("Stop");
    timeFM.start();
    ui.progressBarSM->setRange(0, core->orbitlib->size());
    ui.progressBarSM->setValue(0);
    ui.progressBarSM->show();
    SchwModelRunning=true;
    for(size_t t=0; t<core->CalcManyThr.size(); t++)
        connect(this, SIGNAL(KillFreqMapThreads()), core->CalcManyThr[t], SLOT(stopThread()));
}

void CSmileGUI::SchwOrbitLibraryFinished()    // slot called when core->signalSchwOrbitLibraryFinished() emitted
{
    SchwModelRunning=false;
    ui.ButtonSchwModel->setText("Start");
    ui.progressBarSM->hide();
    ui.groupBox_N_dim->setEnabled(true);
    ui.groupBox_potential->setEnabled(true);
    ui.groupBox_smOptim->setEnabled(true);
    //ui.groupBox_smOptimParams->setEnabled(true);///!!!
    ui.groupBox_smGrid->setEnabled(true);
    ui.ButtonSchwLP->setEnabled(true);
    ui.ButtonSchwQP->setEnabled(true);
    //ui.ButtonSchwLucy->setEnabled(true);
    //ui.ButtonSchwExportGrid->setEnabled(true);
    ui.ButtonSchwNbody->setEnabled(true);
    redrawFreqMap(true);
    info(QString::number(core->orbitlib->numComplete())+" orbits\n" +
        QString::number(core->orbitlib->numComplete()*1000.0/timeFM.elapsed(),'g',3)+" orbits/second");
}

void CSmileGUI::NbodyExportFinished(const QString& result)    // slot called when core->signalSchwExportNbodyFinished() emitted
{
    ExportNbodyRunning=false;
    //ui.tabWidget_sm->setEnabled(true);
    ui.ButtonSchwNbody->setText("Create Nbody IC");
    ui.progressBarSM->hide();
    info(result);
}

void CSmileGUI::on_ButtonSchwNbody_clicked()
{
    if(SchwModelRunning || FreqMapRunning || core->orbitlib==NULL || core->orbitlib->size()==0)
        return;
    if(ExportNbodyRunning)
    {   // stop export (determine which stage we are at now)
        if(core->CalcManyThr.size()>0)
        {
            for(size_t t=0; t<core->CalcManyThr.size(); t++)
                connect(this, SIGNAL(KillFreqMapThreads()), core->CalcManyThr[t], SLOT(stopThread()));
            emit KillFreqMapThreads();
        }
        if(core->NBexportThr!=NULL)  // stop
        {
            QMetaObject::invokeMethod(core->NBexportThr, "stopThread");
        }
        return;
    }
    bool ok;
    int NPoints=QInputDialog::getInteger(this, "Export to N-body", "Enter number of points", 1000000, 0, 1e8, 1, &ok);
    if(!ok) return;
    QString fileName = QFileDialog::getSaveFileName(this, "Export to N-body", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    int refineFactor=QInputDialog::getInteger(this, "Export to N-body", "Mass refine factor (2^N)", 0, 0, 6, 1, &ok);
    if(!ok) refineFactor=0;

    ExportNbodyRunning=true;
    //ui.tabWidget_sm->setEnabled(false);
    ui.ButtonSchwNbody->setText("Cancel");
    ui.progressBarSM->setValue(0);
    ui.progressBarSM->setRange(0, core->orbitlib->size());
    ui.progressBarSM->show();
    core->SchwExportNbody(NPoints, fileName, refineFactor);
}

void CSmileGUI::on_ButtonSchwLP_clicked()
{
    if(core->SchwThr!=NULL || core->model==NULL) return;
    ui.ButtonSchwLP->setEnabled(false);
    ui.ButtonSchwQP->setEnabled(false);
    ui.ButtonSchwLucy->setEnabled(false);
    //ui.groupBox_smOptimParams->setEnabled(false);
    core->SchwStartOptimization(false);
}

void CSmileGUI::on_ButtonSchwQP_clicked()
{
    if(core->SchwThr!=NULL || core->model==NULL) return;
    ui.ButtonSchwLP->setEnabled(false);
    ui.ButtonSchwQP->setEnabled(false);
    ui.ButtonSchwLucy->setEnabled(false);
    //ui.groupBox_smOptimParams->setEnabled(false);
    core->SchwStartOptimization(true);
}

void CSmileGUI::on_ButtonSchwLucy_clicked()
{/*
    if(core->model==NULL) return;
    if(core->SchwThr!=NULL)  // stop
        emit stopSchw();
    else
    {
         ///!!!connect(this, SIGNAL(stopSchw()), core->model, SLOT(stop()));
        ui.radioButton_smviewweights->setChecked(true);
        ui.ButtonSchwLP->setEnabled(false);
        ui.ButtonSchwQP->setEnabled(false);
        ui.ButtonSchwLucy->setText("Stop");
        //ui.groupBox_smOptimParams->setEnabled(false);
        core->SchwStartOptimization(CSchwarzschildThread::SM_LUCY);
    }*/
}

void CSmileGUI::SchwOptimizationFinished(const QString& message)
{
    ui.ButtonSchwLP->setEnabled(true);
    ui.ButtonSchwQP->setEnabled(configCore.sm_useBPMPD);
    //ui.ButtonSchwLucy->setEnabled(true);
    //ui.ButtonSchwLucy->setText("Lucy");
    //ui.groupBox_smOptimParams->setEnabled(true);
    ui.groupBox_smOptim->setEnabled(true);
    SchwModelResult = message;
    redrawSchwModel(true);
}

void CSmileGUI::on_selectEnergyLevel_valueChanged(int)
{
    redrawSchwModel(true);
}

void CSmileGUI::redrawSchwModel(bool updateInfo)
{
    if(!core->model || !core->orbitlib)
        return;
    // init one or more type-specific pointers below if model type is one of known types
    const CBasicShellSchwModel* model_sh = (core->model->ModelType() & CBasicSchwModel::MT_SHELL) ? static_cast<const CBasicShellSchwModel*>(core->model) : NULL;
    const CSchwModelClassic* model_cl = (core->model->ModelType()==CBasicSchwModel::MT_CLASSIC) ? static_cast<const CSchwModelClassic*>(core->model) : NULL;
    const CSchwModelSH* model_sph = (core->model->ModelType() & CBasicSchwModel::MT_SPHERICAL_HARMONIC) ? static_cast<const CSchwModelSH*>(core->model) : NULL;
    const CSchwModelSHGrid* model_shg = (core->model->ModelType()==CBasicSchwModel::MT_SHGRID) ? static_cast<const CSchwModelSHGrid*>(core->model) : NULL;
    int numShell=ui.selectEnergyLevel->value();
    if(model_cl) 
    {
        ui.selectEnergyLevel->setEnabled(true);
    }
    else
    {
        ui.selectEnergyLevel->setEnabled(true);
        numShell=0;
    }
    if(numShell>0 && (core->model->ModelType() & CBasicSchwModel::MT_SHELL))
        ui.lineEdit_E->setText(QString::number(static_cast<const CBasicShellSchwModel*>(core->model)->getShellEnergy(numShell-1)));
    if(ui.radioButton_smviewgrid->isChecked() && model_sh!=NULL)
    {
        vectord CellsFx, CellsFy, CellsIx, CellsIy;
        vectord coefsPenalty(model_sh->calcPenalty(core->orbitlib));
        unsigned int nstart =(numShell==0 || model_cl==NULL ? 0 : (numShell-1)*model_cl->getCellsPerShell());
        unsigned int nend   =(numShell==0 || model_cl==NULL  ? model_sh->getNumCoefsDens() : numShell*model_cl->getCellsPerShell());
        double Xmin=0, Xmax=0, Ymin=0, Ymax=0;
        unsigned int cmin=model_shg!=NULL ? 1 : 0;
        for(unsigned int c=nstart; c<nend; c++)
        {
            double Xp=0, Yp=0;
            if(model_cl)
            {
                double X=0,Y=0,Z=0;
                model_cl->cellCenter(c, X, Y, Z);
                XYZtoXY(X,Y,Z,Xp,Yp);
                if(fabs(Xp)>Xmax) { Xmax=fabs(Xp); Xmin=-Xmax; }
                if(Yp>Ymax) { Ymax=Yp; }
                if(Yp<Ymin) { Ymin=Yp; }
            } 
            else if(model_sph)
            {
                if(c>=cmin)
                {
                    int indAng=static_cast<int>((c-cmin)%model_sph->getNumCoefsAtRadius());
                    int l=static_cast<int>(sqrt(2.*indAng-floor(sqrt(2.*indAng))));
                    int m=indAng-l*(l+1)/2;
                    Xp=l*(l+1)/2 + m*(1.-1./(l+4));  // put the points in row but make somewhat larger spacing between l and l+1 series
                    Yp=static_cast<int>((c-cmin)/model_sph->getNumCoefsAtRadius())+cmin;
                    Xmax=std::max<double>(Xmax,Xp);
                    Ymax=std::max<double>(Ymax,Yp);
                }
            }
            if(coefsPenalty[c]!=0)
            {
                CellsIx.push_back(Xp);
                CellsIy.push_back(Yp);
            }
            else
            {
                CellsFx.push_back(Xp);
                CellsFy.push_back(Yp);
            }
        }
        if(CellsFx.empty())
            CurveSchwCellF->setData(NULL, NULL, 0);
        else
            CurveSchwCellF->setData(&(CellsFx.front()), &(CellsFy.front()), CellsFx.size());
        if(CellsIx.empty())
            CurveSchwCellI->setData(NULL, NULL, 0);
        else
            CurveSchwCellI->setData(&(CellsIx.front()), &(CellsIy.front()), CellsIx.size());
        PlotSchw->setAxisScale(QwtPlot::xBottom, Xmin, Xmax);
        PlotSchw->setAxisScale(QwtPlot::yLeft, Ymin, Ymax);
        PlotSchw->replot();
        PlotSchwHisto->hide();
        PlotSchw->show();
        if(updateInfo)  info(SchwModelResult + "\n" + 
            getOrbitPopulation(numShell, true));
    }
    else if(ui.radioButton_smviewbeta->isChecked())
    {   // plot anisotropy coefficient beta
        if(model_sh)
        {
            unsigned int npoints=(numShell==0 ? model_sh->getNumShellsKinem() : 1);
            unsigned int nstart =(numShell==0 ? 0 : numShell-1);
            unsigned int nend   =(numShell==0 ? model_sh->getNumShellsKinem() : nstart+1);
            vectord val(npoints), np(npoints);
            const CShellOrbitFilteringFnc filter(model_sh, numShell);
            for(unsigned int s=nstart; s<nend; s++)
            {
                double vr2=0, vt2=0;
                for(unsigned int o=0; o<core->orbitlib->size(); o++)
                {
                    const COrbitDesc* orb=core->orbitlib->getOrbitDesc(o);
                    if(orb->getState()==COrbitDesc::OS_DONE && 
                        filter.eval(orb)>=0.5 &&
                        orb->getWeight()>0)
                    {
                        // load grid information data from orbit desc
                        const CSchwInformation* infoGrid=static_cast<const CSchwInformation*>(orb->getInfoByType(CBasicInformation::IT_SCHW));
                        if(infoGrid!=NULL)  ///!!! other types of shell model??
                        {
                            vr2 += infoGrid->getShellVr(s) * infoGrid->getShellTime(s) * orb->getWeight();
                            vt2 += infoGrid->getShellVt(s) * infoGrid->getShellTime(s) * orb->getWeight();
                        }
                    }
                }
                if(vr2>0)
                    val[s-nstart] = 1-vt2/vr2/2;
                else
                    val[s-nstart] = 1;
                np[s-nstart]=s*1.0;
            }
            if(npoints>0)
                CurveSchwHistoR->setData(&(np.front()), &(val.front()), npoints);
            else
                CurveSchwHistoR->setData(NULL, NULL, 0);
            CurveSchwHistoR->setStyle(QwtPlotCurve::Steps);
            CurveSchwHistoC->setData(NULL, NULL, 0);
        } else {
            CurveSchwHistoR->setData(NULL, NULL, 0);
            CurveSchwHistoC->setData(NULL, NULL, 0);
        }
        PlotSchwHisto->setAxisScale(QwtPlot::yLeft, -1, 1);
        PlotSchwHisto->setAxisAutoScale(QwtPlot::xBottom);
        PlotSchwHisto->setAxisScaleEngine(QwtPlot::yLeft, new QwtLinearScaleEngine);
        PlotSchwHisto->replot();
        PlotSchw->hide();
        PlotSchwHisto->show();
        if(updateInfo)  info(SchwModelResult + "\n" + 
            getOrbitPopulation(numShell, true));
    }
    else if(ui.radioButton_smviewweights->isChecked() || ui.radioButton_smviewhisto->isChecked())
    {   // plot histogram of orbit weights
        vectord valuesR, valuesC, npR, npC;
        double minval=1, maxval=0;
        const CShellOrbitFilteringFnc filter(model_sh, numShell);
        const CChaosOrbitFilteringFnc chaos(1, configCore.chaoticMinLambda, configCore.chaoticMinFreqDiff);
        for(size_t o=0; o<core->orbitlib->size(); o++)
        {
            const COrbitDesc* orb=core->orbitlib->getOrbitDesc(o);
            if(orb->getState()==COrbitDesc::OS_DONE && filter.eval(orb)>=0.5)
            {
                if(chaos.eval(orb)>=0.5 && !ui.radioButton_smviewhisto->isChecked())
                {
                    valuesC.push_back(orb->getWeight());
                    npC.push_back(o*1.0);
                }
                else
                {
                    valuesR.push_back(orb->getWeight());
                    npR.push_back(o*1.0);
                }
                if(orb->getWeight()>0) minval = std::min<double>(minval, orb->getWeight());
                maxval = std::max<double>(maxval, orb->getWeight());
            }
        }
        if(ui.radioButton_smviewhisto->isChecked()) 
        {
            std::sort(valuesR.begin(), valuesR.end());
            CurveSchwHistoC->setData(NULL, NULL, 0);
        }
        else
        {
            if(npC.size()==0)
                CurveSchwHistoC->setData(NULL, NULL, 0);
            else
                CurveSchwHistoC->setData(&(npC.front()), &(valuesC.front()), valuesC.size());
        }
        if(npR.size()==0)
            CurveSchwHistoR->setData(NULL, NULL, 0);
        else
            CurveSchwHistoR->setData(&(npR.front()), &(valuesR.front()), valuesR.size());
        CurveSchwHistoR->setStyle(ui.radioButton_smviewhisto->isChecked() ? QwtPlotCurve::Steps : QwtPlotCurve::Dots);
        PlotSchwHisto->setAxisScale(QwtPlot::yLeft, minval, maxval);
        //PlotSchwHisto->setAxisAutoScale(QwtPlot::xBottom);
        PlotSchwHisto->setAxisScale(QwtPlot::xBottom, 0, core->orbitlib->size());
        PlotSchwHisto->setAxisScaleEngine(QwtPlot::yLeft, new QwtLog10ScaleEngine);
        PlotSchwHisto->replot();
        PlotSchw->hide();
        PlotSchwHisto->show();
        if(updateInfo)  info(SchwModelResult + "\n" + 
            getOrbitPopulation(numShell, true));
    }
}

void CSmileGUI::SMclicked(const QPointF& pos)
{
    if(core->orbitlib==NULL || core->orbitlib->size()==0 || core->model==NULL || (core->model->ModelType() & CBasicSchwModel::MT_SHELL) == 0)
        return;
    double px=pos.x();
    double py=pos.y();
    double mindist=1; 

    const CBasicShellSchwModel* model_sh =static_cast<const CBasicShellSchwModel*>(core->model);
    const CSchwModelClassic* model_cl = (core->model->ModelType()==CBasicSchwModel::MT_CLASSIC) ? static_cast<const CSchwModelClassic*>(core->model) : NULL;
    const CSchwModelSH* model_sph = (core->model->ModelType() & CBasicSchwModel::MT_SPHERICAL_HARMONIC) ? static_cast<const CSchwModelSH*>(core->model) : NULL;
    const CSchwModelSHGrid* model_shg = (core->model->ModelType()==CBasicSchwModel::MT_SHGRID) ? static_cast<const CSchwModelSHGrid*>(core->model) : NULL;
    unsigned int cmin=model_shg!=NULL ? 1 : 0;
    int nc=-1;
    for(unsigned int c=0; c<model_sh->getNumCoefsDens(); c++)
    {
        double Xp=0, Yp=0;
        if(model_cl)
        {
            double X=0,Y=0,Z=0;
            model_cl->cellCenter(c, X, Y, Z);
            XYZtoXY(X,Y,Z,Xp,Yp);
        } 
        else if(model_sph)
        {
            if(c>=cmin)
            {
                int indAng=static_cast<int>((c-cmin)%model_sph->getNumCoefsAtRadius());
                int l=static_cast<int>(sqrt(2.*indAng-floor(sqrt(2.*indAng))));
                int m=indAng-l*(l+1)/2;
                Xp=l*(l+1)/2 + m*(1.-1./(l+4));  // put the points in row but make somewhat larger spacing between l and l+1 series
                Yp=static_cast<int>((c-cmin)/model_sph->getNumCoefsAtRadius())+cmin;
            }
        }
        double dist = pow_2(Xp-px) + pow_2(Yp-py);
        if(dist<mindist)
        {
            nc=c;
            mindist=dist;
        }
    }
    if(nc<0) return;
    // nc= index of cell/coefficient that was clicked
    /*int numorbits=0, numused=0;
    double mass=0;
    for(size_t o=0; o<core->orbitlib->OrbitList.size(); o++)
        if(core->orbitlib->OrbitList[o].celltimes[nc]>0)
        {
            numorbits++;
            if(core->orbitlib->OrbitList[o].orbitWeight > SCHW_MIN_ORBIT_WEIGHT_USED/core->orbitlib->OrbitList.size() ) numused++;
            mass += core->orbitlib->OrbitList[o].orbitWeight * core->orbitlib->OrbitList[o].celltimes[nc];
        }
    info(QString::number(numorbits) + " orbits visit cell,\n" +
        QString::number(numused) + " orbits used in model,\n" +
        "Mass difference: " + QString::number((mass/core->model->cellMass[nc]-1) * 100, 'f', 4) + "%\n" +
        getOrbitPopulation(0, nc, nc));
    */
    double coefValue  = model_sh->getCoefsDens().at(nc);
    double coefPenalty= model_sh->calcPenalty(core->orbitlib).at(nc);
    double normFactor = model_sh->getNormFactor(nc);
    info("Coef index="+QString::number(nc) + "\nRequired value="+QString::number(coefValue) + 
        "\nDifference="+QString::number(coefPenalty) + "\nNormalization="+QString::number(normFactor));
}

void CSmileGUI::SMHistoClicked(const QPointF &pos)
{
    if(core->orbitlib==NULL || core->orbitlib->size()==0 || !ui.radioButton_smviewweights->isChecked())
        return;
    int numShell = ui.selectEnergyLevel_fm->value();
    const CBasicShellSchwModel* model_sh = core->model!=NULL && (core->model->ModelType() & CBasicSchwModel::MT_SHELL) ? static_cast<const CBasicShellSchwModel*>(core->model) : NULL;
    const CShellOrbitFilteringFnc filter(model_sh, numShell);
    double px=pos.x();
    double py=pos.y();
    // scale clicked value also to [0..1] interval
    double minyval=1, maxyval=0;
    unsigned int maxxval=core->orbitlib->size(), minxval=0;
    for(unsigned int i=0; i<core->orbitlib->size(); i++)
        if(core->orbitlib->getOrbitDesc(i)->getState()==COrbitDesc::OS_DONE && filter.eval(core->orbitlib->getOrbitDesc(i))>=0.5)
        {
            //if(i<minxval) minxval=i;
            //if(i>maxxval) maxxval=i;
            float weight=core->orbitlib->getOrbitDesc(i)->getWeight();
            if(weight>0)
            {
                minyval=std::min<double>(minyval, weight);
                maxyval=std::max<double>(minyval, weight);
            }
        }
    if(maxxval==0 || maxyval==0) return;
    minyval=log10(minyval);
    maxyval=log10(maxyval);
    px = (px-minxval)*1.0/(maxxval-minxval);
    py = (log10(py)-minyval)/(maxyval-minyval);
    // find nearest point
    double mindist=1; 
    int np=0;
    for(unsigned int i=0; i<core->orbitlib->size(); i++)
    {
        double dist=1;
        if(core->orbitlib->getOrbitDesc(i)->getWeight()>0 && core->orbitlib->getOrbitDesc(i)->getState()==COrbitDesc::OS_DONE && filter.eval(core->orbitlib->getOrbitDesc(i))>=0.5)
        {
            double ox=(i-minxval)*1.0/(maxxval-minxval);
            double oy=(log10(core->orbitlib->getOrbitDesc(i)->getWeight())-minyval)/(maxyval-minyval);
            dist = pow_2(ox-px) + pow_2(oy-py);
        }
        if(dist<mindist)
        {
            np=i;
            mindist=dist;
        }
    }
    if(mindist>=1) return;  // too far
    // load data from np'th orbit
    bool tmpLyapunov=core->InitData.calcLyapunov;
    core->InitData=core->orbitlib->getOrbitDesc(np)->getInitData();
    core->InitData.calcLyapunov=tmpLyapunov;
    // set exact values
    bool tmpUseE=core->useICe;
    core->useICe=false;  // need to put in values for IC, not for E
    core->initIC();
    core->useICe=tmpUseE;
    InternalRecalc=true;
    ui.lineEdit_x->setText(QString::number(core->InitData.initCond.Pos[0]));
    ui.lineEdit_y->setText(QString::number(core->InitData.initCond.Pos[1]));
    ui.lineEdit_z->setText(QString::number(core->InitData.initCond.Pos[2]));
    ui.lineEdit_vx->setText(QString::number(core->InitData.initCond.Vel[0]));
    ui.lineEdit_vy->setText(QString::number(core->InitData.initCond.Vel[1]));
    ui.lineEdit_vz->setText(QString::number(core->InitData.initCond.Vel[2]));
    ui.lineEdit_E->setText(QString::number(core->initCondE));
    ui.lineEdit_Torb->setText(QString::number(core->InitData.timeUnit));
    InternalRecalc=false;
    info(QString::fromStdString(core->orbitlib->getOrbitDesc(np)->toString()));
}

///-------------  export-import -------------///
void CSmileGUI::on_ButtonExportFM_clicked()
{
    if(core->orbitlib==NULL || core->orbitlib->size()==0)
        return;
    QString fileName = QFileDialog::getSaveFileName(this, "Export Frequency map", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);

    if(!core->exportOrbitLib(fileName, false))
        error("Cannot write file "+fileName);
    core->saveSettings(fileName+".ini");
}

void CSmileGUI::on_ButtonImportFM_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Import Frequency map", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    loadApplySettings(fileName+".ini");   // try loading saved potential and other data
    if (!core->importOrbitLib(fileName, false)) {
        error("Cannot read file "+fileName);
        return;
    }
    redrawFreqMap(true);
    info(QString::number(core->orbitlib->size())+" orbits loaded");
}

void CSmileGUI::on_ButtonExportSM_clicked()
{
    if(core->orbitlib==NULL || core->orbitlib->size()==0 || core->model==NULL)
        return;
    QString fileName = QFileDialog::getSaveFileName(this, "Export Schwarzschild model", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);

    if(!core->exportOrbitLib(fileName, true))
        error("Cannot write file "+fileName);
    core->saveSettings(fileName+".ini");
}

void CSmileGUI::on_ButtonImportSM_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Import Schwarzschild model", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    loadApplySettings(fileName+".ini");   // try loading saved potential and other data
    if(core->potential->N_dim!=N_DIM)
    {
        error("N_dim should be equal to 3 for Schwarzschild modelling");
        return;
    }
    if(!core->importOrbitLib(fileName, true)) 
    {
        error("Cannot read file "+fileName);
        return;
    }
    core->initSchwModel();
    ui.groupBox_smOptim->setEnabled(true);
    //ui.groupBox_smOptimParams->setEnabled(true);
    ui.ButtonSchwLP->setEnabled(true);
    ui.ButtonSchwQP->setEnabled(configCore.sm_useBPMPD);
    ui.ButtonSchwLucy->setEnabled(false);  //< disabled
    //ui.ButtonSchwExportGrid->setEnabled(true);
    ui.ButtonSchwNbody->setEnabled(true);
    redrawFreqMap(true);
    on_selectEnergyLevel_valueChanged(ui.selectEnergyLevel->value());  // redraw grid
    SchwModelResult = QString::fromStdString(core->model->getStatistics(core->orbitlib));
    info(QString::number(core->orbitlib->size())+" orbits loaded\n" + SchwModelResult);
}

void CSmileGUI::on_ButtonExportOrbit_clicked()
{
    if(!core->orbit || !core->orbit->getIntTime())
        return;
    QString fileName = QFileDialog::getSaveFileName(this, "Export orbit", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    if(!core->exportOrbit(fileName))
        error("Cannot write file "+fileName);
}

void CSmileGUI::on_ButtonImportOrbit_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Import orbit", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    if(!core->importOrbit(fileName)) {
        error("Cannot read file "+fileName);
        return;
    }
    InternalRecalc=true;
    ui.lineEdit_x->setText(QString::number(core->InitData.initCond.Pos[0]));
    ui.lineEdit_y->setText(QString::number(core->InitData.initCond.Pos[1]));
    ui.lineEdit_z->setText(QString::number(core->InitData.initCond.Pos[2]));
    ui.lineEdit_vx->setText(QString::number(core->InitData.initCond.Vel[0]));
    ui.lineEdit_vy->setText(QString::number(core->InitData.initCond.Vel[1]));
    ui.lineEdit_vz->setText(QString::number(core->InitData.initCond.Vel[2]));
    ui.lineEdit_E->setText(QString::number(core->initCondE));
    ui.lineEdit_intTime->setText(QString::number(configCore.intTimeInPeriods));
    ui.lineEdit_timeStep->setText(QString::number(configCore.intTimeStepsPerPeriod));
    ui.spinBox_intervalCount->setValue(1);
    InternalRecalc=false;
    refreshOrbitAndSpectrum();
}

void CSmileGUI::on_ButtonExportModel_clicked()
{
    if(core->model==NULL || (core->model->ModelType() & CBasicSchwModel::MT_SHELL) == 0)
        return;
    QString fileName = QFileDialog::getSaveFileName(this, "Export Schwarzschild model", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    if(!core->exportSchwModel(fileName))
        error("Cannot write file "+fileName);
}

void CSmileGUI::on_ButtonExportPS_clicked()
{
    if(CurvesPoincare.size()<=1) return;
    QString fileName = QFileDialog::getSaveFileName(this, "Export Schwarzschild model", core->WorkDir);
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    std::ofstream strm(fileName.toStdString().c_str(), std::ios::out);
    if(!strm) { error("Cannot write file "+fileName); return; }
    std::vector< paird > curve;
    for(int c=1; c<CurvesPoincare.size(); c++)
    {
        curve.resize(CurvesPoincare[c]->data().size());
        for(size_t i=0; i<CurvesPoincare[c]->data().size(); i++)
        {
            curve[i].first =CurvesPoincare[c]->data().x(i);
            curve[i].second=CurvesPoincare[c]->data().y(i);
        }
        // sort points in curve in nearest-neighbour order
        for(size_t i=0; i<curve.size()-1; i++)
        {
            double distmin=1e10; 
            size_t indmin=0;
            // find next unused point closest to i-th 
            for(size_t j=i+1; j<curve.size(); j++)
            {
                double dist=pow_2(curve[i].first-curve[j].first)+pow_2(curve[i].second-curve[j].second);
                if(dist<distmin) { distmin=dist; indmin=j; }
            }
            if(indmin>0)
            {
                paird tmp(curve[indmin]);
                curve[indmin]=curve[i+1];
                curve[i+1]=tmp;
            }
        }
        for(size_t i=0; i<curve.size()-1; i++)
            strm << curve[i].first << "\t" << curve[i].second << "\n";
        strm << "\n";
    }
}

void CSmileGUI::on_ButtonPrint_clicked()
{
    /*for(int c=1; c<=ui.spinBox_intervalCount->value(); c++)
    {
        ui.spinBox_IntervalNum->setValue(c);
        QImage P(364,355,QImage::Format_RGB888);
        QwtPlotPrintFilter F=QwtPlotPrintFilter();
        PlotOrbit->print(P, F);
        P.save(QString::number(c+1000)+".png");
    }
    return;*/
    QwtPlot* fig = NULL;  
    if(PlotOrbit->isVisible())
        fig = PlotOrbit;
    if(PlotPoincare->isVisible())
        fig = PlotPoincare;
    if(PlotFreq->isVisible())
        fig = PlotFreq;
    if(PlotLyapunov->isVisible())
        fig = PlotLyapunov;
    if(PlotFreqMap->isVisible())
        fig = PlotFreqMap;
    if(PlotFreqMapHisto->isVisible())
        fig = PlotFreqMapHisto;
    if(PlotFreqMapSS->isVisible())
        fig = PlotFreqMapSS;
    if(PlotSchw->isVisible())
        fig = PlotSchw;
#ifdef USEQWT3D
    bool Plot3dVisible=Plot3d->isVisible();
    // 3d plot is not inherited from QwtPlot, so treat it differently
#else
    bool Plot3dVisible=false;
#endif
    if(fig==NULL && !Plot3dVisible)
        return;
    QString fileName = QFileDialog::getSaveFileName(this, "Print figure to file", core->WorkDir, "PostScript (*.ps);;PDF (*.pdf)");
    if (fileName.isEmpty())
        return;
    core->WorkDir = getDirName(fileName);
    if(fig)
    {
        QPrinter P;
        P.setOutputFileName(fileName);
        fig->print(P);
    }
#ifdef USEQWT3D
    else if(Plot3dVisible)
    {
        QString Format=fileName.endsWith(".pdf", Qt::CaseInsensitive) ? "PDF" : "PS";
        Qwt3D::VectorWriter* printer = (Qwt3D::VectorWriter*)Qwt3D::IO::outputHandler(Format);
        printer->setTextMode(Qwt3D::VectorWriter::NATIVE);
        printer->setSortMode(Qwt3D::VectorWriter::BSPSORT);
        printer->setLandscape(Qwt3D::VectorWriter::OFF);
        Qwt3D::IO::save(Plot3d, fileName, Format);
    }
#endif
}

///------------- Triangulation thread helper class --------------///
#ifdef USE3DFACETS
CTriangThread::CTriangThread(const QString& _appPath) : appPath(_appPath)
{
    finish = restart = false;
}

CTriangThread::~CTriangThread()
{
    mutex.lock();
    finish=restart=true;
    condition.wakeOne();
    mutex.unlock();
}
 
void CTriangThread::run()
{
    while(!finish)
    {
        createTriangulation();
        mutex.lock();
        if(!restart)
            condition.wait(&mutex);
        restart=false;
        mutex.unlock();
    }
}

void CTriangThread::startTriangulation(const std::vector<CPosPoint<double> > &_trajp, double _maxdist, int _slice)
{
    QMutexLocker locker(&mutex);
    // store data 
    nPointsMy = static_cast<int>(_trajp.size());
    trajpMy=_trajp;
    sliceMy = _slice;
    maxdistMy = _maxdist;
    // start or restart thread
    if(!isRunning())
        start(QThread::LowPriority);
    else
    {    
        restart=true;
        condition.wakeOne();
    }
}

bool CTriangThread::toolong(int v1, int v2)
{
    return (pow_2(trjpos[v1].x-trjpos[v2].x) + pow_2(trjpos[v1].y-trjpos[v2].y) + pow_2(trjpos[v1].z-trjpos[v2].z) > pow_2(maxdist));
}
void CTriangThread::addface(int v1, int v2, int v3)
{
    CTriangle T(v1, v2, v3);
    if(Triangles.remove(T)) 
        T.dup=true;
    Triangles << T;
}

void CTriangThread::createTriangulation()
{
    {
        Triangles.clear();
        // before proceeding further, copy data to local variables
        mutex.lock();
        facets.clear();
        trjpos.clear();
        maxdist=maxdistMy;
        int coord = sliceMy-1;
        int nPoints=0;
        for(int i=0; i<nPointsMy; i++)
        {
            if((sliceMy==0) || (trajpMy[i].Pos[coord]>=0))
            {
                trjpos.push_back(Qwt3D::Triple(trajpMy[i].Pos[0], trajpMy[i].Pos[1], trajpMy[i].Pos[2]));
                nPoints++;
            }
        }
        mutex.unlock();

        QProcess process;
        process.setReadChannel(QProcess::StandardOutput);
        process.start(appPath, QStringList("i"));
        if(!process.waitForStarted()) return;
        QTextStream strmI(&process);
        strmI << "3" << endl << nPoints << endl;
        for(int i=0; i<nPoints; i++)
        {
            strmI << trjpos[i].x << " " << trjpos[i].y << " " << trjpos[i].z << endl;
            if(restart) return;
        }
        process.closeWriteChannel();
        int nttot=0, nt=0;
        if(process.waitForFinished())
        {
            QTextStream strmO(&process);
            QString ln = strmO.readLine();
            nttot = ln.toInt();
            while(!strmO.atEnd())
            {
                ln = strmO.readLine();
                QStringList abcd=ln.split(" ");
                if(abcd.size() >= 4)
                {
                    int v1=abcd[0].toInt();
                    int v2=abcd[1].toInt();
                    int v3=abcd[2].toInt();
                    int v4=abcd[3].toInt();
                    if((v1>=nPoints) || (v2>=nPoints) || (v3>=nPoints) || (v4>=nPoints)) continue;
                    if(!(toolong(v1,v2) || toolong(v1,v3) || toolong(v1,v4) ||
                         toolong(v2,v3) || toolong(v2,v4) || toolong(v3,v4)))
                    {
                        addface(v1, v2, v3);
                        addface(v1, v2, v4);
                        addface(v1, v3, v4);
                        addface(v2, v3, v4);
                        nt++;
                    }
                }
                if(restart) return;
            }
        }
        // update faces
        int ntot=0, nout=0;
        QList<CTriangle> TList = Triangles.toList();
        Qwt3D::Cell tricell(3);
        foreach (CTriangle T, TList)
        {
            ntot++;
            if(!T.dup)
            {
                tricell[0] = T.c1;
                tricell[1] = T.c2;
                tricell[2] = T.c3;
              //  if((MyOrbit->tp[0][T.c1]>0) && (MyOrbit->tp[0][T.c2]>0) && (MyOrbit->tp[0][T.c3]>0))
                    facets.push_back(tricell);
                nout++;
            }
            if(restart) return;
        }
        emit finished(true);
    }
}

#endif

QwtText QwtPlotPicker2::trackerText( const QPointF &pos ) const
{
    QString text;

    switch ( rubberBand() )
    {
        case HLineRubberBand:
            text.sprintf( "%.4g", pos.y() );
            break;
        case VLineRubberBand:
            text.sprintf( "%.4g", pos.x() );
            break;
        default:
            text.sprintf( "%.4g, %.4g", pos.x(), pos.y() );
    }
    return QwtText( text );
}

QwtText QwtPlotZoomer2::trackerText( const QPointF &pos ) const
{
    QString text;

    switch ( rubberBand() )
    {
        case HLineRubberBand:
            text.sprintf( "%.4g", pos.y() );
            break;
        case VLineRubberBand:
            text.sprintf( "%.4g", pos.x() );
            break;
        default:
            text.sprintf( "%.4g, %.4g", pos.x(), pos.y() );
    }
    return QwtText( text );
}

}  // namespace