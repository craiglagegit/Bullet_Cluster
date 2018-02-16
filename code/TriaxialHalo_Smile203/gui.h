// this file defines GUI class which handles interactive user input and displays graphical interface. And a helper thread for performing triangulation for showing a 3d orbit as a solid body
#pragma once
#include "common.h"
#include <QMainWindow>
#include <QTime>
#include <QVector>
#include <QThread>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_zoomer.h>
#include "ui_smile.h"

#ifdef USEQWT3D
#include <qwt3d_surfaceplot.h>
#include <qwt3d_io_gl2ps.h>

#ifdef USE3DFACETS
#include <QMutex>
#include <QWaitCondition>
#include <QSet>
#endif

#endif

namespace smile{

#ifdef USE3DFACETS
class CTriangThread;  // forward declaration
#endif

class CSmileCore;

class CSmileGUI : public QMainWindow
{
    Q_OBJECT

public:
    CSmileGUI(QWidget *parent = 0, Qt::WFlags flags = 0);
    ~CSmileGUI();

private:
    Ui::CSmileGUIClass ui;
    QwtPlot *PlotOrbit;
    QwtPlotCurve *CurveOrbit;
    QwtPlotCurve *CurveEquipotential;
    QwtPlotZoomer* ZoomerOrbit;
    QwtPlot *PlotOrbitParam;
    QwtPlotCurve *CurveOrbitParam;
    QwtPlotZoomer* ZoomerOrbitParam;
    QwtPlot *PlotPoincare;
    QwtPlotZoomer *ZoomerPoincare;
    QVector<QwtPlotCurve*> CurvesPoincare;
    QwtPlot *PlotFreq;
    QwtPlotCurve *CurveFreqX;
    QwtPlotCurve *CurveFreqY;
    QwtPlotCurve *CurveFreqZ;
    QwtPlotCurve *CurveFreqPeaksP;
    QwtPlotCurve *CurveFreqPeaks;
    QwtPlotZoomer *ZoomerFreq;
    QwtPlot *PlotLyapunov;
    QwtPlotCurve *CurveDeviation;
    QwtPlotCurve *CurveLyapunov;
    QwtPlot *PlotFreqMap;
    QwtPlotCurve *CurveFreqMapR;
    QwtPlotCurve *CurveFreqMapC;
    QwtPlotZoomer *ZoomerFreqMap;
    QwtPlot *PlotFreqMapHisto;
    QwtPlotCurve *CurveFreqMapHisto;
    QwtPlot *PlotFreqMapSS;
    QwtPlotCurve *CurveFreqMapSSR;
    QwtPlotCurve *CurveFreqMapSSC;
    QwtPlot *PlotSchw;
    QwtPlotCurve *CurveSchwCellF;
    QwtPlotCurve *CurveSchwCellI;
    QwtPlot *PlotSchwHisto;
    QwtPlotCurve *CurveSchwHistoR, *CurveSchwHistoC;
#ifdef USEQWT3D
    Qwt3D::SurfacePlot *Plot3d;
    Qwt3D::TripleField trjpos;  // trajectory, points
    Qwt3D::CellField trjcells;  // trajectory cell (one for whole)
#ifdef USE3DFACETS
    CTriangThread *TriangThr;
#endif
#endif
    CSmileCore* core;                // the instance of core object holding all non-gui data
    double maxdist;                  // max. distance between adjacent points of trajectory (used in triangulation thread)
    bool InternalRecalc;             // indicate that modification of GUI controls is done not by user but internally, so do not react on this
    bool SchwModelRunning;           // state variable
    bool FreqMapRunning;             // yet another state variable
    bool ExportNbodyRunning;         // and yet another..
    double Eprev;                    // prev. energy (if difference is > Eeps, clear Poincare section)
    QString SchwModelResult;         // text message returned by Schwarzschild modelling routine
    int myTimerOrbit;                // timer for single orbit integration
    QTime timeOI;                    // measure time for orbit integration
    double timeOIelapsed;            // keep track of time used for orbit integration (not including analysis)
    QTime timeFM;                    // measure time for frequency map calculation
    virtual void timerEvent(QTimerEvent*);             // update orbit trajectory plot if orbit integration is taking too long
    bool eventFilter(QObject* target, QEvent* Event);  // track 'Enter' keydown events and starts single orbit integration in response
    bool loadApplySettings(QString fileName="");       // load ini file (at start or along with orbit library) and modifies content of GUI controls
    void drawFreqMapResonances();    // draw several most important resonance lines on frequency map plot
    static void XYZtoXY(const double X, const double Y, const double Z, double& x, double& y);   // to display start-space in projection
    QString getOrbitPopulation(int numShell, bool useWeights);         // calls corresponding function from orbit library, if it exists

signals:
    void KillFreqMapThreads();       // raised if user pressed stop during orbit library integration
private slots:
    void on_checkBox_freqLines_stateChanged(int);
    void on_ButtonPrint_clicked();
    void on_checkBox_usePS_stateChanged(int);
    void chooseNbodyFile();
    void typeinNbodyFile();
    void assignNbodyFile(QString fileName);
    void on_spinBox_IntervalNum_valueChanged(int);
    void on_spinBox_intervalCount_valueChanged(int);
    void on_radioButton_ICe_toggled(bool);
    void on_ButtonSchwNbody_clicked();
    void on_ButtonSchwLP_clicked();
    void on_ButtonSchwQP_clicked();
    void on_ButtonSchwLucy_clicked();
    void on_ButtonSchwModel_clicked();
    void on_ButtonFreqMap_clicked();
    void on_checkBox_fmConfigSpace_stateChanged(int);
    void on_selectEnergyLevel_fm_valueChanged(int);
    void on_selectEnergyLevel_valueChanged(int);
    void on_ButtonExportModel_clicked();
    void on_ButtonExportPS_clicked();
    void on_ButtonImportSM_clicked();
    void on_ButtonExportSM_clicked();
    void on_ButtonImportFM_clicked();
    void on_ButtonExportFM_clicked();
    void on_ButtonImportOrbit_clicked();
    void on_ButtonExportOrbit_clicked();
    void on_pushButton_refresh3dmesh_clicked();
    void on_ButtonRandom_clicked();
    void on_radioButton_2d_toggled(bool);
    void on_checkBox_Lyapunov_stateChanged(int);
    void on_ButtonClearPS_clicked();
    void on_ButtonStart_clicked();
    void CalcFinished();
    void FreqMapFinished();
    void SchwOrbitLibraryFinished();
    void SchwOptimizationFinished(const QString &message);
    void NbodyExportFinished(const QString& message);
    void timerEventCore();
    void intTime_changed();
    void fmparams_changed();
    void smparams_changed();
    void ccparams_changed();
    void potentialChanged();
    void on_buttonInitPotential_clicked();
    void ICchanged();
    void PSclicked(const QwtDoublePoint &pos);
    void FMclicked(const QwtDoublePoint &pos);
    void SMclicked(const QwtDoublePoint &pos);
    void SMHistoClicked(const QwtDoublePoint &pos);
    void redrawOrbit(bool);
    void redrawFreqMap(bool);
    void redrawSchwModel(bool);
    void refreshOrbitAndSpectrum();
public slots:   // yeah, black jack and bitches!
    void info(const QString &message);
    void error(const QString &message);
};

/// Helper thread that performs Delaunay triangulation with the help of external program
#ifdef USE3DFACETS
struct CTriangle
{
    int c1, c2, c3;
    bool dup;
    explicit CTriangle(int i1=0, int i2=0, int i3=0) 
    {
        if(i1<i2)
        {
            if(i2<i3)
            {  c1=i1; c2=i2; c3=i3; }
            else if(i1<i3)
            {  c1=i1; c2=i3; c3=i2; }
            else
            {  c1=i3; c2=i1; c3=i2; }
        }
        else
        {
            if(i1<i3)
            {  c1=i2; c2=i1; c3=i3; }
            else if(i2<i3)
            {  c1=i2; c2=i3; c3=i1; }
            else
            {  c1=i3; c2=i2; c3=i1; }
        }
        dup=false;
    }

    inline bool operator==(CTriangle T) const
    {
        return ((T.c1==c1) && (T.c2==c2) && (T.c3==c3));
    }
};
inline unsigned int qHash(const smile::CTriangle T)
{
    return ::qHash(T.c1*T.c2*T.c3);
}
// Thread that performs 3d orbit triangulation
class CTriangThread : public QThread
{
    Q_OBJECT
public:
    volatile bool finish;
    volatile bool restart;
    Qwt3D::TripleField trjpos;
    Qwt3D::CellField facets;      // trajectory boundary mesh facets (triangles)
    CTriangThread(const QString& _appPath);
    ~CTriangThread();
    void startTriangulation(const std::vector<CPosPoint<double> > &_trajp, double _maxdist, int _slice);
signals:
    void finished(bool);
protected:
    void run();
private:
    const QString appPath;
    QMutex mutex;
    QWaitCondition condition;
    // these variables are accessed both from thread and from outside
    int nPointsMy;
    std::vector<CPosPoint<double> > trajpMy;
    int sliceMy;
    double maxdistMy;
    // these are copies that are handled exclusively inside thread (in createTriangulation() and its subroutines)
    std::vector<CPosPoint<double> > trajp;
    double maxdist;
    QSet<CTriangle> Triangles;         // set of triangles for all tetrahedrons of figure

    bool toolong(int v1, int v2);
    void addface(int v1, int v2, int v3);
    void createTriangulation();
    void updateFaces();
};
#endif

/** helper class for indicating position on a 2d plot with arbitrary format string **/
class QwtPlotPicker2: public QwtPlotPicker
{
public:
    QwtPlotPicker2(int xAxis, int yAxis, int selectionFlags, RubberBand rubberBand, DisplayMode trackerMode, QwtPlotCanvas* canvas):
      QwtPlotPicker(xAxis, yAxis, selectionFlags, rubberBand, trackerMode, canvas) {};
protected:
    virtual QwtText trackerText( const QPointF &pos ) const;
};

class QwtPlotZoomer2: public QwtPlotZoomer
{
public:
    QwtPlotZoomer2(int xAxis, int yAxis, int selectionFlags, DisplayMode trackerMode, QwtPlotCanvas* canvas, bool doReplot=true) :
      QwtPlotZoomer(xAxis, yAxis, selectionFlags, trackerMode, canvas, doReplot) {};
protected:
    virtual QwtText trackerText( const QPointF &pos ) const;
};

}  // namespace