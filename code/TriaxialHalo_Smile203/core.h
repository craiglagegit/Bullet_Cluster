// this file defines COrbitCore - the driver class for both console and GUI versions of SMILE, handling ini data, starting and managing computation threads, interpreting scripts, etc. And a number of helper thread classes.
#pragma once
#include "common.h"
#include "potential.h"
#include <QObject>
#include <QThread>
#include <QTime>
#include <QFile>

namespace smile{

///////////////////////////////////// Helper threads that run various tasks in parallel to main program /////

// Thread that runs orbit integration 
class CCalcThread : public QThread   ///!!! rename these thread objects to more meaningful
{
    Q_OBJECT
public:
    virtual void run();
    CCalcThread(COrbit* _orbit, double _time): orbit(_orbit), time(_time) {};
signals:
    void calcFinished();
private:
    COrbit *orbit;   // do not copy orbit - store only pointer
    double time;
};

// Thread that runs orbit library integration (several threads in parallel)
class CCalcManyThread : public QThread
{
    Q_OBJECT
public:
    virtual void run();
    CCalcManyThread(COrbitLibrary* _orbits): orbits(_orbits) {};
signals:
    void calcFinished();
public slots:
    void stopThread();
private:
    COrbitLibrary* orbits;
    bool finished;
};

/** Thread that runs Schw.model optimization solver.
    Initialized with reference to the Schw.model object, an instance of a solver, and an instance of orbit penalty function.
    The latter two are destroyed by the thread upon completion.
    When the optimization is done, emit a signal with text representation of the result from model or the error from the solver
**/
class CSchwarzschildThread : public QThread
{
    Q_OBJECT
public:
    CSchwarzschildThread(const CBasicSchwModel* _model, COrbitLibrary* _orbitlib, const CBasicOptimizationSolver* _solver, const CBasicOrbitFilteringFnc* _orbitPenaltyFnc): 
        model(_model), orbitlib(_orbitlib), solver(_solver), orbitPenaltyFnc(_orbitPenaltyFnc) {};
    ~CSchwarzschildThread() { 
        if(solver!=NULL) delete solver;
        if(orbitPenaltyFnc!=NULL) delete orbitPenaltyFnc; 
    };
    virtual void run();
signals:
    void calcFinished(const QString& result);
private:
    const CBasicSchwModel* model;
    COrbitLibrary* orbitlib;
    const CBasicOptimizationSolver* solver;
    const CBasicOrbitFilteringFnc* orbitPenaltyFnc;
};

/** Thread that calculates statistics after Nbody export (currently, computes virial ratio).
    Initialized with a pointer to point mass set (Nbody model), and an info string.
    Point mass set is destroyed upon thread completion.
    Returns info string + virial ratio computed.
**/
class CPointMassSetHandler;  // defined in orbitlib.h
class CNbodyExportStatThread : public QThread
{
    Q_OBJECT
public:
    CNbodyExportStatThread(const CPointMassSetFloat* _points, const QString& _info) : 
        points(_points), handler(NULL), info(_info) {};
    virtual void run();
signals:
    void exportFinished(const QString& result);
public slots:
    void stopThread();
private:
    const CPointMassSetFloat* points;
    CPointMassSetHandler* handler;
    QString info;
};

struct CConfigCore
{
    double intTimeStepsPerPeriod;    // number of timesteps per unit time described above
    double intTimeInPeriods;         // orbit integration time measured in units of long-axis-orbit periods
    double intTimeMaxAdaptive;       // upper limits for adaptive integration time, expressed in energy-dependent orbital time periods
    bool usePS;                      // whether to plot points on Poincare section
    // data for frequency map
    int fm_numOrbitsStationary;      // number of orbits in stationary start-space
    int fm_numOrbitsPrincipalPlane;  // number of orbits in principal-plane start-space
    int fm_numOrbitsYalpha;          // number of orbits in Y-alpha start space
    int fm_numOrbitsRandom;          // number of orbits in random(ergodic) start-space
    // data for Schwarzschild modelling
    int sm_numOrbitsRandom;
    int sm_numShells;                // grid size: number of radial shells
    int sm_linesPerSegment;          // grid size: number of lines dividing each of 3 segments of each shell into cells
    int sm_numAngularCoefs;          // order of angular spherical-harmonic expansion
    int sm_numRadialCoefs;           // number of radial coefs in BSE 
    double sm_Alpha;                 // Alpha-parameter defining family of basis functions in BSE model
    double sm_innerShellMass;        // fraction of mass in the innermost shell, 0 for default value of 1/(numShells+1)
    double sm_outerShellMass;        // fraction of mass enclosed by the outermost shell, 0 for default value of 1-1/(numShells+1)
    double chaoticMinFreqDiff;       // minimal value of frequency diffusion rate to consider orbit as chaotic
    double chaoticMinLambda;         // minimal value of Lyapunov indicator to consider orbit as chaotic (in OR relation with previous criterion)
    double chaoticWeightFactor;      // weight factor for chaotic (super)orbit in objective function: +1 - disfavor them, -1 - favor, 0 - treat equally as regular
    int sm_numSamplingPoints;        // number of points sampled from trajectory, default NUM_SAMPLING_POINTS
    bool sm_useBPMPD;                // whether to use BPMPD solver (if the executable is present), otherwise use GLPK
};
extern CConfigCore configCore;       // internal smile core/gui parameters (not related to any libsmile modules)

//////////////////////////////////// CSmileCore class holds all non-gui workflow ////

class CSmileCore: public QObject
{
    Q_OBJECT

public:
    // --- configuration data --- //
    COrbitInitData<double> InitData;       // initial conditions for a single orbit
    double initCondE;                      // energy corresponding to initial conditions
    bool useICe;                           // whether to use energy as an initial condition instead of position/velocity (applies to freq.map)
    QString WorkDir, TempDir, AppDir;      // directories: current working directory (for GUI file dialogs); tempdir - for exchanging data with external programs (qdelaunay, bpmpd); appdir - where to locate these programs
    QString NbodyFile;                     // name of file containing discrete points for tree-code or BSE/spline expansion potentials
    // --- end of configuration data --- //

    const CPotential* potential;           // potential currently used in all calculations
    COrbit* orbit;                         // current orbit
    COrbitLibrary* orbitlib;               // library of orbits (for frequency map and Schwarzschild model)
    const CBasicSchwModel* model;          // Schwarzschild model object

    CCalcThread* CalcThr;            // thread that performs orbit integration asynchronously
    std::vector<CCalcManyThread*> CalcManyThr;  // pool of threads for calculating orbits
    CSchwarzschildThread* SchwThr;   // thread that does Schwarzschild modelling
    CNbodyExportStatThread* NBexportThr; // thread that performs nbody export statistics
private:
    int maxNumThread;
    int myTimer;                     // timer object for displaying progress in orbit library integration
    QTime timeElapsed, timeRuntime;  // measure time for orbit integration and for total runtime
    QFile input;                     // input file for scripting
    bool ConsoleMode, ConsoleInput;  // set when running console version; if interactive, set ConsoleInput=true, if processing script, set to false
    vectorRuntimeFncCreators creatorsLib;  // array of runtime function generators for orbit library
    //;   // in exporting Nbody model, this array holds the orbits for which the number of sampling points was insufficient and they need reintegration
    int numPointsNbodyExport;              // internally stored parameters for Nbody export
    QString fileNameNbodyExport; 
    int numBinsRefineNbodyExport;
    bool firstCallNbodyExport;

public:
    CSmileCore();                              // create the instance and set global variable 'core'
    ~CSmileCore();
    bool loadSettings(QString fileName="");    // load configuration data (at start - from predefined file smile.ini, later - when loading orbit library
    void saveSettings(QString fileName="");    // store configuration data (when storing orbit library, or at exit, if demanded in GUI)
    void runScript(QString scriptFileName=""); // in console mode: run script file or process data from console
    // potential
    void initPotential();                      // create an instance of potential specified by core->PotentialType
    void initSchwModel();                      // create an instance of Schwarzschild model of a type specified by config->modeltype
    // orbit integration
    void initIC();                             // set up initial conditions (calc E from X,Y,Z or vice-versa, calc timeUnit, etc.)
    void startOrbit();                         // start thread that runs orbit integration
    // Frequency map and Schwarzschild modelling
    void BuildFreqMapIC(double energy, bool useExistingIC);  // create initial conditions for orbit library for frequency map at given energy
    void BuildSchwIC(bool useExistingIC, CPointCountSetFloat* orbitsForReintegrate=NULL);                    // create initial conditions for orbit library for Schwarzschild model
    void StartOrbitLibrary(const char* finishSignal);        // start several threads to integrate orbits from orbit library
    void SchwStartOptimization(bool quadratic);      // start thread that solves Schw.optimization problem
    void SchwExportNbody(int numPoints, QString fileName, int numBinsRefine);  // export Schw model to Nbody model
    virtual void timerEvent(QTimerEvent*); // to display progress indicator in orbit library integration
    // import/export functions
    bool importOrbit(const QString &fileName);    // reads trajectory from a text file and assigns it to the current orbit
    bool exportOrbit(const QString &fileName);    // export current orbit into a text file
    bool exportPotential(const QString &fileName);      // export potential file (values of potential and forces in certain points of space, for statistics)
    bool importOrbitLib(const QString &fileName, bool withModelData);
    bool exportOrbitLib(const QString &fileName, bool withModelData);
    bool exportSchwModel(const QString &fileName);

private:
    void initPotentialNB();                    // load Nbody coords file, use direct N-body potential or calculate density and potential in basis set / spline expansion
    void initPotentialDefault();               // called when an error occurs during initPotential, since one cannot leave this object uninitialized
    bool exportPotentialAtPoints(const QString &fileName, const vectord points[N_DIM], int blockSize=0);   // actually do export potential, forces and density at given set of points

signals:
    void signalInfo(const QString &info);  // display various information (GUI - in message area, console - print out string)
    void signalTimer();                    // emitted regularly during orbit library integration, to update progress info
    void signalOrbitFinished();            // emitted when orbit integration finished
    void signalFreqMapFinished();          // emitted upon finish of frequency map building
    void signalSchwOrbitLibraryFinished(); // emitted upon finish of building orbit library for Schwarzschild modelling
    void signalSchwOptimizationFinished(const QString& result);  // emitted when Schwarzschild optimization finished
    void signalNbodyExportFinished(const QString& result);   // emitted when export to nbody model is finished
    void scriptNextLine();                 // emitted upon having processed a line in console script
    void scriptDone();                     // emitted when finished processing console script

private slots:
    void coreOrbitFinished();              // called when CalcThread finishes orbit integration
    void coreFreqMapFinished();            // called upon finish of last thread in frequency map calculation; clears unused orbits and emits signalFreqMapFinished
    void coreSchwOrbitLibraryFinished();   // called upon finish of last thread in building orbit library for Schwarzschild modelling; clears unused orbits and emits signalSchwOrbitLibraryFinished
    void coreSchwOptimizationFinished(const QString& result);  // called upon finish of Schwarzschild optimization thread; destroys thread and emits signalSchwOptimizationFinished
    void coreNbodyExport();                // perform export to nbody file, called either from SchwExportNbody or after finishing orbit reintegration
    void coreNbodyExportFinished(const QString& info);        // called nbody export is done
    void scriptProcessLine();              // processes next line in a script 
    void scriptOrbitLibraryFinished();     // displays some statistics (in console mode is connected to signalFreqMapFinished and signalSchwOrbitLibraryFinished)
    void scriptOptimizationFinished(const QString& message);  // displays results of optimization in console mode
    void scriptTimerEvent();               // displays progress in orbit library integration in console mode
    void scriptInfo(const QString& message);  // displays message in console mode
};

QString getDirName(const QString &fileName);  // convenience function, return file path ending with "/"

}