#include "core.h"
#include "orbit.h"
#include "orbitlib.h"
#include "schwarzschild.h"
#include "optimization.h"
#include "fileio.h"
#include <iostream>
#include <cassert>
#include <QCoreApplication>
#include <QSettings>
#include <QStringList>
#include <QTextStream>
#include <QDir>
#include <QMutex>
#include <gsl/gsl_errno.h>

namespace smile{

CConfigCore configCore=       // internal smile core/gui parameters (not related to any libsmile modules)
{
    50,   //double intTimeStepsPerPeriod;    // number of timesteps per unit time described above
    100,  //double intTimeInPeriods;         // orbit integration time measured in units of long-axis-orbit periods
    0,    //double intTimeMaxAdaptive;       // upper limits for adaptive integration time, expressed in energy-dependent orbital time periods
    false,//bool usePS;                      // whether to plot points on Poincare section
    // data for frequency map
    192,  //int fm_numOrbitsStationary;      // number of orbits in stationary start-space
    192,  //int fm_numOrbitsPrincipalPlane;  // number of orbits in principal-plane start-space
    0,    //int fm_numOrbitsYalpha;          // number of orbits in Y-alpha start space
    0,    //int fm_numOrbitsRandom;          // number of orbits in random(ergodic) start-space
    // data for Schwarzschild modelling
    10000,//int sm_numOrbitsRandom;
    25,   //int sm_numShells;                // grid size: number of radial shells
    4,    //int sm_linesPerSegment;          // grid size: number of lines dividing each of 3 segments of each shell into cells
    6,    //int sm_numAngularCoefs;          // order of angular spherical-harmonic expansion
    25,   //int sm_numRadialCoefs;           // number of radial coefs in BSE 
    0,    //double sm_Alpha;                 // Alpha-parameter defining family of basis functions in BSE model, 0=autodetect
    0.01, //double sm_innerShellMass;        // fraction of mass in the innermost shell, 0 for default value of 1/(numShells+1)
    0.99, //double sm_outerShellMass;        // fraction of mass enclosed by the outermost shell, 0 for default value of 1-1/(numShells+1)
    1e-3, //double chaoticMinFreqDiff
    0,    //double chaoticMinLambda
    0,    //double chaoticWeightFactor
    0,    //int sm_numSamplingPoints;        // number of points sampled from trajectory, default NUM_SAMPLING_POINTS
    true  //bool sm_useBPMPD;                // whether to use BPMPD solver (if the executable is present), otherwise use GLPK
};

QMutex mutexErrorInfo;  // mutex for ensuring that error line gets displayed at once
// change default GSL error handler to less nervous one
static void gsl_verbose_error_handler (const char *reason, const char *file, int line, int gsl_errno)
{
    if(gsl_errno == GSL_ETOL || gsl_errno == GSL_EROUND || gsl_errno == GSL_EOVRFLW || gsl_errno == GSL_EUNDRFLW) return;  // do nothing
    mutexErrorInfo.lock();
    std::cerr << "GSL error " << gsl_errno << " in " << file << " at line " << line << ": " << reason << "\n";
    mutexErrorInfo.unlock();
} 

CSmileCore::CSmileCore()
{
    WorkDir = QDir::currentPath(); 
    AppDir  = getDirName(QCoreApplication::applicationFilePath());
    TempDir = "";
    ConsoleMode=false;
    potential=NULL;
    orbit=NULL;
    orbitlib=NULL;
    model=NULL;
    CalcThr=NULL;
    SchwThr=NULL;
    NBexportThr=NULL;
    useICe=false;
    //initPotentialAndSymmetryNameMap();  // map names to enum identifiers
    gsl_set_error_handler(&gsl_verbose_error_handler);
    maxNumThread=std::max<int>(1, QThread::idealThreadCount());  // auto select
}

CSmileCore::~CSmileCore()
{   // commented out because the objects may be destroyed while running multi-threaded orbit library integration, which results in a crash
    /*if(orbit!=NULL)
    {
        delete orbit->getPotential();
        delete orbit;
    }
    if(model!=NULL)
        delete model;
    if(orbitlib!=NULL)
    {
        delete orbitlib->getPotential();
        delete orbitlib;
    }
    if(potential!=NULL)
        delete potential;*/
}

void CSmileCore::initPotential()
{
    if(potential) 
        delete potential;
    switch(configPotential.PotentialType)
    {
    case CDensity::PT_LOG:
        potential=new CPotentialLog(configPotential.N_dim, configPotential.q, configPotential.p, configPotential.Mbh, configPotential.Rc);
        break;
    case CDensity::PT_HARMONIC:
        potential=new CPotentialHarmonic(configPotential.N_dim, configPotential.q, configPotential.p, configPotential.Mbh);
        break;
    case CDensity::PT_DEHNEN:
        potential=new CPotentialDehnen(configPotential.N_dim, configPotential.q, configPotential.p, configPotential.Mbh, configPotential.Gamma);
        break;
    case CDensity::PT_SCALEFREE:
        potential=new CPotentialScaleFree(configPotential.N_dim, configPotential.q, configPotential.p, configPotential.Mbh, configPotential.Gamma);
        break;
    case CDensity::PT_SCALEFREESH:
        potential=new CPotentialScaleFreeSH(configPotential.N_dim, configPotential.q, configPotential.p, configPotential.Mbh, configPotential.Gamma, configPotential.Ncoefs_angular);
        break;
    case CDensity::PT_NB:
        initPotentialNB();
        break;
    case CDensity::PT_SPLINE:
    case CDensity::PT_BSE:   // two variants possible: either init from analytic density model or from Nbody/coef file
        if(configPotential.DensityType==CDensity::PT_NB)
            initPotentialNB();
        else {   // create temporary instance of analytic density profile, pass it to the constructor to compute expansion coefficients
            CDensity* densModel=NULL;
            switch(configPotential.DensityType) 
            {
            case CDensity::PT_DEHNEN: 
                densModel = new CPotentialDehnen(3, configPotential.q, configPotential.p, 0, configPotential.Gamma); 
                break;
            case CDensity::PT_PLUMMER:
                densModel = new CDensityPlummer(configPotential.q, configPotential.p); 
                break;
            case CDensity::PT_ISOCHRONE:
                densModel = new CDensityIsochrone(configPotential.q, configPotential.p);
                break;
            case CDensity::PT_PERFECTELLIPSOID:
                densModel = new CDensityPerfectEllipsoid(configPotential.q, configPotential.p); 
                break;
            case CDensity::PT_NFW:
                double concentration = std::max<double>(configPotential.Rc, 1.0);
                densModel = new CDensityNFW(configPotential.q, configPotential.p, concentration);
                break;
            }
            if(densModel==NULL)  
            {  // error? init in some default way
                if(configPotential.PotentialType==CDensity::PT_BSE) potential = new CPotentialBSE(configPotential.N_dim, configPotential.Mbh, 1, 1, 1); 
                else potential = new CPotentialSpline(configPotential.N_dim, configPotential.Mbh, 1, 1); 
            }
            else
            {
                if(configPotential.PotentialType==CPotential::PT_BSE)
                    potential = new CPotentialBSE(configPotential.N_dim, configPotential.Mbh, configPotential.Alpha, configPotential.Ncoefs_radial, configPotential.Ncoefs_angular, densModel);
                else
                    potential = new CPotentialSpline(configPotential.N_dim, configPotential.Mbh, configPotential.Ncoefs_radial, configPotential.Ncoefs_angular, densModel);
                delete densModel;
            }
        }
        break;
    default:
        initPotentialDefault();
    }
    configPotential.PotentialType = potential->PotentialType();
#ifdef DEBUGPRINT
    //for(double x=1./8; x<10000; x*=2) std::cerr << x << "  " << potential->Mass(x) << "\n";
    std::cerr << "Gamma="<<potential->getGamma() << ", total mass=" << potential->totalMass() << "\n";
#endif
}

void CSmileCore::initPotentialNB()
{
    QFileInfo NbodyFileName(QDir(WorkDir), NbodyFile);
    CFileInputText in(NbodyFileName.absoluteFilePath().toStdString());
    if(!in.ok())
    {
        initPotentialDefault();
        my_error("Error loading potential from file "+NbodyFileName.absoluteFilePath().toStdString());
        return;
    }
    potential=in.readPotential(&configPotential);
    if(!potential)
        initPotentialDefault();
}

void CSmileCore::initPotentialDefault()
{
    potential=new CPotentialLog(configPotential.N_dim, configPotential.q, configPotential.p, configPotential.Mbh, configPotential.Rc);
}

/// Export/import operations ///
bool CSmileCore::exportPotential(const QString& fileName)
{
    if(potential==NULL) return false;
    vectord points[N_DIM];
    // export potential and its derivatives in a set of points: along three principal axes, three diagonals in each plane, and a x=y=z line
    const int Npoints=100;
    const double rmin=0.001, rmax=1000;
    vectord radii;
    createNonuniformGrid(radii, Npoints, rmin, rmax);
    double axiscoefs[7][3] = { {1,0,0}, {0,1,0}, {0,0,1}, {0.7,0.7,0}, {0,0.7,0.7}, {0.7,0,0.7}, {0.58,0.58,0.58} };
    for(int d=0; d<7; d++)
    {
        if((potential->symmetry() & CPotential::ST_PLANESYM) == 0)   // no triaxial symmetry -- output both negative and positive radii
        {
            for(int i=Npoints-1; i>=0; i--)
            {
                double r=-radii[i];
                for(unsigned int c=0; c<N_DIM; c++)
                    points[c].push_back(r*axiscoefs[d][c]);
            }
        }
        for(int i=0; i<Npoints; i++)
        {
            double r=radii[i];
            for(unsigned int c=0; c<N_DIM; c++)
                points[c].push_back(r*axiscoefs[d][c]);
        }
    }
    if(!exportPotentialAtPoints(fileName, points, points[0].size()/7))
        return false;
    // export potential, forces and density at points defined as initial conditions for orbit library (if it exists)
    if(orbitlib && orbitlib->size()>0)
    {
        unsigned int count=orbitlib->size();
        for(unsigned int c=0; c<N_DIM; c++)
            points[c].assign(count, 0);
        for(unsigned int i=0; i<count; i++)
        {
            const COrbitDesc* od=orbitlib->getOrbitDesc(i);
            COrbitInitData<float> initData(od->getInitData());
            for(unsigned int c=0; c<N_DIM; c++)
                points[c][i]=initData.initCond.Pos[c];
        }
        if(!exportPotentialAtPoints(fileName+".points", points))
            return false;
    }
    /*if(potential->PotentialType()==CDensity::PT_NB)
    {
        {QFile filep(fileName+".body");
        if (!filep.open(QFile::WriteOnly | QFile::Text)) 
            return false;
        QTextStream out(&filep);
        out.setRealNumberPrecision(8);
        CPotentialNB* pot = (CPotentialNB*)potential;
        for(int i=0; i<pot->bodytab.size(); i++)
            out << pot->bodytab[i].pos[0] << "\t" << pot->bodytab[i].pos[1] << "\t" << pot->bodytab[i].pos[2] << "\t0\t0\t0\t" << pot->bodytab[i].mass << "\t" << pot->bodytab[i].eps << "\n";
        }
        {QFile filep(fileName+".cell");
        if (!filep.open(QFile::WriteOnly | QFile::Text)) 
            return false;
        QTextStream out(&filep);
        out.setRealNumberPrecision(8);
        CPotentialNB* pot = (CPotentialNB*)potential;
        for(int i=0; i<pot->celltab.size(); i++)
            out << pot->celltab[i].pos[0] << "\t" << pot->celltab[i].pos[1] << "\t" << pot->celltab[i].pos[2] << "\t" << pot->celltab[i].numbody << "\t" << pot->celltab[i].rmax << "\t" << pot->celltab[i].eps << "\t" << 
            (pot->celltab[i].quad[0][0]/pow_2(pot->celltab[i].rmax)/pot->celltab[i].mass) << "\t" << 
            (pot->celltab[i].quad[0][1]/pow_2(pot->celltab[i].rmax)/pot->celltab[i].mass) << "\t" << 
            (pot->celltab[i].quad[0][2]/pow_2(pot->celltab[i].rmax)/pot->celltab[i].mass) << "\t" << 
            (pot->celltab[i].quad[1][1]/pow_2(pot->celltab[i].rmax)/pot->celltab[i].mass) << "\t" << 
            (pot->celltab[i].quad[1][2]/pow_2(pot->celltab[i].rmax)/pot->celltab[i].mass) << "\t" << 
            (pot->celltab[i].quad[2][2]/pow_2(pot->celltab[i].rmax)/pot->celltab[i].mass) << "\t" << 
            "\n";
        }
    }*/
    // export coefficients if it makes sense
    if(potential->PotentialType()==CDensity::PT_BSE || potential->PotentialType()==CDensity::PT_SPLINE || potential->PotentialType()==CDensity::PT_SCALEFREESH)
    {
        CFileOutputText out((fileName+".coef").toStdString());
        return out.writePotential(potential);
    }
    else
        return true;
}

bool CSmileCore::exportPotentialAtPoints(const QString& fileName, const vectord points[N_DIM], int blockSize)
{
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) 
        return false;
    QTextStream out(&file);
    out.setRealNumberPrecision(8);
    out << "x\ty\tz\tPhi\tFx\tFy\tFz\tRho\n";
    for(size_t p=0; p<points[0].size(); p++)
    {
        double v[6]={points[0][p], points[1][p], points[2][p], 0, 0, 0};
        double f[8];
        double Phi = potential->Phi(v[0],v[1],v[2]);
        potential->DiffEq(N_DIM*2, 0, v, f);
        out << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << Phi << "\t" 
            << -f[3] << "\t" << -f[4] << "\t" << -f[5] << "\t" << potential->Rho(v[0],v[1],v[2]) << "\n";
        if(blockSize>0 && (p+1)%blockSize==0) out << "\n";  // separating line
    }
    return true;
}

bool CSmileCore::exportOrbit(const QString& fileName)
{
    if(orbit==NULL || orbit->getIntTime()==0) return false;
    CFileOutputText out(fileName.toStdString());
    out.writeOrbit(orbit);
    return out.ok();
}

bool CSmileCore::importOrbit(const QString &fileName)
{
    if(CalcThr!=NULL) return false;
    CFileInputText in(fileName.toStdString());
    if(!in.ok()) return false;
    CPotential* newPot=potential->clone();
    COrbit* newOrb=in.readOrbit(newPot);
    if(!newOrb)
    {
        delete newPot;
        return false;
    }
    if(orbit!=NULL)
    {
        delete orbit->getPotential();  // manually delete, since it was created by cloning externally to the object
        delete orbit;
    }
    orbit=newOrb;
    InitData = orbit->getInitData();
    initCondE = potential->totalEnergy(InitData.initCond);
    configCore.intTimeInPeriods=orbit->getIntTime()/InitData.timeUnit; 
    configCore.intTimeStepsPerPeriod=InitData.timeUnit/InitData.timeStep;
    return true;
}

bool CSmileCore::exportOrbitLib(const QString &fileName, bool withModelData)
{
    if(orbitlib==NULL || orbitlib->size()==0) return false;
    CFileOutputText out(fileName.toStdString());
    out.writeOrbitLibrary(orbitlib, withModelData);
    return out.ok();
}

bool CSmileCore::importOrbitLib(const QString &fileName, bool withModelData)
{
    assert(potential);
    CFileInputText in(fileName.toStdString());
    if(!in.ok()) return false;
    CPotential* newPot=potential->clone();
    COrbitLibrary* newOrbLib=in.readOrbitLibrary(newPot, withModelData);
    if(!newOrbLib)
    {
        delete newPot;
        return false;
    }
    if(orbitlib!=NULL)
    {
        delete orbitlib->getPotential();  // manually delete, since it was created by cloning externally to the object
        delete orbitlib;
    }
    orbitlib=newOrbLib;
    return true;
}

bool CSmileCore::exportSchwModel(const QString& fileName)
{
    if(model==NULL || (model->ModelType() & CBasicSchwModel::MT_SHELL) == 0)
        return false;
    std::ofstream strm(fileName.toStdString().c_str(), std::ios::out);
    if(!strm) return false;
    const CBasicShellSchwModel* model_sh =static_cast<const CBasicShellSchwModel*>(model);
    const vectord& coefs = model_sh->getCoefsDens();
    vectord coefsPenalty(coefs.size(), 0);
    if(orbitlib!=NULL && orbitlib->size()!=0)
        coefsPenalty = model_sh->calcPenalty(orbitlib);
    const CSchwModelClassic* model_cl = (model->ModelType()==CBasicSchwModel::MT_CLASSIC)? static_cast<const CSchwModelClassic*>(model) : NULL;
    const CSchwModelSHGrid* model_shg = (model->ModelType()==CBasicSchwModel::MT_SHGRID) ? static_cast<const CSchwModelSHGrid*>(model) : NULL;
    const CSchwModelSHBSE* model_shb  = (model->ModelType()==CBasicSchwModel::MT_SHBSE)  ? static_cast<const CSchwModelSHBSE*>(model) : NULL;
    strm << "#index\tRequiredCoef\tActualCoef\tDifference\t";
    if(model_cl) 
        strm << "radius\tx\ty\tz\n";
    else if(model_shg)
        strm << "radius\tl\tm\n";
    else if(model_shb)
        strm << "n\tl\tm\n";
    for(size_t i=0; i<coefs.size(); i++)
    {
        strm << i << "\t" << coefs[i] << "\t" << (coefs[i]+coefsPenalty[i]) << "\t" << coefsPenalty[i] << "\t";
        if(model_cl)
        {
            double X, Y, Z;
            model_cl->cellCenter(i, X, Y, Z);
            strm << sqrt(X*X+Y*Y+Z*Z) << "\t" << X << "\t" << Y << "\t" << Z << "\n";
        }
        else if(model_shg)
        {
            if(i==0)
                strm << "0\t0\t0\n";
            else {
                int indRad= (i-1)/model_shg->getNumCoefsAtRadius();
                int indAng= i-1 - indRad*model_shg->getNumCoefsAtRadius();
                int l=static_cast<int>(sqrt(2.*indAng-floor(sqrt(2.*indAng))));
                int m=indAng-l*(l+1)/2;
                strm << model_sh->getShellRadius(indRad) << "\t" << l << "\t" << m << "\n";
            }
        } 
        else if(model_shb)
        {
            int indRad= i/model_shb->getNumCoefsAtRadius();
            int indAng= i - indRad*model_shb->getNumCoefsAtRadius();
            int l=static_cast<int>(sqrt(2.*indAng-floor(sqrt(2.*indAng))));
            int m=indAng-l*(l+1)/2;
            strm << indRad << "\t" << l << "\t" << m << "\n";
        } 
        else strm << "\n";
    }
    strm << "\n#indexR\tradius\tM(r)\tsigma_r^2\tsigma_t^2\n";

    vectord 
        shellVr(model_sh->getNumShellsKinem(), 0), 
        shellVt(model_sh->getNumShellsKinem(), 0), 
        shellTime(model_sh->getNumShellsKinem(), 0);
    if(orbitlib!=NULL) for(unsigned int o=0; o<orbitlib->size(); o++)
    {
        // load grid information data from orbit desc
        const CSchwInformation* infoGrid=static_cast<const CSchwInformation*> (orbitlib->getOrbitDesc(o)->getInfoByType(CBasicInformation::IT_SCHW));
        if(infoGrid!=NULL && infoGrid->getShellTime().size()==model_sh->getNumShellsKinem())
        {
            double weight=orbitlib->getOrbitDesc(o)->getWeight();
            for(unsigned int s=0; s<model_sh->getNumShellsKinem(); s++)
            {
                shellVr[s] += weight * infoGrid->getShellVr(s) * infoGrid->getShellTime(s);
                shellVt[s] += weight * infoGrid->getShellVt(s) * infoGrid->getShellTime(s);
                shellTime[s] += weight * infoGrid->getShellTime(s);
            }
        }
    }
    for(unsigned int s=0; s<model_sh->getNumShellsKinem(); s++)
    {
        strm << s << "\t" << model_sh->getShellRadius(s) << "\t" << model_sh->getInteriorMass(s) << "\t" << 
            (shellVr[s]/shellTime[s]) << "\t" << (shellVt[s]/shellTime[s]) << "\n";
    }
    return static_cast<bool>(strm);
}

//// ----------- orbit integration ------------ ////
void CSmileCore::initIC()
{
    if(useICe)   // assign from initCondE
    {
        InitData.initCond=CPosVelPoint<double>();
        InitData.initCond.Pos[0]=potential->longaxisradius(initCondE);
    }
    else
    {
        initCondE = potential->totalEnergy(InitData.initCond);
    }
    // calculate T_orb (for x-axis orbit)
    InitData.timeUnit = potential->longaxisperiod(initCondE);
}

void CSmileCore::startOrbit()
{
    if(InitData.timeUnit<=0 || configCore.intTimeInPeriods<=0 || configCore.intTimeStepsPerPeriod<=0) return;  // error
    if(CalcThr!=NULL)
    {
        CalcThr->deleteLater();
    }
    if(orbit!=NULL)
    {
        delete orbit->getPotential();  // manually delete, since it was created by cloning externally to the object
        delete orbit;
    }
    vectorRuntimeFncCreators crts;   // this array of generators is different from core->creatorsLib -- the latter is for all orbits in orbit library, the former is for a single orbit independent of these
    crts.push_back(new COrbitRuntimeTrajectoryCreator());
    if(configCore.usePS)
        crts.push_back(new COrbitRuntimePoincareCreator(0, 1));
    crts.push_back(new COrbitRuntimePericenterCreator());
    orbit = new COrbit(COrbitInitData<double>(potential->clone(), InitData.timeUnit/configCore.intTimeStepsPerPeriod, InitData.timeUnit, InitData.initCond, InitData.calcLyapunov), &crts);
    CalcThr = new CCalcThread(orbit, configCore.intTimeInPeriods*InitData.timeUnit);
    connect(CalcThr, SIGNAL(calcFinished()), SLOT(coreOrbitFinished()));
    CalcThr->start();
}

void CSmileCore::coreOrbitFinished()
{
    emit signalOrbitFinished();   // in GUI - displays various information, called immediately in the same thread before orbit is destroyed

    CalcThr->deleteLater();
    CalcThr=NULL;
}

//// ----------- Frequency map / Schwarzschild modelling ------------- ////

void CSmileCore::StartOrbitLibrary(const char* finishSignal)
{
    for(int t=0; t<maxNumThread; t++)
    {
        CalcManyThr.push_back(new CCalcManyThread(orbitlib));
        connect(CalcManyThr[t], SIGNAL(calcFinished()), finishSignal);
        CalcManyThr[t]->start(QThread::LowPriority);
    }
    emit signalInfo("Starting integration..");
    myTimer = startTimer(2000);  // to refresh progress
}

void CSmileCore::BuildFreqMapIC(double energy, bool useExistingIC)
{
    assert(potential);
    for(vectorRuntimeFncCreators::const_iterator iter=creatorsLib.begin(); iter!=creatorsLib.end(); ++iter)
        delete (*iter);
    creatorsLib.clear();
    creatorsLib.push_back(new COrbitRuntimeTrajectoryCreator());
    creatorsLib.push_back(new COrbitRuntimePericenterCreator());
    if(useExistingIC && orbitlib!=NULL && orbitlib->size()>0)
    {
        delete orbitlib->getPotential();  // manually delete, since it was created by cloning externally to the object
        orbitlib->modifyInitialConditions(potential->clone(), &creatorsLib, configCore.intTimeStepsPerPeriod, configCore.intTimeInPeriods, configCore.intTimeMaxAdaptive, InitData.calcLyapunov, true);
    }
    else
    {
        if(orbitlib)
        {
            delete orbitlib->getPotential();  // manually delete, since it was created by cloning externally to the object
            delete orbitlib;
        }
        orbitlib=new COrbitLibrary(potential->clone(), &creatorsLib);
        emit signalInfo("Preparing initial conditions..");
        orbitlib->prepareInitialConditionsE(configCore.fm_numOrbitsStationary, configCore.fm_numOrbitsPrincipalPlane, configCore.fm_numOrbitsYalpha, configCore.fm_numOrbitsRandom, energy, configCore.intTimeStepsPerPeriod, configCore.intTimeInPeriods, configCore.intTimeMaxAdaptive, InitData.calcLyapunov);
    }
}

void CSmileCore::coreFreqMapFinished()
{
    for(size_t t=0; t<CalcManyThr.size(); t++)
        CalcManyThr[t]->deleteLater();
    CalcManyThr.clear();
    killTimer(myTimer);
    myTimer=0;
    // delete non-used orbits (if integration threads were terminated prematurely, there exist some)
    orbitlib->removeUnused();
    // finally, send signal to the parent
    emit signalFreqMapFinished();
}

void CSmileCore::initSchwModel()
{
    if(model)
        delete model;
    assert(potential->N_dim==N_DIM);
    emit signalInfo("Preparing model..");
    switch(configSchw.ModelType) {
        case CBasicSchwModel::MT_CLASSIC:
            model=new CSchwModelClassic(potential, configCore.sm_numShells, configCore.sm_innerShellMass, configCore.sm_outerShellMass, configCore.sm_linesPerSegment);
            break;
        case CBasicSchwModel::MT_SHGRID:
            model=new CSchwModelSHGrid(potential, configCore.sm_numShells, configCore.sm_innerShellMass, configCore.sm_outerShellMass, configCore.sm_numAngularCoefs);
            break;
        case CBasicSchwModel::MT_SHBSE:
            model=new CSchwModelSHBSE(potential, configCore.sm_numShells, configCore.sm_innerShellMass, configCore.sm_outerShellMass, configCore.sm_numAngularCoefs, configCore.sm_numRadialCoefs, configCore.sm_Alpha);
            break;
        default: my_error("Unrecognized type of Schwarzschild model"); model=NULL;
    }
}
void CSmileCore::BuildSchwIC(bool useExistingIC, CPointCountSetFloat* orbitsForReintegrate)
{
    assert(potential);
    if(potential->N_dim!=3) return;
    // Prepare initial conditions
    if(useExistingIC && orbitlib!=NULL && orbitlib->size()>0)
    {
        delete orbitlib->getPotential();  // manually delete, since it was created by cloning externally to the object
        orbitlib->modifyInitialConditions(potential->clone(), &creatorsLib, configCore.intTimeStepsPerPeriod, configCore.intTimeInPeriods, configCore.intTimeMaxAdaptive, InitData.calcLyapunov, orbitsForReintegrate==NULL);
    }
    else
    {
        if(orbitlib)
        {
            delete orbitlib->getPotential();  // manually delete, since it was created by cloning externally to the object
            delete orbitlib;
        }
        orbitlib=new COrbitLibrary(potential->clone(), &creatorsLib);
        emit signalInfo("Preparing initial conditions..");
        double maxRadiusMass=configCore.sm_outerShellMass>0 && configCore.sm_outerShellMass<1 ? 1-0.25*(1-configCore.sm_outerShellMass) : 1-0.25/configCore.sm_numShells;
        orbitlib->prepareInitialConditions(configCore.sm_numOrbitsRandom, potential->getRadiusByMass(maxRadiusMass*potential->totalMass()), configCore.intTimeStepsPerPeriod, configCore.intTimeInPeriods, configCore.intTimeMaxAdaptive, InitData.calcLyapunov);
    }
    if(orbitsForReintegrate==NULL || model==NULL)
    {   // need to create model 
        initSchwModel();
    }
    // fill the list of runtime function creators
    for(vectorRuntimeFncCreators::const_iterator iter=creatorsLib.begin(); iter!=creatorsLib.end(); ++iter)
        delete (*iter);
    creatorsLib.clear();
    creatorsLib.push_back(new COrbitRuntimeTrajectoryCreator());
    creatorsLib.push_back(new COrbitRuntimePericenterCreator());
    if(configCore.sm_numSamplingPoints>0 || orbitsForReintegrate!=NULL) 
    {
        creatorsLib.push_back(new COrbitRuntimeTrajSampleCreator(configCore.sm_numSamplingPoints, orbitsForReintegrate));
    }
    switch(model->ModelType())
    {   // create grid runtime function associated with the model, for each orbit
    case CBasicSchwModel::MT_CLASSIC:
        creatorsLib.push_back(new COrbitRuntimeSchwClassicCreator(static_cast<const CSchwModelClassic*>(model)));
        break;
    case CBasicSchwModel::MT_SHGRID:
    case CBasicSchwModel::MT_SHBSE:
        creatorsLib.push_back(new COrbitRuntimeSchwSHCreator(static_cast<const CSchwModelSH*>(model)));
        break;
    default:
        my_error("Unknown model type");
    }
}

void CSmileCore::coreSchwOrbitLibraryFinished()
{
    for(size_t t=0; t<CalcManyThr.size(); t++)
        CalcManyThr[t]->deleteLater();
    CalcManyThr.clear();
    killTimer(myTimer);
    myTimer=0;
    // delete non-used orbits (if integration threads were terminated prematurely, there exist some)
    orbitlib->removeUnused();
    // finally, send signal to the parent
    emit signalSchwOrbitLibraryFinished();
}

void CSmileCore::SchwStartOptimization(bool quadratic)
{
    if(SchwThr!=NULL || model==NULL)
        return;      // do not allow to run two threads at once
    const CBasicOptimizationSolver* solver = 
        configCore.sm_useBPMPD ?
        static_cast<const CBasicOptimizationSolver*>(new COptimizationSolverBPMPD(AppDir.toStdString()+"bpmpd", TempDir.toStdString(), quadratic)) :
        static_cast<const CBasicOptimizationSolver*>(new COptimizationSolverGLPK());
    const CBasicOrbitFilteringFnc* orbitPenaltyFnc = 
        configCore.chaoticWeightFactor!=0 ? 
        new CChaosOrbitFilteringFnc(configCore.chaoticWeightFactor, configCore.chaoticMinLambda, configCore.chaoticMinFreqDiff) : NULL;
    SchwThr = new CSchwarzschildThread(model, orbitlib, solver, orbitPenaltyFnc);
    // solver and orbitPenaltyFnc will be automatically deleted upon finishing optimization
    connect(SchwThr, SIGNAL(calcFinished(const QString&)), SLOT(coreSchwOptimizationFinished(const QString&)));
    SchwThr->start();
}

void CSmileCore::coreSchwOptimizationFinished(const QString& result)
{
    emit signalSchwOptimizationFinished(result);
    SchwThr->deleteLater();
    SchwThr=NULL;
}

void CSmileCore::SchwExportNbody(int numPoints, QString fileName, int numBinsRefine)
{
    if(NBexportThr!=NULL || orbitlib==NULL || orbitlib->size()==0)
        return;  // don't allow to run twice
    numPointsNbodyExport=numPoints;
    fileNameNbodyExport=fileName;
    numBinsRefineNbodyExport=numBinsRefine;
    firstCallNbodyExport=true;
    coreNbodyExport();
}

void CSmileCore::coreNbodyExport()
{
    if(fileNameNbodyExport.isEmpty()) return;  // incorrect call
    // if it was called after reintegration, cleanup threads
    if(CalcManyThr.size()>0)
    {
        for(size_t t=0; t<CalcManyThr.size(); t++)
            CalcManyThr[t]->deleteLater();
        CalcManyThr.clear();  // clean up threads used to reintegrate
    }
    if(myTimer!=0){  // to stop showing up progress bar after reintegration is complete
        killTimer(myTimer);
        myTimer=0;
    }
    const CEnergyOrbitFilteringFnc rankFnc;
    const CMassRefinementFnc massRefineFnc(orbitlib, numBinsRefineNbodyExport, &rankFnc);
    CPointMassSetFloat* points = NULL;
    CPointCountSetFloat* orbitsForReintegrate = NULL;
    QString info = QString::fromStdString(
        orbitlib->exportNbody(numPointsNbodyExport, numBinsRefineNbodyExport ? &massRefineFnc : NULL, 
        points, orbitsForReintegrate) );
    if(points==NULL)
    {
        emit signalNbodyExportFinished("Error in creating Nbody model: "+info);
        fileNameNbodyExport.clear();
        return;
    }
    if(orbitsForReintegrate!=NULL && !orbitsForReintegrate->empty())
    {   // requires reintegration
        if(!firstCallNbodyExport)
        {
            emit signalNbodyExportFinished("Export of Nbody model cancelled or was unsuccessful");
            fileNameNbodyExport.clear();
            delete points;
            delete orbitsForReintegrate;
            return;
        }
        firstCallNbodyExport=false;
        BuildSchwIC(true, orbitsForReintegrate);
        StartOrbitLibrary(SLOT(coreNbodyExport()));   // will be called once again
        emit signalInfo("Reintegrating "+QString::number(orbitsForReintegrate->size())+" orbits..");
        delete points;
        delete orbitsForReintegrate;  // a working copy was made by TrajSampleCreator
        return;
    }
    emit signalInfo("Exporting Nbody model..");
    {
        CFileOutputNemo out((fileNameNbodyExport+".nb").toStdString());
        out.writePointMassSet(points);
    }
    {
        CFileOutputText out(fileNameNbodyExport.toStdString());
        out.writePointMassSet(points);
    }
    fileNameNbodyExport.clear();
    emit signalInfo(info+"\nComputing virial ratio...");
    if(NBexportThr!=NULL)
        return;      // do not allow to run two threads at once
    NBexportThr = new CNbodyExportStatThread(points, info);
    connect(NBexportThr, SIGNAL(exportFinished(const QString&)), SLOT(coreNbodyExportFinished(const QString&)));
    NBexportThr->start();   // point set will be destroyed by export thread upon completion
}

void CSmileCore::coreNbodyExportFinished(const QString& info)
{
    NBexportThr->deleteLater();
    NBexportThr=NULL;
    emit signalNbodyExportFinished(info);
    if(ConsoleMode) emit scriptNextLine();
}

//// ----------- load/save settings ------------- ////

bool CSmileCore::loadSettings(QString fileName)
{
    if(fileName=="") 
        fileName=AppDir+DEFAULT_SETTINGS_FILE_NAME;
    else if(!QFile::exists(fileName))
        return false;
    QSettings settings(fileName, QSettings::IniFormat);
    QString WorkDir1 = settings.value("Common/WorkDir", "").toString();
    if(WorkDir1!="")
    {
        QFileInfo fi(WorkDir1);
        if(fi.isRelative()) WorkDir = AppDir+WorkDir1; else WorkDir = WorkDir1;
    }
    QString TempDir1 = settings.value("Common/TempDir", "").toString();
    if(TempDir1!="")
    {
        QFileInfo fi(TempDir1);
        if(fi.isRelative()) TempDir = AppDir+TempDir1; else TempDir = TempDir1;
    } 
    else if(TempDir=="") TempDir=WorkDir;
    if(!WorkDir.endsWith('/')) WorkDir += '/';
    if(!TempDir.endsWith('/')) TempDir += '/';
    int newMaxNumThread=settings.value("Common/MaxThreads", "0").toInt();
    if(newMaxNumThread>=1) maxNumThread=newMaxNumThread;   // if no value exists in ini file, leave current
    // potential section
    configPotential.PotentialType = getPotentialTypeByName(settings.value("Potential/Type", "").toString().toStdString());
    configPotential.DensityType = getDensityTypeByName(settings.value("Potential/DensityModel", "").toString().toStdString());
    configPotential.N_dim = std::max<int>(2, std::min<int>(3, settings.value("Potential/N_dim", N_DIM).toInt()));
    configPotential.q = std::max<double>(0, std::min<double>(1, settings.value("Potential/q", configPotential.q).toDouble())); 
    configPotential.p = std::max<double>(0, std::min<double>(configPotential.q, settings.value("Potential/p", configPotential.p).toDouble())); 
    configPotential.Mbh = std::max<double>(0, settings.value("Potential/Mbh", configPotential.Mbh).toDouble());
    configPotential.Rc = std::max<double>(0, settings.value("Potential/Rc", configPotential.Rc).toDouble());
    configPotential.Gamma = std::max<double>(0, std::min<double>(2, settings.value("Potential/Gamma", configPotential.Gamma).toDouble()));
    configPotential.Alpha = std::max<double>(0, settings.value("Potential/Alpha", configPotential.Alpha).toDouble());  // allowed range is [0.5 - infinity], alpha=0 means try to autodetect..
    configPotential.Omega = settings.value("Potential/Omega", configPotential.Omega).toDouble();
    if(configPotential.Omega<0) configPotential.Omega=-configPotential.Omega;
    configPotential.Ncoefs_radial = std::min<size_t>(MAX_NCOEFS_RADIAL, settings.value("Potential/Ncoefs_radial", configPotential.Ncoefs_radial).toInt());
    configPotential.Ncoefs_angular = std::min<size_t>(MAX_NCOEFS_ANGULAR, settings.value("Potential/Ncoefs_angular", configPotential.Ncoefs_angular).toInt());
    NbodyFile  = settings.value("Potential/NbodyFile").toString();
    configPotential.SymmetryType = getSymmetryTypeByName(settings.value("Potential/Symmetry", "").toString().toStdString());
    configPotential.nbpotEps   = settings.value("Potential/treecodeEps", configPotential.nbpotEps).toDouble();
    configPotential.nbpotTol   = std::min<double>(1, std::max<double>(0, settings.value("Potential/treecodeTheta", configPotential.nbpotTol).toDouble()));
    initPotential();
    // orbit integration section
    configOrbit.accuracyRel= settings.value("Potential/accuracyRelative", configOrbit.accuracyRel).toDouble();
    configOrbit.accuracyAbs= settings.value("Potential/accuracyAbsolute", configOrbit.accuracyAbs).toDouble();
    configOrbit.accTreecode= settings.value("Potential/accuracyTreecode", configOrbit.accTreecode).toDouble();
    configOrbit.treecodeSymmetrizeTimestep= settings.value("Potential/treecodeSymmetrizeTimestep", configOrbit.treecodeSymmetrizeTimestep).toBool();
    configOrbit.adaptiveTimeThreshold = settings.value("Orbit/adaptiveTimeThreshold", configOrbit.adaptiveTimeThreshold).toDouble();
    double icX = settings.value("Orbit/x", 1).toDouble();
    double icY = settings.value("Orbit/y", 0).toDouble();
    double icZ = settings.value("Orbit/z", 0).toDouble();
    double icVx = settings.value("Orbit/vx", 0).toDouble();
    double icVy = settings.value("Orbit/vy", 0).toDouble();
    double icVz = settings.value("Orbit/vz", 0).toDouble();
    initCondE = settings.value("Orbit/E", 0).toDouble();
    useICe = settings.value("Orbit/useE", false).toBool();
    InitData.initCond=CPosVelPoint<double>(icX, icY, icZ, icVx, icVy, icVz);
    InitData.calcLyapunov = settings.value("Orbit/calcLyapunov", false).toBool();
    configCore.intTimeInPeriods = settings.value("Orbit/intTime", configCore.intTimeInPeriods).toDouble();
    configCore.intTimeStepsPerPeriod = settings.value("Orbit/intTimeStep", configCore.intTimeStepsPerPeriod).toDouble();
    configCore.intTimeMaxAdaptive = settings.value("Orbit/intTimeMax", configCore.intTimeMaxAdaptive).toDouble();
    configCore.usePS = settings.value("Orbit/usePS", configCore.usePS).toBool();
    // frequency map
    configCore.fm_numOrbitsStationary = settings.value("Frequency_map/numOrbitsStationary", configCore.fm_numOrbitsStationary).toInt();
    configCore.fm_numOrbitsPrincipalPlane = settings.value("Frequency_map/numOrbitsPrincipalPlane", configCore.fm_numOrbitsPrincipalPlane).toInt();
    configCore.fm_numOrbitsYalpha = settings.value("Frequency_map/numOrbitsYalpha", configCore.fm_numOrbitsYalpha).toInt();
    configCore.fm_numOrbitsRandom = settings.value("Frequency_map/numOrbitsRandom", configCore.fm_numOrbitsRandom).toInt();
    // schwarzschild modelling
    configSchw.ModelType = getSchwModelTypeByName(settings.value("Schwarzschild_model/modelType", "").toString().toStdString());
    configCore.sm_numOrbitsRandom = settings.value("Schwarzschild_model/numOrbitsRandom", configCore.sm_numOrbitsRandom).toInt();
    configCore.sm_linesPerSegment = settings.value("Schwarzschild_model/linesPerSegment", configCore.sm_linesPerSegment).toInt();
    configCore.sm_numAngularCoefs = settings.value("Schwarzschild_model/numAngularCoefs", configCore.sm_numAngularCoefs).toInt();
    configCore.sm_numRadialCoefs = settings.value("Schwarzschild_model/numRadialCoefs", configCore.sm_numRadialCoefs).toInt();
    configCore.sm_Alpha = settings.value("Schwarzschild_model/Alpha", configCore.sm_Alpha).toDouble();
    configCore.sm_numShells = settings.value("Schwarzschild_model/numShells", configCore.sm_numShells).toInt();
    configCore.sm_innerShellMass = settings.value("Schwarzschild_model/innerShellMass", configCore.sm_innerShellMass).toDouble();
    configCore.sm_outerShellMass = settings.value("Schwarzschild_model/outerShellMass", configCore.sm_outerShellMass).toDouble();
    configCore.sm_numSamplingPoints = settings.value("Schwarzschild_model/numSamplingPoints", configCore.sm_numSamplingPoints).toInt();
    configCore.chaoticMinFreqDiff = settings.value("Schwarzschild_model/chaoticMinFreqDiff", configCore.chaoticMinFreqDiff).toDouble();
    configCore.chaoticMinLambda = settings.value("Schwarzschild_model/chaoticMinLambda", configCore.chaoticMinLambda).toDouble();
    configCore.chaoticWeightFactor = settings.value("Schwarzschild_model/chaoticWeightFactor", configCore.chaoticWeightFactor).toDouble();
    configSchw.sm_maxWeight = settings.value("Schwarzschild_model/maxWeight", configSchw.sm_maxWeight).toDouble();
    configSchw.sm_constrainBeta = settings.value("Schwarzschild_model/constrainBeta", configSchw.sm_constrainBeta).toBool();
    configSchw.sm_betaIn = settings.value("Schwarzschild_model/betaIn", configSchw.sm_betaIn).toDouble();
    configSchw.sm_betaOut= settings.value("Schwarzschild_model/betaOut",configSchw.sm_betaOut).toDouble();
    configCore.sm_useBPMPD = settings.value("Schwarzschild_model/useBPMPD", true).toBool();
    if(!QFile::exists(AppDir+"bpmpd.exe")) configCore.sm_useBPMPD=false;  // by default, use fast BPMPD solver, but if it does not exist, resort to GLPK

    return true;
}

void CSmileCore::saveSettings(QString fileName)
{
    if(fileName=="") fileName=AppDir+DEFAULT_SETTINGS_FILE_NAME;
    QSettings settings(fileName, QSettings::IniFormat);
    // potential section
    settings.setValue("Potential/Type", QString::fromStdString(PotentialNames[configPotential.PotentialType]));
    settings.setValue("Potential/DensityModel", QString::fromStdString(DensityNames[configPotential.DensityType]));
    settings.setValue("Potential/N_dim", configPotential.N_dim);
    settings.setValue("Potential/q", configPotential.q);
    settings.setValue("Potential/p", configPotential.p);
    settings.setValue("Potential/Mbh", configPotential.Mbh);
    settings.setValue("Potential/Rc", configPotential.Rc);
    settings.setValue("Potential/Gamma", configPotential.Gamma);
    settings.setValue("Potential/Alpha", configPotential.Alpha);
    settings.setValue("Potential/Omega", configPotential.Omega);
    settings.setValue("Potential/Ncoefs_radial", configPotential.Ncoefs_radial);
    settings.setValue("Potential/Ncoefs_angular", configPotential.Ncoefs_angular);
    settings.setValue("Potential/Symmetry", QString::fromStdString(SymmetryNames[configPotential.SymmetryType]));
    settings.setValue("Potential/NbodyFile", NbodyFile);
    settings.setValue("Potential/treecodeEps", configPotential.nbpotEps);
    settings.setValue("Potential/treecodeTheta", configPotential.nbpotTol);
    settings.setValue("Potential/accuracyRelative", configOrbit.accuracyRel);
    settings.setValue("Potential/accuracyAbsolute", configOrbit.accuracyAbs);
    settings.setValue("Potential/accuracyTreecode", configOrbit.accTreecode);
    settings.setValue("Potential/treecodeSymmetrizeTimestep", configOrbit.treecodeSymmetrizeTimestep);
    // orbit integration section
    settings.setValue("Orbit/x", InitData.initCond.Pos[0]);
    settings.setValue("Orbit/y", InitData.initCond.Pos[1]);
    settings.setValue("Orbit/z", InitData.initCond.Pos[2]);
    settings.setValue("Orbit/vx", InitData.initCond.Vel[0]);
    settings.setValue("Orbit/vy", InitData.initCond.Vel[1]);
    settings.setValue("Orbit/vz", InitData.initCond.Vel[2]);
    settings.setValue("Orbit/E", initCondE);
    settings.setValue("Orbit/useE", useICe);
    settings.setValue("Orbit/intTime", configCore.intTimeInPeriods);
    settings.setValue("Orbit/intTimeStep", configCore.intTimeStepsPerPeriod);
    settings.setValue("Orbit/intTimeMax", configCore.intTimeMaxAdaptive);
    settings.setValue("Orbit/adaptiveTimeThreshold", configOrbit.adaptiveTimeThreshold);
    settings.setValue("Orbit/calcLyapunov", InitData.calcLyapunov);
    settings.setValue("Orbit/usePS", configCore.usePS);
    // frequency map
    settings.setValue("Frequency_map/numOrbitsStationary", configCore.fm_numOrbitsStationary);
    settings.setValue("Frequency_map/numOrbitsPrincipalPlane", configCore.fm_numOrbitsPrincipalPlane);
    settings.setValue("Frequency_map/numOrbitsYalpha", configCore.fm_numOrbitsYalpha);
    settings.setValue("Frequency_map/numOrbitsRandom", configCore.fm_numOrbitsRandom);
    // schwarzschild modelling
    settings.setValue("Schwarzschild_model/modelType", QString::fromStdString(getSchwModelNameByType(configSchw.ModelType)));
    settings.setValue("Schwarzschild_model/numOrbitsRandom", configCore.sm_numOrbitsRandom);
    settings.setValue("Schwarzschild_model/linesPerSegment", configCore.sm_linesPerSegment);
    settings.setValue("Schwarzschild_model/numAngularCoefs", configCore.sm_numAngularCoefs);
    settings.setValue("Schwarzschild_model/numRadialCoefs", configCore.sm_numRadialCoefs);
    settings.setValue("Schwarzschild_model/Alpha", configCore.sm_Alpha);
    settings.setValue("Schwarzschild_model/numShells", configCore.sm_numShells);
    settings.setValue("Schwarzschild_model/innerShellMass", configCore.sm_innerShellMass);
    settings.setValue("Schwarzschild_model/outerShellMass", configCore.sm_outerShellMass);
    settings.setValue("Schwarzschild_model/numSamplingPoints", configCore.sm_numSamplingPoints);
    settings.setValue("Schwarzschild_model/chaoticMinFreqDiff", configCore.chaoticMinFreqDiff);
    settings.setValue("Schwarzschild_model/chaoticMinLambda", configCore.chaoticMinLambda);
    settings.setValue("Schwarzschild_model/chaoticWeightFactor", configCore.chaoticWeightFactor);
    settings.setValue("Schwarzschild_model/constrainBeta", configSchw.sm_constrainBeta);
    settings.setValue("Schwarzschild_model/betaIn",  configSchw.sm_betaIn);
    settings.setValue("Schwarzschild_model/betaOut", configSchw.sm_betaOut);
    settings.setValue("Schwarzschild_model/maxWeight", configSchw.sm_maxWeight);
}

QString getDirName(const QString &fileName)  // return file path ending with "/"
{
    QFileInfo fi(fileName);
    QString fn = fi.absolutePath();
    if(!fn.endsWith('/')) fn += '/';
    return fn;
}

//// ----------- console scripting ------------ //// 

void CSmileCore::runScript(QString scriptFileName)
{
    ConsoleMode=true;
    std::cout << "SMILE console version "VERSION"  build "__DATE__ << std::endl;
    loadSettings();
    initIC();
    if(scriptFileName!="")
    {
        input.setFileName(scriptFileName);
        ConsoleInput=false;
        if(!input.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            std::cout << "Error! cannot read script file " << scriptFileName.toStdString() << std::endl;
            emit scriptDone();
        }
    }
    else
    {
        input.open(stdin, QIODevice::ReadOnly | QIODevice::Text);
        ConsoleInput=true;
    }
    connect(this, SIGNAL(scriptNextLine()), SLOT(scriptProcessLine()));
    connect(this, SIGNAL(signalFreqMapFinished()), SLOT(scriptOrbitLibraryFinished()));
    connect(this, SIGNAL(signalSchwOrbitLibraryFinished()), SLOT(scriptOrbitLibraryFinished()));
    connect(this, SIGNAL(signalSchwOptimizationFinished(const QString&)), SLOT(scriptOptimizationFinished(const QString&)));
    connect(this, SIGNAL(signalNbodyExportFinished(const QString&)), SLOT(scriptInfo(const QString&)));
    connect(this, SIGNAL(signalInfo(const QString&)), SLOT(scriptInfo(const QString&)));
    connect(this, SIGNAL(signalTimer()), SLOT(scriptTimerEvent()));
    timeRuntime.start();
    myTimer=0;
    emit scriptNextLine();
}

void CSmileCore::scriptInfo(const QString& message)
{
    std::cout << message.toStdString() << std::endl;
}

void CSmileCore::timerEvent(QTimerEvent *)
{
    emit signalTimer();
}

void CSmileCore::scriptTimerEvent()
{
    if(!orbitlib) return;
    std::cout << orbitlib->numComplete() << " out of " << orbitlib->size() << " orbits done\015" << std::flush;
}

void CSmileCore::scriptOrbitLibraryFinished()
{
    std::cout << orbitlib->numComplete() << " orbits total, " << 
        QString::number(orbitlib->numComplete()*1000.0/timeElapsed.elapsed(),'g',3).toStdString() << " orbits/s" << std::endl;
    if(myTimer!=0){
        killTimer(myTimer);
        myTimer=0;
    }
    emit scriptNextLine();
}

void CSmileCore::scriptOptimizationFinished(const QString& message)
{
    emit scriptInfo(message);
    emit scriptNextLine();
}

void CSmileCore::scriptProcessLine()
{
    QString Line;
    std::cout << "> ";
    if(!input.atEnd())
        Line = input.readLine().trimmed();
    else
        Line = "Exit";
    if(Line=="" || Line.startsWith("#"))
    {
        emit scriptNextLine();
        return;
    }
    if(Line.startsWith("Exit", Qt::CaseInsensitive) || Line.startsWith("Quit", Qt::CaseInsensitive) || Line.startsWith("Bye", Qt::CaseInsensitive))
    {
        std::cout << "Done!  Total runtime: " << timeRuntime.elapsed()/1000 << " seconds" << std::endl;
        emit scriptDone();
        return;
    }
    if(!ConsoleInput) 
        std::cout << (Line.toStdString()) << std::endl;
    
    // parse input string and perform desired action
    if(Line.startsWith("ReadIni", Qt::CaseInsensitive) || Line.startsWith("ReadConfig", Qt::CaseInsensitive) || Line.startsWith("LoadIni", Qt::CaseInsensitive) || Line.startsWith("LoadConfig", Qt::CaseInsensitive))
    {
        int i1=Line.indexOf("(""");
        int i2=Line.indexOf(""")");
        QString FileName = (i1>=2 && i2>i1+3) ? Line.mid(i1+2, i2-i1-3) : (i1==0&&i2==0 ? DEFAULT_SETTINGS_FILE_NAME : "");
        if(!FileName.isEmpty())
        {
            if(loadSettings(FileName))
            {
                std::cout << "Loaded config from " << FileName.toStdString() << std::endl;
               // delete orbits;
               // orbits=NULL;
                initIC();
            }
            else
                std::cout << "Error! Cannot load config from " << FileName.toStdString() << std::endl;
        }
        else std::cout << "Error! No valid filename specified" << std::endl;
    }
    else if(Line.startsWith("ExportPotential", Qt::CaseInsensitive))
    {
        int i1=Line.indexOf("(""");
        int i2=Line.indexOf(""")");
        QString FileName = Line.mid(i1+2, i2-i1-3);
        if(i1>0 && i2>i1)
        {
            if(exportPotential(FileName))
                std::cout << "Potential exported to " << FileName.toStdString() << std::endl;
            else
                std::cout << "Error! Cannot export potential to " << FileName.toStdString() << std::endl;
        }
        else std::cout << "Error! No filename specified" << std::endl;
    }
    else if(Line.startsWith("BuildFreqMap", Qt::CaseInsensitive))
    {
        BuildFreqMapIC(initCondE, Line.contains("exist", Qt::CaseInsensitive));
        StartOrbitLibrary(SLOT(coreFreqMapFinished()));
        timeElapsed.start();
        return;
    }
    else if(Line.startsWith("SchwBuildOrbitLibrary", Qt::CaseInsensitive))
    {
        BuildSchwIC(Line.contains("exist", Qt::CaseInsensitive));
        StartOrbitLibrary(SLOT(coreSchwOrbitLibraryFinished()));
        timeElapsed.start();
        return;
    }
    else if(Line.startsWith("SchwLinearOptimization", Qt::CaseInsensitive) || Line.startsWith("SchwQuadraticOptimization", Qt::CaseInsensitive) || Line.startsWith("SchwLucyOptimization", Qt::CaseInsensitive))
    {
        if(model!=NULL && orbitlib!=NULL)
        {
            //if(Line.startsWith("SchwLucyOptimization")) operation=CSchwarzschildThread::SM_LUCY;
            SchwStartOptimization(Line.startsWith("SchwQuadraticOptimization"));
            return;
        }
        else
            std::cout << "No Schwarzschild model or orbit library defined" << std::endl;
    }
    //else if(Line.startsWith("", Qt::CaseInsensitive))
    else if(Line.startsWith("ExportOrbits", Qt::CaseInsensitive))
    {
        if(orbitlib!=NULL)
        {
            int i1=Line.indexOf("(""");
            int i2=Line.indexOf(""")");
            QString FileName = Line.mid(i1+2, i2-i1-3);
            if(i1>0 && i2>i1)
            {
                bool withmodel=Line.startsWith("ExportOrbitsCell", Qt::CaseInsensitive);
                if(!exportOrbitLib(FileName, withmodel))
                    std::cout << "Error writing file: " << FileName.toStdString() << std::endl;
            }
            else std::cout << "No filename specified" << std::endl;
        }
        else
            std::cout << "No orbit library defined" << std::endl;
    }
    else if(Line.startsWith("ImportOrbits", Qt::CaseInsensitive))
    {
        int i1=Line.indexOf("(""");
        int i2=Line.indexOf(""")");
        QString FileName = Line.mid(i1+2, i2-i1-3);
        if(i1>0 && i2>i1)
        {
            bool withmodel=Line.startsWith("ImportOrbitsCell", Qt::CaseInsensitive);
            if(!importOrbitLib(FileName, withmodel))
                std::cout << "Error reading file: " << FileName.toStdString() << std::endl;
            else {
                std::cout << orbitlib->numComplete() << " orbits loaded" << std::endl;
                if(withmodel) initSchwModel();
            }
        }
        else std::cout << "No filename specified" << std::endl;
    }
    else if(Line.startsWith("ExportNbody", Qt::CaseInsensitive))
    {
        if(orbitlib)
        {
            int i1=Line.indexOf("(");
            int i2=Line.lastIndexOf(")");
            QString paramstr = Line.mid(i1+1, i2-i1-1);
            QStringList params=paramstr.split(QRegExp("[,\\s]+"));
            if(params.count()<2)
                std::cout << "Syntax: ExportNbody(NumBodies, \"FileName\"[, RefineFactor=0])" << std::endl;
            else{
                int NPoints=params[0].toInt();
                i1=params[1].indexOf("\"");
                i2=params[1].lastIndexOf("\"");
                QString fileName=params[1].mid(i1+1, i2-i1-1);
                int refineFactor=0;
                if(params.count()>2)
                    refineFactor=params[2].toInt();
                if(NPoints==0)
                    std::cout << "Incorrect first parameter: NumBodies" << std::endl;
                else
                if(fileName!="")
                {
                    SchwExportNbody(NPoints, fileName, refineFactor);
                    return;
                }
                else std::cout << "No filename specified" << std::endl;
            }
        }
        else
            std::cout << "No orbit library defined" << std::endl;
    }
    else 
        std::cout << "Invalid command" << std::endl;
    emit scriptNextLine();
}

///////////////////////////////////// Helper threads that run various tasks in parallel to main program /////

// CalcThread - runs single orbit integration
void CCalcThread::run()
{
    orbit->integrateToTime(time);
    orbit->finishIntegration();
    emit calcFinished();
}

// CalcManyThr - runs orbit library integration in parallel
volatile int nthreads=0;
QMutex mutexL;
QMutex mutexC;

void CCalcManyThread::run()
{
    mutexC.lock();
    ++nthreads;
    mutexC.unlock();
    finished=false;
    while(!finished)
    {
        mutexL.lock();
        int no = orbits->findUnstartedOrbit();
        mutexL.unlock();
        if(no>=0)
            orbits->runOrbit(no);
        else finished=true;
    }
    mutexC.lock();
    --nthreads;
    if(nthreads==0)
        emit calcFinished();
    mutexC.unlock();
}

void CCalcManyThread::stopThread()
{
    finished=true;
    orbits->halt();  // will be called "nthread" times, but it doesn't hurt
}

void CSchwarzschildThread::run()
{ 
    if(model==NULL || solver==NULL || orbitlib==NULL) 
        emit calcFinished( "Error!" );
    int result = model->solveOptimizationProblem(orbitlib, solver, orbitPenaltyFnc);
    emit calcFinished( result >= 0 ? 
        QString::fromStdString(model->getStatistics(orbitlib)) :
        QString::fromStdString(solver->errorDescription(result))
    );
}

void CNbodyExportStatThread::run()
{
    handler = new CPointMassSetHandler(points);
    double T=0, W=0;
    handler->computeVirialRatio(&T, &W);
    delete handler;
    handler=NULL;
    delete points;
    emit exportFinished(info + 
        "\nT="+QString::number(T, 'g', 5) + 
        ", W="+QString::number(W, 'g', 5) + 
        "\nVirial ratio="+QString::number(2*T/W, 'g', 5));
}

void CNbodyExportStatThread::stopThread()
{
    if(handler!=NULL) handler->halt();
}

}