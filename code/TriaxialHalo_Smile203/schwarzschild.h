#pragma once
#include "common.h"

namespace smile{

/** Class that keeps density and velocity data for all shell-based variants of Schwarzschild model.
    Velocity data is kept in three arrays of length "numShells", containing fraction of time spent in shell, 
    radial and tangential velocity dispersion.
    Other (model-specific) density data is kept in a single array "densityData", whose meaning depends on the variant of model.
**/
class CSchwInformation: public CBasicInformation
{
public:
    template<typename NumT> CSchwInformation(const std::vector<NumT> &_densityData, const std::vector<NumT> &_shellTime, const std::vector<NumT> &_shellVr, const std::vector<NumT> &_shellVt)
    {
        densityData.reserve(_densityData.size());
        for(size_t c=0; c<_densityData.size(); c++) densityData.push_back(static_cast<smnumtype>(_densityData[c]));
        shellTime.reserve(_shellTime.size());
        for(size_t c=0; c<_shellTime.size(); c++) shellTime.push_back(static_cast<smnumtype>(_shellTime[c]));
        shellVr.reserve(_shellVr.size());
        for(size_t s=0; s<_shellVr.size(); s++) shellVr.push_back(static_cast<smnumtype>(_shellVr[s]));
        shellVt.reserve(_shellVt.size());
        for(size_t s=0; s<_shellVt.size(); s++) shellVt.push_back(static_cast<smnumtype>(_shellVt[s]));
    };
    virtual INFOTYPE infoType() const { return IT_SCHW; };
    virtual std::string toString() const { return std::string(); };
    // get functions
    const vectorn& getDensityData() const { return densityData; }
    const vectorn& getShellTime() const { return shellTime; }
    const vectorn& getShellVr() const { return shellVr; }
    const vectorn& getShellVt() const { return shellVt; }
    // overloaded functions for convenience access
    smnumtype getDensityData(size_t i) const { return densityData.at(i); }
    smnumtype getShellVr(size_t i) const { return shellVr.at(i); }
    smnumtype getShellVt(size_t i) const { return shellVt.at(i); }
    smnumtype getShellTime(size_t i) const { return shellTime.at(i); }
private:
    vectorn densityData;
    vectorn shellTime;
    vectorn shellVr, shellVt;
};

/** Base class for performing Schwarzschild modelling, i.e. choosing the weights of orbits in orbitlibrary 
    so that to satisfy some constraints. These constraints may be e.g. mass in cells of grid in space, or 
    coefficients of spherical-harmonic expansion, or/and kinematical data, etc.
    Orbit library must be created before model, the base class has one virtual method which invokes linear/quadratic solver
    (which should be a descendant of base solver class) and returns the result code from the solver.
    Each derived class has its corresponding orbit runtime function, which has pointers to both an orbit being integrated 
    and the model, and collects the necessary data during orbit integration. The data is then stored in the corresponding 
    descendant class of information container, and the model obtains these data from COrbitInformation items from 
    the orbit library.
    Pointers to potential and orbits must remain valid throughout the lifetime of this object.
**/
class CBasicSchwModel
{
public:
enum MODELTYPE {
    MT_SHELL=16,
    MT_SPHERICAL_HARMONIC=8,
    MT_CLASSIC=1 | MT_SHELL,
    MT_SHGRID=2 | MT_SHELL | MT_SPHERICAL_HARMONIC,
    MT_SHBSE=3 | MT_SHELL | MT_SPHERICAL_HARMONIC
};
    CBasicSchwModel() {};
    virtual ~CBasicSchwModel() {};
    virtual MODELTYPE ModelType() const = 0;
    virtual std::string ModelName() const = 0;
    virtual int solveOptimizationProblem(COrbitLibrary* orbits, const CBasicOptimizationSolver* solver, const CBasicOrbitFilteringFnc* orbitPenaltyFnc) const = 0;
    virtual std::string getStatistics(const COrbitLibrary* orbits) const = 0;   // compute deviation from solution, report other statistical info
};

/** Base class for models in which radial and tangential velocity dispersion are recorded in spherical shells,
    and density data is somehow recorded in a 1-d array of coefficients (implemented differently in derived classes).
    Corresponding orbit runtime data also contains the same information
    The optimization procedure is common for all derived classes, the only difference is in assigning weights to 
    coefficients.
**/
class CBasicShellSchwModel: public CBasicSchwModel
{
public:
    // create Schw.model for given potential model (it is not used after constructor call)
    CBasicShellSchwModel(const CPotential* potential, unsigned int _numShells, double innerShellMass, double outerShellMass);
    virtual int solveOptimizationProblem(COrbitLibrary* orbits, const CBasicOptimizationSolver* solver, const CBasicOrbitFilteringFnc* orbitPenaltyFnc) const;
    unsigned int whichShellR(double radius) const;  // returns index of shell containing given radius
    unsigned int whichShellE(double energy) const;  // returns index of shell containing given energy
    // get functions
    double getShellRadius(unsigned int index) const { return shellRadius.at(index); }  // return index of shell in which a given radius lies, or shellRadius.size() if it is beyond the last shell
    double getShellEnergy(unsigned int index) const { return shellEnergy.at(index); }
    double getInteriorMass(unsigned int index) const { return massInteriorToShell.at(index); }
    unsigned int getNumShellsKinem() const { return static_cast<unsigned int>(shellRadius.size()); }
    unsigned int getNumCoefsDens() const { return static_cast<unsigned int>(densityData.size()); }
    virtual std::string getStatistics(const COrbitLibrary* orbits) const;   // compute deviation from solution, report other statistical info
    vectord calcPenalty(const COrbitLibrary* orbits) const;
    const vectord& getCoefsDens() const { return densityData; }
    virtual double getNormFactor(unsigned int /*index*/) const { return 1; }
protected:
    double totalMass;      // mass of entire model, approximate
    vectord shellRadius;   // shells for kinematic data
    vectord shellEnergy;
    vectord massInteriorToShell;  // for each radius of a shell, records mass within all shells up to and including this one (approximate)
    vectord densityData;   // implementation-specific density data
};

/** Classic grid-based Schwarzschild model, 
    with division of the positive octant of space into three segments and their further division into NxN cells.
**/
class CSchwModelClassic: public CBasicShellSchwModel
{
public:
    CSchwModelClassic(const CPotential* potential, unsigned int _numShells, double innerShellMass, double outerShellMass, unsigned int _linesPerSegment);
    virtual MODELTYPE ModelType() const { return MT_CLASSIC; };
    virtual std::string ModelName() const { return myName(); };
    static std::string myName() { return "Classic"; };
    unsigned int whichCell(const double X, const double Y, const double Z) const;   // return index of cell containing given coordinates, called from COrbitRuntimeSchwClassic to compute time spent in each cell
    void cellCenter(const unsigned int cell, double& X, double& Y, double& Z) const; // return coordinates of geometric center of a cell
    // get functions
    unsigned int getCellsPerShell() const { return cellsPerShell; }
    virtual double getNormFactor(unsigned int index) const;
private:
    const unsigned int linesPerSegment;    // Schwarzschild's scheme for dividing octant into three segments, each of them divided into linesPerSegment^2 squares
    const unsigned int cellsPerSegment;
    const unsigned int cellsPerShell;
    const unsigned int numCells;
    vectord cellBoundary;   // angles dividing segment into linesPerSegment parts
    void initCellMass(const CPotential* potential);  // called from constructor
};

/** Base class for variants of Schw.model based on spherical-harmonic expansion of density of both target model and orbits.
    Two derived classes share several virtual methods, and the runtime function computing coefs of expansion for an orbit 
    is the same for both variants.
**/
class CSchwModelSH: public CBasicShellSchwModel
{
public:
    CSchwModelSH(const CPotential* potential, unsigned int _numShells, double innerShellMass, double outerShellMass, unsigned int _Nang) :
      CBasicShellSchwModel(potential, _numShells, innerShellMass, outerShellMass), Nang(_Nang) {};
    virtual vectord computeDensityData(const CPointMassSetDouble& traj) const=0;  // evaluate coefs for orbit trajectory data, called from COrbitRuntimeSchwSH, implemented differently in derived classes
    // get functions
    unsigned int getNumCoefsAtRadius() const { return (Nang/2+1)*(Nang/2+2)/2; }
protected:
    const unsigned int Nang;     // order of expansion in angular harmonics (lmax)
};

/** Variant of Schw.model based on spherical-harmonic expansion of potential coefficients at radial grid points.
    Only even non-negative harmonics in l,m are used because of triaxial symmetry. 
    Nang is the expansion order (lmax), and the number of coefficients at each radius is (Nang/2+1)*(Nang/2+2)/2
    Note: the potential in which the model is evolved, needs not to be itself a spline expansion,
    but if it is, then the model coefs are initialized from the potential rather than by computing integrals over density
    (which are expensive). 
**/
class CSchwModelSHGrid: public CSchwModelSH
{
public:
    CSchwModelSHGrid(const CPotential* potential, unsigned int _numShells, double innerShellMass, double outerShellMass, unsigned int _Nang);
    virtual MODELTYPE ModelType() const { return MT_SHGRID; };
    virtual std::string ModelName() const { return myName(); };
    static std::string myName() { return "SHGrid"; };
    virtual vectord computeDensityData(const CPointMassSetDouble& traj) const;  // evaluate coefs for orbit trajectory data, called from COrbitRuntimeSchwSH
    // get functions
    virtual double getNormFactor(unsigned int index) const;
private:
    const unsigned int numGrid;  // number of grid points for computing expansion coefficients (may be different from numShells used in recording velocity dispersion data)
    vectord GridRadii;           // expansion coefficients are evaluated at given radii
};

/** Variant of Schw.model based on basis-set expansion in terms of Zhao(1995) basis set for density of both required profile and of orbits.
    Nrad+1 is the number of radial functions and Nang is the order of angular expansion in spherical harmonics, the number of angular coefs 
    for each radial function is (Nang/2+1)*(Nang/2+2)/2. Only even non-negative harmonics in l,m are used because of triaxial symmetry. 
    The density profile is thus fitted in 'nonlocal' sense, but the kinematic data (radial and tangential velocity dispersions) 
    are still evaluated at a radial grid, so the number of radial shells must also be provided (which is unrelated to Nrad).
    Alpha is the parameter controlling the family of basis functions, same as in CPotentialBSE.
    Note: the potential in which the model is evolved, needs not to be itself a BSE expansion,
    but if it is, then the model coefs are initialized from the potential rather than by computing integrals over density
    (which are expensive). 
**/
class CSchwModelSHBSE: public CSchwModelSH
{
public:
    CSchwModelSHBSE(const CPotential* potential, unsigned int _numShells, double innerShellMass, double outerShellMass, unsigned int _Nang, unsigned int _Nrad, double _Alpha);
    virtual MODELTYPE ModelType() const { return MT_SHBSE; };
    virtual std::string ModelName() const { return myName(); };
    static std::string myName() { return "BSE"; };
    virtual vectord computeDensityData(const CPointMassSetDouble& traj) const;
    virtual double getNormFactor(unsigned int index) const;
private:
    const unsigned int Nrad;     // number of radial basis functions
    double Alpha;                // parameter controlling the form of Zhao(1995) basis functions
};

/** base class for runtime function recording spatial and kinematic data for 
    several variants of shell-based Schwarzschild modelling approaches.
    This class only records kinematic data (time spent in each shell and radial/tangential velocity dispersion in a shell),
    derived classes should implement recording of density data.
**/
class COrbitRuntimeBasicShellSchw: public CBasicOrbitRuntimeFnc
{
public:
    COrbitRuntimeBasicShellSchw(const COrbit* _orbit, const CBasicShellSchwModel* _model);
    virtual FNCTYPE FncType() const { return FT_SCHW; };
    virtual void Timestep(const double told, const double t, const double y[]);
    virtual void Finish();
protected:
    const CBasicShellSchwModel* model;   // pointer to Schw.model
    double timeAccuracy;
    vectord shellTime, shellVr, shellVt; // time spent in each shell and corresponding radial/tangential velocity dispersions
private:
    unsigned int shellprev;
    void calcShellTime(double tlow, double tupp, unsigned int slow, unsigned int supp, double vr2, double vt2);  // record time-in-shell with fraction-of-period accuracy
};

/** Runtime function corresponding to classic Schwarzschild model **/
class COrbitRuntimeSchwClassic: public COrbitRuntimeBasicShellSchw
{
public:
    COrbitRuntimeSchwClassic(const COrbit* _orbit, const CSchwModelClassic* _model);
    virtual void Timestep(const double told, const double t, const double y[]);
    virtual void Finish();
    virtual CBasicInformation* getData() const { return new CSchwInformation(cellTime, shellTime, shellVr, shellVt); }
private:
    const CSchwModelClassic* model;  // points to the same object as COrbitRuntimeBasicSchw::model, but with a correct type cast
    vectord cellTime;                // time spent in each cell of Schwarzschild's model
    unsigned int cellprev;           // cell index at previous timestep (for recording time-spent-in-cell)
    void calcCellTime(double tlow, double tupp, unsigned int clow, unsigned int cupp);  // record time-in-cell with fraction-of-period accuracy
};

/** corresponding creator class for classic Schw runtime function **/
class COrbitRuntimeSchwClassicCreator: public CBasicOrbitRuntimeFncCreator
{
public:
    COrbitRuntimeSchwClassicCreator(const CSchwModelClassic* _model) : model(_model) {};
    virtual CBasicOrbitRuntimeFnc* createRuntimeFnc(const COrbit* orbit) const { return new COrbitRuntimeSchwClassic(orbit, model); };
private:
    const CSchwModelClassic* model;
};

/** Runtime function corresponding to of Schw.model based on spherical-harmonic expansion of potential coefficients 
    at radial grid points or for some radial basis functions (serves both variants)
**/
class COrbitRuntimeSchwSH: public COrbitRuntimeBasicShellSchw
{
public:
    COrbitRuntimeSchwSH(const COrbit* _orbit, const CSchwModelSH* _model);
    virtual void Timestep(const double told, const double t, const double y[]);
    virtual void Finish();  // everything is done here
    virtual CBasicInformation* getData() const { return new CSchwInformation(SHcoefs, shellTime, shellVr, shellVt); }
private:
    const CSchwModelSH* model;       // points to the same object as COrbitRuntimeBasicSchw::model, but with a correct type cast
    CPointMassSetDouble traj;        // keeps trajectory data to compute SH coefficients when finished
    vectord SHcoefs;                 // coefficients of sph-harm expansion
    double timeStep;                 // defined maximum interval between stored trajectory points (need not be stored equidistantly in time)
};

/** corresponding creator class for grid-SH Schw runtime function **/
class COrbitRuntimeSchwSHCreator: public CBasicOrbitRuntimeFncCreator
{
public:
    COrbitRuntimeSchwSHCreator(const CSchwModelSH* _model) : model(_model) {};
    virtual CBasicOrbitRuntimeFnc* createRuntimeFnc(const COrbit* orbit) const { return new COrbitRuntimeSchwSH(orbit, model); };
private:
    const CSchwModelSH* model;
};

/** Orbit filtering function that evaluates chaotic properties of an orbit, based on threshold in frequency diffusion rate 
    and in Lyapunov exponent value, multiplied by some predefined factor.
    May be used for telling regular from chaotic orbits in e.g. computing orbit population of a model, 
    or for assigning positive or negative bias for penalizing or favouring orbits of particular type in Schwarzschild model.
    Positive factor favours regular orbits, negative - chaotic.
**/
class CChaosOrbitFilteringFnc: public CBasicOrbitFilteringFnc
{
public:
    CChaosOrbitFilteringFnc(double _chaoticWeightFactor, double _chaoticMinLambda, double _chaoticMinFreqDiff) :
      chaoticWeightFactor(_chaoticWeightFactor), chaoticMinLambda(_chaoticMinLambda), chaoticMinFreqDiff(_chaoticMinFreqDiff) {};
    virtual double eval(const COrbitDesc* orbit) const;
private:
    double chaoticWeightFactor;
    double chaoticMinLambda;
    double chaoticMinFreqDiff;
};

/** Orbit filtering function based on whether an orbit belongs to a given energy shell.
    if model==NULL or numShell==0, any orbit is acceptable, otherwise only if orbit energy is within numShell-1 energy shell.
**/
class CShellOrbitFilteringFnc: public CBasicOrbitFilteringFnc
{
public:
    CShellOrbitFilteringFnc(const CBasicShellSchwModel* _model, unsigned int _numShell) :
      model(_model), numShell(_numShell) {};
    virtual double eval(const COrbitDesc* orbit) const;
private:
    const CBasicShellSchwModel* model;
    unsigned int numShell;
};

/** structure that contains tunable parameters for schw.modelling **/
struct CConfigSchw 
{
    CBasicSchwModel::MODELTYPE ModelType;
    double sm_maxWeight;             // upper bound on orbit weight
    bool sm_constrainBeta;           // whether to constrain anisotropy
    double sm_betaIn, sm_betaOut;    // anisotropy coefficient requested (values for the inner and outer shells, linearly interpolated in between)
};
extern CConfigSchw configSchw;

/** Correspondence between Schw.model names and types **/
std::string getSchwModelNameByType(CBasicSchwModel::MODELTYPE ModelType);
CBasicSchwModel::MODELTYPE getSchwModelTypeByName(std::string ModelName);

}