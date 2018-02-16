/// this file defines COrbit class which performs actual orbit integration and classification, and some auxiliary classes
#pragma once
#include "common.h"
#include <complex>

namespace smile{

class CBasicOrbitIntegrator;

/// general parameters for orbit integration and classification (common for all orbits)
struct CConfigOrbit
{
    // integration accuracy
    double accuracyRel;              // relative accuracy parameter to 8th order Runge-Kutta integrator
    double accuracyAbs;              // absolute accuracy parameter 
    double accTreecode;              // factor in leap-frog timestep selection for treecode potential
    bool treecodeSymmetrizeTimestep; // whether to use time-symmetrized leapfrog (Hut, Makino & McMillan 1995) - doubles the number of evaluations but better conserves energy
    // orbit integration and classification data
    double adaptiveTimeThreshold;    // if adaptive int.time is used (with Pfenniger's method), controls max difference in cell times for subsequent intervals
    bool freqStrictOrder;            // whether to require leading frequencies to be sorted so that w_x <= w_y <= w_z  (with small tolerance described below)
    double minFreqX;                 // search for leading frequency in X coordinate which is not less than MIN_FREQ_X
    double minFreqYZ;                // search for leading freq. in Y and Z which are not less than MIN_FREQ_YZ * (leading frequency in X and Y correspondingly)
};
extern CConfigOrbit configOrbit;

/// defines orbit with given initial conditions, with routines to do integration in the given potential, analysis, etc.
class COrbit
{
public:
enum ORBITSTATE{
    OS_INITIALIZED,
    OS_RUNNING,
    OS_NEEDTOTERMINATE,
    OS_DONE
};
    template<typename NumT> explicit COrbit(const COrbitInitData<NumT> &_InitData, const vectorRuntimeFncCreators* creators=NULL);
    template<typename NumT> COrbit(const CPotential* _potential, const double _timeStep, const double _timeUnit, const std::vector< CPosVelPoint<NumT> > &_traj);
    ~COrbit();
private:
    const COrbitInitData<double> InitData;
    const CPotential* potential;     // copied from InitData for convenience
    const unsigned int N_dim;        // copied from potential->N_dim, for convenience
    volatile ORBITSTATE state;
    CPosVelPoint<double> endCond;    // final coordinates/velocities (after integration)
    double intTime;                  // integration time 
    vectorRuntimeFncs RuntimeFnc;  // holds all attached 'timestep functions' (if any). They are destroyed along with COrbit!!!
    CPosVelPoint<double> devVector;  // holds value of deviation vector between integration periods
    //std::string DescriptionLFCCdiff; // description based on rate of LFCC diffusion

public:
    void integrateToTime(double time);   // start or continue integration up to given time
    void integrateAdaptiveTime(double minTime, double maxTime);  // perform integration for a variable time interval, using Pfenniger's method for checking convergence of orbit cell occupation
    void finishIntegration();        // call Finish() for all attached timestep functions and set state to DONE
    void halt() { if(state==OS_RUNNING) state=OS_NEEDTOTERMINATE; };
    double getInterpolatedTrajectory(unsigned int c, double t) const;   // returns c-th coordinate interpolated at time t within currect timestep interval
    // get functions
    ORBITSTATE getState() const { return state; };
    const COrbitInitData<double> getInitData() const { return InitData; };
    const CPotential* getPotential() const { return potential; };
    double getIntTime() const { return state==OS_RUNNING ? xout : intTime; };
    const CBasicOrbitRuntimeFnc* getRuntimeFnc(unsigned int index) const { return RuntimeFnc.at(index); };
    const CBasicInformation * getNewInfo(unsigned int index) const { return RuntimeFnc.at(index)->getData(); };
    unsigned int getInfoNum() const { return static_cast<unsigned int>(RuntimeFnc.size()); };
    std::string toString(bool withRuntimeFncInfo=true) const;
private:
    //vectord latestcelltimes;         // used in Pfenniger method: holds celltimes on the latest interval of integration ///!!!
    //double latestIntTime;            // duration of latest integration interval
    //double diffCellTimes();          // calculate difference in celltimes between latest interval and whole integration time
    unsigned int nfcn, nstep, naccpt, nrejct;   // statistics of force evaluations, timesteps, etc.
    int drdtsign;                    // direction of motion at prev timestep
    bool needrenorm, allowrenorm;    // flag set when renormalization of deviation vector required (usually close to pericenter), and when it is allowed to do (in apocenter)
    double lnwadd;                   // addition term - required for renormalization so that |w| remains not too large
    // dop853.c encapsulated variables and functions - 7/8 order Runge-Kutta integrator
    double    hout, xold, xout;
    unsigned int nrds;
    unsigned int* indir;
    double    *rcont1, *rcont2, *rcont3, *rcont4;
    double    *rcont5, *rcont6, *rcont7, *rcont8;

    void fcn(unsigned int n, double x, double *y, double *f);   // rhs of the differential equation
    void solout(long nr, double xold, double x, double y[], unsigned int n);  // performs action at each timestep

    int dop853
     (unsigned n,      /* dimension of the system <= UINT_MAX-1*/
      //  FcnEqDiff fcn,   /* function computing the value of f(x,y) */
      double x,        /* initial x-value */
      double* y,       /* initial values for y */
      double xend,     /* final x-value (xend-x may be positive or negative) */
      double* rtoler,  /* relative error tolerance */
      double* atoler,  /* absolute error tolerance */
      int itoler,      /* switch for rtoler and atoler */
      //  SolTrait solout, /* function providing the numerical solution during integration */
      int iout,        /* switch for calling solout */
      //  FILE* fileout,   /* messages stream */
      double uround,   /* rounding unit */
      double safe,     /* safety factor */
      double fac1,     /* parameters for step size selection */
      double fac2,
      double beta,     /* for stabilized step size control */
      double hmax,     /* maximal step size */
      double h,        /* initial step size */
      long nmax,       /* maximal number of allowed steps */
      int meth,        /* switch for the choice of the coefficients */
      long nstiff,     /* test for stiffness */
      unsigned int nrdens, /* number of components for which dense outpout is required */
      unsigned int* icont, /* indexes of components for which dense output is required, >= nrdens */
      unsigned int licont  /* declared length of icon */
     );
    double contd8      /* returns interpolated coordinate at time x within current timestep */
     (unsigned int ii, /* index of desired component */
      double t         /* approximation at t */
    ) const;
    // internal integration initialization routine
    double hinit (unsigned int n, /*FcnEqDiff fcn,*/ double x, double* y,
      double posneg, double* f0, double* f1, double* yy1, int iord,
      double hmax, double* atoler, double* rtoler, int itoler);
    
    // internal integration routine
    int dopcor (unsigned int n, /*FcnEqDiff fcn,*/ double x, double* y, double xend,
      double hmax, double h, double* rtoler, double* atoler,
      int itoler, /*FILE* fileout, SolTrait solout,*/ int iout,
      unsigned int nmax, double uround, int meth, long nstiff, double safe,
      double beta, double fac1, double fac2, unsigned* icont);

    // low-order leapfrog method provided for treecode, which has discontinuous potential/forces and is poorly handled by RK integrator
    bool leapfrog;
    double yprev[4*N_DIM], aprev[2*N_DIM], anew[N_DIM];
    int intlf(unsigned int n, double x, double* y, double xend);
    double contlf      /* returns interpolated coordinate at time x within current timestep */
     (unsigned int ii, /* index of desired component */
      double x         /* approximation at x */
    ) const;

};

/// info classes ///
/// template class for orbit analysis results
template<typename NumT> class COrbitInformation: public CBasicInformation
{
public:
    explicit COrbitInformation(const std::string &_Description) :
      Description(_Description), Einit(0), Ediff(0), lfccdiff(0), lambda(0) 
      { for(unsigned int d=0; d<N_DIM; d++) { lf[d]=0; inertia[d]=0; Lavg[d]=0; Lvar[d]=0; } };
    COrbitInformation(const std::string &_Description, const NumT _Einit, const NumT _Ediff, const NumT _lf[N_DIM], const NumT _lfccdiff, const NumT _lambda, const NumT _inertia[N_DIM], const NumT _Lavg[N_DIM], const NumT _Lvar[N_DIM]) :
      Description(_Description), Einit(_Einit), Ediff(_Ediff), lfccdiff(_lfccdiff), lambda(_lambda) 
      { for(unsigned int d=0; d<N_DIM; d++) { lf[d]=_lf[d]; inertia[d]=_inertia[d]; Lavg[d]=_Lavg[d]; Lvar[d]=_Lvar[d]; } };
    COrbitInformation(const COrbitInformation<float> &src) :
      Description(src.Description), Einit(static_cast<NumT>(src.Einit)), Ediff(static_cast<NumT>(src.Ediff)), lfccdiff(static_cast<NumT>(src.lfccdiff)), lambda(static_cast<NumT>(src.lambda)) 
      { for(unsigned int d=0; d<N_DIM; d++) { lf[d]=static_cast<NumT>(src.lf[d]); inertia[d]=static_cast<NumT>(src.inertia[d]); Lavg[d]=static_cast<NumT>(src.Lavg[d]); Lvar[d]=static_cast<NumT>(src.Lvar[d]); } };
    COrbitInformation(const COrbitInformation<double> &src) :
      Description(src.Description), Einit(static_cast<NumT>(src.Einit)), Ediff(static_cast<NumT>(src.Ediff)), lfccdiff(static_cast<NumT>(src.lfccdiff)), lambda(static_cast<NumT>(src.lambda)) 
      { for(unsigned int d=0; d<N_DIM; d++) { lf[d]=static_cast<NumT>(src.lf[d]); inertia[d]=static_cast<NumT>(src.inertia[d]); Lavg[d]=static_cast<NumT>(src.Lavg[d]); Lvar[d]=static_cast<NumT>(src.Lvar[d]); } };
    virtual INFOTYPE infoType() const { return sizeof(NumT)==sizeof(float) ? IT_TRAJ_ANALYSIS : IT_UNKNOWN; };  // streams support only float instantiations of this class
    virtual std::string toString() const;
    // get functions
    std::string getDescription() const { return Description; }
    NumT getEinit() const { return Einit; }
    NumT getEdiff() const { return Ediff; }
    NumT getlf(unsigned int index) const { if(index<N_DIM) return lf[index]; else return -1; }
    NumT getlfccdiff() const { return lfccdiff; }
    NumT getlambda() const { return lambda; }
    NumT getinertia(unsigned int index) const { if(index<N_DIM) return inertia[index]; else return -1; }
    NumT getLavg(unsigned int index) const { if(index<N_DIM) return Lavg[index]; else return -1; }
    NumT getLvar(unsigned int index) const { if(index<N_DIM) return Lvar[index]; else return -1; }
private:
    std::string Description;     // description of orbit type
    NumT Einit, Ediff; 
    NumT lf[N_DIM];              // leading frequencies
    NumT lfccdiff;               // rate of LFCC diffusion
    NumT lambda;                 // Lyapunov exponent
    NumT inertia[N_DIM];         // inertia tensor components (only diagonal in Cartesian coordinates)
    NumT Lavg[N_DIM];            // average components of angular momentum
    NumT Lvar[N_DIM];            // variation in angular momentum components
};
// explicit instantiations of float<->double conversion initializers
/*template COrbitInformation<float >::COrbitInformation(const smile::COrbitInformation<float > &src);
template COrbitInformation<float >::COrbitInformation(const smile::COrbitInformation<double> &src);
template COrbitInformation<double>::COrbitInformation(const smile::COrbitInformation<float > &src);
template COrbitInformation<double>::COrbitInformation(const smile::COrbitInformation<double> &src);
*/
template<typename NumT> class CPoincareInformation: public CBasicInformation
{
public:
    CPoincareInformation(const std::vector< std::pair<NumT, NumT> > &_ps) : ps(_ps) {};
    virtual INFOTYPE infoType() const { return IT_POINCARE; };
    virtual std::string toString() const;
    // get functions
    unsigned int size() const { return static_cast<unsigned int>(ps.size()); }
    const std::pair<NumT, NumT>& at(unsigned int index) const { return ps.at(index); }
private:
    std::vector< std::pair<NumT, NumT> > ps;
};

template<typename NumT> class CPericenterInformation: public CBasicInformation
{
public:
    CPericenterInformation(const NumT _fitIntercept, const NumT _fitSlope, const NumT _fitScatter, const NumT _fitSignificance, const NumT _Lcirc2) : 
      fitIntercept(_fitIntercept), fitSlope(_fitSlope), fitScatter(_fitScatter), fitSignificance(_fitSignificance), Lcirc2(_Lcirc2) {};
    virtual INFOTYPE infoType() const { return sizeof(NumT)==sizeof(float) ? IT_PERICENTER : IT_UNKNOWN; };
    virtual std::string toString() const;
    // get functions
    void getParams(NumT* _fitIntercept, NumT* _fitSlope, NumT* _fitScatter, NumT* _fitSignificance, NumT* _Lcirc2) const 
      { *_fitIntercept=fitIntercept; *_fitSlope=fitSlope; *_fitScatter=fitScatter; *_fitSignificance=fitSignificance; *_Lcirc2=Lcirc2; };
private:
    NumT fitIntercept, fitSlope, fitScatter, fitSignificance;  // fitting coefficients for distribution of squared angular momentum at pericenter
    NumT Lcirc2;            // approximate squared ang.mom. of a circular orbit
};
// explicit instantiations of a template function
//template std::string CPericenterInformation<double>::toString() const;
//template std::string CPericenterInformation<float>::toString() const;

template<typename NumT> class CTrajSampleInformation: public CBasicInformation
{
public:
    explicit CTrajSampleInformation(const std::vector< CPosVelPoint<NumT> > &src) : trajSample(src) {};
    virtual INFOTYPE infoType() const { return sizeof(NumT)==sizeof(float) ? IT_TRAJ_SAMPLE : IT_UNKNOWN; };
    virtual std::string toString() const;
    // get functions
    const CPosVelDataFloat& getTraj() const { return trajSample; }
private:
    CPosVelDataFloat trajSample;
};

// a couple of examples of runtime functions

struct CSpectralLine
{
    double freq;  // frequency
    double ampl;  // amplitude
    double phase; // phase (complex amplitude = ampl * exp(I*phase)
    bool precise; // whether it was obtained using more precise Hunter's method
    CSpectralLine(double Freq=0, double Ampl=0, double Phase=0, bool Precise=true)
    {   freq=Freq; ampl=Ampl; phase=Phase; precise=Precise; }
    inline bool operator<(CSpectralLine S) const
    {   return (ampl>S.ampl); }    // sort in reverse order of amplitude
    inline operator double() const { return freq; }
};
typedef std::vector<CSpectralLine> CSpectrum;
typedef std::complex<double> complexd;

/// Class for recording orbit trajectory and doing various analysis (frequency, Lyapunov exponent, ...)
class COrbitRuntimeTrajectory: public CBasicOrbitRuntimeFnc
{
public:
    COrbitRuntimeTrajectory(const COrbit* _orbit);
    template<typename NumT> COrbitRuntimeTrajectory(const COrbit* _orbit, const std::vector< CPosVelPoint<NumT> > &_traj);   // create from existing trajectory
    virtual FNCTYPE FncType() const {return FT_TRAJ_ANALYSIS; };
    virtual void Timestep(const double told, const double t, const double y[]);
    virtual CBasicInformation* getData() const;

    template<typename NumT> COrbitInformation<NumT>* analyzeOrbit(size_t start, size_t count, CSpectrum *spectrum=NULL, vectord *fa=NULL) const;
    // get functions
    unsigned int getTrajSize() const { return static_cast<unsigned int>(traj.size()); };
    const CPosVelPoint<double>& getTraj(unsigned int index) const { return traj[index]; }
    double getLnDevVec(unsigned int index) const { return lnw.at(index); }
    double getMaxDist(size_t start=0, size_t count=0) const;  // return max distance between adjacent points in given interval
private:
    const COrbitInitData<double> InitData;
    const CPotential* potential;     // copied from InitData for convenience
    const unsigned int N_dim;        // copied from potential->N_dim, for convenience

    CPosVelDataDouble traj;
    vectord lnw;                     // ln(|w|) - norm of deviation vector, lambda - finite-time Lyapunov exponent
    double Einit, Ediff;

    double calcLambda() const;             // estimate Lyapunov exponent
    std::string performFrequencyAnalysis(size_t start, size_t count, CSpectrum sl[N_DIM], vectord fa[N_DIM], double Lavg[N_DIM], double Lvar[N_DIM]) const;   // extract frequencies and perform orbit classification; by default on whole orbit, but may use part of orbit as well. Output spectral lines in sl, if necessary - output whole spectrum in fa 
    void findLines(size_t start, size_t count, CSpectrum sl[N_DIM], vectord fa[N_DIM]) const;  // perform Fourier transform of trajectory, find spectral lines and store in sl[N_dim]; if fa!=NULL, store amplitude of spectrum in fa[N_dim]
    void calcLFCCoverInterval(size_t start, size_t count, double leadfreq[N_DIM]) const;   // calculate LFCCs over a part of orbit: output sli=double[N_DIM] - array of LFCCs
    double calcLFCCdiffusion(size_t start, size_t count) const;             // evaluate rate of diffusion of fundamental frequencies over given interval (whole orbit by default)
    void calcInertiaOverInterval(size_t start, size_t count, double inertia[N_DIM]) const;   // calculate LFCCs over a part of orbit: output sli=double[N_DIM] - array of LFCCs
};

/// corresponding creator class for Trajectory analysis timestep fnc
class COrbitRuntimeTrajectoryCreator: public CBasicOrbitRuntimeFncCreator
{
public:
    virtual CBasicOrbitRuntimeFnc* createRuntimeFnc(const COrbit* orbit) const { return new COrbitRuntimeTrajectory(orbit); };
};

// Trajectory sampler
class COrbitRuntimeTrajSample: public CBasicOrbitRuntimeFnc
{
public:
    COrbitRuntimeTrajSample(const COrbit* _orbit, size_t _numSamplingPoints);
    virtual FNCTYPE FncType() const {return FT_TRAJ_SAMPLE; };
    virtual void Timestep(const double told, const double t, const double y[]);
    virtual void Finish();   // draw sampling points from trajectory
    virtual CBasicInformation* getData() const { return new CTrajSampleInformation<float>(trajSample); };
private:
    size_t numSamplingPoints;
    const double timeStep;
    const unsigned int N_dim;
    CPosVelDataFloat trajEntire;
    CPosVelDataFloat trajSample;
};

/// corresponding creator class for traj.sampletimestep fnc
class COrbitRuntimeTrajSampleCreator: public CBasicOrbitRuntimeFncCreator
{
public:
    explicit COrbitRuntimeTrajSampleCreator(unsigned int _numSamplingPointsDefault, const CPointCountSetFloat* _numSamplingPointsSpecial) : 
      numSamplingPointsDefault(_numSamplingPointsDefault)
    { 
        if(numSamplingPointsDefault==0) my_error("Trajectory sampler: numSamplingPoints=0"); 
        if(_numSamplingPointsSpecial!=NULL) numSamplingPointsSpecial.assign(_numSamplingPointsSpecial->begin(), _numSamplingPointsSpecial->end());
    };
    virtual CBasicOrbitRuntimeFnc* createRuntimeFnc(const COrbit* orbit) const;
private:
    const unsigned int numSamplingPointsDefault;
    CPointCountSetFloat numSamplingPointsSpecial;
};

// records points in Poincare section  given by condition: when coord2=0 and coord2'>0, record a pair (coord1, coord1')
class COrbitRuntimePoincare: public CBasicOrbitRuntimeFnc
{
public:
    COrbitRuntimePoincare(const COrbit* _orbit, unsigned int _PS1, unsigned int _PS2): CBasicOrbitRuntimeFnc(_orbit), PS1(_PS1), PS2(_PS2) {};
    virtual FNCTYPE FncType() const {return FT_POINCARE; };
    virtual void Timestep(const double told, const double t, const double y[]);
    virtual CBasicInformation* getData() const { return new CPoincareInformation<float>(ps); };
private:
    const unsigned int PS1, PS2;              // index of coordinates to use in Poincare section
    std::vector< std::pair<float, float> > ps;                     // Poincare section  (x - xdot plane)
};

/// corresponding creator class for Poincare timestep fnc
class COrbitRuntimePoincareCreator: public CBasicOrbitRuntimeFncCreator
{
public:
    COrbitRuntimePoincareCreator(unsigned int _PS1, unsigned int _PS2) : PS1(_PS1), PS2(_PS2) {};
    virtual CBasicOrbitRuntimeFnc* createRuntimeFnc(const COrbit* orbit) const { return new COrbitRuntimePoincare(orbit, PS1, PS2); };
private:
    const unsigned int PS1, PS2;
};
        
// records pericenter passages
class COrbitRuntimePericenter: public CBasicOrbitRuntimeFnc
{
public:
    COrbitRuntimePericenter(const COrbit* _orbit);
    template<typename NumT> COrbitRuntimePericenter(const COrbit* _orbit, const std::vector< CPosVelPoint<NumT> > &_traj);   // create from existing trajectory
    virtual FNCTYPE FncType() const {return FT_PERICENTER; };
    virtual void Timestep(const double told, const double t, const double y[]);
    virtual void Finish();   // fit pericenter distribution
    virtual CBasicInformation* getData() const { return new CPericenterInformation<float>((float)fitIntercept, (float)fitSlope, (float)fitScatter, (float)fitSignificance, (float)Lcirc2); };
    size_t size() const { return pericentert.size(); };
    double gett(size_t index) const { return pericentert[index]; };
    double getr(size_t index) const { return pericenterr[index]; };
    double getl2(size_t index)const { return pericenterl2[index]; };
private:
    int drdtsign;  // whether approaching or receding from center on the previous timestep
    vectord pericentert, pericenterr;  // times and radii of consequent pericenter passages
    vectord pericenterl2;              // L2 - squared magnitude of angular momentum
    double fitIntercept, fitSlope;     // fitting coefficients for L2 distribution
    double fitScatter, fitSignificance;// scatter in linear fit (less=better, >0.5-poor fit) and significance of non-zero value of intercept (in sigmas)
    double Lcirc2;                     // approximate squared ang.mom. of a circular orbit with this energy
    double Lavg[N_DIM], L2avg[N_DIM];  // accumulator for average value of i-th component of angular momentum and its square (to determine if there is a conserved component)
};

class COrbitRuntimePericenterCreator: public CBasicOrbitRuntimeFncCreator
{
public:
    virtual CBasicOrbitRuntimeFnc* createRuntimeFnc(const COrbit* orbit) const { return new COrbitRuntimePericenter(orbit); };
};

/**** utility function to ensure that n does not have large prime divisors, if needed, decrease it ****/
void tuneNonPrime(size_t &n);

/**** Utility function to find a 1:1 periodic orbit in a plane specified by coord1:coord2 ****/
bool findPeriodicOrbit(const CPotential* potential, double E, int coord1, int coord2, double& cross1);    // find 1:1 periodic orbit at given energy in coord1-coord2 plane, return starting position in 1st coordinate

/**** Utility function to find approximate L^2 of a circular orbit with given energy ****/
double findLcirc2(const CPotential* potential, double E);

}  // namespace