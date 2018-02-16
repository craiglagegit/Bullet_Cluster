/** \file    common.h
    \brief   provides forward class declarations, several abstract data classes, a couple of global functions, and global constants
    \author  EV
    \date    2009-2012

*/
#pragma once
#include <utility>
#include <vector>
#include <string>

namespace smile{

/*------- define several convenience types -------*/
typedef std::vector<double> vectord;
typedef std::vector<float> vectorf;
typedef std::pair<double, double> paird;
typedef std::vector< paird > vectorpd;
typedef std::pair<float, float> pairf;
typedef std::vector< pairf > vectorpf;
/*------- Declare several global convenience functions that are defined in common.cpp -------*/

/// a shorthand for squaring a number
inline double pow_2(double x) { return x*x; }
inline float pow_2(float x) { return x*x; }
inline unsigned int pow_2(unsigned int x) { return x*x; }

/// generates a grid with exponentially growing spacing, so that 1st element is at xmin and last at xmax; 
/// if zeroelem=true, 0th node is at zero (otherwise at xmin)
void createNonuniformGrid(vectord &grid, int nnodes, double xmin, double xmax, bool zeroelem=false);

//*------- typedef and pointer to message and error reporting routines -------*//
//* (implementation is program-specific, may dump text to stderr or do something else *//
typedef void show_message_type(const std::string &message);

/// global variable to the routine that displays errors (if it is NULL then nothing is done)
extern show_message_type* my_error_ptr;

/// global variable to the routine that shows information messages 
/// (emitted by various time-consuming routines to display progress), if NULL then no messaging is done
extern show_message_type* my_message_ptr;

/// the interface routine for error reporting (redirects the call to my_error_ptr if it is defined)
inline void my_error(const std::string &message)
{ if(my_error_ptr!=NULL) my_error_ptr(message); }

/// the interface routine for message reporting (redirects the call to my_message_ptr if it is defined)
inline void my_message(const std::string &message)
{ if(my_message_ptr!=NULL) my_message_ptr(message); }

/*------- Global constants not changeable in runtime -------*/
/// to dump out some internal debug information
//#define DEBUGPRINT

/// to use standard math constants such as M_PI defined in cmath
#define _USE_MATH_DEFINES

/// number of significant digits for initial conditions
#define NUM_DIGITS_ROUND 6

/// SMILE version identifier
#define VERSION "2.0"

/*----------- Parameters for potential and force evaluation -------------*/

// Dimensions of the real world. DO NOT CHANGE (unless you are God)!
const unsigned int N_DIM=3;

// Black hole smoothing radius (for calculation of potential and forces)
const double BH_SMOOTH=1e-8;

// max number of basis-function expansion members (radial and angular). 
const unsigned int MAX_NCOEFS_RADIAL=100;
const unsigned int MAX_NCOEFS_ANGULAR=20;

// relative accuracy of density computation
const double EPSREL_DENSITY_INT=1e-4;

// relative accuracy of potential computation (integration tolerance parameter)
const double EPSREL_POTENTIAL_INT=1e-6;

// inner and outer cutoff radii appearing in calculations of M(r), r(E) and related
const double MIN_RADIUS_CUTOFF=1e-6;
const double MAX_RADIUS_CUTOFF=1e+6;

/*------------- Orbit integration parameters ----------*/

/// maximal number of integration steps (in case something went wrong with the timestep, do not dead-lock)
const unsigned int NUMSTEP_MAX=((unsigned int)1e7);

/// with black hole, if radius is less than Mbh*{const below}, increase accuracy parameter in Runge-Kutta integrator (to make timestep smaller)
const double BH_RADIUS_INCREASE_ACCURACY=1e-2;

/// for Nbody treecode, another safety measure is to decrease timestep near BH by this factor (in leap-frog integrator)
const double TREECODE_TIMESTEP_BH_FACTOR=0.2;

/// whether to use quadrupole moments in tree expansion of Nbody potential (no reason to switch it off)
#define TREECODE_QUADRUPOLE

/// if defined, use compact Epanechnikov kernel for potential/force smoothing, otherwise use Plummer softening (not recommended?)
#define TREECODE_EPS_COMPACT_KERNEL

/// if defined, use Dehnen's (2000) multipole acceptance criterion based on r_max - max distance between cell's c.o.m. and its child nodes, otherwise use standard Barnes&Hut criterion based on cell size
//#define TREECODE_DEHNEN_MAC

/// density estimate in treecode by using distance to k-th neighbour. If using spline, put more (>~20), since this kernel is more concentrated
const unsigned int TREECODE_DENSITY_KTH_NEIGHBOUR=32;

/// method to smooth density: standard SPH spline, or Epanechnikov kernel (both have finite support equal to distance to Kth neighbour)
//#define TREECODE_DENSITY_SPLINE

/// if LYAPUNOVVAREQ is defined, orbit integration uses variational equation for deviation vector.
/// if not, a nearby trajectory is integrated. Note: variational equation is currently not implemented for scale-free potential
//#define LYAPUNOVVAREQ

/// initial value of orbit separation (if variational equation is not used and two nearby trajectories are followed)
const double LYAP_DEV_INIT=1e-10;

/// if distance between nearby trajectories reaches the following threshold, renormalize it back to LYAP_DEV_INIT
const double LYAP_DEV_MAX=1e-6;

/// if variational equation is used, perform renormalization of deviation vector if greater than threshold (only for numerical convenience, since the var.equation is anyway linear_
const double LYAP_VAREQ_DEV_MAX=1e5;

/*------------- Orbit spectral analysis and classification options ----------*/
/// orbit spectrum is calculated by fast fourier transform from GSL.
/// to avoid worst-case behavior of FFT, we ensure that its length never contains primes larger that MAX_FFT_PRIME (here equal to 42 minus unity)
const unsigned int MAX_FFT_PRIME=41;

/// Use Hunter's method for precise determination of frequency of a given line if the spectrum neat its maximum is reasonably similar to that of a single line. Otherwise use less precise method without Hanning filtering. Controls threshold for selection.
const double FREQ_PRECISE_THRESHOLD=5.0;

/// maximal number of spectral lines in each coordinate to search for
const unsigned int MAX_SPECTRAL_LINES=15;

/// do not use lines which amplitude is X times lower than amplitude of the leading line 
const double MIN_SPECTRAL_LINE_AMPLITUDE=1e-2;

/// method to estimate freq.diffusion rate: 0 - difference in the largest of three frequencies (Valluri&Merritt 98), 
/// 1 - average of differences in all three freq, 2 - average over those coords for which |f1-f2|/(f1+f2)<0.5 (more robust to spurious errors in lf determination)
#define FREQ_DIFF_METHOD 2

/// accuracy of frequency comparison (to determine resonances), in units of Nyquist frequency
const double FREQ_ACCURACY=0.1;

/// minimum value of freq.diff threshold separating regular and chaotic orbits (actual threshold depends on integration interval)
const double FREQ_DIFF_REG_MIN=1e-6;

/// maximum order of m:n frequency resonance (or linear dependence of three frequencies)
const int FREQ_RES_2=20;

/// maximum order of l*a+m*b+n*c frequency commensurability (thin orbit)
const int FREQ_RES_3=10;

/// threshold value of |<L>| / <|L|> to classify orbit as tube
const double OC_TUBE_AVGL=0.9;

/// threshold value of |<z>| / rmax to classify orbit as pyramid
const double OC_PYRAMID_AVGZ=0.1;

/// max number of points in Poincare section
const unsigned int MAX_POINCARE_POINTS=10000;

/// fraction of closest pericenter passages which are used in fitting
const double PERICENTER_FIT_FRACTION=0.1;

/// minimum number of points used in pericenter fitting
const int PERICENTER_FIT_MIN_POINTS=10;

/// defines max points in pericenter fitting
const int PERICENTER_FIT_MAX_POINTS=50;

/*-------------  Schwarzschild modelling options ---------------*/

/// if using Schwarzschild stationary start space (3d), put points on a square grid in 3 sectors adjacent to X, Y and Z axes. (same as in Schwarzschild 1993)
/// otherwise, put points along lines of constant Z (=constant theta), number of points being proportional to sin(theta). (same as in Terzic 2002)
#define SCHWARZSCHILD_START_SPACE

/// if checked, use principal-plane start-space confined to annulus with inner radius of closed loop orbit; otherwise put points in the whole sector
#define ANNULUS_PP_START_SPACE

/// calculation of cell mass is done by refining cell into AxAxR smaller cells and computing density in each node. Below specify number of subdivisions
///!!! replace naive uniform subdivision by something fancier! (Note: integration in radius is done by computing power-law slope and integrating it analytically)
#define NUM_SUBDIV_ANGULAR 5
#define NUM_SUBDIV_RADIAL 10

/// max. relative difference in cell mass to consider the cell as feasible: abs(1-Mcell/Mcellrequired)
const double SCHW_MAX_CELLMASS_REL_DIFFERENCE=1e-2;

/// min. orbit weight to use the orbit after optimization problem solved or iterations finished (divide by number of orbits!)
const double SCHW_MIN_ORBIT_WEIGHT_USED=1e-4;

/// weight of quadratic terms in objective function
const double SCHW_CONSTRAINT_QUAD_WEIGHT=0.02;

/// weight of anisotropy constraints in objective function, (weight of orbit-cell-mass constraint is 1.0)
const double SCHW_CONSTRAINT_BETA_WEIGHT=0.1;

/*/// min. and max. number of Lucy iterations
#define SCHW_MIN_ITERATIONS 10
#define SCHW_MAX_ITERATIONS 50000
/// number of Lucy iterations during which the residual is allowed to grow (if more, then stop and use the best values for weights)
#define SCHW_MAX_GROWING_ITERATIONS 20
/// number of iterations to update weights
#define SCHW_ITERATIONS_UPDATE 100
*/
/// Nbody export: relative importance of nullifying total angular momentum at the expense of pericenter velocity * pericenter offset
const double NBEXPORT_XV_FACTOR=10.0;

/*-------------- Options for GUI and high-level operations ------------*/

/// number of points in curve representing equipotential surface (actually *4)
const unsigned int NUM_POINTS_PLOT_EQUIPOTENTIAL=100;

/// number of points in outer bounding curve in Poincare section (actually *4)
const unsigned int NUM_POINTS_PLOT_POINCARE=100;

/// min. weight of orbit to be displayed in GUI
const double SCHW_MIN_ORBIT_WEIGHT_DISPLAYED=1e-4;

/// width of frequency spectrum displayed (in units of orbital frequency)
const double PLOT_FREQ_MAX_F=5.0;

/// whether to use 3d plotting library Qwt3D in GUI
#define USEQWT3D

#ifdef USEQWT3D
/// whether to perform triangulation and display 'envelope' of an orbit (useless without Qwt3D; requires external program QDelaunay)
#define USE3DFACETS
#endif

/// auto-load settings from file
#define DEFAULT_SETTINGS_FILE_NAME "smile.ini"

/*------- Declare several classes implemented in common.cpp and various other modules --------*/

class CDensity;
class CPotential;
class COrbit;
class COrbitDesc;
class COrbitLibrary;
class CBasicOptimizationSolver;
class CBasicOrbitFilteringFnc;
class CBasicSchwModel;
class CBasicOrbitRuntimeFnc;
class CBasicOrbitRuntimeFncCreator;
class CBasicInformation;
template<typename NumT> struct CPosPoint;
template<typename NumT> struct CPosVelPoint;
template<typename NumT> struct COrbitInitData;
typedef float smnumtype;  // how to store data for schw.modelling (float or double)
typedef std::vector<smnumtype> vectorn;

/*----------- definition of some basic (abstract) classes to be implemented in various modules -------------*/

/** solves the optimization problem defined as  
   M sol = rhs, where M is N_v * N_c matrix, sol is N_v vector of variables to be found, rhs is N_c vector of constraints to be satisfied;
   tries to achieve exact solution but allows for deviation from it, which is penalized in the cost function;
   denote deviation dev[c] = | \Sum_o M[v][c] sol[v] - rhs[c] |
   with cost function = (quad.coef.)*\Sum_v sol[v]^2 + \Sum_v solWeight[v]*sol[v] + \Sum_c rhsWeight[c]*dev[c]
   if solved by linear programming, the first terms is omitted
   Return value is 0 for success, <0 for error code
*/
class CBasicOptimizationSolver
{
public:
    virtual ~CBasicOptimizationSolver() {};
    virtual int callSolver(
        const std::vector< vectorn > &linearMatrix, // matrix M of linear equations for the optimization problem  "M w = rhs"
        const vectorn &rhs,         // constraints to be satisfied
        const vectorn &rhsWeight,   // constraint importance
        const vectorn &solWeight,   // penalties for elements in sol vector
        const double maxSolValue,   // if nonzero, specifies upper limit on sol[v] in addition to lower limit of zero
        vectorn *sol) const=0;      // returns solution in this vector
    virtual std::string errorDescription(int result) const=0;  // reports text information for given error code
};

/** Parent class for orbit filtering functions which evaluate whether an orbit should be used (in some calculation).
    Return value is a float number (for example, between 0 and 1 for telling how well an orbit satisfies some criterion, but may be >1 or negative).
    Example of the usage of these functions is in COrbitDesc::getOrbitPopulation routine which computes statistics 
    based on acceptance criteria returned by such a function, or in assigning penalty to orbits of specified sort in Schwarzschild modelling.
**/
class CBasicOrbitFilteringFnc
{
public:
    virtual ~CBasicOrbitFilteringFnc() {};
    virtual double eval(const COrbitDesc* orbit) const = 0;
};

/** defines a routine which may be called after each integration timestep, to perform user-specific data collection 
    (such as computing cell occupation times in Schw model or recording Poincare section points)
**/
class CBasicOrbitRuntimeFnc
{
public:
enum FNCTYPE {
    FT_TRAJ_ANALYSIS,
    FT_TRAJ_SAMPLE,
    FT_POINCARE,
    FT_PERICENTER,
    FT_SCHW
};
    explicit CBasicOrbitRuntimeFnc(const COrbit* _orbit): orbit(_orbit) {} ;   // constructs instance and initializes some internal data
    virtual ~CBasicOrbitRuntimeFnc() {};    // destroys all necessary internal data
    virtual FNCTYPE FncType() const=0;  // derived classes return their type
    virtual void Timestep(const double told, const double t, const double y[])=0;  // called after each timestep (t_old -> t), with y[] holding current values of x,v (and possibly deviation vectors if Lyapunov exp.is calculated). This procedure may also call orbit->getInterpolatedTrajectory() to get interpolated coordinates on the interval t_old, t
    virtual void Finish() {};   // called when integration is finished, performs cleanup/postprocessing but doesn't destroy collected data
    virtual CBasicInformation* getData() const = 0;
protected:
    const COrbit* orbit;
};

/** defines a basic class that creates some sort of timestep function for any orbit passed as a parameter.
    Derived classes may receive additional data in constructor and emit specific types of timestep functions for a given orbit.
**/
class CBasicOrbitRuntimeFncCreator
{
public:
    virtual ~CBasicOrbitRuntimeFncCreator() {};
    virtual CBasicOrbitRuntimeFnc* createRuntimeFnc(const COrbit* orbit) const=0;
};

/// Data-only class containing 3 coordinates and 3 velocities. 
/// Template parameter NumT may be float or double
/// Copy constructor provided for both instantiations so that each one may be seamlessly converted to another one by (implicit) type cast
template<typename NumT> struct CPosVelPoint
{
public:
    /// default constructor
    CPosVelPoint() { for(unsigned int i=0; i<N_DIM; i++) { Pos[i]=0; Vel[i]=0; } }
    /// construct from six float/double numbers
    CPosVelPoint( const NumT Pos0, const NumT Pos1, const NumT Pos2, const NumT Vel0, const NumT Vel1, const NumT Vel2)
    { Pos[0]=Pos0; Pos[1]=Pos1; Pos[2]=Pos2; Vel[0]=Vel0; Vel[1]=Vel1; Vel[2]=Vel2; };
    /// copy constructor from CPosVelPoint<double>
    CPosVelPoint( const CPosVelPoint<double> &src ) { for(unsigned int i=0; i<N_DIM; i++) 
    { Pos[i]=static_cast<NumT>(src.Pos[i]); Vel[i]=static_cast<NumT>(src.Vel[i]); } };
    /// copy constructor from CPosVelPoint<float>
    CPosVelPoint( const CPosVelPoint<float> &src ) { for(unsigned int i=0; i<N_DIM; i++) 
    { Pos[i]=static_cast<NumT>(src.Pos[i]); Vel[i]=static_cast<NumT>(src.Vel[i]); } };
    NumT Pos[N_DIM];   //< Components of position vector. Warning: no range check performed when setting Pos or Vel[n>=N_dim] !
    NumT Vel[N_DIM];   //< Components of velocity vector
};

/// Data-only class containing 3 coordinates. 
/// Template parameter NumT may be float or double
/// Copy constructor provided for both instantiations so that each one may be seamlessly converted to another one by (implicit) type cast
template<typename NumT> struct CPosPoint
{
public:
    /// default constructor
    CPosPoint() { for(unsigned int i=0; i<N_DIM; i++) { Pos[i]=0; } }
    /// construct from three float/double numbers
    CPosPoint( const NumT Pos0, const NumT Pos1, const NumT Pos2)
    { Pos[0]=Pos0; Pos[1]=Pos1; Pos[2]=Pos2; };
    /// copy constructor from CPosPoint<double>
    CPosPoint( const CPosPoint<double> &src ) { for(unsigned int i=0; i<N_DIM; i++) 
    { Pos[i]=static_cast<NumT>(src.Pos[i]); } };
    /// copy constructor from CPosPoint<float>
    CPosPoint( const CPosPoint<float> &src ) { for(unsigned int i=0; i<N_DIM; i++) 
    { Pos[i]=static_cast<NumT>(src.Pos[i]); } };
    /// copy constructor from CPosVelPoint<double>
    CPosPoint( const CPosVelPoint<double> &src ) { for(unsigned int i=0; i<N_DIM; i++) 
    { Pos[i]=static_cast<NumT>(src.Pos[i]); } };
    /// copy constructor from CPosVelPoint<float>
    CPosPoint( const CPosVelPoint<float> &src ) { for(unsigned int i=0; i<N_DIM; i++) 
    { Pos[i]=static_cast<NumT>(src.Pos[i]); } };
    NumT Pos[N_DIM];   //< Components of position vector. Warning: no range check performed when setting Pos or Vel[n>=N_dim] !
};

/** template class for specifying orbit initial conditions.
    Template parameter NumT may be float or double
    Initial conditions include: pointer to potential, triplet of coordinates/velocities, 
    time unit used for dealing with dimensionless orbit time (normally is equal to the period of X-axis orbit with given energy),
    time step for outputting trajectory (has nothing to do with the actual integration timestep),
    flag signalling whether to use twice as many integration variables to estimate Lyapunov exponent (variation equation or nearby orbit)
*/
template<typename NumT> struct COrbitInitData
{
    COrbitInitData() {potential=NULL; timeStep=0; timeUnit=0; calcLyapunov=false; };
    COrbitInitData(const CPotential* _potential, const NumT _timeStep, const NumT _timeUnit, const CPosVelPoint<NumT> &_initCond, const bool _calcLyapunov) : 
      potential(_potential), 
      timeStep(_timeStep), 
      timeUnit(_timeUnit),
      initCond(_initCond), 
      calcLyapunov(_calcLyapunov) 
    {};
    template<typename NumT1> COrbitInitData(const COrbitInitData<NumT1> &src) :
      potential(src.potential), 
      timeStep(static_cast<NumT>(src.timeStep)),
      timeUnit(static_cast<NumT>(src.timeUnit)),
      initCond(src.initCond),
      calcLyapunov(src.calcLyapunov) 
    {};
    const CPotential* potential;  //< pointer for potential; must remain valid throughout orbit integration and analysis
    NumT timeStep;           //< constant timestep of output trajectory
    NumT timeUnit;           //< period of x-axis orbit with given energy (unit of time used in dimensionless frequencies)
    CPosVelPoint<NumT> initCond; //< initial position and velocity
    bool calcLyapunov;       //< whether to calculate Lyapunov exponent
};
/** helper function for initializing missing timestep/timeunit in incomplete initial conditions data **/
template<typename NumT> COrbitInitData<NumT> completeInitData(const COrbitInitData<NumT>& InitData);

/// root data-container class for exchange between various objects
class CBasicInformation
{
public:
enum INFOTYPE {
    IT_UNKNOWN,
    IT_TRAJ_ANALYSIS,
    IT_TRAJ_SAMPLE,
    IT_POINCARE,
    IT_PERICENTER,
    IT_SCHW
};
    virtual ~CBasicInformation() {};
    virtual INFOTYPE infoType() const=0;
    virtual std::string toString() const = 0;
};

/// convenience definition of vectors of timestep functions and function creators
typedef std::vector<CBasicOrbitRuntimeFnc*> vectorRuntimeFncs;
typedef std::vector<const CBasicOrbitRuntimeFncCreator*> vectorRuntimeFncCreators;
typedef std::vector<const CBasicInformation*> vectorInformation;

/// convenience definitions of templated vectors containing particle sets with/without masses
typedef std::vector< CPosVelPoint<double> > CPosVelDataDouble;
typedef std::vector< CPosVelPoint<float > > CPosVelDataFloat;
typedef std::vector< std::pair< CPosVelPoint<double>, double> > CPointMassSetDouble;
typedef std::vector< std::pair< CPosVelPoint<float >, float > > CPointMassSetFloat;
typedef std::vector< std::pair< CPosVelPoint<float >, unsigned int> > CPointCountSetFloat;
}