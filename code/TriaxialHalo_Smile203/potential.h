// This file holds expressions for potential, forces and density for various models.
#pragma once
#include "common.h"
#include <gsl/gsl_spline.h>
#include <map>

namespace smile{

/** Base class defining density model (without corresponding potential!). May be used for initializing density for general-purpose potential approximations (BSE and Spline) **/
class CDensity
{
public:
enum POTENTIALTYPE       // list of all existing types of density or density/potential models, each of them implemented in its own class
{
    PT_UNKNOWN,
    PT_NB,
    PT_BSE,
    PT_SPLINE,
    PT_LOG,
    PT_HARMONIC,
    PT_SCALEFREE,
    PT_SCALEFREESH,
    PT_DEHNEN,
    PT_PLUMMER,
    PT_ISOCHRONE,
    PT_PERFECTELLIPSOID,
    PT_NFW
};
enum SYMMETRYTYPE{     // Type of symmetry (used mainly in spherical-harmonic expansion potentials to disregard some coefficients based on the symmetry properties of the underlying potential)
    ST_NONE = 0,       // all coefficients should be considered
    ST_REFLECTION = 1, // reflection about origin (change of sign of all coordinates simultaneously) - terms with odd l are zero
    ST_PLANESYM = 2,   // reflection about principal planes (change of sign of any coordinate) - no terms with sin(m phi) or odd m
    ST_ZROTSYM = 4,    // axisymmetric model (all terms with m>0 are zero)
    ST_SPHSYM = 8,     // spherical symmetry (all terms with l>0 are zero)
    ST_TRIAXIAL = ST_REFLECTION | ST_PLANESYM,
    ST_AXISYMMETRIC = ST_TRIAXIAL | ST_ZROTSYM,
    ST_SPHERICAL = ST_AXISYMMETRIC | ST_SPHSYM
};
    CDensity() {};
    virtual ~CDensity() {};
    virtual CDensity* clone() const=0;              // derivative classes would return pointers to cloned instances of themselves
    virtual POTENTIALTYPE PotentialType() const=0;  // enumerable potential type
    virtual std::string PotentialName() const=0;    // string representation of potential type
    virtual SYMMETRYTYPE symmetry() const=0;        // returns symmetry type of this potential
    virtual double Rho(double X, double Y, double Z, double t=0) const=0;    // returns density at given coordinates, this should obviously be overriden in derivative classes
    virtual double Mass(const double r) const;      // returns mass inside given radius (approximately! not necessary to integrate density over sphere, just a rough estimate used e.g. in choosing radial nodes of Schwarzschild grid). Excluding black hole!
    double totalMass() const;                       // returns estimated M(r=infinity) or -1 if mass is infinite
    double getRadiusByMass(const double m) const;   // solves for Mass(r)=m
    bool checkMassMonotonic() const;                // safety measure: check (roughly) that mass is increasing with radius
    bool checkDensityNonzero() const;               // another safety measure: check that density doesn't drop to zero along any of three axes (important to assess spherical-harmonic approximation quality)
    virtual double getGamma() const;                // returns inner density slope estimate
private:
    CDensity& operator=(const CDensity& );          // assignment forbidden, use clone()
};

class CDensityPlummer: public CDensity
{
public:
    CDensityPlummer(double _q, double _p): CDensity(), q(_q), p(_p) {};
    virtual CDensity* clone() const { return new CDensityPlummer(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_PLUMMER; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Plummer"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double getGamma() const { return 0; };
private:
    const double q, p;          // axis ratio (y/x and z/x) of equidensity surfaces
};

class CDensityIsochrone: public CDensity
{
public:
    CDensityIsochrone(double _q, double _p): CDensity(), q(_q), p(_p) {};
    virtual CDensity* clone() const { return new CDensityIsochrone(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_ISOCHRONE; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Isochrone"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double getGamma() const { return 0; };
private:
    const double q, p;          // axis ratio (y/x and z/x) of equidensity surfaces
};

class CDensityPerfectEllipsoid: public CDensity
{
public:
    CDensityPerfectEllipsoid(double _q, double _p): CDensity(), q(_q), p(_p) {};
    virtual CDensity* clone() const { return new CDensityPerfectEllipsoid(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_PERFECTELLIPSOID; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Perfect Ellipsoid"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double getGamma() const { return 0; };
private:
    const double q, p;          // axis ratio (y/x and z/x) of equidensity surfaces
};

class CDensityNFW: public CDensity
{
public:
    CDensityNFW(double _q, double _p, double c);  // c == concentration parameter
    virtual CDensity* clone() const { return new CDensityNFW(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_NFW; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "NFW"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double getGamma() const { return 1.0; };
private:
    const double q, p; // axis ratio (y/x and z/x) of equidensity surfaces
    double rcutoff;    // cutoff radius
    double norm;       // normalization factor for density
};

/** base class for potential-density model including potential derivatives (=forces)  **/
class CPotential: public CDensity
{
public:
    const unsigned int N_dim; // dimensions (2 or 3)
    double Mbh;         // black hole mass
    ///!!!nonpublic
    const double Omega;       // rotation speed (around z axis)
    
    //CPotential(): CDensity(), N_dim(3), Mbh(0), Omega(0) {};   // overriden; 
    CPotential(unsigned int _N_dim, double _Mbh): CDensity(), N_dim(std::min<unsigned int>(3,std::max<unsigned int>(2,_N_dim))), Mbh(_Mbh), Omega(0) {};
    virtual CPotential* clone() const=0;    // override definition from CDensity because of different return type
    virtual double Phi(double X, double Y, double Z, double t=0) const=0;     // derived classes: return potential at given point
    virtual void DiffEq(unsigned n, double t, double *v, double *f) const=0;  // derived classes: evaluate forces (and possibly force derivatives) at given point
    virtual std::string info() const;       // provides one-line textual representation of potential parameters 
    double longaxisradius(double E) const;  // find maximal elongation of long-axis orbit of given energy
    double longaxisperiod(double E, double Xmax=-1) const;  // find period of oscillations for long-axis orbit of given energy
    double findintersection(double E, double X, double Y, double Z) const;       // finds k so that Phi(k*X, k*Y, k*Z) = E
    virtual bool finite() const=0;   // derived classes: report whether potential -> 0 at infinity or not
    double totalEnergy(const CPosVelPoint<double> &point) const { return Phi(point.Pos[0], point.Pos[1], point.Pos[2]) + (pow_2(point.Vel[0])+pow_2(point.Vel[1])+pow_2(point.Vel[2]))/2; };   // convenience function to return total energy for given position/velocity point
    ///nonpublic
private:
    double corotationradius() const;
};

/** Logarithmic potential, with possible core **/
class CPotentialLog: public CPotential
{
public:
    CPotentialLog(unsigned int _N_dim, double _q, double _p, double _Mbh, double _Rc): CPotential(_N_dim, _Mbh), q(_q), p(_p), Rc(_Rc) {};
    virtual CPotential* clone() const { return new CPotentialLog(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_LOG; };
    virtual std::string PotentialName() const  { return myName(); };
    static std::string myName() { return "Logarithmic"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double Phi(double X, double Y, double Z, double t=0) const;
    virtual void DiffEq(unsigned n, double t, double *v, double *f) const;
    virtual bool finite() const { return false; };
private:
    const double q, p;    // axis ratio (y/x and z/x) of equipotential surfaces
    const double Rc;      // core radius of log potential
};

/** Potential of a constant-density (harmonic) core **/
class CPotentialHarmonic: public CPotential 
{
public:
    CPotentialHarmonic(unsigned int _N_dim, double _q, double _p, double _Mbh): CPotential(_N_dim, _Mbh), q(_q), p(_p) {};
    virtual CPotential* clone() const { return new CPotentialHarmonic(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_HARMONIC; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Harmonic"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double Phi(double X, double Y, double Z, double t=0) const;
    virtual void DiffEq(unsigned n, double t, double *v, double *f) const;
    virtual bool finite() const { return false; };
private:
    const double q, p;          // axis ratio (y/x and z/x) of equidensity surfaces
};

/** Dehnen(1993) double power-law model **/
class CPotentialDehnen: public CPotential
{
public:
    CPotentialDehnen(unsigned int _N_dim, double _q, double _p, double _Mbh, double _Gamma);
    virtual CPotential* clone() const { return new CPotentialDehnen(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_DEHNEN; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Dehnen"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double Phi(double X, double Y, double Z, double t=0) const;
    virtual void DiffEq(unsigned n, double t, double *v, double *f) const;
    virtual double getGamma() const {return Gamma; }
    virtual bool finite() const { return true; };
private:
    const double q, p;          // axis ratio (y/x and z/x) of equidensity surfaces
    const double Gamma;      // cusp exponent for Dehnen potential
};

/** parent class for all potential expansions based on spherical harmonics for angular variables **/
class CPotentialSH: public CPotential
{
public:
    CPotentialSH(unsigned int _N_dim, double _Mbh, unsigned int _Ncoefs_angular): 
        CPotential(_N_dim, _Mbh), Ncoefs_angular(std::min<unsigned int>(MAX_NCOEFS_ANGULAR, _Ncoefs_angular)), mysymmetry(ST_TRIAXIAL) {};
    virtual SYMMETRYTYPE symmetry() const { return mysymmetry; };
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double Phi(double X, double Y, double Z, double t=0) const;
    virtual void DiffEq(unsigned n, double t, double *v, double *f) const;  // common function for all derivative classes
    virtual bool finite() const { return true; };
    unsigned int getNcoefs_angular() const { return static_cast<unsigned int>(Ncoefs_angular); }
protected:
    unsigned int Ncoefs_angular;
    SYMMETRYTYPE mysymmetry;    // may have different type of symmetry
    virtual void computeSHCoefs(const double r, double coefsF[], double coefsdFdr[], double coefsd2Fdr2[]) const = 0;   
    // overriden in derived classes, computes spherical-harmonic coefficients for potential and its radial (first/second) derivative at given radius, 
    // used in evaluation of potential and forces; unnecessary coefficients are indicated by passing NULL for coefs** and should not be computed
};

/** basis-set expansion on Zhao96 basis set (alpha models) **/
class CPotentialBSE: public CPotentialSH
{
public:
    CPotentialBSE(unsigned int _N_dim, double _Mbh, double _Alpha, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
        const CPointMassSetDouble &points=CPointMassSetDouble(), SYMMETRYTYPE _sym=ST_TRIAXIAL);  // init coefficients from N point masses
    CPotentialBSE(unsigned int _N_dim, double _Mbh, double _Alpha, 
        const std::vector< vectord > &coefs);  // init from stored coefficients
    CPotentialBSE(unsigned int _N_dim, double _Mbh, double _Alpha, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
        const CDensity* density);    // init potential from analytic mass model
    virtual CPotential* clone() const { return new CPotentialBSE(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_BSE; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "BSE"; };
    virtual double Mass(const double r) const;  // faster estimate of M(r)
    // get functions
    void getCoefs(std::vector< vectord > *coefsArray) const { if(coefsArray!=NULL) *coefsArray=SHcoefs; };  // return BSE coefficients array
    double getAlpha() const { return Alpha; };
    unsigned int getNcoefs_radial() const { return Ncoefs_radial; }
private:
    unsigned int Ncoefs_radial;
    std::vector<vectord> SHcoefs;
    double Alpha;      // model parameter controlling inner and outer slopes
    virtual void computeSHCoefs(const double r, double coefsF[], double coefsdFdr[], double coefsd2Fdr2[]) const;
    void prepareCoefsDiscrete(const CPointMassSetDouble &points);
    void prepareCoefsAnalytic(const CDensity* density);
    void initDefault();    // called as a default initialization when everything else fails
    void checkSymmetry();  // assigns symmetry class if some coefficients are (near-)zero
};

/** spherical-harmonic expansion as in Aguilar&Merritt 1990 (or possibly earlier refs) **/
class CPotentialSpline: public CPotentialSH
{
public:
    // init potential from N point masses with given assumed symmetry type
    CPotentialSpline(unsigned int _N_dim, double _Mbh, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
        const CPointMassSetDouble &points=CPointMassSetDouble(), SYMMETRYTYPE _sym=ST_TRIAXIAL);
    // init potential from stored SHE coefficients at given radii
    CPotentialSpline(unsigned int _N_dim, double _Mbh, 
        const vectord &_SHradii, const std::vector< vectord > &_coefs);
    // init potential from analytic mass model, may also supply desired grid radii (if not given, assign automatically)
    CPotentialSpline(unsigned int _N_dim, double _Mbh, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
        const CDensity* density, const vectord *radii=NULL);
    // not really a constructor, but a way to compute expansion coefficients 
    // for a set of discrete point masses at given radii, and output coefs in coefsArray. 
    // Splines are not initialized and the potential is not usable after this call!
    CPotentialSpline(unsigned int _Ncoefs_angular, const CPointMassSetDouble &points, const vectord *radii, std::vector< vectord > *coefsArray);
    ~CPotentialSpline();
    virtual CPotential* clone() const;
    virtual POTENTIALTYPE PotentialType() const { return PT_SPLINE; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Spline"; };
    virtual double Mass(const double r) const;  // faster estimate of M(r)
    // get functions
    unsigned int getNcoefs_radial() const { return Ncoefs_radial; }
    void getCoefs(vectord *radii, std::vector< vectord > *coefsArray, bool useNodes=true) const;  // return SHE coefficients (calculated by some internal rule) at given radii (internal spline nodes if useNodes=true, otherwise at supplied radii)
private:
    unsigned int Ncoefs_radial;
    vectord SHradii;
    double minr, maxr;                   // definition range of splines; extrapolation beyond this radius 
    double ascale;                       // scaling radius for l>0 coefs: they are functions of xi=ln(a+r)
    double gammain,  coefin;             // slope and coef. for extrapolating potential inside minr (spherically-symmetric part, l=0)
    double gammaout, coefout;            // slope and coef. for extrapolating potential outside maxr (spherically-symmetric part, l=0)
    double potcenter, potmax, potminr;   // (abs.value) potential in the center (for transformation of l=0 spline), at the outermost spline node, and at 1st spline node
    std::vector<gsl_spline*> splines;    // spline coefficients at each harmonic
    vectord slopein, slopeout;           // slope of coefs for l>0 for extrapolating inside rmin/outside rmax

    virtual void computeSHCoefs(const double r, double coefsF[], double coefsdFdr[], double coefsd2Fdr2[]) const;
    void checkSymmetry(const std::vector< vectord > &coefsArray);       // assigns symmetry class if some coefficients are (near-)zero, called on an intermediate step after computing coefs but before initializing splines
    void initDefault();    // called as a default initialization when everything else fails
    void initSpline(const vectord &radii, const std::vector< vectord > &coefsArray);  // initialize spline coefficients; radii should have "Ncoefs_radial" elements and coefsArray - "Ncoefs_radial * (Ncoefs_angular+1)^2" elements
    void prepareCoefsDiscrete(const CPointMassSetDouble &points, const vectord *radii=NULL, std::vector< vectord > *coefsArray=NULL);
    void prepareCoefsAnalytic(const CDensity* density, const vectord *radii);
    void coef0(double r, double *val, double *der, double *der2) const;     // evaluate value and optionally up to two derivatives of l=0 coefficient, taking into account extrapolation and argument scaling
    void coeflm(size_t lm, double r, double xi, double *val, double *der, double *der2, 
        double c0val, double c0der=0, double c0der2=0) const;  // evaluate value, and optionally first and second derivative of l>0 coefficients (lm is the combined index of angular harmonic >0); corresponding values for 0th coef must be known
};

/** power-law potential **/
class CPotentialScaleFree: public CPotential
{
public:
    CPotentialScaleFree(unsigned int _N_dim, double _q, double _p, double _Mbh, double _Gamma): CPotential(_N_dim, _Mbh), q(_q), p(_p), Gamma(_Gamma) {};
    virtual CPotential* clone() const { return new CPotentialScaleFree(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_SCALEFREE; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Scale-free"; };
    virtual SYMMETRYTYPE symmetry() const { return (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL); };   // returns symmetry type of this potential based on p, q being unity or not
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual double Phi(double X, double Y, double Z, double t=0) const;
    virtual void DiffEq(unsigned n, double t, double *v, double *f) const;
    virtual double getGamma() const {return Gamma; }
    virtual bool finite() const { return false; };
private:
    const double q, p;          // axis ratio (y/x and z/x) of equidensity surfaces
    const double Gamma;      // scale-free profile exponent
};

/** angular expansion of scale-free potential in spherical harmonics **/
class CPotentialScaleFreeSH: public CPotentialSH
{
public:
    CPotentialScaleFreeSH(unsigned int _N_dim, double _q, double _p, double _Mbh, double _Gamma, unsigned int _Ncoefs_angular);
    CPotentialScaleFreeSH(unsigned int _N_dim, double _Mbh, double _Gamma, const vectord &_coefs);    // init potential from stored SHE coefficients
    virtual CPotential* clone() const { return new CPotentialScaleFreeSH(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_SCALEFREESH; };
    virtual std::string PotentialName() const { return myName(); };
    static std::string myName() { return "Scale-free SH"; };
    virtual bool finite() const { return false; };
    // get functions
    void getCoefs(vectord *coefsArray) const { if(coefsArray!=NULL) *coefsArray=SHcoefs; };  // return spherical-harmonic coefficients array
private:
    const double Gamma;      // scale-free profile exponent
    vectord SHcoefs;
    virtual void computeSHCoefs(const double r, double coefsF[], double coefsdFdr[], double coefsd2Fdr2[]) const;
    void prepareCoefs(double q, double p);
};

/** frozen N-body potential calculated by tree-code algorithm (based on hackcode1.c) **/
class CPotentialNB: public CPotential   // frozen N-body potential
{
private:

enum NODETYPE{
    NODE,
    BODY,
    CELL
};

// NODE: data common to BODY and CELL structures.
typedef double vec3[N_DIM], matrix3[N_DIM][N_DIM];
struct node {
    NODETYPE type;    /* code for node type */
    double mass;      /* total mass of node */
    vec3 pos;         /* position of node */
    double eps;       /* softening length associated with this node */
};
typedef node *nodeptr;

// BODY: data structure used to represent particles.
typedef node body, *bodyptr;

// CELL: structure used to represent internal nodes of tree. First few fields coincide with the definition of body structure
#define IMAX (1 << (8 * sizeof(int) - 2))       /* max integer coordinate */
#define NSUB (1 << N_DIM)  /* subcells per cell */
struct cell{
    NODETYPE type;
    double mass;           /* total mass of cell */
    vec3 pos;              /* cm. position of cell */
    double eps;            /* softening length associated with this node */
//#ifdef TREECODE_DEHNEN_MAC
    double rmax;           /* max distance between cell's c.o.m. and any of child nodes */
//#endif
    unsigned int numbody;           /* total number of bodies inside this cell */
    nodeptr subp[NSUB];    /* descendents of cell */
#ifdef TREECODE_QUADRUPOLE
    matrix3 quad;          /* quad. moment of cell - sum of m*dx_i*dx_j */
#endif
    cell() { type=CELL; for(int s=0; s<NSUB; s++) subp[s]=NULL; mass=-1; eps=rmax=0; numbody=0; }
};
typedef cell *cellptr;

enum WALKOP{   /* operation to perform during tree-walk: quantity to calculate */
    CALCPHI,   /* potential */
    CALCACC    /* acceleration(3 dim) */
};  
public:
    
    CPotentialNB(unsigned int _N_dim, double _Mbh);  // default (empty) initialization
    template<typename NumT> CPotentialNB(unsigned int _N_dim, double _Mbh, double _eps, double _tol, 
        const std::vector< std::pair< CPosVelPoint<NumT>, NumT> > &points);  // NumT=float or double
    virtual CPotential* clone() const { return new CPotentialNB(*this); };
    virtual POTENTIALTYPE PotentialType() const { return PT_NB; };
    virtual std::string PotentialName() const { return myName(); };
    virtual SYMMETRYTYPE symmetry() const { return ST_NONE; };  // nothing assumed about the position of points, right?
    static std::string myName() { return "Nbody"; };
    virtual double Phi(double X, double Y, double Z, double t=0) const;
    virtual double Rho(double X, double Y, double Z, double t=0) const;
    virtual void DiffEq(unsigned n, double t, double *v, double *f) const;
    virtual double Mass(const double r) const;   // override method by counting bodies inside given radius
    virtual bool finite() const { return true; };
    // get functions
    size_t bodyNum() const { return bodytab.size(); };
    template<typename NumT> CPosPoint<NumT> bodyPos(size_t index) const
    { return CPosPoint<NumT>(static_cast<NumT>(bodytab.at(index).pos[0]), static_cast<NumT>(bodytab.at(index).pos[1]), static_cast<NumT>(bodytab.at(index).pos[2])); };
    double bodyMass(size_t index) const { return bodytab.at(index).mass; };
private:
    const double eps;         /* potential softening parameter; negative value means position-dependent softening proportional to local density^{-1/3} */
    const double tol;         /* accuracy parameter: 0.0 => exact N^2 force calculation (don't ever try!) */
    std::vector<body> bodytab;
    std::vector<cell> celltab;
    
    nodeptr troot;      /* tree root */
    vec3 rmin;          /* lower-left corner of coord. box */
    double rsize;       /* side-length of int. coord. box */
    // tree construction
    void maketree();                   // initialize tree for the constructor
    bool expandbox(const size_t b);             // expand bounding box (aka root cell) if particle does not fit into existing one
    bool appendtree(const size_t b);              // append particle to the tree 
    bool intcoord(int xp[N_DIM], const vec3 rp) const; // calculate integerized coords
    int subindex(const int x[N_DIM], const int l) const;     // determine subcell that a given position belongs to (l is cell level)
    void centerofmass(const nodeptr q, const vec3 corner, const double cellsize, int l);      // recursive calculation of cell center-of-mass
    void assigneps(const nodeptr q, const double epsparent, const double cellsize);   // recursive initialization of individual smoothing radii of cells (epsparent is smoothing radius of parent cell), from the topmost cell down to every leaf
    void propagateeps(const nodeptr q);      // recursively assign epsilon to tree cells from underlying nodes
    // tree traversal
    void walktree(const nodeptr p, const double dsq, const vec3 pos0, const WALKOP operation, double data[]) const;
    void walktreedens(const nodeptr p, const double dsq, const double searchradsq, const vec3 pos0, const int xp[N_DIM], const int lev, vectorpd* data) const;
};

/** Correspondence between potential/density names and corresponding classes **/
typedef std::map<CDensity::POTENTIALTYPE, std::string> PotentialNameMapType;   // lists all 'true' potentials, i.e. those providing a complete density-potential(-force) pair
typedef std::map<CDensity::POTENTIALTYPE, std::string> DensityNameMapType;     // lists all analytic density profiles (including those that don't have corresponding potential, but excluding general-purpose expansions)
typedef std::map<CDensity::SYMMETRYTYPE, std::string > SymmetryNameMapType;    // lists available symmetry types
extern PotentialNameMapType PotentialNames;
extern DensityNameMapType DensityNames;
extern SymmetryNameMapType SymmetryNames;
CDensity::POTENTIALTYPE getPotentialTypeByName(std::string PotentialName);  // perform reverse lookup in PotentialNames map
CDensity::POTENTIALTYPE getDensityTypeByName(std::string DensityName);      // same for DensityNames
CDensity::SYMMETRYTYPE getSymmetryTypeByName(std::string SymmetryName);     // same for symmetry

/** structure that contains parameters for all possible potentials **/
struct CConfigPotential
{
    unsigned int N_dim;                    // dimensions of space (2 or 3)
    double q, p, Mbh, Rc, Gamma, Alpha;    // these values will be used to initialize proper potential type
    double Omega;                          // rotation frequency of potential figure (unused)
    unsigned int Ncoefs_radial, Ncoefs_angular;  // number of radial and angular coefficients in spherical-harmonic expansion
    double nbpotEps;                       // treecode smooothing length (negative means adaptive smoothing based on local density, absolute value is the proportionality coefficient between eps and mean interparticle distance)
    double nbpotTol;                       // tree cell opening angle
    CDensity::POTENTIALTYPE PotentialType; // currently selected potential type
    CDensity::POTENTIALTYPE DensityType;   // if pot.type == BSE or Spline, this gives the underlying density profile approximated by these expansions or flags that an Nbody file should be used
    CDensity::SYMMETRYTYPE SymmetryType;   // if using Nbody file with the above two potential expansions, may assume certain symmetry on the coefficients (do not compute them but just assign to zero)
};
extern CConfigPotential configPotential;   // default parameters

}