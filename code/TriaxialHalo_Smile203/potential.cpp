#include "potential.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include "stringconv.h"

#ifdef DEBUGPRINT
#include <iostream>
#endif

namespace smile{

//----------------------------------------------------------------------------//
// base class for density-only model, without corresponding potential: CDensity

// helper functions for integrating density over volume
struct CDensityParam{
    const CDensity* P;
    double r, rcostheta, rsintheta;
};
double getDensityPhi(double phi, void* params)
{
    return ((CDensityParam*)params)->P->Rho(((CDensityParam*)params)->rsintheta*cos(phi), ((CDensityParam*)params)->rsintheta*sin(phi), ((CDensityParam*)params)->rcostheta);
}
double getDensityTheta(double theta, void* params)
{
    gsl_function F;
    F.function=&getDensityPhi;
    F.params=params;
    ((CDensityParam*)params)->rcostheta=((CDensityParam*)params)->r*cos(theta);
    ((CDensityParam*)params)->rsintheta=((CDensityParam*)params)->r*sin(theta);
    double result, error;
    size_t neval;
    gsl_integration_qng(&F, 0, M_PI/2, 0, EPSREL_DENSITY_INT, &result, &error, &neval);
    return result*sin(theta);
}
double getDensityR(double r, void* params)
{
    gsl_function F;
    F.function=&getDensityTheta;
    F.params=params;
    ((CDensityParam*)params)->r=r;
    double result, error;
    size_t neval;
    gsl_integration_qng(&F, 0, M_PI/2, 0, EPSREL_DENSITY_INT, &result, &error, &neval);
    return result*8*r*r;
}
double CDensity::Mass(const double r) const
{
    if(r<=0) return 0;
    // by default, integrate over density inside given radius; may be replaced by cheaper and more approximate evaluation for child classes
    gsl_function F;
    CDensityParam params;
    F.function=&getDensityR;
    F.params=&params;
    params.P=this;
    double result, error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    gsl_integration_qag(&F, 0, r, 0, EPSREL_DENSITY_INT, 1000, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

double CDensity::totalMass() const
{
    double rad=32;
    double mass1, mass2=Mass(rad), mass3=Mass(rad*2), massest=0, massestprev;
    int num_neg=0;
    const int max_num_neg=4;
    do{
        rad*=2;
        mass1=mass2; mass2=mass3; mass3=Mass(rad*2);
        if(mass2==mass3) return mass3;  // mass doesn't seem to grow with raduis anymore
        massestprev=massest>0 ? massest : mass3;
        massest=(mass2*mass2-mass1*mass3)/(2*mass2-mass1-mass3);
        if(gsl_isnan(massest) || massest<0) num_neg++;  // increase counter of 'bad' attempts (negative means that mass is growing at least logarithmically with radius)
    } while(rad<MAX_RADIUS_CUTOFF && mass2!=mass3 && (massest<0 || fabs(massestprev-massest)/massest>EPSREL_DENSITY_INT) && num_neg<max_num_neg);
    if(fabs((massestprev-massest)/massest)>EPSREL_DENSITY_INT) return -1;   // total mass seems to be infinite
    else return massest;
}

double findRadiusByMass(double r, void* params)
{
    return ((CDensityParam*)params)->P->Mass(r) - ((CDensityParam*)params)->r;  // params->r holds mass fraction
}

double CDensity::getRadiusByMass(const double m) const
{
    double totMass=totalMass();
    if(totMass<0 || m>=totMass) return -1;
    double Xlow=0, Xupp=32, Xroot;
    while(Mass(Xupp)<m && Xupp<MAX_RADIUS_CUTOFF) Xupp*=2;
    if(Xupp>=MAX_RADIUS_CUTOFF) return MAX_RADIUS_CUTOFF; // error, but not a fatal one...
    gsl_function F;
    CDensityParam PP;
    PP.r=m;
    PP.P=this;
    F.function=&findRadiusByMass;
    F.params=&PP;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set (s, &F, Xlow, Xupp);
    int status=0, iter=0;
    do{
        iter++;
        status = gsl_root_fsolver_iterate (s);
        Xroot= gsl_root_fsolver_root (s);
        Xlow = gsl_root_fsolver_x_lower (s);
        Xupp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (Xlow, Xupp, 0, 1e-8);
    }
    while (status == GSL_CONTINUE && iter < 100);
    gsl_root_fsolver_free (s);
    return Xroot;
}

bool CDensity::checkMassMonotonic() const
{
    double massprev=Mass(0);
    if(massprev!=massprev || massprev!=0) return false;  // mass should be zero at origin (not including Mbh)
    for(double rad=0.01; rad<std::min(MAX_RADIUS_CUTOFF,10000.0); rad*=1.25)
    {
        double mass=Mass(rad);
        if(mass<massprev*(1-1e-10)) return false;
        massprev=mass;
    }
    return true;  // no apparent evidence for weird behaviour..
}

bool CDensity::checkDensityNonzero() const
{
    for(double rad=1e-4; rad<std::min(MAX_RADIUS_CUTOFF,10000.0); rad*=1.5)
    {
        if(Rho(rad, 0, 0) <=0 || Rho(0, rad, 0) <=0 || Rho(0, 0, rad) <=0)
            return false;
    }
    return true;  // no apparent evidence for weird behaviour..
}

double CDensity::getGamma() const
{
    double mass1, mass2, mass3;
    double rad=1./1024;
    do{
        mass2=Mass(rad);
        if(mass2<=0) rad*=2;
    }while(rad<1 && mass2==0);
    if(mass2<=0) return -1; // apparent error
    mass3=Mass(rad*2); // hopefully >0!
    double alpha1, alpha2=log(mass3/mass2)/log(2.), gamma1=-1, gamma2=3-alpha2;
    do{
        rad/=2;
        mass1=Mass(rad);
        if(mass1<=0) return gamma2;
        alpha1=log(mass2/mass1)/log(2.);
        gamma2=gamma1<0 ? 3-alpha1 : gamma1;  // rough estimate
        gamma1=3- (2*alpha1-alpha2); // extrapolated estimate
        alpha2=alpha1;
        mass3=mass2; mass2=mass1;
    } while(fabs(gamma1-gamma2)>0.001);
    return gamma1;
}
//----------------------------------------------------------------------------//
// base class for potential/density model: CPotential

// helper functions for long axis radius/period evaluation
struct CPotentialParam1{
    const CPotential* P;
    double E;
};
struct CPotentialParam4{
    const CPotential* P;
    double E, X, Y, Z;
};

// general functions for any potential
double twoV2(double x, void* params)
{
    double E=((CPotentialParam1*)params)->E;
    return E-((CPotentialParam1*)params)->P->Phi(x,0,0) -0.5*pow_2(((CPotentialParam1*)params)->P->Omega*x);
}
double Torb(double x, void* params)
{
    if(x==0) return 0;
    double E=((CPotentialParam1*)params)->E;
    double V2=2*(E-((CPotentialParam1*)params)->P->Phi(x,0,0) -0.5*pow_2(((CPotentialParam1*)params)->P->Omega*x));
    if(V2<=0) return 0; else return 1/sqrt(V2);
}
double k_intersect(double k, void* params)
{
    double E= ((CPotentialParam4*)params)->E;
    double X= ((CPotentialParam4*)params)->X;
    double Y= ((CPotentialParam4*)params)->Y;
    double Z= ((CPotentialParam4*)params)->Z;
    return E - ((CPotentialParam4*)params)->P->Phi(X*k,Y*k,Z*k) -0.5*pow_2(((CPotentialParam4*)params)->P->Omega*k)*(X*X+Y*Y);
}
double minusPhi(double x, void* params)
{
    if(((CPotentialParam4*)params)->E)
        return -((CPotentialParam1*)params)->P->Phi(0,x,0);
    else
        return -((CPotentialParam1*)params)->P->Phi(x,0,0);
}

double CPotential::longaxisradius(double E) const
{
    // find maximal x-axis elongation for given energy
    double Xlow=1, Xupp=1, Xmax;
    while ((Phi(Xlow,0,0)+0.5*pow_2(Omega*Xlow)>E) && (Xlow>MIN_RADIUS_CUTOFF)) Xlow/=2;
    ///nonpublic
    if(Xlow<MIN_RADIUS_CUTOFF) 
    {
        /*if(Omega>0) // it may happen that we're outside corotation radius and the effective potential is declining with radius
        {
            Xlow=1;
            while((Phi(Xupp,0,0)>E) && (Xupp<1e4)) Xupp*=2;
            if(Xupp>1e4)
                return -1;  // still something erroneous happened
        }
        else*/ return -1;  // signal of error (E lower than minimum of the potential)
    }
    else
    {
        while ((Phi(Xupp,0,0)+0.5*pow_2(Omega*Xupp)<E) && (Xupp<MAX_RADIUS_CUTOFF)) Xupp*=2;
        if(Xupp>MAX_RADIUS_CUTOFF) 
        {
            /*if(Omega>0) // we might catch potential maximum between Xlow and Xupp
            {
                int niter=0;
                double Eupp, Elow=Phi(Xlow,0,0);
                Xupp=2*Xlow;
                while(niter<100)
                {
                    Eupp=Phi(Xupp,0,0);
                    if(Eupp>E) niter=100;
                    else if(Eupp>Elow) Xupp*=1.25; else Xupp=(Xupp+Xlow)/2;
                    niter++;
                }
                if(Eupp<E) return -1;  // still unable to enclose root
            }
            else*/ return -1;  // signal of error (e.g. E>0 and potential is everywhere negative)
        }
    }
    gsl_function F;
    CPotentialParam1 PP;
    PP.E=E;
    PP.P=this;
    F.function=&twoV2;
    F.params=&PP;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set (s, &F, Xlow, Xupp);
    int status=0, iter=0;
    do{
        iter++;
        status = gsl_root_fsolver_iterate (s);
        Xmax = gsl_root_fsolver_root (s);
        Xlow = gsl_root_fsolver_x_lower (s);
        Xupp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (Xlow, Xupp, 0, 1e-9);
    }
    while (status == GSL_CONTINUE && iter < 100);
    gsl_root_fsolver_free (s);
    return Xmax;    
}

double CPotential::longaxisperiod(double E, double Xmax) const
{
    // calculate T_orb (for x-axis orbit)
    if(Xmax<0)
        Xmax=longaxisradius(E);
    if(Xmax<0) return -1;  // error
    double result, error;
    gsl_function F;
    CPotentialParam1 PP;
    PP.E=E;
    PP.P=this;
    F.function=&Torb;
    F.params=&PP;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    gsl_integration_qags (&F, 0, Xmax, 0, 1e-7, 1000, w, &result, &error); 
    gsl_integration_workspace_free (w);
    return 4*result;
}

double CPotential::findintersection(double E, double X, double Y, double Z) const
{
    double klow=1, kupp=1, k;
    while ((Phi(klow*X,klow*Y,klow*Z)>E) && (klow>MIN_RADIUS_CUTOFF)) klow/=2;
    if(klow<MIN_RADIUS_CUTOFF) 
        return -1;  // signal of error (E lower than minimum of the potential)
    while ((Phi(kupp*X,kupp*Y,kupp*Z)<E) && (kupp<MAX_RADIUS_CUTOFF)) kupp*=2;
    if(kupp>MAX_RADIUS_CUTOFF) 
        return -1;  // signal of error (e.g. E>0 and potential is everywhere negative)
    gsl_function F;
    CPotentialParam4 PP;
    PP.E=E;
    PP.P=this;
    PP.X=X;
    PP.Y=Y;
    PP.Z=Z;
    F.function=&k_intersect;
    F.params=&PP;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set (s, &F, klow, kupp);
    int status=0, iter=0;
    do{
        iter++;
        status = gsl_root_fsolver_iterate (s);
        k = gsl_root_fsolver_root (s);
        klow = gsl_root_fsolver_x_lower (s);
        kupp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (klow, kupp, 0, 1e-9);
    }
    while (status == GSL_CONTINUE && iter < 100);
    gsl_root_fsolver_free (s);
    return k;
}

///!!! nonpublic?
double CPotential::corotationradius() const
{
    if(Omega==0) return 0;  // actually 1/0

    double Xlow=1, Xupp=3, Xmax, Ymax;
    gsl_function F;
    CPotentialParam1 PP;
    PP.E=0;
    PP.P=this;
    F.function=&minusPhi;
    F.params=&PP;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    // X axis
    while ((Phi(Xlow,0,0)>Phi(Xlow*2,0,0)) && (Xlow>MIN_RADIUS_CUTOFF)) Xlow/=2;
    if(Xlow<MIN_RADIUS_CUTOFF) return -1;  // signal of error
    while ((Phi(Xupp,0,0)>Phi(Xlow*2,0,0)) && (Xupp<MAX_RADIUS_CUTOFF)) Xupp*=2;
    if(Xupp>MAX_RADIUS_CUTOFF) return -1;  
    gsl_min_fminimizer_set (s, &F, Xlow*2, Xlow, Xupp);
    int status=0, iter=0;
    do{
        iter++;
        status = gsl_min_fminimizer_iterate (s);
        Xmax = gsl_min_fminimizer_x_minimum (s);
        Xlow = gsl_min_fminimizer_x_lower (s);
        Xupp = gsl_min_fminimizer_x_upper (s);
        status = gsl_root_test_interval (Xlow, Xupp, 0, 1e-9);
    }
    while (status == GSL_CONTINUE && iter < 100);
    // Y axis
    Xlow=1; Xupp=3;
    while ((Phi(0,Xlow,0)>Phi(0,Xlow*2,0)) && (Xlow>MIN_RADIUS_CUTOFF)) Xlow/=2;
    if(Xlow<MIN_RADIUS_CUTOFF) return -1;  // signal of error
    while ((Phi(0,Xupp,0)>Phi(0,Xlow*2,0)) && (Xupp<MAX_RADIUS_CUTOFF)) Xupp*=2;
    if(Xupp>MAX_RADIUS_CUTOFF) return -1;  
    PP.E=1;
    gsl_min_fminimizer_set (s, &F, Xlow*2, Xlow, Xupp);
    status=0, iter=0;
    do{
        iter++;
        status = gsl_min_fminimizer_iterate (s);
        Ymax = gsl_min_fminimizer_x_minimum (s);
        Xlow = gsl_min_fminimizer_x_lower (s);
        Xupp = gsl_min_fminimizer_x_upper (s);
        status = gsl_root_test_interval (Xlow, Xupp, 0, 1e-9);
    }
    while (status == GSL_CONTINUE && iter < 100);
    gsl_min_fminimizer_free (s);
    return (Xmax+Ymax)/2;
}

std::string CPotential::info() const
{
    return PotentialName()+//" y/x="+convertToString(q)+" z/x="+convertToString(p)+
        (Mbh>0 ? " Mbh="+convertToString(Mbh) : "");
}

//----------------------------------------------------------------------------//
// Logarithmic potential

double CPotentialLog::Rho(double X, double Y, double Z, double /*t*/) const
{
    double m2=pow_2(Rc) + pow_2(X) + pow_2(Y)/pow_2(q) + pow_2(Z)/pow_2(p);
    return (0.5*(1+1/pow_2(q)+1/pow_2(p))/m2 - (pow_2(X) + pow_2(Y/pow_2(q)) + pow_2(Z/pow_2(p)))/pow_2(m2))/M_PI;
}

double CPotentialLog::Phi(double X, double Y, double Z, double /*t*/) const
{
    return log(pow_2(Rc) + pow_2(X) + pow_2(Y)/pow_2(q) + pow_2(Z)/pow_2(p)) - ((Mbh>0)?(Mbh/sqrt(sqrt(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z))):0);
}

// r.h.s. of the differential equation: for given time 't' and variables 'v', return their difference 'f'
void CPotentialLog::DiffEq(unsigned n, double /*t*/, double *v, double *f) const
{
    bool calcLyapunov = n>2*N_DIM;
    double X=v[0];
    double Y=v[1];
    double Z=(N_dim==3)?v[2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
    f[0] = v[N_dim];
    f[1] = v[N_dim+1];
    if(N_dim==3)
        f[2] = v[N_dim+2];
    double m2=(pow_2(X) + pow_2(Y)/pow_2(q) + pow_2(Z)/pow_2(p) + pow_2(Rc));
    double r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
    f[N_dim]   = -2*X/m2    - ((Mbh>0)?(Mbh*X/r3):0);
    f[N_dim+1] = -2*Y/m2/pow_2(q) - ((Mbh>0)?(Mbh*Y/r3):0);
    if(N_dim==3)
        f[N_dim+2] = -2*Z/m2/pow_2(p) - ((Mbh>0)?(Mbh*Z/r3):0);
    if(calcLyapunov)
    {
#ifdef LYAPUNOVVAREQ
        // variational equation
        double r5=pow(pow_2(X) + pow_2(Y) + pow_2(Z), 2.5);
        double wx = v[N_dim*2];     // variation vector
        double wy = v[N_dim*2+1];
        double wz = (N_dim==3)? v[N_dim*2+2] : 0;
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2]=v[N_dim*3+2];
        double d2Vdxdy = -X*Y * (4/pow_2(q*m2) + ((Mbh>0)?(3*Mbh/r5):0));
        double d2Vdxdz = (N_dim==3)?( -X*Z * (4/pow_2(p*m2)   + ((Mbh>0)?(3*Mbh/r5):0)) ):0;
        double d2Vdydz = (N_dim==3)?( -Y*Z * (4/pow_2(q*p*m2) + ((Mbh>0)?(3*Mbh/r5):0)) ):0;
        f[N_dim*3]   = -wx * (2/m2   -4*pow_2(X/m2   ) - 3*Mbh*pow_2(X)/r5 + Mbh/r3) - wy*d2Vdxdy - wz*d2Vdxdz;
        f[N_dim*3+1] = -wy * (2/m2/pow_2(q)-4*pow_2(Y/m2/pow_2(q)) - 3*Mbh*pow_2(Y)/r5 + Mbh/r3) - wx*d2Vdxdy - wz*d2Vdydz;
        if(N_dim==3)
            f[N_dim*3+2] = -wz * (2/m2/pow_2(p)-4*pow_2(Z/m2/pow_2(p)) - 3*Mbh*pow_2(Z)/r5 + Mbh/r3) - wx*d2Vdxdz - wy*d2Vdydz;
#else
        // integrate nearby orbit
        X=v[N_dim*2];
        Y=v[N_dim*2+1];
        Z=(N_dim==3)?v[N_dim*2+2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2] = v[N_dim*3+2];
        m2=(pow_2(X) + pow_2(Y)/pow_2(q) + pow_2(Z)/pow_2(p) + pow_2(Rc));
        r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
        f[N_dim*3]   = -2*X/m2    - ((Mbh>0)?(Mbh*X/r3):0);
        f[N_dim*3+1] = -2*Y/m2/pow_2(q) - ((Mbh>0)?(Mbh*Y/r3):0);
        if(N_dim==3)
            f[N_dim*3+2] = -2*Z/m2/pow_2(p) - ((Mbh>0)?(Mbh*Z/r3):0);
#endif
    }
}

//----------------------------------------------------------------------------//
// Dehnen potential
CPotentialDehnen::CPotentialDehnen(unsigned int _N_dim, double _q, double _p, double _Mbh, double _Gamma): 
    CPotential(_N_dim, _Mbh), q(_q), p(_p),
    Gamma(std::max<double>(0, std::min<double>(_Gamma, 2.0))) {}

double CPotentialDehnen::Rho(double X, double Y, double Z, double /*t*/) const
{
    double m = sqrt(pow_2(X) + pow_2(Y/q) + pow_2(Z/p));
    return (3-Gamma)/(4*M_PI*p*q) * pow(m, -Gamma) * pow(1+m, Gamma-4);
}

struct CPotentialParamD{
    double X, Y, Z;
    double wx, wy, wz;
    double p, q, Gamma;
};

double intPhiD(double s, void* params)
{
    double X= ((CPotentialParamD*)params)->X;
    double Y= ((CPotentialParamD*)params)->Y;
    double Z= ((CPotentialParamD*)params)->Z;
    double p= ((CPotentialParamD*)params)->p;
    double q= ((CPotentialParamD*)params)->q;
    double Gamma= ((CPotentialParamD*)params)->Gamma;
    double m = s * sqrt( pow_2(X) + pow_2(Y)/(1-(1-pow_2(q))*pow_2(s)) + pow_2(Z)/(1-(1-pow_2(p))*pow_2(s)));
    double numerator = (Gamma==2)? (log((1+m)*s/m) - 1/(1+m) - log(s)) :
        (1 - (3-Gamma)*pow(m/(m+1), 2-Gamma) + (2-Gamma)*pow(m/(m+1), 3-Gamma))/(2-Gamma);
    return numerator / sqrt( (1-(1-pow_2(q))*pow_2(s)) * (1-(1-pow_2(p))*pow_2(s)));
}

// test runs: accurate adaptive integration of potential and forces (muuuuch slower, used only in high-accuracy runs)
//#define DEHNENPOTENTIAL_ACCURATE_QUADRATURE

double CPotentialDehnen::Phi(double X, double Y, double Z, double /*t*/) const
{
    double result, error;
    CPotentialParamD PP;
    PP.X=X; PP.Y=Y; PP.Z=Z;
    PP.p=p; PP.q=q; PP.Gamma=Gamma;
    gsl_function F;
    F.function=&intPhiD;
    F.params=&PP;
    size_t neval;
    gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
#ifdef DEHNENPOTENTIAL_ACCURATE_QUADRATURE
    if(fabs(error/result)>1e-3) {
        gsl_integration_workspace * ws = gsl_integration_workspace_alloc (1000);
        gsl_integration_qag (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, 1000, 3,ws, &result, &error); 
        gsl_integration_workspace_free(ws);
    }
#endif
    return -result - ((Mbh>0)?(Mbh/sqrt(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z))):0);
}

double intForceDx(double s, void* params)
{
    double X= ((CPotentialParamD*)params)->X;
    double Y= ((CPotentialParamD*)params)->Y;
    double Z= ((CPotentialParamD*)params)->Z;
    double p= ((CPotentialParamD*)params)->p;
    double q= ((CPotentialParamD*)params)->q;
    double Gamma= ((CPotentialParamD*)params)->Gamma;
    double m = s * sqrt( pow_2(X) + pow_2(Y)/(1+(pow_2(q)-1)*pow_2(s)) + pow_2(Z)/(1+(pow_2(p)-1)*pow_2(s)));
    return pow_2(s) * pow(m, -Gamma) * pow(1+m, Gamma-4) /
        sqrt( (1+(pow_2(q)-1)*pow_2(s)) * (1+(pow_2(p)-1)*pow_2(s)));
}
double intForceDy(double s, void* params)
{
    double X= ((CPotentialParamD*)params)->X;
    double Y= ((CPotentialParamD*)params)->Y;
    double Z= ((CPotentialParamD*)params)->Z;
    double p= ((CPotentialParamD*)params)->p;
    double q= ((CPotentialParamD*)params)->q;
    double Gamma= ((CPotentialParamD*)params)->Gamma;
    double m = s * sqrt( pow_2(X)/(pow_2(q)+(1-pow_2(q))*pow_2(s)) + pow_2(Y)/pow_2(q) + pow_2(Z)/(pow_2(q)+(pow_2(p)-pow_2(q))*pow_2(s)));
    return pow_2(s) * pow(m, -Gamma) * pow(1+m, Gamma-4) /
        sqrt( (pow_2(q)+(pow_2(p)-pow_2(q))*pow_2(s)) * (pow_2(q)+(1-pow_2(q))*pow_2(s)));
}
double intForceDz(double s, void* params)
{
    double X= ((CPotentialParamD*)params)->X;
    double Y= ((CPotentialParamD*)params)->Y;
    double Z= ((CPotentialParamD*)params)->Z;
    double p= ((CPotentialParamD*)params)->p;
    double q= ((CPotentialParamD*)params)->q;
    double Gamma= ((CPotentialParamD*)params)->Gamma;
    double m = s * sqrt( pow_2(X)/(pow_2(p)+(1-pow_2(p))*pow_2(s)) + pow_2(Y)/(pow_2(p)+(pow_2(q)-pow_2(p))*pow_2(s)) + pow_2(Z)/pow_2(p));
    return pow_2(s) * pow(m, -Gamma) * pow(1+m, Gamma-4) /
        sqrt( (pow_2(p)+(1-pow_2(p))*pow_2(s)) * (pow_2(p)+(pow_2(q)-pow_2(p))*pow_2(s)));
}

double intForceD2x(double s, void* params)
{
    double X= ((CPotentialParamD*)params)->X;
    double Y= ((CPotentialParamD*)params)->Y;
    double Z= ((CPotentialParamD*)params)->Z;
    double wx=((CPotentialParamD*)params)->wx;
    double wy=((CPotentialParamD*)params)->wy;
    double wz=((CPotentialParamD*)params)->wz;
    double p= ((CPotentialParamD*)params)->p;
    double q= ((CPotentialParamD*)params)->q;
    double Gamma= ((CPotentialParamD*)params)->Gamma;
    double m = s * sqrt( pow_2(X) + pow_2(Y)/(1+(pow_2(q)-1)*pow_2(s)) + pow_2(Z)/(1+(pow_2(p)-1)*pow_2(s)));
    return pow_2(pow_2(s)) * pow(m, -Gamma-2) * pow(1+m, Gamma-5) * (Gamma+4*m) * 
      ( X * wx / sqrt( (1+(pow_2(q)-1)*pow_2(s)) * (1+(pow_2(p)-1)*pow_2(s))) +
        Y * wy * pow( 1 + (pow_2(q)-1)*pow_2(s), -1.5) / sqrt(1 + (pow_2(p)-1)*pow_2(s)) +
        Z * wz * pow( 1 + (pow_2(p)-1)*pow_2(s), -1.5) / sqrt(1 + (pow_2(q)-1)*pow_2(s)) );
}
double intForceD2y(double s, void* params)
{
    double X= ((CPotentialParamD*)params)->X;
    double Y= ((CPotentialParamD*)params)->Y;
    double Z= ((CPotentialParamD*)params)->Z;
    double wx=((CPotentialParamD*)params)->wx;
    double wy=((CPotentialParamD*)params)->wy;
    double wz=((CPotentialParamD*)params)->wz;
    double p= ((CPotentialParamD*)params)->p;
    double q= ((CPotentialParamD*)params)->q;
    double Gamma= ((CPotentialParamD*)params)->Gamma;
    double m = s * sqrt( pow_2(X)/(pow_2(q)+(1-pow_2(q))*pow_2(s)) + pow_2(Y)/pow_2(q) + pow_2(Z)/(pow_2(q)+(pow_2(p)-pow_2(q))*pow_2(s)));
    return pow_2(pow_2(s)) * pow(m, -Gamma-2) * pow(1+m, Gamma-5) * (Gamma+4*m) * 
      ( Y * wy / pow_2(q) / sqrt( (pow_2(q)+(pow_2(p)-pow_2(q))*pow_2(s)) * (pow_2(q)+(1-pow_2(q))*pow_2(s))) +
        Z * wz * pow( pow_2(q) + (pow_2(p)-pow_2(q))*pow_2(s), -1.5) / sqrt(pow_2(q) + (1-pow_2(q))*pow_2(s)) +
        X * wx * pow( pow_2(q) + (1-pow_2(q))*pow_2(s), -1.5) / sqrt(pow_2(q) + (pow_2(p)-pow_2(q))*pow_2(s)) );
}
double intForceD2z(double s, void* params)
{
    double X= ((CPotentialParamD*)params)->X;
    double Y= ((CPotentialParamD*)params)->Y;
    double Z= ((CPotentialParamD*)params)->Z;
    double wx=((CPotentialParamD*)params)->wx;
    double wy=((CPotentialParamD*)params)->wy;
    double wz=((CPotentialParamD*)params)->wz;
    double p= ((CPotentialParamD*)params)->p;
    double q= ((CPotentialParamD*)params)->q;
    double Gamma= ((CPotentialParamD*)params)->Gamma;
    double m = s * sqrt( pow_2(X)/(pow_2(p)+(1-pow_2(p))*pow_2(s)) + pow_2(Y)/(pow_2(p)+(pow_2(q)-pow_2(p))*pow_2(s)) + pow_2(Z)/pow_2(p));
    return pow_2(pow_2(s)) * pow(m, -Gamma-2) * pow(1+m, Gamma-5) * (Gamma+4*m) * 
      ( Z * wz / pow_2(p) / sqrt( (pow_2(p)+(1-pow_2(p))*pow_2(s)) * (pow_2(p)+(pow_2(q)-pow_2(p))*pow_2(s))) +
        X * wx * pow( pow_2(p) + (1-pow_2(p))*pow_2(s), -1.5) / sqrt(pow_2(p) + (pow_2(q)-pow_2(p))*pow_2(s)) +
        Y * wy * pow( pow_2(p) + (pow_2(q)-pow_2(p))*pow_2(s), -1.5) / sqrt(pow_2(p) + (1-pow_2(p))*pow_2(s)) );
}

void CPotentialDehnen::DiffEq(unsigned n, double /*t*/, double *v, double *f) const
{
    bool calcLyapunov=n>2*N_DIM;
    double X=v[0];
    double Y=v[1];
    double Z=(N_dim==3)?v[2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
    f[0] = v[N_dim];
    f[1] = v[N_dim+1];
    if(N_dim==3)
        f[2] = v[N_dim+2];
    double r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
    // force calculation requires integration
    double resultx, resulty, resultz, error; 
    size_t neval;
    CPotentialParamD PP;
    PP.X=X; PP.Y=Y; PP.Z=Z;
    PP.p=p; PP.q=q; PP.Gamma=Gamma;
    gsl_function F;
    F.params=&PP;
    F.function=&intForceDx;
    gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &resultx, &error, &neval); 
#ifdef DEHNENPOTENTIAL_ACCURATE_QUADRATURE
    gsl_integration_workspace * ws = NULL;
    if(fabs(error/resultx)>0.001) {
        if(!ws) ws=gsl_integration_workspace_alloc (1000);
        gsl_integration_qag (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, 1000, 3,ws, &resultx, &error); 
    }
#endif
    f[N_dim]   = -(3-Gamma)*X * resultx   - ((Mbh>0)?(Mbh*X/r3):0);
    F.function=&intForceDy;
    gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &resulty, &error, &neval); 
#ifdef DEHNENPOTENTIAL_ACCURATE_QUADRATURE
    if(fabs(error/resulty)>0.001) {
        if(!ws) ws=gsl_integration_workspace_alloc (1000);
        gsl_integration_qag (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, 1000, 3,ws, &resulty, &error); 
    }
#endif
    f[N_dim+1] = -(3-Gamma)*Y/q * resulty - ((Mbh>0)?(Mbh*Y/r3):0);
    if(N_dim==3)
    {
        F.function=&intForceDz;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &resultz, &error, &neval); 
#ifdef DEHNENPOTENTIAL_ACCURATE_QUADRATURE
        if(fabs(error/resultz)>0.001) {
            if(!ws) ws=gsl_integration_workspace_alloc (1000);
            gsl_integration_qag (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, 1000, 3,ws, &resultz, &error); 
        }
#endif
        f[N_dim+2] = -(3-Gamma)*Z/p * resultz - ((Mbh>0)?(Mbh*Z/r3):0);
    }
#ifdef DEHNENPOTENTIAL_ACCURATE_QUADRATURE
    gsl_integration_workspace_free(ws);
#endif
    if(calcLyapunov)
    {
#ifdef LYAPUNOVVAREQ
        // variational equation
        double r5=pow(pow_2(X) + pow_2(Y) + pow_2(Z), 2.5);
        double wx = v[N_dim*2];     // variation vector
        double wy = v[N_dim*2+1];
        double wz = (N_dim==3)? v[N_dim*2+2] : 0;
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2]=v[N_dim*3+2];
        CPotentialParamD2 P2;
        P2.X=X; P2.Y=Y; P2.Z=Z;
        P2.wx=wx; P2.wy=wy; P2.wz=wz;
        P2.P=this;
        double result;
        F.params=&P2;
        F.function=&intForceD2x;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
        f[N_dim*3]   = -(3-Gamma) * (wx*resultx - X*result) + ((Mbh>0)?Mbh*( (3*pow_2(X)/r5-1/r3)*wx + 3*X*Y/r5*wy + 3*X*Z/r5*wz ):0);
        F.function=&intForceD2y;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
        f[N_dim*3+1] = -(3-Gamma)/q * (wy*resulty - Y*result) + ((Mbh>0)?Mbh*( (3*pow_2(Y)/r5-1/r3)*wy + 3*Y*X/r5*wx + 3*Y*Z/r5*wz ):0);
        if(N_dim==3)
        {
            F.function=&intForceD2z;
            gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
            f[N_dim*3+2] = -(3-Gamma)/p * (wz*resultz - Z*result) + ((Mbh>0)?Mbh*( (3*pow_2(Z)/r5-1/r3)*wz + 3*Z*X/r5*wx + 3*Z*Y/r5*wy ):0);
        }
        /*double d2Vdxdy = -X*Y * (4/pow_2(q*m2) + ((Mbh>0)?(3*Mbh/r5):0));
        double d2Vdxdz = (N_dim==3)?( -X*Z * (4/pow_2(p*m2)   + ((Mbh>0)?(3*Mbh/r5):0)) ):0;
        double d2Vdydz = (N_dim==3)?( -Y*Z * (4/pow_2(q*p*m2) + ((Mbh>0)?(3*Mbh/r5):0)) ):0;
        f[N_dim*3]   = -wx * (2/m2   -4*pow_2(X/m2   ) - 3*Mbh*pow_2(X)/r5 + Mbh/r3) - wy*d2Vdxdy - wz*d2Vdxdz;
        f[N_dim*3+1] = -wy * (2/m2/pow_2(q)-4*pow_2(Y/m2/pow_2(q)) - 3*Mbh*pow_2(Y)/r5 + Mbh/r3) - wx*d2Vdxdy - wz*d2Vdydz;
        if(N_dim==3)
            f[N_dim*3+2] = -wz * (2/m2/pow_2(p)-4*pow_2(Z/m2/pow_2(p)) - 3*Mbh*pow_2(Z)/r5 + Mbh/r3) - wx*d2Vdxdz - wy*d2Vdydz;*/
#else
        // integrate nearby orbit
        X=v[N_dim*2];
        Y=v[N_dim*2+1];
        Z=(N_dim==3)?v[N_dim*2+2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2] = v[N_dim*3+2];
        r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
        // force calculation requires integration
        PP.X=X; PP.Y=Y; PP.Z=Z;
        F.function=&intForceDx;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &resultx, &error, &neval); 
        f[N_dim*3]   = -(3-Gamma)*X * resultx    - ((Mbh>0)?(Mbh*X/r3):0);
        F.function=&intForceDy;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &resulty, &error, &neval); 
        f[N_dim*3+1] = -(3-Gamma)*Y/q * resulty    - ((Mbh>0)?(Mbh*Y/r3):0);
        if(N_dim==3)
        {
            F.function=&intForceDz;
            gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &resultz, &error, &neval); 
            f[N_dim*3+2] = -(3-Gamma)*Z/p * resultz    - ((Mbh>0)?(Mbh*Z/r3):0);
        }
#endif
    }
}

//----------------------------------------------------------------------------//
// CPotentialSH -- parent class for all potentials using angular expansion in spherical harmonics

// helper functions for integration over angles (common for BSE and Spline potentials), used in computing coefficients from analytic density model
struct CPotentialParamSH{
    const CDensity* P;
    int sh_n, sh_l, sh_m;
    double costheta, sintheta, r, Alpha;
    double theta_max, phi_max;    // Pi/2, Pi or even 2*Pi, depending on potential symmetry
};
double intSH_phi(double phi, void* params)
{
    double costheta= ((CPotentialParamSH*)params)->costheta;
    double sintheta= ((CPotentialParamSH*)params)->sintheta;
    double r= ((CPotentialParamSH*)params)->r;
    const CDensity* P=((CPotentialParamSH*)params)->P;
    int sh_l=((CPotentialParamSH*)params)->sh_l;
    int sh_m=((CPotentialParamSH*)params)->sh_m;
    return P->Rho(r*sintheta*cos(phi), r*sintheta*sin(phi), r*costheta)
        * gsl_sf_legendre_sphPlm(sh_l, abs(sh_m), costheta) * (sh_m>=0 ? cos(sh_m*phi) : sin(-sh_m*phi));
}
double intSH_theta(double theta, void* params)
{
    double result, error;
    size_t neval;
    ((CPotentialParamSH*)params)->costheta=cos(theta);
    ((CPotentialParamSH*)params)->sintheta=sin(theta);
    gsl_function F;
    F.function=&intSH_phi;
    F.params=params;
    gsl_integration_qng(&F, 0, ((CPotentialParamSH*)params)->phi_max, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval);
    return sin(theta)*result;
}

double CPotentialSH::Rho(double X, double Y, double Z, double /*t*/) const
{
    double costheta=(X==0&&Y==0&&Z==0)?0:(Z/sqrt(pow_2(X)+pow_2(Y)+pow_2(Z)));
    double phi=atan2(Y, X);
    double r=sqrt(X*X+Y*Y+Z*Z);
    double result=0;
    // arrays where angular expansion coefficients will be accumulated by calling computeSHcoefs() for derived classes
    double coefsF[MAX_NCOEFS_ANGULAR*MAX_NCOEFS_ANGULAR];      // F(theta,phi)
    double coefsdFdr[MAX_NCOEFS_ANGULAR*MAX_NCOEFS_ANGULAR];   // dF(theta,phi)/dr
    double coefsd2Fdr2[MAX_NCOEFS_ANGULAR*MAX_NCOEFS_ANGULAR]; // d2F(theta,phi)/dr2
    computeSHCoefs(r, coefsF, coefsdFdr, coefsd2Fdr2);         // implemented in derived classes
    double legendre_array[MAX_NCOEFS_ANGULAR-1];
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : Ncoefs_angular; // if axisymmetric model, use only m=0 terms
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    bool msine = mysymmetry & ST_PLANESYM ? false : true;      // if triaxial symmetry, do not use sine terms
    for(int m=0; m<=mmax; m+=mstep)
    {
        gsl_sf_legendre_sphPlm_array(lmax, m, costheta, legendre_array);
        int lmin = lstep==2 ? (m+1)/2*2 : m;   // if lstep is even and m is odd, start from next even number greater than m
        double cosmphi = cos(m*phi);
        for(int l=lmin; l<=lmax; l+=lstep)
        {
            result += (coefsd2Fdr2[l*(l+1)+m] + 2/r*coefsdFdr[l*(l+1)+m] - l*(l+1)/(r*r)*coefsF[l*(l+1)+m]) * legendre_array[l-m] * cosmphi;
        }
        if(msine && m>0)
        {
            double sinmphi = sin(m*phi);
            for(int l=lmin; l<=lmax; l+=mstep)
                result += (coefsd2Fdr2[l*(l+1)-m] + 2/r*coefsdFdr[l*(l+1)-m] - l*(l+1)/(r*r)*coefsF[l*(l+1)-m]) * legendre_array[l-m] * sinmphi;
        }
    }
    result *= -1/(4*M_PI);
    return std::max<double> (result,0);
}

double CPotentialSH::Phi(double X, double Y, double Z, double /*t*/) const
{
    double costheta=(X==0&&Y==0&&Z==0)?0:(Z/sqrt(pow_2(X)+pow_2(Y)+pow_2(Z)));
    double phi=atan2(Y, X);
    double r=sqrt(X*X+Y*Y+Z*Z);
    double res=0;
    // arrays where angular expansion coefficients will be accumulated by calling computeSHcoefs() for derived classes
    double coefsF[MAX_NCOEFS_ANGULAR*MAX_NCOEFS_ANGULAR];      // F(theta,phi)
    computeSHCoefs(r, coefsF, NULL, NULL);  // implemented in derived classes
    double legendre_array[MAX_NCOEFS_ANGULAR-1];
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : Ncoefs_angular;   // if axisymmetric model, use only m=0 terms
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    bool msine = mysymmetry & ST_PLANESYM ? false : true;      // if triaxial symmetry, do not use sine terms
    for(int m=0; m<=mmax; m+=mstep)
    {
        gsl_sf_legendre_sphPlm_array(lmax, m, costheta, legendre_array);
        double cosmphi = cos(m*phi);
        int lmin = lstep==2 ? (m+1)/2*2 : m;   // if lstep is even and m is odd, start from next even number greater than m
        for(int l=lmin; l<=lmax; l+=lstep)
            res +=  coefsF[l*(l+1)+m] * legendre_array[l-m] * cosmphi;
        if(msine && m>0)
        {
            double sinmphi = sin(m*phi);
            for(int l=lmin; l<=lmax; l+=lstep)
                res +=  coefsF[l*(l+1)-m] * legendre_array[l-m] * sinmphi;
        }
    }
    return -res - ((Mbh>0)?(Mbh/sqrt(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z))):0) - ((Omega>0)?pow_2(Omega)*(pow_2(X)+pow_2(Y))/2:0);
}

void CPotentialSH::DiffEq(unsigned n, double /*t*/, double *v, double *f) const
{
    bool calcLyapunov=n>2*N_DIM;
    double X=v[0];
    double Y=v[1];
    double Z=(N_dim==3)?v[2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
    f[0] = v[N_dim];
    f[1] = v[N_dim+1];
    if(N_dim==3)
        f[2] = v[N_dim+2];

    double dFdR=0, dFdTheta=0, dFdPhi=0;   
#ifdef LYAPUNOVVAREQ
    double d2FdR2=0, d2FdRdTheta=0, d2FdRdPhi=0, d2FdTheta2=0, d2FdThetadPhi=0, d2FdPhi2=0;
#endif
    double r2 = pow_2(X) + pow_2(Y) + pow_2(Z);
    double r  = sqrt(r2);
    double rxy2=pow_2(X)+pow_2(Y);
    double rxy =sqrt(rxy2);
    double theta=atan2(rxy, Z); 
    double phi  =atan2(Y, X);
    double costheta=(Z==0)?0:cos(theta);
    double sintheta=sin(theta);
    double cosphi=(X==0)?0:cos(phi);
    double sinphi=(Y==0)?0:sin(phi);
    // arrays where angular expansion coefficients will be accumulated by calling computeSHcoefs() for derived classes
    double coefsF[MAX_NCOEFS_ANGULAR*MAX_NCOEFS_ANGULAR];      // F(theta,phi)
    double coefsdFdr[MAX_NCOEFS_ANGULAR*MAX_NCOEFS_ANGULAR];   // dF(theta,phi)/dr
#ifdef LYAPUNOVVAREQ
    double coefsd2Fdr2[MAX_NCOEFS_ANGULAR*MAX_NCOEFS_ANGULAR]; // d2F(theta,phi)/dr2
#else
    double* coefsd2Fdr2=NULL;  // not needed to calculate
#endif
    computeSHCoefs(r, coefsF, coefsdFdr, coefsd2Fdr2);  // implemented in derived classes
    double legendre_array[MAX_NCOEFS_ANGULAR-1];
    double legendre_deriv_array[MAX_NCOEFS_ANGULAR-1];
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : Ncoefs_angular;   // if axisymmetric model, use only m=0 terms
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    bool msine = mysymmetry & ST_PLANESYM ? false : true;      // if triaxial symmetry, do not use sine terms
    for(int m=0; m<=mmax; m+=mstep)
    {
        gsl_sf_legendre_sphPlm_deriv_array(lmax, m, costheta, legendre_array, legendre_deriv_array);
        double cosmphi = cos(m*phi);
        double sinmphi = sin(m*phi);
        int lmin = lstep==2 ? (m+1)/2*2 : m;   // if lstep is even and m is odd, start from next even number greater than m
        for(int l=lmin; l<=lmax; l+=lstep)
        {
            int indx=l*(l+1)+m;
            dFdR +=  coefsdFdr[indx] * legendre_array[l-m] * cosmphi;
            dFdTheta += coefsF[indx] * legendre_deriv_array[l-m] * (-sintheta) * cosmphi;
            dFdPhi   += coefsF[indx] * legendre_array[l-m] * (-m)*sinmphi;
#ifdef LYAPUNOVVAREQ
            if(calcLyapunov)
            {
                d2FdR2 +=  coefsd2Fdr2[indx] * legendre_array[l-m] * cosmphi;
                d2FdRdTheta+=coefsdFdr[indx] * legendre_deriv_array[l-m] * (-sintheta) * cosmphi;
                d2FdRdPhi += coefsdFdr[indx] * legendre_array[l-m] * (-m)*sinmphi;
                d2FdTheta2   += coefsF[indx] * (costheta * legendre_deriv_array[l-m] - (l*(l+1)-pow_2(m/sintheta)) * legendre_array[l-m]) * cosmphi;
                d2FdThetadPhi+= coefsF[indx] * legendre_deriv_array[l-m] * (-sintheta) * (-m)*sinmphi;
                d2FdPhi2     += coefsF[indx] * legendre_array[l-m] * cosmphi * -pow_2(m);
            }
#endif
            if(msine && m>0)
            {
                indx=l*(l+1)-m;
                dFdR +=  coefsdFdr[indx] * legendre_array[l-m] * sinmphi;
                dFdTheta += coefsF[indx] * legendre_deriv_array[l-m] * (-sintheta) * sinmphi;
                dFdPhi   += coefsF[indx] * legendre_array[l-m] * m*cosmphi;
#ifdef LYAPUNOVVAREQ
                if(calcLyapunov)
                {
                    d2FdR2 +=  coefsd2Fdr2[indx] * legendre_array[l-m] * sinmphi;
                    d2FdRdTheta+=coefsdFdr[indx] * legendre_deriv_array[l-m] * (-sintheta) * sinmphi;
                    d2FdRdPhi += coefsdFdr[indx] * legendre_array[l-m] * m*cosmphi;
                    d2FdTheta2   += coefsF[indx] * (costheta * legendre_deriv_array[l-m] - (l*(l+1)-pow_2(m/sintheta)) * legendre_array[l-m]) * sinmphi;
                    d2FdThetadPhi+= coefsF[indx] * legendre_deriv_array[l-m] * (-sintheta) * m*cosmphi;
                    d2FdPhi2     += coefsF[indx] * legendre_array[l-m] * sinmphi * -pow_2(m);
                }
#endif
            }
        }
    }
        
    //if((X==0)||(Y==0)) dFdPhi=0; // to avoid singularity
    //if(Z==0) dFdTheta=0;
    if(rxy==0) rxy=BH_SMOOTH; 
    double r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
    double J[3][3]   = { {sintheta*cosphi, sintheta*sinphi, costheta}, {costheta*cosphi/r, costheta*sinphi/r, -sintheta/r}, {-sinphi/rxy, cosphi/rxy, 0} };
    f[N_dim]   = ( dFdR * J[0][0] + dFdTheta * J[1][0] + dFdPhi * J[2][0]) - ((Mbh>0)?(Mbh*X/r3):0) + ((Omega>0)?  2*Omega*v[N_dim+1]+pow_2(Omega)*X : 0) ;
    f[N_dim+1] = ( dFdR * J[0][1] + dFdTheta * J[1][1] + dFdPhi * J[2][1]) - ((Mbh>0)?(Mbh*Y/r3):0) + ((Omega>0)? -2*Omega*v[N_dim]  +pow_2(Omega)*Y : 0) ;
    if((N_dim==3))
        f[N_dim+2] = ( dFdR * J[0][2] + dFdTheta * J[1][2]) - ((Mbh>0)?(Mbh*Z/r3):0);
#ifdef LYAPUNOVVAREQ
    if(calcLyapunov)
    {
        double X2=X*X, Y2=Y*Y, Z2=Z*Z;
        double d2F[3][3] = { { d2FdR2, d2FdRdTheta, d2FdRdPhi }, {d2FdRdTheta, d2FdTheta2, d2FdThetadPhi}, {d2FdRdPhi, d2FdThetadPhi, d2FdPhi2} };
        double m=1/(r*r2);
        double d2R[3][3] = { {(Y2+Z2)*m, -X*Y*m, -X*Z*m}, {-X*Y*m, (X2+Z2)*m, -Y*Z*m}, {-X*Z*m, -Y*Z*m, (X2+Y2)*m} };
        m=1/(rxy*rxy2*r2*r2);
        double d2Theta[3][3] = { 
            {-Z*(2*X2*X2+Y2*X2-Y2*Y2-Z2*Y2)*m, -X*Y*Z*(3*X2+3*Y2+Z2)*m, X*(X2+Y2-Z2)*rxy2*m},
            {-X*Y*Z*(3*X2+3*Y2+Z2)*m, -Z*(2*Y2*Y2+Y2*X2-X2*X2-Z2*X2)*m, Y*(X2+Y2-Z2)*rxy2*m},
            {X*(X2+Y2-Z2)*rxy2*m, Y*(X2+Y2-Z2)*rxy2*m, 2*Z*rxy2*rxy2*m} };
        m=1/(rxy2*rxy2);
        double d2Phi[3][3] = { {2*X*Y*m, (Y2-X2)*m, 0}, {(Y2-X2)*m, -2*X*Y*m, 0}, {0, 0, 0} };
        // convert to Cartesian coordinates
        double r5=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 2.5);
        double k0 = (J[0][0]*d2F[0][0]+J[1][0]*d2F[0][1]+J[2][0]*d2F[0][2]);
        double k1 = (J[0][0]*d2F[0][1]+J[1][0]*d2F[1][1]+J[2][0]*d2F[1][2]);
        double k2 = (J[0][0]*d2F[0][2]+J[1][0]*d2F[1][2]+J[2][0]*d2F[2][2]);
        double d2Fdx2 = (k0*J[0][0] + k1*J[1][0] + k2*J[2][0] + dFdR*d2R[0][0] + dFdTheta*d2Theta[0][0] + dFdPhi*d2Phi[0][0]) +Mbh*(3*X2/r5-1/r3);
        double d2Fdxdy= (k0*J[0][1] + k1*J[1][1] + k2*J[2][1] + dFdR*d2R[0][1] + dFdTheta*d2Theta[0][1] + dFdPhi*d2Phi[0][1]) +Mbh*(3*X*Y/r5);
        double d2Fdxdz= (k0*J[0][2] + k1*J[1][2] + k2*J[2][2] + dFdR*d2R[0][2] + dFdTheta*d2Theta[0][2] + dFdPhi*d2Phi[0][2]) +Mbh*(3*X*Z/r5);
        k0 = (J[0][1]*d2F[0][0]+J[1][1]*d2F[0][1]+J[2][1]*d2F[0][2]);
        k1 = (J[0][1]*d2F[0][1]+J[1][1]*d2F[1][1]+J[2][1]*d2F[1][2]);
        k2 = (J[0][1]*d2F[0][2]+J[1][1]*d2F[1][2]+J[2][1]*d2F[2][2]);
        double d2Fdydx= (k0*J[0][0] + k1*J[1][0] + k2*J[2][0] + dFdR*d2R[1][0] + dFdTheta*d2Theta[1][0] + dFdPhi*d2Phi[1][0]) +Mbh*(3*Y*X/r5);
        double d2Fdy2 = (k0*J[0][1] + k1*J[1][1] + k2*J[2][1] + dFdR*d2R[1][1] + dFdTheta*d2Theta[1][1] + dFdPhi*d2Phi[1][1]) +Mbh*(3*Y2/r5-1/r3);
        double d2Fdydz= (k0*J[0][2] + k1*J[1][2] + k2*J[2][2] + dFdR*d2R[1][2] + dFdTheta*d2Theta[1][2] + dFdPhi*d2Phi[1][2]) +Mbh*(3*Y*Z/r5);
        k0 = (J[0][2]*d2F[0][0]+J[1][2]*d2F[0][1]+J[2][2]*d2F[0][2]);
        k1 = (J[0][2]*d2F[0][1]+J[1][2]*d2F[1][1]+J[2][2]*d2F[1][2]);
        k2 = (J[0][2]*d2F[0][2]+J[1][2]*d2F[1][2]+J[2][2]*d2F[2][2]);
        double d2Fdzdx= (k0*J[0][0] + k1*J[1][0] + k2*J[2][0] + dFdR*d2R[2][0] + dFdTheta*d2Theta[2][0] + dFdPhi*d2Phi[2][0]) +Mbh*(3*Z*X/r5);
        double d2Fdzdy= (k0*J[0][1] + k1*J[1][1] + k2*J[2][1] + dFdR*d2R[2][1] + dFdTheta*d2Theta[2][1] + dFdPhi*d2Phi[2][1]) +Mbh*(3*Z*Y/r5);
        double d2Fdz2 = (k0*J[0][2] + k1*J[1][2] + k2*J[2][2] + dFdR*d2R[2][2] + dFdTheta*d2Theta[2][2] + dFdPhi*d2Phi[2][2]) +Mbh*(3*Z2/r5-1/r3);

        double wx = v[N_dim*2];     // variation vector
        double wy = v[N_dim*2+1];
        double wz = (N_dim==3)? v[N_dim*2+2] : 0;
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2]=v[N_dim*3+2];
        f[N_dim*3]   = (wx*d2Fdx2  + wy*d2Fdxdy + wz*d2Fdxdz);
        f[N_dim*3+1] = (wx*d2Fdydx + wy*d2Fdy2  + wz*d2Fdydz);
        if(N_dim==3)
            f[N_dim*3+2] = (wx*d2Fdzdx + wy*d2Fdzdy + wz*d2Fdz2);
    }
#else
    if(calcLyapunov)
    {
        // integrate nearby orbit
        X=v[N_dim*2];
        Y=v[N_dim*2+1];
        Z=(N_dim==3)?v[N_dim*2+2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2] = v[N_dim*3+2];
        r2 = pow_2(X) + pow_2(Y) + pow_2(Z);
        r  = sqrt(r2);
        rxy2=pow_2(X)+pow_2(Y);
        rxy =sqrt(rxy2);
        theta=atan2(rxy, Z); 
        phi  =atan2(Y, X);
        costheta=(Z==0)?0:cos(theta);
        sintheta=sin(theta);
        cosphi=(X==0)?0:cos(phi);
        sinphi=(Y==0)?0:sin(phi);
        dFdR=0, dFdTheta=0, dFdPhi=0;   

        computeSHCoefs(r, coefsF, coefsdFdr, coefsd2Fdr2);  // implemented in derived classes

        for(int m=0; m<=mmax; m+=mstep)
        {
            gsl_sf_legendre_sphPlm_deriv_array(lmax, m, costheta, legendre_array, legendre_deriv_array);
            double cosmphi = cos(m*phi);
            double sinmphi = sin(m*phi);
            int lmin = lstep==2 ? (m+1)/2*2 : m;   // if lstep is even and m is odd, start from next even number greater than m
            for(int l=lmin; l<=lmax; l+=lstep)
            {
                int indx=l*(l+1)+m;
                dFdR +=  coefsdFdr[indx] * legendre_array[l-m] * cosmphi;
                dFdTheta += coefsF[indx] * legendre_deriv_array[l-m] * (-sintheta) * cosmphi;
                dFdPhi   += coefsF[indx] * legendre_array[l-m] * (-m)*sinmphi;
                if(msine && m>0)
                {
                    int indx=l*(l+1)-m;
                    dFdR +=  coefsdFdr[indx] * legendre_array[l-m] * sinmphi;
                    dFdTheta += coefsF[indx] * legendre_deriv_array[l-m] * (-sintheta) * sinmphi;
                    dFdPhi   += coefsF[indx] * legendre_array[l-m] * m*cosmphi;
                }
            }
        }
        //if((X==0)||(Y==0)) dFdPhi=0; // to avoid singularity
        //if(Z==0) dFdTheta=0;
        if(rxy==0) rxy=BH_SMOOTH; 
        double r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
        double J[3][3]   = { {sintheta*cosphi, sintheta*sinphi, costheta}, {costheta*cosphi/r, costheta*sinphi/r, -sintheta/r}, {-sinphi/rxy, cosphi/rxy, 0} };
        f[N_dim*3]   = ( dFdR * J[0][0] + dFdTheta * J[1][0] + dFdPhi * J[2][0]) - ((Mbh>0)?(Mbh*X/r3):0) + ((Omega>0)?  2*Omega*v[N_dim*3+1]+pow_2(Omega)*X : 0) ;
        f[N_dim*3+1] = ( dFdR * J[0][1] + dFdTheta * J[1][1] + dFdPhi * J[2][1]) - ((Mbh>0)?(Mbh*Y/r3):0) + ((Omega>0)? -2*Omega*v[N_dim*3]  +pow_2(Omega)*Y : 0) ;
        if((N_dim==3))
            f[N_dim*3+2] = ( dFdR * J[0][2] + dFdTheta * J[1][2]) - ((Mbh>0)?(Mbh*Z/r3):0);

    }
#endif
}

//----------------------------------------------------------------------------//
// Basis-set expansion for arbitrary potential (using Zhao(1995) basis set)

CPotentialBSE::CPotentialBSE(unsigned int _N_dim, double _Mbh, double _Alpha, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
    const CPointMassSetDouble &points, SYMMETRYTYPE _sym):
    CPotentialSH(_N_dim, _Mbh, _Ncoefs_angular), 
    Ncoefs_radial(std::min<unsigned int>(MAX_NCOEFS_RADIAL, _Ncoefs_radial)),
    Alpha(_Alpha)
{
    mysymmetry=_sym;
    if(points.empty()) { initDefault(); return; }
    prepareCoefsDiscrete(points);
    checkSymmetry();
}

CPotentialBSE::CPotentialBSE(unsigned int _N_dim, double _Mbh, double _Alpha, const std::vector< vectord > &coefs):
    CPotentialSH(_N_dim, _Mbh, (assert(coefs.size()>0), static_cast<unsigned int>(sqrt(coefs[0].size()*1.0)-1))), 
    Ncoefs_radial(std::min<unsigned int>(MAX_NCOEFS_RADIAL, static_cast<unsigned int>(coefs.size()-1))),
    Alpha(std::max<double>(0.5,_Alpha))  // here Alpha!=0 - no autodetect
{
    if(_Alpha<0.5) { my_error("Error, invalid parameter Alpha="+convertToString(_Alpha)); initDefault(); return; }
    assert(coefs[0].size()==pow_2(Ncoefs_angular+1));
    SHcoefs = coefs;  // do not check array size -- assumed to be done beforehand (in load routine)...
    checkSymmetry();
}

CPotentialBSE::CPotentialBSE(unsigned int _N_dim, double _Mbh, double _Alpha, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
    const CDensity* density):    // init potential from analytic mass model
    CPotentialSH(_N_dim, _Mbh, _Ncoefs_angular), 
    Ncoefs_radial(std::min<unsigned int>(MAX_NCOEFS_RADIAL, _Ncoefs_radial)),
    Alpha(_Alpha)
{
    mysymmetry=density->symmetry();
    prepareCoefsAnalytic(density);
    checkSymmetry();
}
void CPotentialBSE::initDefault()
{
    Ncoefs_radial=0; 
    Ncoefs_angular=0;
    Alpha=1;
    SHcoefs.resize(1);
    SHcoefs[0].assign(1,1);
}

void CPotentialBSE::checkSymmetry()
{   
    SYMMETRYTYPE sym=ST_SPHERICAL;  // too optimistic:))
    const double MINCOEF=1e-8;
    for(size_t n=0; n<=Ncoefs_radial; n++)
    {
        for(int l=0; l<=(int)Ncoefs_angular; l++)
            for(int m=-l; m<=l; m++)
                if(fabs(SHcoefs[n][l*(l+1)+m])>MINCOEF)  
                {   // nonzero coef.: check if that breaks any symmetry
                    if(l%2==1)  sym = (SYMMETRYTYPE)(sym & ~ST_REFLECTION);
                    if(m<0 || m%2==1)  sym = (SYMMETRYTYPE)(sym & ~ST_PLANESYM);
                    if(m!=0) sym = (SYMMETRYTYPE)(sym & ~ST_ZROTSYM);
                    if(l>0) sym = (SYMMETRYTYPE)(sym & ~ST_SPHSYM);
                }
    }
    // now set all coefs excluded by the inferred symmetry  to zero
    for(size_t n=0; n<=Ncoefs_radial; n++)
    {
        for(int l=0; l<=(int)Ncoefs_angular; l++)
            for(int m=-l; m<=l; m++)
                if( (l>0 && (sym & ST_SPHSYM)) ||
                    (m!=0 && (sym & ST_ZROTSYM)) ||
                    ((m<0 || m%2==1) && (sym & ST_PLANESYM)) ||
                    (l%2==1 && (sym & ST_REFLECTION)) )  
                        SHcoefs[n][l*(l+1)+m] = 0;
    }
    mysymmetry = sym; 
    bool densityNonzero = checkDensityNonzero();
    bool massMonotonic  = checkMassMonotonic();
    if(!massMonotonic || !densityNonzero) 
        my_message(std::string("Warning, ") + 
        (!massMonotonic ? "mass does not monotonically increase with radius" : "") +
        (!massMonotonic && !densityNonzero ? " and " : "") + 
        (!densityNonzero ? "density drops to zero at a finite radius" : "") + "!");
}

double intBSE_xi(double xi, void* params)   // integrate over scaled radial variable
{
    if(xi>=1.0 || xi<=-1.0) return 0;
    double Alpha=((CPotentialParamSH*)params)->Alpha;
    double r=pow((1+xi)/(1-xi), Alpha);
    int n=((CPotentialParamSH*)params)->sh_n;
    int l=((CPotentialParamSH*)params)->sh_l;
    double gp=gsl_sf_gegenpoly_n(n, (2*l+1)*Alpha+0.5, xi);
    ((CPotentialParamSH*)params)->r=r;
    gsl_function F;
    F.function=&intSH_theta;
    F.params=params;
    double result, error;
    size_t neval;
    gsl_integration_qng(&F, 0, ((CPotentialParamSH*)params)->theta_max, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval);
    return result * gp * pow(r, 3-1/Alpha+l) * pow(1+pow(r,1/Alpha), 2-(2*l+1)*Alpha) * Alpha * 0.5 * 2*M_SQRTPI;
}

void CPotentialBSE::prepareCoefsAnalytic(const CDensity* density)
{
    if(Alpha<0.5)
    {  // find best-suited value depending on inner density slope
        double Gamma=density->getGamma();
        Alpha = (Gamma<=1? 1.0+Gamma :
        (Gamma<=1.5? 2.0 :
        (Gamma<1.75? 1./(2-Gamma) : 4.0)));
        if(density->PotentialType()==CDensity::PT_PLUMMER) Alpha=0.5;  // the default value for plummer potential only; however even in this case alpha=1 might perform better..
    }
    SHcoefs.resize(Ncoefs_radial+1);
    for(size_t n=0; n<=Ncoefs_radial; n++)
        SHcoefs[n].assign(pow_2(Ncoefs_angular+1), 0);
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    CPotentialParamSH PP;
    PP.P=density;
    PP.Alpha=Alpha;
    PP.theta_max= mysymmetry & ST_REFLECTION ? M_PI_2 : M_PI;  // if symmetries exist, no need to integrate over whole space
    PP.phi_max= mysymmetry & ST_PLANESYM ? M_PI_2 : 2*M_PI;
    int multfactor = (mysymmetry & ST_PLANESYM ? 4 : 1) * (mysymmetry & ST_REFLECTION ? 2 : 1);  // compensates integration of only half- or 1/8-space
    gsl_integration_workspace * ws = gsl_integration_workspace_alloc (1000);
    double interval[2]={-1.0, 1.0};
    gsl_function F;
    F.function=&intBSE_xi;
    F.params=&PP;
    for(unsigned int n=0; n<=Ncoefs_radial; n++)
        for(int l=0; l<=lmax; l+=lstep)
        {
            PP.sh_n=n;  PP.sh_l=l;
            double w=(2*l+1)*Alpha+0.5;
            double Knl = (4*pow_2(n+w)-1)/8/pow_2(Alpha);
            double Inl = Knl * 4*M_PI*Alpha * exp( gsl_sf_lngamma(n+2*w) - 2*gsl_sf_lngamma(w) - gsl_sf_lnfact(n) - 4*w*log(2.0)) / (n+w);
            for(int m=l*mmin; m<=l*mmax; m+=mstep)
            {
                PP.sh_m=m;
                double result, error;
                gsl_integration_qagp(&F, interval, 2, 0, EPSREL_POTENTIAL_INT, 1000, ws, &result, &error);
                SHcoefs[n][l*(l+1)+m] = result * multfactor * (1+(m!=0)) / Inl;
            }
        }
    gsl_integration_workspace_free (ws);
}

void CPotentialBSE::prepareCoefsDiscrete(const CPointMassSetDouble &points)
{
    SHcoefs.resize(Ncoefs_radial+1);
    for(size_t n=0; n<=Ncoefs_radial; n++)
        SHcoefs[n].assign(pow_2(1+Ncoefs_angular), 0);
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    size_t npoints=points.size();
    if(Alpha<0.5)  // autodetect.. the heuristic below is not terribly good, so if you know the inner density slope beforehand, it's better to set alpha=1/(2-gamma) manually
    {
        vectorpd radmass(npoints);
        for(size_t i=0; i<npoints; i++) 
            radmass[i] = paird(sqrt(pow_2(points[i].first.Pos[0])+pow_2(points[i].first.Pos[1])+pow_2(points[i].first.Pos[2])), points[i].second);
        std::sort(radmass.begin(), radmass.end());
        unsigned int n_est=std::min<unsigned int>(std::max<unsigned int>( static_cast<unsigned int>(npoints)/50,1000),static_cast<unsigned int>(npoints)-1); 
        unsigned int n_min=n_est/100;
        double inm=0, r_n_est=radmass[n_est].first;
        while(n_min<n_est && radmass[n_min].first<r_n_est*0.03) n_min++;  // neglect some first entries
        gsl_matrix* X=gsl_matrix_alloc(n_est-n_min, 3);  // matrix of coefficients
        gsl_vector* Y=gsl_vector_alloc(n_est-n_min);     // vector of r.h.s. values (mass(r))
        gsl_vector* W=gsl_vector_alloc(n_est-n_min);     // weights 
        for(unsigned int i=0; i<n_est; i++)
        {
            inm+=radmass[i].second;
            if(i>=n_min){
            gsl_matrix_set(X, i-n_min, 0, 1);  // constant term
            gsl_matrix_set(X, i-n_min, 1, log(radmass[i].first+BH_SMOOTH));  // to avoid log zero values
            gsl_matrix_set(X, i-n_min, 2, radmass[i].first);
            gsl_vector_set(Y, i-n_min, log(inm));
            gsl_vector_set(W, i-n_min, 1/*radmass[i].first/r_n_est*(1-radmass[i].first/r_n_est)*/);  // simple weighting function disfavouring lower and higher radii
            }
        }
        if(n_est-n_min>100)  { // meaningful fit
            // approximate m(r) by power law function:  m(r) = C * r^(3-Gamma) * (1 + k*r), three parameters: C, Gamma, k
            gsl_vector* fit=gsl_vector_alloc(3);
            gsl_matrix* cov=gsl_matrix_alloc(3,3);
            double chisq;
            gsl_multifit_linear_workspace* ws=gsl_multifit_linear_alloc(n_est-n_min, 3);
            gsl_multifit_wlinear(X, W, Y, fit, cov, &chisq, ws);
            gsl_multifit_linear_free(ws);
            double Gamma = 3-gsl_vector_get(fit, 1);   // mass ~ r^(3-gamma)
            Alpha = (Gamma<=1? 1.0+Gamma :
              (Gamma<=1.5? 2.0 :
              (Gamma<1.75? 1./(2-Gamma) : 4.0)));
            gsl_vector_free(fit);
            gsl_matrix_free(cov);
        }  else Alpha=1.;  // default value
        gsl_matrix_free(X);
        gsl_vector_free(Y);
        gsl_vector_free(W);
    }
    double legendre_array[MAX_NCOEFS_ANGULAR][MAX_NCOEFS_ANGULAR-1];
    double gegenpoly_array[MAX_NCOEFS_RADIAL];
    double Inl[MAX_NCOEFS_RADIAL][MAX_NCOEFS_ANGULAR];
    for(int l=0; l<=lmax; l+=lstep)
    {   // pre-compute coefficients
        double w=(2*l+1)*Alpha+0.5;
        for(unsigned n=0; n<=Ncoefs_radial; n++)
            Inl[n][l] = 2*M_SQRTPI*Alpha * exp( gsl_sf_lngamma(n+2*w) - 2*gsl_sf_lngamma(w) - gsl_sf_lnfact(n) - 4*w*log(2.0)) / (n+w) * (4*(n+w)*(n+w)-1)/(8*Alpha*Alpha);
    }
    for(size_t i=0; i<npoints; i++)
    {
        double massi=points[i].second;
        double r=sqrt(pow_2(points[i].first.Pos[0])+pow_2(points[i].first.Pos[1])+pow_2(points[i].first.Pos[2]));
        double costheta= (r==0) ? 0 : (points[i].first.Pos[2]/r);
        double phi=atan2(points[i].first.Pos[1], points[i].first.Pos[0]);
        double ralpha=pow(r, 1/Alpha);
        double xi=(ralpha-1)/(ralpha+1);
        for(int m=0; m<=lmax; m+=mstep)
            gsl_sf_legendre_sphPlm_array(lmax, m, costheta, legendre_array[m]);

        for(int l=0; l<=lmax; l+=lstep)
        {
            double w=(2*l+1)*Alpha+0.5;
            double phil=pow(r, l) * pow(1+ralpha, -(2*l+1)*Alpha);
            gsl_sf_gegenpoly_array(Ncoefs_radial, w, xi, gegenpoly_array);
            for(unsigned n=0; n<=Ncoefs_radial; n++)
                {
                    //double Inl = 2*M_SQRTPI*Alpha * exp( gsl_sf_lngamma(n+2*w) - 2*gsl_sf_lngamma(w) - gsl_sf_lnfact(n) - 4*w*log(2.0)) / (n+w) * (4*(n+w)*(n+w)-1)/(8*Alpha*Alpha);
                    double mult= massi * gegenpoly_array[n] * phil / Inl[n][l];
                    for(int m=0; m<=l*mmax; m+=mstep)
                        SHcoefs[n][l*(l+1)+m] += mult * legendre_array[m][l-m] * cos(m*phi) * (1+(m!=0));
                    if(mmin)
                        for(int m=mmin*l; m<0; m+=mstep)
                            SHcoefs[n][l*(l+1)+m] += mult * legendre_array[-m][l+m] * sin(-m*phi) * 2;
                }
        }
    }
}

double CPotentialBSE::Mass(const double r) const
{
    if(r<=0) return 0;
    double ralpha=pow(r, 1/Alpha);
    double xi=(ralpha-1)/(ralpha+1);
    double gegenpoly_array[MAX_NCOEFS_RADIAL];
    gsl_sf_gegenpoly_array(Ncoefs_radial, Alpha+0.5, xi, gegenpoly_array);
    double multr = pow(1+ralpha, -Alpha);
    double multdr= -ralpha/((ralpha+1)*r);
    if(multdr!=multdr) multdr=0;  // safety measure to avoid NaN
    double result=0;
    for(int n=0; n<=static_cast<int>(Ncoefs_radial); n++)
    {
        double dGdr=(n>0 ? (-n*xi*gegenpoly_array[n] + (n+2*Alpha)*gegenpoly_array[n-1])/(2*Alpha*r) : 0);
        result += SHcoefs[n][0] * multr * (multdr * gegenpoly_array[n] + dGdr);
    }
    return -result * r*r;   // d Phi(r)/d r = - G M(r) /r^2
}

void CPotentialBSE::computeSHCoefs(const double r, double coefsF[], double coefsdFdr[], double coefsd2Fdr2[]) const
{
    double ralpha=pow(r, 1/Alpha);
    double xi=(ralpha-1)/(ralpha+1);
    double gegenpoly_array[MAX_NCOEFS_RADIAL];
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    if(coefsF)      for(unsigned int k=0; k<pow_2(Ncoefs_angular+1); k++) coefsF     [k] = 0;
    if(coefsdFdr)   for(unsigned int k=0; k<pow_2(Ncoefs_angular+1); k++) coefsdFdr  [k] = 0;
    if(coefsd2Fdr2) for(unsigned int k=0; k<pow_2(Ncoefs_angular+1); k++) coefsd2Fdr2[k] = 0;
    for(int l=0; l<=lmax; l+=lstep)
    {
        double w=(2*l+1)*Alpha+0.5;
        gsl_sf_gegenpoly_array(Ncoefs_radial, w, xi, gegenpoly_array);
        double multr = pow(r, l) * pow(1+ralpha, -(2*l+1)*Alpha) * 2*M_SQRTPI;
        double multdr= (l-(l+1)*ralpha)/((ralpha+1)*r);
        if(multdr!=multdr) multdr=0;  // safety measure to avoid NaN
        for(int n=0; n<=(int)Ncoefs_radial; n++)
        {
            double multdFdr=0, multd2Fdr2=0, dGdr=0;
            if(coefsdFdr!=NULL) {
                dGdr=(n>0 ? (-n*xi*gegenpoly_array[n] + (n+2*w-1)*gegenpoly_array[n-1])/(2*Alpha*r) : 0);
                multdFdr= multdr * gegenpoly_array[n] + dGdr;
                if(coefsd2Fdr2!=NULL)
                    multd2Fdr2 = ((l+1)*(l+2)*pow_2(ralpha) + ((1-2*l*(l+1)) - (2*n+1)*(2*l+1)/Alpha - n*(n+1)/pow_2(Alpha))*ralpha + l*(l-1))/pow_2((1+ralpha)*r)*gegenpoly_array[n] - dGdr*2/r;
            }
            for(int m=l*mmin; m<=l*mmax; m+=mstep)
            {
                int indx=l*(l+1)+m;
                double coef = SHcoefs[n][indx] * multr;
                if(coefsF)      coefsF     [indx] += coef * gegenpoly_array[n];
                if(coefsdFdr)   coefsdFdr  [indx] += coef * multdFdr;
                if(coefsd2Fdr2) coefsd2Fdr2[indx] += coef * multd2Fdr2;
            }
        }
    }
}

//----------------------------------------------------------------------------//
// Scale-free potential 

double CPotentialScaleFree::Rho(double X, double Y, double Z, double /*t*/) const
{
    return pow(sqrt(pow_2(X) + pow_2(Y/q) + pow_2(Z/p)), -Gamma);
}

struct CPotentialParamSF2{
    double theta, phi;
    double p, q, Gamma;
};
struct CPotentialParamSF3{
    const CPotentialScaleFree* P;
    double X, Y, Z;
    double p, q, Gamma;
};

double intPhiSF_angular(double s, void* params)
{
    double theta= ((CPotentialParamSF2*)params)->theta;
    double phi  = ((CPotentialParamSF2*)params)->phi;
    double p= ((CPotentialParamSF2*)params)->p;
    double q= ((CPotentialParamSF2*)params)->q;
    double Gamma= ((CPotentialParamSF2*)params)->Gamma;
    return 2*pow(pow_2(s)*( pow_2(sin(theta))*(pow_2(cos(phi)) + pow_2(sin(phi))/(1-(1-pow_2(q))*pow_2(s))) + pow_2(cos(theta))/(1-(1-pow_2(p))*pow_2(s)) ), 1-Gamma/2) /
        sqrt( (1-(1-pow_2(q))*pow_2(s)) * (1-(1-pow_2(p))*pow_2(s)));
}

double intForceSFx(double s, void* params)
{
    double X= ((CPotentialParamSF3*)params)->X;
    double Y= ((CPotentialParamSF3*)params)->Y;
    double Z= ((CPotentialParamSF3*)params)->Z;
    double p= ((CPotentialParamSF3*)params)->p;
    double q= ((CPotentialParamSF3*)params)->q;
    double Gamma= ((CPotentialParamSF3*)params)->Gamma;
    return pow(s, 2-Gamma) * pow( pow_2(X) + pow_2(Y)/(1-(1-pow_2(q))*pow_2(s)) + pow_2(Z)/(1-(1-pow_2(p))*pow_2(s)), -Gamma/2) /
        sqrt( (1-(1-pow_2(q))*pow_2(s)) * (1-(1-pow_2(p))*pow_2(s)));
}
double intForceSFy(double s, void* params)
{
    double X= ((CPotentialParamSF3*)params)->X;
    double Y= ((CPotentialParamSF3*)params)->Y;
    double Z= ((CPotentialParamSF3*)params)->Z;
    double p= ((CPotentialParamSF3*)params)->p;
    double q= ((CPotentialParamSF3*)params)->q;
    double Gamma= ((CPotentialParamSF3*)params)->Gamma;
    return pow(s, 2-Gamma) * pow( pow_2(X) + pow_2(Y)/(1-(1-pow_2(q))*pow_2(s)) + pow_2(Z)/(1-(1-pow_2(p))*pow_2(s)), -Gamma/2) /
        sqrt( (1-(1-pow_2(q))*pow_2(s)) * (1-(1-pow_2(p))*pow_2(s))) / (1-(1-pow_2(q))*pow_2(s));
}
double intForceSFz(double s, void* params)
{
    double X= ((CPotentialParamSF3*)params)->X;
    double Y= ((CPotentialParamSF3*)params)->Y;
    double Z= ((CPotentialParamSF3*)params)->Z;
    double p= ((CPotentialParamSF3*)params)->p;
    double q= ((CPotentialParamSF3*)params)->q;
    double Gamma= ((CPotentialParamSF3*)params)->Gamma;
    return pow(s, 2-Gamma) * pow( pow_2(X) + pow_2(Y)/(1-(1-pow_2(q))*pow_2(s)) + pow_2(Z)/(1-(1-pow_2(p))*pow_2(s)), -Gamma/2) /
        sqrt( (1-(1-pow_2(q))*pow_2(s)) * (1-(1-pow_2(p))*pow_2(s))) / (1-(1-pow_2(p))*pow_2(s));
}

double CPotentialScaleFree::Phi(double X, double Y, double Z, double /*t*/) const
{
    double result, error;
    size_t neval;
    CPotentialParamSF2 PP;
    PP.theta=atan2(sqrt(pow_2(X)+pow_2(Y)), Z); 
    PP.phi=atan2(Y, X);
    PP.p=p; PP.q=q; PP.Gamma=Gamma;
    gsl_function F;
    F.function=&intPhiSF_angular;
    F.params=&PP;
    gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
    return result * 2*M_PI/(2-Gamma)*q*p * pow(pow_2(X)+pow_2(Y)+pow_2(Z),1-Gamma/2) - ((Mbh>0)?(Mbh/sqrt(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z))):0);
}

void CPotentialScaleFree::DiffEq(unsigned n, double /*t*/, double *v, double *f) const
{
    bool calcLyapunov=n>2*N_DIM;
    double X=v[0];
    double Y=v[1];
    double Z=(N_dim==3)?v[2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
    f[0] = v[N_dim];
    f[1] = v[N_dim+1];
    if(N_dim==3)
        f[2] = v[N_dim+2];
    double mult=-4*M_PI*q*p;
    double r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
    double result, error; 
    size_t neval;
    CPotentialParamSF3 PP;
    PP.X=X; PP.Y=Y; PP.Z=Z;
    PP.p=p; PP.q=q; PP.Gamma=Gamma;
    gsl_function F;
    F.params=&PP;
    F.function=&intForceSFx;
    gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
    f[N_dim]   = mult * X * result - ((Mbh>0)?(Mbh*X/r3):0);
    F.function=&intForceSFy;
    gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
    f[N_dim+1] = mult * Y * result - ((Mbh>0)?(Mbh*Y/r3):0);
    if((N_dim==3))
    {
        F.function=&intForceSFz;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
        f[N_dim+2] = mult * Z * result - ((Mbh>0)?(Mbh*Z/r3):0);
    }
    if(calcLyapunov)
    {
#ifdef LYAPUNOVVAREQ
        // variational equation - not implemented in scale-free model
        f[N_dim*2]   = 0;
        f[N_dim*2+1] = 0;
        if(N_dim==3) f[N_dim*2+2]=0;
        f[N_dim*3]   = 0;
        f[N_dim*3+1] = 0;
        if(N_dim==3) f[N_dim*3+2]=0;
#else
        // integrate nearby orbit
        X=v[N_dim*2];
        Y=v[N_dim*2+1];
        Z=(N_dim==3)?v[N_dim*2+2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2] = v[N_dim*3+2];
        r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
        PP.X=X; PP.Y=Y; PP.Z=Z;
        F.function=&intForceSFx;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
        f[N_dim*3]   = mult * X * result - ((Mbh>0)?(Mbh*X/r3):0);
        F.function=&intForceSFy;
        gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
        f[N_dim*3+1] = mult * Y * result - ((Mbh>0)?(Mbh*Y/r3):0);
        if((N_dim==3))
        {
            F.function=&intForceSFz;
            gsl_integration_qng (&F, 0, 1, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval); 
            f[N_dim*3+2] = mult * Z * result - ((Mbh>0)?(Mbh*Z/r3):0);
        }
#endif
    }
}

//----------------------------------------------------------------------------//
// scale-free potential - spherical harmonic expansion
CPotentialScaleFreeSH::CPotentialScaleFreeSH(unsigned int _N_dim, double q, double p, double _Mbh, double _Gamma, unsigned int _Ncoefs_angular): 
CPotentialSH(_N_dim, _Mbh, _Ncoefs_angular), Gamma(_Gamma)
{
    prepareCoefs(q, p);
    mysymmetry = (q==1 ? (p==1 ? ST_SPHERICAL : ST_AXISYMMETRIC) : ST_TRIAXIAL);
}

CPotentialScaleFreeSH::CPotentialScaleFreeSH(unsigned int _N_dim, double _Mbh, double _Gamma, const vectord &_coefs): 
CPotentialSH(_N_dim, _Mbh, static_cast<unsigned int>(sqrt(_coefs.size()*1.0)-1)), Gamma(_Gamma)
{
    assert(_coefs.size()==pow_2(Ncoefs_angular+1));
    SHcoefs=_coefs;  // no array size check performed!! 
    mysymmetry=ST_TRIAXIAL;  // no check performed if it is the case
}

struct CPotentialParamSFSH{
    int sh_l, sh_m;
    double costheta;
    double p, q, Gamma;
};
double intSF_phi(double phi, void* params)
{
    double costheta= ((CPotentialParamSFSH*)params)->costheta;
    int sh_l=((CPotentialParamSFSH*)params)->sh_l;
    int sh_m=((CPotentialParamSFSH*)params)->sh_m;
    double q=((CPotentialParamSFSH*)params)->q;
    double p=((CPotentialParamSFSH*)params)->p;
    double Gamma=((CPotentialParamSFSH*)params)->Gamma;
    return pow((1-pow_2(costheta))*(1+(1/pow_2(q)-1)*pow_2(sin(phi))) + pow_2(costheta/p), -Gamma/2)
        * gsl_sf_legendre_sphPlm(sh_l*2, sh_m*2, costheta) * cos(sh_m*2*phi);
}
double intSF_theta(double theta, void* params)
{
    double result, error;
    ((CPotentialParamSFSH*)params)->costheta=cos(theta);
    gsl_function F;
    F.function=&intSF_phi;
    F.params=params;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag(&F, 0, M_PI_2, 0, EPSREL_POTENTIAL_INT, 1000, 3, w, &result, &error);
    gsl_integration_workspace_free (w);
    return sin(theta)*result;
}
void CPotentialScaleFreeSH::prepareCoefs(double q, double p)
{   // this function assumes that potential in question has triaxial symmetry, unlike the more general function for BSE potential
    SHcoefs.assign(pow_2(Ncoefs_angular+1), 0);
    CPotentialParamSFSH PP;
    PP.q=q; PP.p=p;
    PP.Gamma=std::max<double>(Gamma,1e-5);  // HACK to avoid singularity in gamma=0 coefficient computation
    //gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    F.function=&intSF_theta;
    F.params=&PP;
    for(PP.sh_l=0; PP.sh_l<=static_cast<int>(Ncoefs_angular)/2; PP.sh_l++)
        for(PP.sh_m=0; PP.sh_m<=PP.sh_l; PP.sh_m++)
        {
            double result, error;
            size_t neval;
            gsl_integration_qng(&F, 0, M_PI_2, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval);
            int l=2*PP.sh_l, m=2*PP.sh_m;
            SHcoefs[l*(l+1)+m]= -8*result * (1.0/(l+3-Gamma) + 1.0/(l-2+Gamma)) * 1.0/(2*l+1)  * (1+(m!=0)) * 4*M_PI;
        }
    //gsl_integration_workspace_free (w);
}


void CPotentialScaleFreeSH::computeSHCoefs(const double r, double coefsF[], double coefsdFdr[], double coefsd2Fdr2[]) const
{
    double multr=pow(r, 2-Gamma);
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    for(int l=0; l<=lmax; l+=lstep)
    {
        for(int m=l*mmin; m<=l*mmax; m+=mstep)
        {
            int indx=l*(l+1)+m;
            double coef = -SHcoefs[indx] * multr;
            if(coefsF)      coefsF     [indx] = coef;
            if(coefsdFdr)   coefsdFdr  [indx] = coef * (2-Gamma)/r;
            if(coefsd2Fdr2) coefsd2Fdr2[indx] = coef * (2-Gamma)*(1-Gamma)/(r*r);
        }
    }
}

//----------------------------------------------------------------------------//
// simple harmonic oscillator potential

double CPotentialHarmonic::Rho(double /*X*/, double /*Y*/, double /*Z*/, double /*t*/) const
{
    return 2*(1+1/q/q+1/p/p);
}

double CPotentialHarmonic::Phi(double X, double Y, double Z, double /*t*/) const
{
    return X*X + Y*Y/q/q + Z*Z/p/p - Mbh/sqrt(pow_2(BH_SMOOTH)+ pow_2(X) + pow_2(Y) + pow_2(Z)) - ((Omega>0)?pow_2(Omega)*(pow_2(X)+pow_2(Y))/2:0);
}

void CPotentialHarmonic::DiffEq(unsigned n, double /*t*/, double *v, double *f) const
{
    bool calcLyapunov=n>2*N_DIM;
    double X=v[0];
    double Y=v[1];
    double Z=(N_dim==3)?v[2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
    f[0] = v[N_dim];
    f[1] = v[N_dim+1];
    if(N_dim==3)
        f[2] = v[N_dim+2];
    double gradr= Mbh/pow( pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);

    f[N_dim]   = -X*(gradr + 2)     + ((Omega>0)?  2*Omega*v[N_dim+1]+pow_2(Omega)*X : 0);
    f[N_dim+1] = -Y*(gradr + 2/q/q) + ((Omega>0)? -2*Omega*v[N_dim]  +pow_2(Omega)*Y : 0);
    if((N_dim==3))
        f[N_dim+2] = -Z*(gradr + 2/p/p);
    if(calcLyapunov)
    {
#ifdef LYAPUNOVVAREQ
        // variational equation
        double Mr3=Mbh/pow(pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
        double Mr5=Mbh/pow(pow_2(X) + pow_2(Y) + pow_2(Z), 2.5);
        double wx = v[N_dim*2];     // variation vector
        double wy = v[N_dim*2+1];
        double wz = (N_dim==3)? v[N_dim*2+2] : 0;
        f[N_dim*2]   = v[N_dim*3];
        f[N_dim*2+1] = v[N_dim*3+1];
        if(N_dim==3)
            f[N_dim*2+2]=v[N_dim*3+2];
        f[N_dim*3]   = -wx * (2     - 3*Mr5*pow_2(X) + Mr3) + wy*X*Y*3*Mr5 + wz*X*Z*3*Mr5;
        f[N_dim*3+1] = -wy * (2/q/q - 3*Mr5*pow_2(Y) + Mr3) + wx*Y*X*3*Mr5 + wz*Y*Z*3*Mr5;
        if(N_dim==3)
            f[N_dim*3+2] = -wz * (2/p/p - 3*Mr5*pow_2(Z) + Mr3) - wx*Z*X*3*Mr5 + wy*Z*Y*3*Mr5;
#else
        // integrate nearby orbit
        X=v[N_dim*2];
        Y=v[N_dim*2+1];
        Z=(N_dim==3)?v[N_dim*2+2]:0;    // for universality: if only 2 dimensions, Z=0 will add nothing
        f[N_dim*2]   = v[N_dim*3]   + ((Omega>0)?  2*Omega*v[N_dim*3+1]+pow_2(Omega)*X : 0);
        f[N_dim*2+1] = v[N_dim*3+1] + ((Omega>0)? -2*Omega*v[N_dim*3]  +pow_2(Omega)*Y : 0);
        if(N_dim==3)
            f[N_dim*2+2] = v[N_dim*3+2];
        gradr= Mbh/pow( pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
        f[N_dim*3]   = -X*(gradr + 2);
        f[N_dim*3+1] = -Y*(gradr + 2/q/q);
        if((N_dim==3))
            f[N_dim*3+2] = -Z*(gradr + 2/p/p);
#endif
    }
}

//----------------------------------------------------------------------------//
// Frozen-N-body potential (adapted from Barnes&Hut treecode from 'hackforce')
CPotentialNB::CPotentialNB(unsigned int _N_dim, double _Mbh) :
    CPotential(_N_dim, _Mbh), eps(0), tol(0)
{ 
    maketree();
}

template<typename NumT> CPotentialNB::CPotentialNB(
unsigned int _N_dim, double _Mbh, double _eps, double _tol, const std::vector< std::pair< CPosVelPoint<NumT>, NumT> > &points):
    CPotential(_N_dim, _Mbh), eps(_eps), tol(_tol)
{
    if(N_dim==N_DIM)
    {
        bodytab.resize(points.size());
        for(size_t i=0; i<points.size(); i++)
        {
            bodytab[i].type=BODY;
            bodytab[i].mass=points[i].second;
            bodytab[i].pos[0]=points[i].first.Pos[0];
            bodytab[i].pos[1]=points[i].first.Pos[1];
            bodytab[i].pos[2]=points[i].first.Pos[2];
            bodytab[i].eps=eps;
        }
    }
    maketree();  // should be called anyway
}
// explicit instantiations of constructor for two variants of template argument
template CPotentialNB::CPotentialNB(unsigned int _N_dim, double _Mbh, double _eps, double _tol, const std::vector< std::pair< CPosVelPoint<float >, float > > &points);
template CPotentialNB::CPotentialNB(unsigned int _N_dim, double _Mbh, double _eps, double _tol, const std::vector< std::pair< CPosVelPoint<double>, double> > &points);

double CPotentialNB::Phi(double X, double Y, double Z, double /*t*/) const
{
    double phi=0;
    vec3 point={X,Y,Z};
    walktree(troot, pow_2(rsize), point, CALCPHI, &phi);
    return -phi - ((Mbh>0)?(Mbh/sqrt(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z))):0) - ((Omega>0)?pow_2(Omega)*(pow_2(X)+pow_2(Y))/2:0);
}

bool comparepd(paird elem1, paird elem2)
{
    return(elem1.first<elem2.first);
}

double CPotentialNB::Rho(double X, double Y, double Z, double /*t*/) const
{   
    vectorpd vec;
    vec3 point={X,Y,Z};
    // first determine which cell does this point belong to, and its size (which will be used as an estimate of neighbour search radius)
    double searchrad=rsize;
    nodeptr *qptr=(nodeptr*)&troot;
    int xp[N_DIM];
    if(!intcoord(xp, point)) return 0;   // if point is outside root cell, density is zero
    int l = IMAX >> 1;                /* start with top bit       */
    while (*qptr != NULL && l>0) {  /* loop descending tree     */
        if ((*qptr)->type == BODY) {/*   reached a "leaf"?      */
            l=0;
        }
        else{
            qptr = &(((cellptr)(*qptr))->subp[subindex(xp, l)]);    /*   move down one level    */
            l = l >> 1;                /*   and test next bit      */
            searchrad/=2;
        }
    }
    searchrad *= 4;  // to have somewhat more than one neighbour!
    // now look for neighbours inside searchrad
    size_t numneighbour=0;
    do{
        vec.clear();
        walktreedens(troot, rsize*rsize, searchrad*searchrad, point, xp, IMAX >> 1, &vec);
        numneighbour=vec.size();
        searchrad*=2;
    } while(numneighbour<TREECODE_DENSITY_KTH_NEIGHBOUR);
    std::sort(vec.begin(), vec.end(), comparepd);
    unsigned int kth=std::min<unsigned int>(static_cast<unsigned int>(vec.size()), TREECODE_DENSITY_KTH_NEIGHBOUR);
    if(kth<3) return 0; // failure to get density
    double mass=0;
    double smoothrad=vec[kth-1].first;
    for(size_t k=0; k<vec.size(); k++)
    {
        double radrel=vec[k].first/smoothrad;  // smoothing radius is the distance to Kth neighbour
        double weight=0;  
#ifdef TREECODE_DENSITY_SPLINE
        if(radrel<0.5) weight=8-48*radrel*radrel+48*radrel*radrel*radrel;  // use spline weighting function with compact support
        else if(radrel<1) weight=16*(1-radrel)*(1-radrel)*(1-radrel);
#else
        if(radrel<1) weight=15.0/8*(1-radrel*radrel);  // use simpler Ferrer F1 kernel
#endif
        mass += vec[k].second*weight;
        if(weight==0) k=vec.size(); // all other neighbours are farther than spline support
    }
    mass /= (M_PI*pow(smoothrad, 3.0));  // for spline
    return mass;
}

void CPotentialNB::DiffEq(unsigned n, double /*t*/, double *v, double *f) const
{
    double X=v[0];
    double Y=v[1];
    double Z=v[2];
    f[0] = v[3];
    f[1] = v[4];
    f[2] = v[5];
    f[3]=f[4]=f[5]=0;
    f[N_DIM*2]=rsize*rsize;  // init min distance (squared)
    if(n<2*N_DIM) return;   // do not live in Flatland...
    vec3 point={X,Y,Z};
    walktree(troot, pow_2(rsize), point, CALCACC, f+N_DIM);
    f[2*N_DIM]=sqrt(f[2*N_DIM]);  // this holds minimum distance (unsoftened) to nearest particle (un-squared:)
    if(Mbh>0 || Omega>0)
    {
        double r3=pow(pow_2(BH_SMOOTH) + pow_2(X) + pow_2(Y) + pow_2(Z), 1.5);
        f[N_DIM]   += - ((Mbh>0)?(Mbh*X/r3):0) - ((Omega>0)?  2*Omega*v[N_dim+1]+pow_2(Omega)*X : 0) ;
        f[N_DIM+1] += - ((Mbh>0)?(Mbh*Y/r3):0) - ((Omega>0)? -2*Omega*v[N_dim]  +pow_2(Omega)*Y : 0) ;
        f[N_DIM+2] += - ((Mbh>0)?(Mbh*Z/r3):0);
    }
}

double CPotentialNB::Mass(const double r) const
{
    double m=0, r2=r*r;
    for(size_t i=0; i<bodytab.size(); i++)
        if(pow_2(bodytab[i].pos[0])+pow_2(bodytab[i].pos[1])+pow_2(bodytab[i].pos[2]) < r2)
            m+=bodytab[i].mass;
    return m;
}

void CPotentialNB::maketree()    // initialize tree from set of bodies
{
    int celltabsize = (int) (100+bodytab.size()*0.52);  // should be enough
    bool result=true; 
    do{
        troot=NULL;
        for(unsigned int d=0; d<N_DIM; d++) rmin[d]=-2.0;
        rsize=4.0;  // init box size
        celltab.clear();
        celltab.reserve(celltabsize); 
        result=true; 
        size_t b=0;
        for(; b<bodytab.size() && result; b++)
        {
            result = expandbox(b);   /* expand root to fit */
            result&= appendtree(b);  /* insert into tree   */
        }
        if(!result)
        {
            if(celltab.size()>=celltab.capacity() && b>bodytab.size()/2)  // error due to insufficient allocation of celltab
            {
                celltabsize = (int) (celltab.capacity() *1.05* bodytab.size()/b); // increase it and repeat
                my_error("Warning, restarting tree initialization because of insufficient initial allocation of celltab");
            }
            else
            {
                troot=NULL;
                celltab.clear();
                my_error("Error (unknown) in tree construction, perhaps pathological mass distribution");
                return;  // unknown error...
            }
        }
    } while(!result);

    centerofmass(troot, rmin, rsize,0);                /* find c-of-m coordinates  */

    if(eps<0)  // need to initialize each particle's smoothing radius according to local density
        for(int iter=0; iter<10; iter++)   // iterations to reduce 'binning' effect when eps ~ size of the smallest cell containing this one particle
            assigneps(troot, 0, rsize);
    /*{    // this is more exact method, based on local density, but much more costly
        for(int i=0; i<bodytab.size(); i++)
        {
            double dens=Rho(bodytab[i].pos[0], bodytab[i].pos[1], bodytab[i].pos[2]);
            bodytab[i].eps = fabs(eps)*pow(dens/bodytab[i].mass, -1.0/3);
        }
    }*/
    propagateeps(troot);   // initialize cells' epsilon by averaging from underlying nodes
    
}

bool CPotentialNB::expandbox(const size_t b)  // expand bounding box if needed to accomodate for a new body
{
    int xtmp[N_DIM];
    while (! intcoord(xtmp, bodytab[b].pos)) {  /* expand box (rarely)      */
       vec3 rmid;
       for(unsigned int d=0; d<N_DIM; d++)
            rmid[d]=rmin[d]+0.5*rsize;          /*   find box midpoint      */
        for(unsigned int d=0; d<N_DIM; d++)              /*   loop over dimensions   */
            if (bodytab[b].pos[d] < rmid[d])    /*     is p left of mid?    */
                rmin[d] -= rsize;               /*       extend to left     */
        rsize = 2.0 * rsize;                    /*   double length of box   */
        if (troot != NULL) {                    /*   repot existing tree?   */
            cell c;                             /*     create new root cell */
            int xmid[N_DIM];
            if(!intcoord(xmid, rmid)) return false;   /*     locate old root cell */
            int s = subindex(xmid, IMAX >> 1);  /*     find old tree index  */
            c.subp[s] = troot;                  /*     graft old on new     */
            if(celltab.size()>=celltab.capacity()) return false;  /* cannot reallocate vector, otherwise all pointers will be confused */
            celltab.push_back(c);
            troot = (nodeptr) (&(celltab.back()));             /*     plant new tree       */
        }
    }
    return true;
}

bool CPotentialNB::appendtree(const size_t b)            // append a body to the tree, creating new cell if necessary
{
    int xp[N_DIM], xq[N_DIM];
    nodeptr *qptr;

    intcoord(xp, bodytab[b].pos);   /* form integer coords      */
    int l = IMAX >> 1;                /* start with top bit       */
    qptr = &troot;                  /* start with tree root     */
    while (*qptr != NULL) {         /* loop descending tree     */
        if(l==0) {                  /*   dont run out of bits   */
            (*qptr)->mass += bodytab[b].mass;  /* if two particles coincide, just sum up their mass */
            return true;
        }
        if ((*qptr)->type == BODY) {/*   reached a "leaf"?      */
            cell c;                 /*     alloc a new cell     */
            if(!intcoord(xq, (*qptr)->pos)) return false;   /*     get integer coords   */
            c.subp[subindex(xq, l)] = *qptr;    /*     put body in cell     */
            if(celltab.size()>=celltab.capacity()) return false;  /* cannot reallocate vector, otherwise all pointers will be confused */
            celltab.push_back(c);
            *qptr = (nodeptr) (&(celltab.back()));        /*     link cell in tree    */
        }
        qptr = &(((cellptr)(*qptr))->subp[subindex(xp, l)]);    /*   move down one level    */
        l = l >> 1;                /*   and test next bit      */
    }
    *qptr = (nodeptr) (&(bodytab[b]));            /* found place, store p     */
    return true;
}

bool CPotentialNB::intcoord(    // compute integer coords, return false if point is out of bounds for current root cell
            int xp[N_DIM],   /* integerized coordinate vec3 [0,IMAX) */
            const vec3 rp)   const   /* double coordinate vec3 (system coords) */
{
    bool inb = true;                    /* use to check bounds      */
    for (unsigned int k = 0; k < N_DIM; k++) {        /* loop over dimensions     */
        double xsc = (rp[k] - rmin[k]) / rsize; /*   scale to range [0,1)   */
        if (0.0 <= xsc && xsc < 1.0)            /*   within unit interval?  */
            xp[k] = (int) floor(IMAX * xsc);    /*     then integerize      */
        else                                    /*   out of range           */
            inb = false;                        /*     then remember that   */
    }
    return inb;
}

int CPotentialNB::subindex(   // determine the subcell that a given position belongs to
           const int x[N_DIM],        /* integerized coordinates of particle */
           const int l)     const      /* current level of tree */
{
    int i = 0;                                  /* sum index in i           */
    for (unsigned int k = 0; k < N_DIM; k++)             /* check each dimension     */
        if (x[k] & l)                           /*   if beyond midpoint     */
            i += NSUB >> (k + 1);               /*     skip over subcells   */
    return (i);
}

void CPotentialNB::centerofmass(const nodeptr q, const vec3 corner, const double cellsize, int l)      // recursive center-of-mass finder, also initializes rmax and quadrupole moment
{
    if(q==NULL) return;
    if(q->type == CELL) {                      /* is this a cell?          */
        q->mass = 0.0;                          /*   init total mass        */
        for(unsigned int d=0; d<N_DIM; d++) q->pos[d]=0; /*   and c. of m.           */
        for(unsigned int i=0; i< NSUB; i++) {            /*   loop over subcells     */
            vec3 subcellcorner;
            nodeptr r = ((cellptr)q)->subp[i];
            if(r != NULL) {                    /*     does subcell exist?  */
                if(r->type==BODY) 
                    ((cellptr)q)->numbody++;   
                else
                {
                    double halfsize=cellsize/2;
                    for(unsigned int d=0; d<N_DIM; d++)
                        subcellcorner[d]=corner[d] + halfsize*((i & (NSUB >> (d+1))) > 0);
                    centerofmass(r, subcellcorner, halfsize, l+1);  /*       find subcell cm    */
                    ((cellptr)q)->numbody+=((cellptr)r)->numbody;   /*  count total number of bodies */
                }
                q->mass += r->mass;             /*       sum total mass     */
                for(unsigned int d=0; d<N_DIM; d++)      /*   sum c.o.m. coordinates */
                    q->pos[d] += r->mass * r->pos[d];
            }
        }
        if(q->mass>0)
            for(unsigned int d=0; d<N_DIM; d++)
                q->pos[d] /= q->mass;
        // now we have center of mass coordinates, may compute quadrupole moment and other stuff
#ifdef TREECODE_QUADRUPOLE
        for(unsigned int d1=0; d1<N_DIM; d1++)
            for(unsigned int d2=0; d2<N_DIM; d2++)
                ((cellptr)q)->quad[d1][d2]=0;
        for(unsigned int i = 0; i < NSUB; i++) {        /*   loop over subnodes     */
            nodeptr r = ((cellptr)q)->subp[i];
            if (r != NULL) {                    /*     does subnode exist?  */
                vec3 dr={r->pos[0]-q->pos[0], r->pos[1]-q->pos[1], r->pos[2]-q->pos[2]};  /* displacement vect. */
                for(unsigned int d1=0; d1<N_DIM; d1++)
                    for(unsigned int d2=0; d2<N_DIM; d2++)
                    {
                        ((cellptr)q)->quad[d1][d2] += dr[d1]*dr[d2] * r->mass;
                        if(r->type==CELL)       /*       if subnode is cell */
                            ((cellptr)q)->quad[d1][d2] += ((cellptr)r)->quad[d1][d2];  /*         use its moment   */
                    }
            }
        }
#endif
#ifdef TREECODE_DEHNEN_MAC
        double rmax1=0, rmax2=0;
        // loop over sub-nodes
        for (unsigned int i = 0; i < NSUB; i++) {        /*   loop over subnodes     */
            nodeptr r = ((cellptr)q)->subp[i];
            if (r != NULL) {            /*     does subnode exist?  */
                vec3 dr={r->pos[0]-q->pos[0], r->pos[1]-q->pos[1], r->pos[2]-q->pos[2]};  /* displacement vect. */
                double drmag=sqrt(pow_2(dr[0])+pow_2(dr[1])+pow_2(dr[2]));
                rmax1 = std::max<double>(rmax1, (r->type==BODY ? 0 : ((cellptr)r)->rmax) + drmag );
            }
        }
        // find distance from this cell's c.o.m. to its most distant corner
        for(unsigned int c=0; c<NSUB; c++)  // loop over corners
        {
            double dist2corner=0;
            for(unsigned int d=0; d<N_DIM; d++)
                dist2corner += pow_2(corner[d] + ((c & (NSUB >> (d+1))) > 0)*cellsize - q->pos[d]);
            rmax2 = std::max<double>(rmax2, sqrt(dist2corner));
        }
        // assign cell's rmax to lower of these two values (Dehnen 2000)
        ((cellptr)q)->rmax=std::min<double>(rmax1, rmax2);
#ifdef TREECODE_QUADRUPOLE
        // check that quadrupole terms are not too large, otherwise increase rmax
        //double quadfact = (((cellptr)q)->quad[0][0]+((cellptr)q)->quad[1][1]+((cellptr)q)->quad[2][2])/(q->mass*pow_2(((cellptr)q)->rmax));
        //if(quadfact>0.3)  ((cellptr)q)->rmax *= sqrt(quadfact/0.3);
#endif
#else
        ((cellptr)q)->rmax=cellsize;
#endif
    }
}

void CPotentialNB::assigneps(const nodeptr q, const double epsparent, const double cellsize)      // recursive assign of smoothing length
{
    if(q==NULL) return;
    double epsavg=0; int numsub=0;
    double epscur = q->eps>0 ? q->eps : fabs(eps) * cellsize/(q->type==BODY ? 1.0 : pow(((cellptr)q)->numbody*1.0, 1.0/3) );
    // weigth mix with the parent cell's eps
    if(epsparent!=0)
        epscur = sqrt(epsparent*epscur); //pow( (pow(epsparent, -0.5) + pow(epscur, -0.5))/2, -2.0);
    q->eps=epscur;   // assign some preliminary guess
    if(q->type==BODY) return;  // no descendants
    for(unsigned int i=0; i< NSUB; i++) {
        nodeptr r = ((cellptr)q)->subp[i];
        if (r != NULL) {
            assigneps(r, ((cellptr)q)->eps, cellsize/2);
            epsavg += pow(r->eps, -1./3);//* (r->type==BODY ? 1 : ((cellptr)r)->numbody);
            numsub++;
        }
    }
    if(numsub>0)
        q->eps = pow(epsavg/numsub, -3.);   // update guess
}

void CPotentialNB::propagateeps(const nodeptr q)
{
    if(q==NULL) return;
    double epsavg=0; int numsub=0;
    for(unsigned int i=0; i< NSUB; i++) {    //assign eps to tree cells by averaging the values in sub-cells/leafs
        nodeptr r = ((cellptr)q)->subp[i];
        if (r != NULL) {
            if(r->type==CELL) propagateeps(r);
            epsavg += pow(r->eps, -2.);
            numsub++;
        }
    }
    if(numsub>0)
        q->eps = pow(epsavg/numsub, -0.5);
}

void CPotentialNB::walktree(const nodeptr p,   // pointer to interacting node or cell 
           const double dsq,                   // cell size squared
           const vec3 pos0,                    // position to calculate quantity in question
           const WALKOP operation,             // operation to perform during tree-walk
           double data[])  const               // data accumulated according to the operation
{
    if(p==NULL) return;
    vec3 dr={p->pos[0]-pos0[0], p->pos[1]-pos0[1], p->pos[2]-pos0[2]};
    double distsq=pow_2(dr[0])+pow_2(dr[1])+pow_2(dr[2]);
    // test whether cell should be opened
    if(p->type==CELL &&
//#ifdef TREECODE_DEHNEN_MAC
        tol*tol*distsq < pow_2(((cellptr)p)->rmax)
//#else
//        tol*tol*distsq < dsq
//#endif
        ) {
        for(unsigned int c=0; c<NSUB; c++)
            if(((cellptr)p)->subp[c] != NULL)  /* walk into sub-cells */
                walktree(((cellptr)p)->subp[c], dsq/4.0, pos0, operation, data);
    }
    else if(distsq>0) 
    {   // perform requested operation on leaf or unopened cell
//        if( p->type==CELL && ((cellptr)p)->rmax<sqrt(dsq)*0.3)
  //          distsq+=sin(0.);
        // assign smoothed distance, based either on Plummer softening or on F1 compact kernel;  notation g1-g3 is taken from Springel et al.2001
#ifdef TREECODE_EPS_COMPACT_KERNEL
        double eps2=p->eps*p->eps, eps3=eps2*p->eps;    /* use Epanechnikov kernel with compact support */
        double dist_sq_over_eps = p->eps==0?2 : distsq/eps2;
        double dist_inv_eps = dist_sq_over_eps>=1 ? 1/sqrt(distsq) : (15./8 - 5./4*dist_sq_over_eps + 3./8*pow_2(dist_sq_over_eps) )/p->eps;
        double g1 = dist_sq_over_eps>=1 ? dist_inv_eps/distsq : (2.5-1.5*dist_sq_over_eps)/eps3;
        double g2 = dist_sq_over_eps>=1 ? -3*g1/distsq : -3/(eps2*eps3);
        double g3 = dist_sq_over_eps>=1 ? -5*g2/distsq : 0;
#else
        double dist_sq_eps = distsq + pow_2(p->eps);    /* use Plummer softening */
        double dist_inv_eps = 1/sqrt(dist_sq_eps);
        double g1 = dist_inv_eps/dist_sq_eps;
        double g2 = -3*g1/dist_sq_eps;
        double g3 = -5*g2/dist_sq_eps;
#endif
#ifdef TREECODE_QUADRUPOLE
        vec3 quaddr;
        double drquaddr=0, quadtrace=0;
        if(p->type == CELL){        /* if cell, add quad. term  */
            quadtrace=( ((cellptr)p)->quad[0][0]+((cellptr)p)->quad[1][1]+((cellptr)p)->quad[2][2] );
            for(unsigned int d1=0; d1<N_DIM; d1++)
            {
                quaddr[d1] = 
                    ((cellptr)p)->quad[d1][0] * dr[0] + 
                    ((cellptr)p)->quad[d1][1] * dr[1] + 
                    ((cellptr)p)->quad[d1][2] * dr[2];
                drquaddr += quaddr[d1]*dr[d1];
            }
        }
#endif
        switch(operation){
            case CALCPHI: 
                data[0] += p->mass*dist_inv_eps; 
#ifdef TREECODE_QUADRUPOLE
                data[0] -= 0.5*(g2*drquaddr + g1*quadtrace);
#endif
                break;
            case CALCACC: 
                for(unsigned int d=0; d<N_DIM; d++)
                    data[d] += p->mass*dr[d] * g1;
#ifdef TREECODE_QUADRUPOLE
                if(p->type == CELL) {        /* if cell, add quad. term  */
                    for(unsigned int d=0; d<N_DIM; d++)
                        data[d] += g2*quaddr[d] + 0.5*(g3*drquaddr+g2*quadtrace)*dr[d];
                }
#endif
                if(distsq<data[N_DIM]) data[N_DIM]=distsq;   // update minimum distance
                //data[N_DIM]=std::min<double>(data[N_DIM], 1/pow_2(dist_inv_eps));
                break;
        }
    }
}

void CPotentialNB::walktreedens(const nodeptr p,   // pointer to interacting node or cell 
           const double dsq,      // cell size squared
           const double searchradsq, // search radius squared
           const vec3 pos0,       // position to calculate quantity in question
           const int xp[N_DIM],   // integerized coordinates of the point in question
           const int lev,         // current grid level
           vectorpd* data)  const // array to store neighbours
{
    if(p==NULL) return;
    vec3 dr={p->pos[0]-pos0[0], p->pos[1]-pos0[1], p->pos[2]-pos0[2]};
    double distsq=pow_2(dr[0])+pow_2(dr[1])+pow_2(dr[2]);
    if(p->type==CELL) {  
        if(distsq < dsq*3 + searchradsq) {  // worst case - CoM may lie on the other side of diagonal
            for(unsigned int c=0; c<NSUB; c++)
                if(((cellptr)p)->subp[c] != NULL)  /* walk into sub-cells */
                    walktreedens(((cellptr)p)->subp[c], dsq/4, searchradsq, pos0, xp, lev>>1, data);
        }
        else {  // point is far from our CoM, but lies within one of child cells, which still needs to be opened
            int subcell=subindex(xp, lev);
            if(((cellptr)p)->subp[subcell] != NULL)  /* walk into this sub-cell */
                walktreedens(((cellptr)p)->subp[subcell], dsq/4, searchradsq, pos0, xp, lev>>1, data);
        }
    }
    else if (distsq < searchradsq)   // this is a leaf within search radius
        data->push_back(paird(sqrt(distsq),p->mass));
}

//----------------------------------------------------------------------------//
// Spherical-harmonic expansion of arbitrary potential, radial part is spline interpolated on a grid

CPotentialSpline::CPotentialSpline(unsigned int _N_dim, double _Mbh, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
    const CPointMassSetDouble &points, SYMMETRYTYPE _sym):
    CPotentialSH(_N_dim, _Mbh, _Ncoefs_angular),
    Ncoefs_radial(std::max<unsigned int>(4,_Ncoefs_radial))
{
    mysymmetry=_sym;
    if(points.empty()) { initDefault(); return; }
    prepareCoefsDiscrete(points);
}

CPotentialSpline::CPotentialSpline(unsigned int _N_dim, double _Mbh, 
    const vectord &_SHradii, const std::vector< vectord > &_coefs):
    CPotentialSH(_N_dim, _Mbh, (assert(_coefs.size()>0), static_cast<unsigned int>(sqrt(_coefs[0].size()*1.0)-1))), 
    Ncoefs_radial(std::min<unsigned int>(MAX_NCOEFS_RADIAL, static_cast<unsigned int>(_coefs.size()-1)))
{
    assert(_coefs[0].size()==pow_2(Ncoefs_angular+1));
    initSpline(_SHradii, _coefs);
}

CPotentialSpline::CPotentialSpline(unsigned int _N_dim, double _Mbh, unsigned int _Ncoefs_radial, unsigned int _Ncoefs_angular, 
    const CDensity* density, const vectord *radii):    // init potential from analytic mass model
    CPotentialSH(_N_dim, _Mbh, _Ncoefs_angular), 
    Ncoefs_radial(std::max<unsigned int>(4,_Ncoefs_radial))
{
    mysymmetry=density->symmetry();
    prepareCoefsAnalytic(density, radii);
}

// not a constructor, but a method to compute coefs at given radii
CPotentialSpline::CPotentialSpline(unsigned int _Ncoefs_angular, 
        const CPointMassSetDouble &points, const vectord *radii, std::vector< vectord > *coefsArray):
    CPotentialSH(N_DIM, 0, _Ncoefs_angular),
    Ncoefs_radial(static_cast<unsigned int>(radii->size()-1))
{
    mysymmetry=ST_TRIAXIAL;
    prepareCoefsDiscrete(points, radii, coefsArray);
}

void CPotentialSpline::initDefault()
{
    Ncoefs_radial=4; Ncoefs_angular=0;
    vectord nodes(Ncoefs_radial+1);
    nodes[0]=0; nodes[1]=1; nodes[2]=2; nodes[3]=4; nodes[4]=8;
    std::vector< vectord > values(Ncoefs_radial+1);
    values[0]=vectord(1, 1);
    values[1]=vectord(1, 0.5);
    values[2]=vectord(1, 0.333);
    values[3]=vectord(1, 0.21);
    values[4]=vectord(1, 0.1);
    initSpline(nodes, values);
}

CPotentialSpline::~CPotentialSpline()
{
    for(size_t i=0; i<splines.size(); i++)
            if(splines[i]!=NULL)  gsl_spline_free(splines[i]);
}

CPotential* CPotentialSpline::clone() const 
{   // need to duplicate spline arrays explicitly; everything else is copied automatically
    CPotentialSpline* newpot = new CPotentialSpline(*this); 
    for(size_t i=0; i<splines.size(); i++)
        if(splines[i]!=NULL)
        {
            gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, splines[i]->size);
            gsl_spline_init (spline, splines[i]->x, splines[i]->y, splines[i]->size);
            newpot->splines[i]=spline;
        }
        else newpot->splines[i]=NULL;
    return newpot;
}

#ifdef SPLINE_OLD_METHOD
double intSpline_r(double s, void* params)
{
    double r0=((CPotentialParamSH*)params)->r;
    int l=((CPotentialParamSH*)params)->sh_l;
    double result, error;
    size_t neval;
    ((CPotentialParamSH*)params)->r=s;
    gsl_function F;
    F.function=&intSH_theta;
    F.params=params;
    gsl_integration_qng(&F, 0, ((CPotentialParamSH*)params)->theta_max, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval);
    ((CPotentialParamSH*)params)->r=r0;
    if(r0==0) return (l>0 ? 0 : result*s);
    return result * (s>r0 ? s*pow(s/r0, -l) : s*pow(s/r0, l+1));
}
#else
double intSpline_r(double r, void* params)
{
    if(r==0) return 0;
    int n=((CPotentialParamSH*)params)->sh_n;  // power index for s^n, either l+2 or 1-l
    double result, error;
    size_t neval;
    ((CPotentialParamSH*)params)->r=r;
    gsl_function F;
    F.function=&intSH_theta;
    F.params=params;
    /*if(((CPotentialParamSH*)params)->sh_l>=6)
    {
        double result1;
        gsl_integration_qng(&F, 0, ((CPotentialParamSH*)params)->theta_max/2, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval);
        gsl_integration_qng(&F, ((CPotentialParamSH*)params)->theta_max/2, ((CPotentialParamSH*)params)->theta_max, 0, EPSREL_POTENTIAL_INT, &result1, &error, &neval);
        result+=result1;
    } else*/
        gsl_integration_qng(&F, 0, ((CPotentialParamSH*)params)->theta_max, 0, EPSREL_POTENTIAL_INT, &result, &error, &neval);
    //if(fabs(result)<1e-14) result=0;
    return result * pow(r, n*1.0);
}
#endif
void CPotentialSpline::prepareCoefsAnalytic(const CDensity* density, const vectord *srcradii)
{
    std::vector< vectord > coefsArray(Ncoefs_radial+1);  // SHE coefficients to pass to initspline routine
    vectord radii(Ncoefs_radial+1);  // true radii to pass to initspline routine
    for(size_t i=0; i<=Ncoefs_radial; i++)
        coefsArray[i].assign(pow_2(1+Ncoefs_angular), 0);
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    bool initUserRadii = (srcradii!=NULL && srcradii->size()==Ncoefs_radial+1 && srcradii->front()==0);
    if(srcradii!=NULL && !initUserRadii)  // something went wrong with manually supplied radii, initializing by default
        my_error("Warning, invalid call to constructor of SHGrid Schwarzschild model");
    if(initUserRadii)
        radii= *srcradii;
    else
    {
        // find inner/outermost radius
        double totalMass=density->totalMass();
        double epsout = 0.1/sqrt(pow_2(Ncoefs_radial)+0.01*pow(Ncoefs_radial*1.0,4.0));      // how far should be the outer node (leave out this fraction of mass)
        double epsin = 5.0/pow(Ncoefs_radial*1.0,3.0);                                       // how close can we get to zero, in terms of innermost grid node
        double rout = totalMass<0 ? 10 : density->getRadiusByMass(totalMass*(1-epsout));     // totalMass<0 should not really happen! means that mass is infinite
        double rin = std::min<double>(epsin, density->getRadiusByMass(totalMass*epsin*0.1)); // somewhat arbitrary choice for min/max radii, but probably reasonable
        createNonuniformGrid(radii, Ncoefs_radial+1, rin, rout, true); 
#ifdef DEBUGPRINT
        std::cerr << "Mtotal="<<totalMass << ", rin="<<radii[1]<<", rout="<<rout<<"\n";
#endif
    }
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    CPotentialParamSH PP;
    PP.P=density;
    PP.theta_max= mysymmetry & ST_REFLECTION ? M_PI_2 : M_PI;  // if symmetries exist, no need to integrate over whole space
    PP.phi_max= mysymmetry & ST_PLANESYM ? M_PI_2 : 2*M_PI;
    int multfactor = (mysymmetry & ST_PLANESYM ? 4 : 1) * (mysymmetry & ST_REFLECTION ? 2 : 1);  // compensates integration of only half- or 1/8-space
    gsl_function F;
    F.function=&intSpline_r;
    F.params=&PP;
#ifdef SPLINE_OLD_METHOD
    for(int l=0; l<=lmax; l+=lstep)
    {
        PP.sh_l=l;
        for(int m=l*mmin; m<=l*mmax; m+=mstep)
        {
            PP.sh_m=m;
            for(size_t c=0; c<=Ncoefs_radial; c++)
            {   // compute integrals on subintervals
                PP.r=radii[c];
                double result, error, result1;
                if(c>0 && (radii[c]<0.027 || radii[c]>38.0)) {   // more accurate integration near the ends of interval
                    gsl_integration_qags(&F, 0, radii[c], 1e-20, EPSREL_POTENTIAL_INT, 1000, w, &result1, &error);
                    gsl_integration_qagiu(&F, radii[c], 1e-20, EPSREL_POTENTIAL_INT, 1000, w, &result, &error);
                    result+=result1; 
                } else {   // the reason we give absolute error threshold is that for flat-density-core profiles the high-l,m coefs at small radii are extremely tiny and their exact calculation is impractically slow
                    gsl_integration_qagiu(&F, 0, 1e-20, EPSREL_POTENTIAL_INT, 1000, w, &result, &error);
                }
                coefsArray[c][l*(l+1) + m] = multfactor*result * sqrt(4*M_PI)/(2*l+1) * (1+(m!=0));
            }
        }
    }
#else
    vectord coefsInner, coefsOuter;
    radii.front()=BH_SMOOTH*radii[1];  // to prevent log divergence for gamma=2 potentials
    for(int l=0; l<=lmax; l+=lstep)
    {
        PP.sh_l=l;
        for(int m=l*mmin; m<=l*mmax; m+=mstep)
        {
            PP.sh_m=m;
            // first precompute inner and outer density integrals at each radial grid point, summing contributions from each interval of radial grid
            coefsInner.assign(Ncoefs_radial+1, 0);
            coefsOuter.assign(Ncoefs_radial+1, 0);
            // loop over inner intervals
            double result, error;
            PP.sh_n = l+2;
            for(size_t c=0; c<Ncoefs_radial; c++)
            {   // the reason we give absolute error threshold is that for flat-density-core profiles the high-l,m coefs at small radii are extremely tiny and their exact calculation is impractically slow
                gsl_integration_qags(&F, radii[c], radii[c+1], 1e-20, EPSREL_POTENTIAL_INT, 1000, w, &result, &error);
                coefsInner[c+1] = result + coefsInner[c];
            }
            // loop over outer intervals, starting from infinity backwards
            PP.sh_n = 1-l;
            gsl_integration_qagiu(&F, radii.back(), 1e-20, EPSREL_POTENTIAL_INT, 1000, w, &result, &error);
            coefsOuter.back() = result;
            for(size_t c=Ncoefs_radial; c>static_cast<size_t>(l==0 ? 0 : 1); c--)
            {
                gsl_integration_qags(&F, radii[c-1], radii[c], EPSREL_POTENTIAL_INT*fabs(coefsOuter[c]), EPSREL_POTENTIAL_INT, 1000, w, &result, &error);
                coefsOuter[c-1] = result + coefsOuter[c];
            }
            // now compute the coefs of potential expansion themselves
            for(size_t c=0; c<=Ncoefs_radial; c++)
            {
                coefsArray[c][l*(l+1) + m] = ((c>0 ? coefsInner[c]*pow(radii[c], -l-1.0) : 0) + coefsOuter[c]*pow(radii[c], l*1.0)) *
                    multfactor * sqrt(4*M_PI)/(2*l+1) * (1+(m!=0));
            }
        }
    }
    radii.front()=0;
#endif
    gsl_integration_workspace_free (w);
    initSpline(radii, coefsArray);
}

struct CParticle{
    double r, theta, phi, mass;
    CParticle() {r=theta=phi=mass=0;}
    CParticle(const std::pair< CPosVelPoint<double>, double> &src) {
        r=sqrt(pow_2(src.first.Pos[0])+pow_2(src.first.Pos[1])+pow_2(src.first.Pos[2])); 
        theta=r>0?acos(src.first.Pos[2]/r):0; 
        phi=atan2(src.first.Pos[1],src.first.Pos[0]); 
        mass=src.second; 
    }
    inline bool operator<(CParticle p) const
    { return (r<p.r); }
};

void CPotentialSpline::prepareCoefsDiscrete(const CPointMassSetDouble &srcpoints, const vectord *srcradii, std::vector< vectord > *outputCoefsArray)
{
    try{
    bool initUserRadii = srcradii!=NULL;

    std::vector< vectord > coefsArray(Ncoefs_radial+1);  // SHE coefficients to pass to initspline routine
    for(size_t i=0; i<=Ncoefs_radial; i++)
        coefsArray[i].assign(pow_2(Ncoefs_angular+1), 0);
    double totalMass=0;
    size_t npoints=srcpoints.size();
    std::vector<CParticle> points;
    points.reserve(npoints);
    std::vector<double> SHradii(npoints);
    std::vector< vectord > SHcoefs(npoints);  // coefficients of expansion at each point
    for(size_t n=0; n<SHcoefs.size(); n++)
    {
        SHcoefs[n].assign(pow_2(Ncoefs_angular+1), 0);
    }
    for(size_t i=0; i<npoints; i++)
    {
        if(srcpoints[i].second>0)  // don't consider zero-mass points
            points.push_back(CParticle(srcpoints[i]));
        if(srcpoints[i].second<0)
        {
            my_error("Some negative mass creeped in, unphysical and against world harmony");
            initDefault();
        }
        totalMass+=srcpoints[i].second;
    }
    std::sort(points.begin(), points.end());
    npoints=points.size();
    if(npoints<=Ncoefs_radial*10)
    {
        if(!initUserRadii) my_error("Number of points is too small");
        initDefault();
    }
        
    // having sorted particles in radius, may now initialize coefs
    vectord 
        CoefsOuter(pow_2(Ncoefs_angular+1)), 
        CoefsInner(pow_2(Ncoefs_angular+1)), 
        CoefsOuterMax(pow_2(Ncoefs_angular+1));
    double legendre_array[MAX_NCOEFS_ANGULAR][MAX_NCOEFS_ANGULAR-1];
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    unsigned int indUserRadii=0;
    size_t ind=0, numRestart=0;
    do{ // may require several passes because of round-off error accumulation. Remember the last initialized point in "ind"
        double rpivot=points[ind].r;
        // first compute coefs at r_pivot
        CoefsOuter.assign(pow_2(Ncoefs_angular+1), 0);
        CoefsInner.assign(pow_2(Ncoefs_angular+1), 0);
        CoefsOuterMax.assign(pow_2(Ncoefs_angular+1), 0);
        for(size_t i=0; i<npoints; i++)
        {
            double costheta=cos(points[i].theta);
            for(int m=0; m<=lmax; m+=mstep)
                gsl_sf_legendre_sphPlm_array(lmax, m, costheta, legendre_array[m]);
            for(int l=0; l<=lmax; l+=lstep)
                for(int m=l*mmin; m<=l*mmax; m+=mstep)
                {
                    int absm=abs(m);  // negative m correspond to sine, positive - to cosine
                    double mult= sqrt(4*M_PI)/(2*l+1) * (1+(m!=0)) * points[i].mass *
                        legendre_array[absm][l-absm] * (m>=0 ? cos(m*points[i].phi) : sin(m*points[i].phi));
                    if(points[i].r>=rpivot) {
                        CoefsOuter[l*(l+1)+m] += mult * (points[i].r>0 ? pow(points[i].r, -(1+l)) : 0);   // no particle with nonzero mass at r=0 is allowed (checked above)
                        CoefsOuterMax[l*(l+1)+m]=std::max<double>(CoefsOuterMax[l*(l+1)+m], fabs(CoefsOuter[l*(l+1)+m]));  // update max value
                    }
                    else 
                        CoefsInner[l*(l+1)+m] += mult * pow(points[i].r, l);
                }
        }
        if(initUserRadii && ind==0)  // init coefficients at grid radii that are below 1st point (there are only l=0 coefs)
        {
            while(indUserRadii<=Ncoefs_radial && srcradii->at(indUserRadii)<points[0].r)
            {
                coefsArray[indUserRadii][0] = CoefsOuter[0];
                indUserRadii++;
            }
        }
        // then iterate over points subtracting and adding outer/inner expansion coefs, unless a loss of precision is suspected
        bool ok=true;
        for(; ind<npoints && ok; ind++)
        {
            SHradii[ind]=points[ind].r;
            // if using manually supplied srcradii, find out indices of particles at which these grid nodes lie
            unsigned int countIndUserRadii=0;
            if(initUserRadii)
            {
                while(indUserRadii+countIndUserRadii<=Ncoefs_radial && 
                    (ind==npoints-1 || srcradii->at(indUserRadii+countIndUserRadii)<points[ind+1].r))
                    countIndUserRadii++;     // counting how many user grid points lie between previous and current particle radii
            }
            double costheta=cos(points[ind].theta);
            for(int m=0; m<=lmax; m+=mstep)
                gsl_sf_legendre_sphPlm_array(lmax, m, costheta, legendre_array[m]);
            for(int l=0; l<=lmax; l+=lstep)
                for(int m=l*mmin; m<=l*mmax; m+=mstep)
                {
                    int absm=abs(m);
                    int coefind=l*(l+1)+m;
                    double mult= sqrt(4*M_PI)/(2*l+1) * (1+(m!=0)) * points[ind].mass *
                        legendre_array[absm][l-absm] * (m>=0 ? cos(m*points[ind].phi) : sin(m*points[ind].phi));
                    CoefsOuter[coefind] -= mult * (points[ind].r>0 ? pow(points[ind].r, -(1+l)) : 0);
                    CoefsInner[coefind] += mult * pow(points[ind].r, l);
                    if(fabs(CoefsOuter[coefind]) < 1e-10*CoefsOuterMax[coefind]) 
                        ok=false;  // need to re-init coefs at current radius because of possible loss of precision
                    if(ind==npoints-1) CoefsOuter[coefind]=0;  // exact value beyond all points
                    if(initUserRadii)  // need to compute coefs at a given set of radii
                    {
                        for(unsigned int indGrid=indUserRadii; indGrid<indUserRadii+countIndUserRadii; indGrid++)
                            coefsArray[indGrid][coefind] = CoefsInner[coefind] * pow(srcradii->at(indGrid), -(1+l)) + CoefsOuter[coefind] * pow(srcradii->at(indGrid), l);
                    }
                    else  // compute coefs at every particle radius
                        SHcoefs[ind][coefind] = CoefsInner[coefind] * (points[ind].r>0 ? pow(points[ind].r, -(1+l)) : 0) + CoefsOuter[coefind] * pow(points[ind].r, l);
                }
            indUserRadii+=countIndUserRadii;
        }
        numRestart++;
        if(numRestart>100) {   // bad convergence -- abandon initialization unless an explicit symmetry has been specified which avoids computation of near-zero coefficients
            my_error("Warning, in Spline potential: slow convergence due to roundoff errors. Likely the input Nbody snapshot had some degree of symmetry not explicitly set up in constructor. Initializing by default."); 
            initDefault(); 
            return; 
        }
    } while(ind<npoints);

    if(initUserRadii)
    {
        outputCoefsArray->assign(coefsArray.begin(), coefsArray.end());
        return;  // no spline initialization, potential not usable!
    }

    // now construct basis splines to approximate the computed coefffficients on a much sparser grid
    vectord radii(Ncoefs_radial+1);  // true radii to pass to initspline routine
    // choose the radial grid parameters: innermost cell contains minBinMass and outermost radial node should encompass cutoffMass
    // spline definition region extends up to outerRadius which is ~several times larger than outermost radial node, however coefficients at that radius are not used in the potential computation later
    double minBinMass=totalMass* std::min<double>(10./npoints, 0.1/Ncoefs_radial);
    double cutoffMass=totalMass* (1-0.1/Ncoefs_radial);  // total mass within spline definition region 
    double cutoffRadius=0, innerBinRadius=0, outerRadius=0;
    size_t npointsCutoff=0, npointsInner=0, npointsOuter=0;
    double innerMass=0;
    for(size_t i=0; i<npoints; i++)
    {
        innerMass+=points[i].mass;
        if(innerMass>= minBinMass && innerBinRadius==0) { innerBinRadius=points[i].r; npointsInner=i; }
        if(innerMass>= cutoffMass && cutoffRadius==0)   { cutoffRadius=points[i].r; npointsCutoff=i; }
    }
    createNonuniformGrid(radii, Ncoefs_radial+1, innerBinRadius, cutoffRadius, true);
    // find index of the outermost point which is used at all
    outerRadius = pow_2(points.back().r)/points[points.size()-2].r;  // roughly equally logarithmically spaced from the last two nodes
    for(npointsOuter=npointsCutoff+1; npointsOuter<npoints && points[npointsOuter].r<outerRadius; npointsOuter++);

    size_t numBSplineKnots = Ncoefs_radial+2;  // including zero and outermost point; only interior nodes are actually used for computing best-fit coefs (except l=0 coef, for which r=0 is also used)
    gsl_vector* knots = gsl_vector_alloc(numBSplineKnots);  // log-scaled radii to use in the least-square approximator
    gsl_bspline_workspace* splinewkspace = gsl_bspline_alloc(4, numBSplineKnots);
    /* least-square problem:  
       (x_i, y_i) - original values, i=0..numPointsUsed-1
       B_p(x) - basis functions (b-splines), p=0..numBasisFnc-1
       C_ip - basis matrix, its values = B_p(x_i)
       w_p - weight coefficients to be found
       Y(x) = \sum_p=0^{numBasisFnc-1} w_p B_p(x) - interpolated function
       task: minimize  |y_i-Y(x_i)|^2
    */
    npointsInner=2;  ///!!??
    size_t numPointsUsed = npointsOuter-npointsInner;  // outer and inner points are ignored
    size_t numBasisFnc = Ncoefs_radial+4;  
    gsl_matrix* splinematr = gsl_matrix_alloc(numPointsUsed, numBasisFnc);  // matrix C_ip
    gsl_vector* splinevals = gsl_vector_alloc(numBSplineKnots+2);   // values of B_p(x)
    gsl_vector* originvals = gsl_vector_alloc(numPointsUsed);       // y_i
    gsl_vector* splinecoefs= gsl_vector_alloc(numBasisFnc);         // w_p
    gsl_matrix* covarmatr  = gsl_matrix_alloc(numBasisFnc, numBasisFnc);
    gsl_multifit_linear_workspace* fitwkspace = gsl_multifit_linear_alloc(numPointsUsed, numBasisFnc);
    if(fitwkspace==NULL)
        throw std::bad_alloc();

    double chisq;

    // first construct spline for zero-order term (radial dependence)
    potcenter=SHcoefs[0][0];  // value of potential at origin (times 1/2\sqrt{\pi} )
    for(size_t p=1; p<=Ncoefs_radial; p++)
        gsl_vector_set(knots, p, log(radii[p]));
    gsl_vector_set(knots, 0, 2*log(radii[1])-log(radii[2]));   // inner and outer B-spline grid points
    gsl_vector_set(knots, Ncoefs_radial+1, log(outerRadius));  // are not used in approximation, only points within the interior b-spline nodes are used
    gsl_vector_set(knots, 0, log(points[npointsInner].r));
    gsl_bspline_knots(knots, splinewkspace);
    for(size_t i=0; i<numPointsUsed; i++)
    {
        gsl_bspline_eval(log(SHradii[i+npointsInner]), splinevals, splinewkspace);
        for(size_t p=0; p<numBasisFnc; p++)
            gsl_matrix_set(splinematr, i, p, gsl_vector_get(splinevals, p+0));
        gsl_vector_set(originvals, i, log(1/(1/SHcoefs[i+npointsInner][0] - 1/potcenter)));
    }
    gsl_multifit_linear(splinematr, originvals, splinecoefs, covarmatr, &chisq, fitwkspace);
    // now store fitted values in coefsArray to pass to initspline routine
    coefsArray[0][0] = SHcoefs[0][0];
    for(size_t c=1; c<=Ncoefs_radial; c++)
    {
        gsl_bspline_eval(gsl_vector_get(knots, c), splinevals, splinewkspace);
        double value=0;
        for(size_t p=0; p<numBasisFnc; p++)
            value += gsl_vector_get(splinecoefs, p)*gsl_vector_get(splinevals, p+0);
        coefsArray[c][0] = 1./(exp(-value)+1/potcenter);
    }

    // construct splines for all spherical-harmonic terms separately
    for(size_t p=0; p<=Ncoefs_radial; p++)
        gsl_vector_set(knots, p, log(1+radii[p]));
    gsl_vector_set(knots, Ncoefs_radial+1, log(1+outerRadius));
    gsl_bspline_knots(knots, splinewkspace);
    for(int l=lstep; l<=lmax; l+=lstep)
    {
        for(int m=l*mmin; m<=l*mmax; m+=mstep)
        {
            // init matrix of values to fit
            for(size_t i=0; i<numPointsUsed; i++)
            {
                gsl_bspline_eval(log(1+SHradii[i+npointsInner]), splinevals, splinewkspace);
                for(size_t p=0; p<numBasisFnc; p++)
                    gsl_matrix_set(splinematr, i, p, gsl_vector_get(splinevals, p+0));
                gsl_vector_set(originvals, i, SHcoefs[i+npointsInner][l*(l+1)+m]/SHcoefs[i+npointsInner][0] * std::min<double>(points[i].r/innerBinRadius, 1));
            }
            gsl_multifit_linear(splinematr, originvals, splinecoefs, covarmatr, &chisq, fitwkspace);
            // now store fitted values in coefsArray to pass to initspline routine
            coefsArray[0][l*(l+1) + m]=0;  // not used
            for(size_t c=1; c<=Ncoefs_radial; c++)
            {
                gsl_bspline_eval(gsl_vector_get(knots, c), splinevals, splinewkspace);
                double value=0;
                for(size_t p=0; p<numBasisFnc; p++)
                    value += gsl_vector_get(splinecoefs, p)*gsl_vector_get(splinevals, p+0);
                value *= coefsArray[c][0];  // scale back (multiply by l=0,m=0 coefficient)
                coefsArray[c][l*(l+1) + m] = value;
            }
        }
    }
    initSpline(radii, coefsArray);
    
    gsl_multifit_linear_free(fitwkspace);
    gsl_bspline_free(splinewkspace);
    gsl_matrix_free(splinematr);
    gsl_vector_free(splinevals);
    gsl_vector_free(originvals);
    gsl_vector_free(splinecoefs);
    gsl_matrix_free(covarmatr);
    }  //< try
    catch(const std::bad_alloc&) {
        my_error("Not enough memory in initializing Spline potential");
        initDefault();
    }
}

void CPotentialSpline::checkSymmetry(const std::vector< vectord > &coefsArray)
{   
    SYMMETRYTYPE sym=ST_SPHERICAL;  // too optimistic:))
    const double MINCOEF=1e-8;   // if ALL coefs of a certain subset of indices are below this value, assume some symmetry
    for(size_t n=0; n<=Ncoefs_radial; n++)
    {
        for(int l=0; l<=(int)Ncoefs_angular; l++)
            for(int m=-l; m<=l; m++)
                if(fabs(coefsArray[n][l*(l+1)+m])>MINCOEF)  
                {   // nonzero coef.: check if that breaks any symmetry
                    if(l%2==1)  sym = (SYMMETRYTYPE)(sym & ~ST_REFLECTION);
                    if(m<0 || m%2==1)  sym = (SYMMETRYTYPE)(sym & ~ST_PLANESYM);
                    if(m!=0) sym = (SYMMETRYTYPE)(sym & ~ST_ZROTSYM);
                    if(l>0) sym = (SYMMETRYTYPE)(sym & ~ST_SPHSYM);
                }
    }
    mysymmetry = sym; 
}

// A special, 'unnatural' type of cubic spline, with given values of first derivative at endpoints. 
// Hacking into GSL spline machinery by defining a new spline type which differs from natural cubic spline only in initialization
typedef struct
{
  double * c;
  double * g;
  double * diag;
  double * offdiag;
} cspline_state_t;

// 1st derivatives at left- and rightmost points are given in ya[n] & ya[n+1]
static int mycspline_init (void * vstate, const double xa[], const double ya[], size_t size)
{
  cspline_state_t *state = (cspline_state_t *) vstate;

  size_t i;
  size_t num_points = size;
  size_t max_index = num_points - 1;  /* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index - 1;    /* linear system is sys_size x sys_size */

  double der1=ya[size], der2=ya[size+1];  // numbers must be supplied!
  for (i = 0; i < sys_size; i++)
    {
      const double h_i   = xa[i + 1] - xa[i];
      const double h_ip1 = xa[i + 2] - xa[i + 1];
      const double ydiff_i   = ya[i + 1] - ya[i];
      const double ydiff_ip1 = ya[i + 2] - ya[i + 1];
      const double g_i = (h_i != 0.0) ? 1.0 / h_i : 0.0;
      const double g_ip1 = (h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
      state->offdiag[i] = h_ip1;
      state->diag[i] = 2.0 * (h_ip1 + h_i);
      state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 -  ydiff_i * g_i);
      if(i == 0) {
        state->diag[i] = 1.5 * h_i + 2.0 * h_ip1;
        state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 - 1.5 * ydiff_i * g_i + 0.5 * der1);
      }
      if(i == sys_size-1) {
        state->diag[i] = 1.5 * h_ip1 + 2.0 * h_i;
        state->g[i] = 3.0 * (1.5 * ydiff_ip1 * g_ip1 - 0.5 * der2 - ydiff_i * g_i);
      }
    }

  if (sys_size == 1)
    {
      state->c[1] = state->g[0] / state->diag[0];
      return GSL_SUCCESS;
    }
  else
    {
      gsl_vector_view g_vec = gsl_vector_view_array(state->g, sys_size);
      gsl_vector_view diag_vec = gsl_vector_view_array(state->diag, sys_size);
      gsl_vector_view offdiag_vec = gsl_vector_view_array(state->offdiag, sys_size - 1);
      gsl_vector_view solution_vec = gsl_vector_view_array ((state->c) + 1, sys_size);
      
      int status = gsl_linalg_solve_symm_tridiag(&diag_vec.vector, &offdiag_vec.vector, &g_vec.vector, &solution_vec.vector);
      state->c[0] = ( 3.0*(ya[1]-ya[0])/(xa[1]>xa[0] ? xa[1]-xa[0] : 1) - 3.0*der1 - state->c[1]*(xa[1]-xa[0]) )*0.5/(xa[1]>xa[0] ? xa[1]-xa[0] : 1);
      state->c[max_index] = -( 3*(ya[max_index]-ya[max_index-1])/(xa[max_index]-xa[max_index-1]) - 3*der2 + state->c[max_index-1]*(xa[max_index]-xa[max_index-1]) )*0.5/(xa[max_index]-xa[max_index-1]);
      return status;
    }
}

static const gsl_interp_type mycspline_type = 
{
  "mycspline", 
  3,
  gsl_interp_cspline->alloc,
  &mycspline_init,
  gsl_interp_cspline->eval,
  gsl_interp_cspline->eval_deriv,
  gsl_interp_cspline->eval_deriv2,
  gsl_interp_cspline->eval_integ,
  gsl_interp_cspline->free
};

const gsl_interp_type * gsl_interp_mycspline = &mycspline_type;

// helper function for finding outer density slope
struct CGammaOutParam{
    double r1,r2,r3,K;
};
double findGammaOut(double y, void* params)
{
    return( pow(((CGammaOutParam*)params)->r2, 3-y) - pow(((CGammaOutParam*)params)->r1, 3-y))/( pow(((CGammaOutParam*)params)->r3, 3-y) - pow(((CGammaOutParam*)params)->r2, 3-y)) - ((CGammaOutParam*)params)->K;
}

void CPotentialSpline::initSpline(const vectord &radii, const std::vector< vectord > &coefsArray)
{
    if(radii[0]!=0)  my_error("Warning, in Spline potential: radii[0]!=0 (assumed to be zero)");
    if(radii.size()!=Ncoefs_radial+1 || coefsArray.size()!=Ncoefs_radial+1) my_error("Warning, in Spline potential: Ncoefs_radial+1 != array size");  // likely to get segmentation fault shortly soon
    checkSymmetry(coefsArray);   // assign nontrivial symmetry class if some of coefs are equal or close to zero
    SHradii=radii;    // copy real radii
    minr=SHradii[1];
    maxr=SHradii.back();
    std::vector<double> spnodes(Ncoefs_radial);  // scaled radii
    std::vector<double> spvalues(Ncoefs_radial+2);
    splines.assign(pow_2(Ncoefs_angular+1), (gsl_spline*)NULL);
    slopein. assign(pow_2(Ncoefs_angular+1), 1.);
    slopeout.assign(pow_2(Ncoefs_angular+1), -1.);

    // estimate outermost slope  (needed for accurate extrapolation beyond last grid point)
    gammaout=4; coefout=0;  // default
    if(Ncoefs_radial>=5)
    {
        gsl_function F;
        CGammaOutParam PP;
        PP.r1=radii[Ncoefs_radial];
        PP.r2=radii[Ncoefs_radial-1];
        PP.r3=radii[Ncoefs_radial-2];
        PP.K = (coefsArray[Ncoefs_radial][0]*PP.r1 - coefsArray[Ncoefs_radial-1][0]*PP.r2)/(coefsArray[Ncoefs_radial-1][0]*PP.r2 - coefsArray[Ncoefs_radial-2][0]*PP.r3);
        F.function=&findGammaOut;
        F.params=&PP;
        double Ylow=3.01, Yupp=10;
        if(findGammaOut(Ylow, &PP) * findGammaOut(Yupp, &PP) <0)
        {
            gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
            gsl_root_fsolver_set (s, &F, Ylow, Yupp);
            int status=0, iter=0;
            do{
                iter++;
                status = gsl_root_fsolver_iterate (s);
                gammaout= gsl_root_fsolver_root (s);
                Ylow = gsl_root_fsolver_x_lower (s);
                Yupp = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (Ylow, Yupp, 0, 1e-4);
            }
            while (status == GSL_CONTINUE && iter < 100);
            gsl_root_fsolver_free (s);
        }  // else leave gamma=4
        coefout = (1 - coefsArray[Ncoefs_radial-1][0]*PP.r2/(coefsArray[Ncoefs_radial][0]*PP.r1)) / (pow(PP.r2/PP.r1, 3-gammaout) - 1);
        if(coefout<0) coefout=0;
    }
    else { gammaout=4; coefout=0; }  // assume zero density beyond cutoff
    // first init l=0 spline which has radial scaling "log(r)" and nontrivial transformation 1/(1/phi-1/phi0)
    potcenter=coefsArray[0][0];    // here assumed that SHradii[0]=0, i.e. innermost node is at origin
    potmax=coefsArray.back()[0];
    potminr=coefsArray[1][0];
    for(size_t i=0; i<Ncoefs_radial; i++)
    {
        spnodes[i] = log(SHradii[i+1]);
        spvalues[i]= log(1/ (1/coefsArray[i+1][0] - 1/potcenter));
    }
    gammain = 2-log((coefsArray[1][0]-potcenter)/(coefsArray[2][0]-potcenter))/log(SHradii[1]/SHradii[2]);
    if(gammain<0) gammain=0; 
    if(gammain>2) gammain=2;
    coefin = (1-coefsArray[1][0]/potcenter)/pow(SHradii[1],2-gammain);
    spvalues[Ncoefs_radial] = -(2-gammain)*potcenter/coefsArray[1][0];   // derivative at leftmost node
    spvalues[Ncoefs_radial+1] = - (1+coefout*(3-gammaout))/(1 - potmax/potcenter);  // derivative at rightmost node
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_mycspline, Ncoefs_radial);
    gsl_spline_init (spline, &(spnodes.front()), &(spvalues.front()), Ncoefs_radial);
    splines[0]=spline;
#ifdef DEBUGPRINT
    std::cerr << "#gammain="<<gammain<<";  gammaout="<<gammaout<<"\n";
#endif
    ascale=1;
    // next init all higher-order splines which have radial scaling log(a+r) and value scaled to l=0,m=0 coefficient
    for(size_t i=0; i<Ncoefs_radial; i++)
        spnodes[i]=log(ascale+SHradii[i+1]);
    double C00val, C00der;
    coef0(minr, &C00val, &C00der, NULL);
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    for(int l=lstep; l<=lmax; l+=lstep)
    {
        for(int m=l*mmin; m<=l*mmax; m+=mstep)
        {
            int coefind=l*(l+1)+m;
            for(size_t i=0; i<Ncoefs_radial; i++)
                spvalues[i] = coefsArray[i+1][coefind]/coefsArray[i+1][0];
            slopein[coefind] = log(coefsArray[2][coefind]/coefsArray[1][coefind]) / log(SHradii[2]/SHradii[1]);   // estimate power-law slope of Clm(r) at r->0
            if(gsl_isnan(slopein[coefind])) slopein[coefind]=1.0;  // default
            slopein[coefind] = std::max<double>(slopein[coefind], std::min<double>(l, 2-gammain));  // the asymptotic power-law behaviour of the coefficient expected for power-law density profile
            spvalues[Ncoefs_radial] = spvalues[0] * (1+ascale/minr) * (slopein[coefind] - minr*C00der/C00val);   // derivative at innermost node
            slopeout[coefind] = log(coefsArray[Ncoefs_radial][coefind]/coefsArray[Ncoefs_radial-1][coefind]) / log(SHradii[Ncoefs_radial]/SHradii[Ncoefs_radial-1]) + 1;   // estimate slope of Clm(r)/C00(r) at r->infinity (+1 is added because C00(r) ~ 1/r at large r)
            if(gsl_isnan(slopeout[coefind])) slopeout[coefind]=-1.0;  // default
            slopeout[coefind] = std::min<double>(slopeout[coefind], std::max<double>(-l, 3-gammaout));
            spvalues[Ncoefs_radial+1] = spvalues[Ncoefs_radial-1] * (1+ascale/maxr) * slopeout[coefind];   // derivative at outermost node
            spline = gsl_spline_alloc (gsl_interp_mycspline, Ncoefs_radial);
            gsl_spline_init (spline, &(spnodes.front()), &(spvalues.front()), Ncoefs_radial);
            splines[coefind]=spline;
#ifdef DEBUGPRINT
            std::cerr << "#l="<<l<<",m="<<m<<" - inner="<<slopein[coefind] << ", outer="<<slopeout[coefind] <<"\n";
#endif
        }
    }
    bool densityNonzero = checkDensityNonzero();
    bool massMonotonic  = checkMassMonotonic();
    if(!massMonotonic || !densityNonzero) 
        my_message(std::string("Warning, ") + 
        (!massMonotonic ? "mass does not monotonically increase with radius" : "") +
        (!massMonotonic && !densityNonzero ? " and " : "") + 
        (!densityNonzero ? "density drops to zero at a finite radius" : "") + "!");
}

void CPotentialSpline::getCoefs(vectord *radii, std::vector< vectord > *coefsArray, bool useNodes) const
{
    if(radii==NULL || coefsArray==NULL) return;
    if(useNodes)
    {
        radii->resize(Ncoefs_radial+1);
        for(size_t i=0; i<=Ncoefs_radial; i++)
            (*radii)[i] = SHradii[i];
    }
    size_t numrad=radii->size();
    coefsArray->resize(numrad);
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    for(size_t i=0; i<numrad; i++)
    {
        double rad=(*radii)[i];
        double xi=log(ascale+rad);
        double Coef00;
        coef0(rad, &Coef00, NULL, NULL);
        (*coefsArray)[i].assign(pow_2(Ncoefs_angular+1), 0);
        (*coefsArray)[i][0] = Coef00;
        for(int l=lstep; l<=lmax; l+=lstep)
            for(int m=l*mmin; m<=l*mmax; m+=mstep)
            {
                int coefind=l*(l+1)+m;
                coeflm(coefind, rad, xi, &((*coefsArray)[i][l*(l+1)+m]), NULL, NULL, Coef00);
            }
    }
}

void CPotentialSpline::coeflm(size_t lm, double r, double xi, double *val, double *der, double *der2, double c0val, double c0der, double c0der2) const  // works only for l2>0
{
    double cval=0, cder=0, cder2=0;   // value and derivatives of \tilde Clm = Clm(r)/C00(r)
    if(r < maxr)
    {
      if(r > minr)  // normal interpolation
      {
        cval= gsl_spline_eval(splines[lm], xi, NULL);
        if(der!=NULL){ cder = gsl_spline_eval_deriv(splines[lm], xi, NULL)/(r+ascale);    // to save time, evaluate these only if needed
        if(der2!=NULL) cder2=(gsl_spline_eval_deriv2(splines[lm], xi, NULL)/(r+ascale) - cder)/(r+ascale); } // obviously, if der2!=NULL then assuming that der!=NULL
      }
      else  // power-law asymptotics at r<minr
      {
        cval = splines[lm]->y[0] * potminr;
        if(val!=NULL)  *val = cval * pow(r/minr, slopein[lm]);
        if(der!=NULL){ *der = (*val) * slopein[lm]/r;
        if(der2!=NULL) *der2= (*der) * (slopein[lm]-1)/r; }
        return;   // for r<minr, Clm is not scaled by C00
      }
    }
    else  // power-law asymptotics for r>maxr
    {
        cval = splines[lm]->y[splines[lm]->interp->size-1] * pow(r/maxr, slopeout[lm]);
        cder = cval * slopeout[lm]/r;
        cder2= cder * (slopeout[lm]-1)/r;
    }
    // scale by C00
    if(val!=NULL)  *val = cval*c0val;
    if(der!=NULL)  *der = cder*c0val + cval*c0der;
    if(der2!=NULL) *der2= cder2*c0val + 2*cder*c0der + cval*c0der2;
}

void CPotentialSpline::coef0(double r, double *val, double *der, double *der2) const  // works only for l=0
{
    if(r<=maxr)
    {
        double logr=log(r);
        double sval, sder, sder2;
        if(r<minr)
        {
            double ratio = 1-coefin*pow(r, 2-gammain);  // C00(r)/C00(0)
            sval = log(potcenter/coefin) - (2-gammain)*logr + log(ratio);
            sder = -(2-gammain)/ratio;
            sder2= -pow_2(sder)*(1-ratio);
        }
        else
        {
            sval = gsl_spline_eval(splines[0], logr, NULL);
            sder = gsl_spline_eval_deriv(splines[0], logr, NULL);
            sder2= gsl_spline_eval_deriv2(splines[0], logr, NULL);
        }
        double sexp = (r>0? exp(-sval) : 0);
        double vval = 1./(sexp+1/potcenter);
        if(val!=NULL)  *val = vval;
        if(der!=NULL)  *der = vval*vval*sexp/r * sder;  // this would not work for r=0 anyway...
        if(der2!=NULL) *der2= pow_2(vval/r)*sexp * (sder2 - sder + sder*sder*(2*vval*sexp-1) );
    }
    else
    {
        double r_over_maxr_g=pow(r/maxr, 3-gammaout);
        if(val!=NULL)  *val =  potmax*maxr/r * (1 - coefout*(r_over_maxr_g-1));
        if(der!=NULL)  *der = -potmax*maxr/r/r * (1 - coefout*(r_over_maxr_g*(gammaout-2) - 1));
        if(der2!=NULL) *der2=2*potmax*maxr/pow(r,3) * (1 - coefout*(r_over_maxr_g*(gammaout-1)*(gammaout-2)/2 - 1));
    }
}

double CPotentialSpline::Mass(const double r) const
{
    if(r<=0) return 0;
    double der;
    coef0(r, NULL, &der, NULL);
    return -der * r*r;   // d Phi(r)/d r = - G M(r) /r^2
}

void CPotentialSpline::computeSHCoefs(const double r, double coefsF[], double coefsdFdr[], double coefsd2Fdr2[]) const
{
    double xi = log(r+ascale);
    double val00, der00, der200;
    coef0(r, &val00, &der00, &der200);  // compute value and two derivatives of l=0,m=0 spline
    if(coefsF)    coefsF[0]    = val00*2*M_SQRTPI;
    if(coefsdFdr) coefsdFdr[0] = der00*2*M_SQRTPI;
    if(coefsd2Fdr2) coefsd2Fdr2[0] = der200*2*M_SQRTPI;
    int lmax = mysymmetry & ST_SPHSYM ? 0 : Ncoefs_angular;    // if spherical model, use only l=0,m=0 term
    int lstep= mysymmetry & ST_REFLECTION ? 2 : 1;             // if reflection symmetry, use only even l
    int mmax = mysymmetry & ST_ZROTSYM ? 0 : 1;                // if axisymmetric model, use only m=0 terms, otherwise all terms up to l (1 is the multiplying factor)
    int mmin = mysymmetry & ST_PLANESYM ? 0 : -1;              // if triaxial symmetry, do not use sine terms which correspond to m<0
    int mstep= mysymmetry & ST_PLANESYM ? 2 : 1;               // if triaxial symmetry, use only even m
    for(int l=lstep; l<=lmax; l+=lstep)
    {
        for(int m=l*mmin; m<=l*mmax; m+=mstep)
        {
            int coefind=l*(l+1)+m;
            coeflm(coefind, r, xi, 
                coefsF!=NULL ? coefsF+coefind : NULL, 
                coefsdFdr!=NULL ? coefsdFdr+coefind : NULL, 
                coefsd2Fdr2!=NULL ? coefsd2Fdr2+coefind : NULL, 
                val00, der00, der200);
            if(coefsF) coefsF[coefind]*=2*M_SQRTPI;
            if(coefsdFdr) coefsdFdr[coefind]*=2*M_SQRTPI;
            if(coefsd2Fdr2) coefsd2Fdr2[coefind]*=2*M_SQRTPI;
        }
    }
}

//----------------------------------------------------------------------------//
// simple density models, without potential
double CDensityPlummer::Rho(double X, double Y, double Z, double /*t*/) const
{
    double m = sqrt(pow_2(X) + pow_2(Y/q) + pow_2(Z/p));
    return 3/(4*M_PI*p*q) * pow(1+m*m, -2.5);
}

double CDensityIsochrone::Rho(double X, double Y, double Z, double /*t*/) const
{
    double m = sqrt(pow_2(X) + pow_2(Y/q) + pow_2(Z/p));
    double a=sqrt(m*m+1);
    return (3*(1+a)*a*a-m*m*(1+3*a))/(4*M_PI*p*q*pow(a*(1+a),3));
}

double CDensityPerfectEllipsoid::Rho(double X, double Y, double Z, double /*t*/) const
{
    double m = sqrt(pow_2(X) + pow_2(Y/q) + pow_2(Z/p));
    return 1/(M_PI*M_PI*p*q*pow_2(1+m*m));
}

// some more machinery for NFW potential with a gradual cutoff in radius
double findNFWrcutoff(double r, void* params)
{
    return 2*M_PI*(r*(-2.0+3.0*r-6.0*r*r)+6.0*gsl_pow_4(r)*log(1.0+1.0/r)+2.0*log(1.0+r)) - *((double*)params);
}

CDensityNFW::CDensityNFW(double _q, double _p, double c):
  q(_q), p(_p)
{
    if(c<1.0) c=1.0;
    double totalMass=(4*M_PI*(log(1+c)-c/(c+1)));
    // "c" is the standard NFW concentration, defined as the ratio of virial radius to scale radius, and assuming that the density drops instantly to zero beyond that radius.
    // we introduce another radius "rcutoff" such that (1) density drops more gradually (as r^-4) beyond rcutoff, and (2) total mass is the same as for truncated NFW model with given concentration
    // now the task is to find rcutoff
    gsl_function F;
    F.function=&findNFWrcutoff;
    F.params=&totalMass;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    double Xlow=0.25, Xupp=c;
    rcutoff=1.0;
    gsl_root_fsolver_set (s, &F, Xlow, Xupp);
    int status=0, iter=0;
    do{
        iter++;
        status = gsl_root_fsolver_iterate (s);
        rcutoff= gsl_root_fsolver_root (s);
        Xlow = gsl_root_fsolver_x_lower (s);
        Xupp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (Xlow, Xupp, 0, 1e-8);
    }
    while (status == GSL_CONTINUE && iter < 100);
    gsl_root_fsolver_free (s);
#ifdef DEBUGPRINT
    std::cerr << "rcutoff="<<rcutoff<<"\n";
#endif
    rcutoff = c; // Ignoring above calculations
    totalMass=(4*M_PI*(log(1+c)-c/(c+1))); // Ignoring above calculations
    norm = 1.0/(totalMass*p*q);
}

double CDensityNFW::Rho(double X, double Y, double Z, double /*t*/) const
{
    double m = sqrt(pow_2(X) + pow_2(Y/q) + pow_2(Z/p));
    double cutoff = (m < rcutoff ? 1.0 : gsl_pow_4(rcutoff / m));
    return norm / (m * pow_2(1.0 + m)) * cutoff;
}

//----------------------------------------------------------------------------//
// 'class factory':)  (makes correspondence between enum potential and symmetry types and string names)
PotentialNameMapType PotentialNames;
DensityNameMapType DensityNames;
SymmetryNameMapType SymmetryNames;

bool initPotentialAndSymmetryNameMap()     // obviously, should be called before using the map
{
    PotentialNames.clear();
    PotentialNames[CDensity::PT_LOG] = CPotentialLog::myName();
    PotentialNames[CDensity::PT_HARMONIC] = CPotentialHarmonic::myName();
    PotentialNames[CDensity::PT_DEHNEN] = CPotentialDehnen::myName();
    PotentialNames[CDensity::PT_SCALEFREE] = CPotentialScaleFree::myName();
    PotentialNames[CDensity::PT_SCALEFREESH] = CPotentialScaleFreeSH::myName();
    PotentialNames[CDensity::PT_BSE] = CPotentialBSE::myName();
    PotentialNames[CDensity::PT_SPLINE] = CPotentialSpline::myName();
    PotentialNames[CDensity::PT_NB] = CPotentialNB::myName();

    // list of density models available for BSE and Spline approximation
    DensityNames.clear();
    DensityNames[CDensity::PT_NB] = CPotentialNB::myName();   // denotes not the actual tree-code potential, but rather the fact that the density model comes from discrete points in Nbody file
    DensityNames[CDensity::PT_DEHNEN] = CPotentialDehnen::myName();
    DensityNames[CDensity::PT_PLUMMER] = CDensityPlummer::myName();
    DensityNames[CDensity::PT_PERFECTELLIPSOID] = CDensityPerfectEllipsoid::myName();
    DensityNames[CDensity::PT_ISOCHRONE] = CDensityIsochrone::myName();
    DensityNames[CDensity::PT_NFW] = CDensityNFW::myName();

    SymmetryNames[CDensity::ST_NONE]="None";
    SymmetryNames[CDensity::ST_REFLECTION]="Reflection";
    SymmetryNames[CDensity::ST_TRIAXIAL]="Triaxial";
    SymmetryNames[CDensity::ST_AXISYMMETRIC]="Axisymmetric";
    SymmetryNames[CDensity::ST_SPHERICAL]="Spherical";
    
    return rand()==42;
}
bool mapinitializer = initPotentialAndSymmetryNameMap();  // this weird move ensures that the function be implicitly called upon program start

CDensity::POTENTIALTYPE getPotentialTypeByName(std::string PotentialName)
{
    for(PotentialNameMapType::const_iterator iter=PotentialNames.begin(); iter!=PotentialNames.end(); iter++)
        if(iter->second==PotentialName) return iter->first;
    return CDensity::PT_UNKNOWN;
}

CDensity::POTENTIALTYPE getDensityTypeByName(std::string DensityName)
{
    for(DensityNameMapType::const_iterator iter=DensityNames.begin(); iter!=DensityNames.end(); iter++)
        if(iter->second==DensityName) return iter->first;
    return CDensity::PT_UNKNOWN;
}

CDensity::SYMMETRYTYPE getSymmetryTypeByName(std::string SymmetryName)
{
    for(SymmetryNameMapType::const_iterator iter=SymmetryNames.begin(); iter!=SymmetryNames.end(); iter++)
        if(iter->second==SymmetryName) return iter->first;
    return CDensity::ST_TRIAXIAL;
}

/** default parameters for all possible potentials **/
CConfigPotential configPotential = {
    3,  // N_dim
    1,  // q
    1,  // p
    0,  // Mbh
    0,  // Rc
    1,  // Gamma
    1,  // Alpha
    0,  // Omega
    10, // Ncoefs_radial
    6,  // Ncoefs_angular
    -2, // nbpotEps
    0.5,// nbpotTol
    CDensity::PT_UNKNOWN, // PotentialType
    CDensity::PT_DEHNEN,  // DensityType
    CDensity::ST_TRIAXIAL // SymmetryType
};

}
