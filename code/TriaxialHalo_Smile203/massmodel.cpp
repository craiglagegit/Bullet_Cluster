#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include "massmodel.h"

namespace smile{

//----------------------------------------------------------------------------//
// modified functions for spline evaluation that extrapolate last segment linearly beyond the definition region
double gsl_spline_eval1(const gsl_spline *spline, double x, gsl_interp_accel *a )
{
    if(x<spline->interp->xmin) 
        return spline->y[0] + (x-spline->x[0])/(spline->x[1]-spline->x[0])*(spline->y[1]-spline->y[0]);
    else if(x>spline->interp->xmax)
        return spline->y[spline->size-1] + (x-spline->x[spline->size-1])/(spline->x[spline->size-2]-spline->x[spline->size-1])*(spline->y[spline->size-2]-spline->y[spline->size-1]);
    else
        return gsl_spline_eval(spline, x, a);
}
double gsl_spline_eval_deriv1(const gsl_spline *spline, double x, gsl_interp_accel *a )
{
    if(x<spline->interp->xmin) 
        return (spline->y[1]-spline->y[0])/(spline->x[1]-spline->x[0]);
    else if(x>spline->interp->xmax)
        return (spline->y[spline->size-2]-spline->y[spline->size-1])/(spline->x[spline->size-2]-spline->x[spline->size-1]);
    else
        return gsl_spline_eval_deriv(spline, x, a);
}
inline double pow_2(const double x) { return x*x; }

// interpolated mass inside radius
double CMassModel::mass(double r)
{
    if(r==0) return masscenter;
    double mu = exp(gsl_spline_eval1(spline_loginmass, log(r), acc_loginmass));
    return (mu*totalmass + masscenter)/(mu+1);
}
// interpolated radius as a function of enclosed mass
double CMassModel::rad(double m)
{
    if(m<=masscenter) return 0;
    if(m>=totalmass) return 1/sin(0.0);
    return exp(gsl_spline_eval1(spline_logrmass, log( (m-masscenter)/(totalmass-m) ), acc_logrmass));
}
// interpolated density at radius
double CMassModel::dens(double r)
{
    if(densconst && r<exp(logradii[1])) return rho0;
    double inm=mass(r);
    return gsl_spline_eval_deriv1(spline_loginmass, log(r), acc_loginmass) * (inm-masscenter) * (totalmass-inm) / (totalmass-masscenter) / (12.56637*pow(r,3.0));
}
// interpolated potential at radius
double CMassModel::pot(double r)
{
    if(r==0) return potfinite ? potentialcenter : -1/sin(0.0);
    double phi = -exp(gsl_spline_eval1(spline_logpotential, log(r), acc_logpotential));
    if(potfinite) return 1/(1/phi+1/potentialcenter); else return phi;
}
// interpolated d(rho)/d(phi) at given phi
double CMassModel::drhodphi(double phi)
{
    if(densconst && phi<potential[2]) return -coefdrhodphi*pow(phi-potentialcenter, alpha/2-1);
    double logpot = potfinite ? log(-1/(1/phi-1/potentialcenter)) : log(-phi);
    double rho = exp(gsl_spline_eval1(spline_logdenspot, logpot, acc_logdenspot));
    return gsl_spline_eval_deriv1(spline_logdenspot, logpot, acc_logdenspot) * rho / phi / (potfinite ? 1-phi/potentialcenter : 1);
}
// interpolated distribution function of energy
double CMassModel::distrf(double phi)
{
    if(phi>=0 || (minEdfpositive && phi<=minEdfpositive)) return 0;
    double logpot = potfinite ? log(-1/(1/phi-1/potentialcenter)) : log(-phi);
    double distrfint = -exp(gsl_spline_eval1(spline_logdistrf, logpot, acc_logdistrf));
    return 0.035822448 * gsl_spline_eval_deriv1(spline_logdistrf, logpot, acc_logdistrf) * distrfint / phi / (potfinite ? 1-phi/potentialcenter : 1);
}
// interpolated differential energy distribution
double CMassModel::ded(double phi)
{
    if(phi>=0) return 0;
    double logpot = potfinite ? log(-1/(1/phi-1/potentialcenter)) : log(-phi);
    return exp(gsl_spline_eval1(spline_logded, logpot, acc_logded));
}
// interpolated radial period at energy
double CMassModel::trad(double phi)
{
    double logpot = potfinite ? log(-1/(1/phi-1/potentialcenter)) : log(-phi);
    return exp(gsl_spline_eval1(spline_logtrad, logpot, acc_logtrad));
}
// interpolated radius of circular orbit at energy
double CMassModel::rcirc(double phi)
{
    double logpot = potfinite ? log(-1/(1/phi-1/potentialcenter)) : log(-phi);
    return exp(gsl_spline_eval1(spline_logrcirc, logpot, acc_logrcirc));
}
// interpolated angular momentum of circular orbit at energy
double CMassModel::lcirc(double phi)
{
    double rc=rcirc(phi);
    return sqrt(2*(phi-pot(rc)))*rc;
}
// interpolated radius of radial orbit at energy (rmax)
double CMassModel::rmax(double phi)
{
    double logpot = potfinite ? log(-1/(1/phi-1/potentialcenter)) : log(-phi);
    return exp(gsl_spline_eval1(spline_logrmax, logpot, acc_logrmax));
}
// tell whether model has constant-density core (and no BH). if yes, give its parameters
bool CMassModel::hasDensityCore(double &rho0_, double &drhodr_, double &alpha_)
{
    if(densconst) { rho0_=rho0; drhodr_=coefdrhodr; alpha_=alpha; }
    return densconst;
}

///--- helper functions for integration and root-finding ---///
// helper structure keeping params for GSL call-back functions
struct CHelperParam{
    CMassModel* model;
    double E;
    double M;
    double Phi;
    double L2;
    double kprime_over_k;
    double n;
    double p2,p3,r1,r2,r3;  // for constant-density core estimator
};
// integrand for potential
double int_pot(double r, void* params)
{
    if(r==0) return 0;
    return (((CHelperParam*)params)->model->mass(r) - ((CHelperParam*)params)->model->mass(0))/r/r;
}
// integrand for Eddington inversion formula
double int_distr(double Phi, void* params)
{
    return ((CHelperParam*)params)->model->drhodphi(-Phi)/sqrt(((CHelperParam*)params)->E-Phi);
}
// integrand for finding radial period of purely radial orbit
double int_trad(double r, void* params)
{
    return 1/sqrt(2*( ((CHelperParam*)params)->E - ((CHelperParam*)params)->model->pot(r)));
}
// integrand for finding precession period of nearly circular orbit
/*double int_tprec(double r, void* params)
{
    return 1/sqrt(2*( ((CHelperParam*)params)->E - ((CHelperParam*)params)->model->pot(r)));
}*/
// integrand for differential energy distribution
double int_ded(double r, void* params)
{
    if(r==0) return 0;
    return r*r*sqrt(std::max<double>(2*( ((CHelperParam*)params)->E - ((CHelperParam*)params)->model->pot(r)),0));
}
// integrand for f*g
double int_fg(double E, void* params)
{
    return ((CHelperParam*)params)->model->ded(E) * ((CHelperParam*)params)->model->distrf(E);
}
// integrand for sigma^2(r)
double int_sigma(double v, void* params)
{
    return ((CHelperParam*)params)->model->distrf(((CHelperParam*)params)->Phi+v*v/2) *4*M_PI/3*pow(v,4);
}
// integrand for projected density
double int_densproj(double z, void* params)
{
    return ((CHelperParam*)params)->model->dens(sqrt(pow_2(z)+pow_2(((CHelperParam*)params)->r1)));
}
// helper root-finder for a constant-density core
double find_alpha(double a, void* params)
{
    return (1-((CHelperParam*)params)->p2)*(pow(((CHelperParam*)params)->r3,a)-((CHelperParam*)params)->p3*pow(((CHelperParam*)params)->r1,a)) 
         - (1-((CHelperParam*)params)->p3)*(pow(((CHelperParam*)params)->r2,a)-((CHelperParam*)params)->p2*pow(((CHelperParam*)params)->r1,a));
}
// helper root-finder for circular radius at energy
double find_rcirc(double r, void* params)
{
    return ((CHelperParam*)params)->model->pot(r) + ((CHelperParam*)params)->model->mass(r)/2/r - ((CHelperParam*)params)->E;
}
// helper root-finder for radius at mass
double find_radius_by_mass(double r, void* params)
{
    return ((CHelperParam*)params)->model->mass(r) - ((CHelperParam*)params)->M;
}
// helper root-finder for maximum of f()*v^2
double find_max_fv(double v, void* params)
{
    return -((CHelperParam*)params)->model->distrf(v*v/2 - ((CHelperParam*)params)->Phi)*v*v;
}

// compute velocity dispersion at given radius
double CMassModel::veldisp(double rad)
{
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(1000);
    CHelperParam params;
    params.Phi=pot(rad);
    params.model=this;
    gsl_function F;
    F.function=&int_sigma;
    F.params=&params;
    double result, error, vmax=sqrt(-2*params.Phi);
    gsl_integration_qags(&F, 0, vmax, 0, 1e-7, 1000, ws, &result, &error);
    gsl_integration_workspace_free(ws);
    return sqrt(result/dens(rad));
}
// compute projected density at given radius
double CMassModel::densproj(double rad)
{
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(1000);
    CHelperParam params;
    params.r1=rad;
    params.model=this;
    gsl_function F;
    F.function=&int_densproj;
    F.params=&params;
    double result, error;
    gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, ws, &result, &error);
    gsl_integration_workspace_free(ws);
    return 2*result;
}
// compute mass of particles having energy <E
double CMassModel::masse(double E)
{
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(1000);
    double elow=potfinite?potential[0]:potential[1]*100;
    double result, error;
    gsl_function F;
    CHelperParam params;
    F.function=&int_fg;
    F.params=&params;
    params.model=this;
    gsl_integration_qags(&F, elow, E, 0, 1e-7, 1000, ws, &result, &error);
    gsl_integration_workspace_free(ws);
    return result;
}
// generate particle from f(E) at given r
double CMassModel::samplev(double radius)
{
    if(!gsl_finite(radius)) return 0;  // infinity
    gsl_min_fminimizer *min_fv = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_function F;
    CHelperParam params;
    F.function=&find_max_fv;
    F.params=&params;
    params.model=this;
    params.Phi=-pot(radius);
    // determine maximum of v^2*f(Psi-v^2/2)
    double vesc=sqrt(2*params.Phi), fmax, vlow, vupp;
    gsl_min_fminimizer_set(min_fv, &F, vesc/2, 1e-8, vesc-1e-8);
    int status, iter=0;
    do{
        iter++;
        status = gsl_min_fminimizer_iterate (min_fv);
        fmax = -gsl_min_fminimizer_f_minimum (min_fv);
        vlow = gsl_min_fminimizer_x_lower (min_fv);
        vupp = gsl_min_fminimizer_x_upper (min_fv);
        status = gsl_min_test_interval (vlow, vupp, 0, 1e-4);
    }
    while (status == GSL_CONTINUE && iter < 100);
    fmax*=1.2;
    double v, f1;
    do{
        v=vesc*rand()/RAND_MAX;//gsl_rng_uniform(randn);
        f1=fmax*rand()/RAND_MAX;//gsl_rng_uniform(randn);
    } while(f1>distrf(-params.Phi+v*v/2)*v*v);
    gsl_min_fminimizer_free(min_fv);
    return v;
}

///--- Calculation of diffusion coefficients ---///

// integrand for I0
double int_I0(double Ef, void* params)
{
    return ((CHelperParam*)params)->model->distrf(-Ef);
}
// integrand for In
double int_In(double Ef, void* params)
{
    return ((CHelperParam*)params)->model->distrf(-Ef) * pow( ((CHelperParam*)params)->Phi-Ef, ((CHelperParam*)params)->n);
}
// root-finder for v_r^2=0
double find_vr2(double r, void* params)
{
    return 2*(-((CHelperParam*)params)->model->pot(r) - ((CHelperParam*)params)->E) - ((CHelperParam*)params)->L2/(r*r);
}
// coefficients appearing in delta v_par, delta v_per
double coef_I0(void* params)
{
    double result, error; size_t niter;
    gsl_function F;
    F.function=&int_I0;
    F.params=params;
    gsl_integration_qng(&F, 0, ((CHelperParam*)params)->E, 0, 1e-4, &result, &error, &niter);
    return result;
}
double coef_In(void* params, double n)
{
    double result, error; size_t niter;
    gsl_function F;
    F.function=&int_In;
    F.params=params;
    ((CHelperParam*)params)->n=n;
    gsl_integration_qng(&F, ((CHelperParam*)params)->E, ((CHelperParam*)params)->Phi, 0, 1e-4, &result, &error, &niter);
    return result / pow(((CHelperParam*)params)->Phi-((CHelperParam*)params)->E, n);
}
// local diffusion coefs (E > 0 is binding energy, R is squared relative ang.mom.)
double Dif1(double r, void* params)
{
    double vr2=find_vr2(r, params);
    if(vr2<=0) return 0;
    return 1 / sqrt(vr2);
}
double DifEElocal(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    double vr2=find_vr2(r, params);
    if(vr2<=0) return 0;
    ((CHelperParam*)params)->Phi=Phi; 
    return 2*(Phi-E) * 2.0/3 * (coef_I0(params) + coef_In(params, 1.5)) / sqrt(vr2);
}
double DifRRlocal(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    //double L2=((CHelperParam*)params)->L2;
    double vr2=find_vr2(r, params);
    if(vr2<=0) return 0;
    ((CHelperParam*)params)->Phi=Phi; 
    double v2=2*(Phi-E);
    double mulv2par=pow_2(2+((CHelperParam*)params)->kprime_over_k*v2);
    double mulv2per=2*vr2/(v2-vr2);
    return 2.0/3 * (coef_I0(params)*(mulv2par+2*mulv2per) + 3*coef_In(params, 0.5)*mulv2per + coef_In(params, 1.5)*(mulv2par-mulv2per)) / v2 / sqrt(vr2);
}
double DifLLlocal(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    //double L2=((CHelperParam*)params)->L2;
    double vr2=find_vr2(r, params);
    if(vr2<=0) return 0;
    ((CHelperParam*)params)->Phi=Phi; 
    double v2=2*(Phi-E);
    double mulv2par=1;
    double mulv2per=0.5*vr2/(v2-vr2);
    return 2.0/3 * (coef_I0(params)*(mulv2par+2*mulv2per) + 3*coef_In(params, 0.5)*mulv2per + coef_In(params, 1.5)*(mulv2par-mulv2per)) / v2 / sqrt(vr2);
}
// L-averaged diffusion coefs (E > 0 is binding energy)
double Dif1E(double r, void* params)
{
    double v2=2*(-((CHelperParam*)params)->model->pot(r) - ((CHelperParam*)params)->E);
    if(v2<=0) return 0;
    return sqrt(v2)*r*r;
}
double DifEElocE(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    double v2=2*(Phi-E);
    if(v2<=0) return 0;
    ((CHelperParam*)params)->Phi=Phi; 
    return pow(v2,1.5)*r*r * 2.0/3 * (coef_I0(params) + coef_In(params, 1.5));
}
double DifLLlocE(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    double v2=2*(Phi-E);
    if(v2<=0) return 0;
    ((CHelperParam*)params)->Phi=Phi; 
    return sqrt(v2)*r*r * r*r/3 * (2*coef_I0(params) + coef_In(params, 0.5) + coef_In(params, 1.5));
}
double DifELlocE(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    double v2=2*(Phi-E);
    if(v2<=0) return 0;
    ((CHelperParam*)params)->Phi=Phi; 
    return v2*pow(r,3) *M_PI/6 *(coef_I0(params) + coef_In(params, 1.5));
}
double DifRRlocE(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    double v2=2*(Phi-E);
    if(v2<=0) return 0;
    ((CHelperParam*)params)->Phi=Phi; 
    double Q2=pow_2(1+((CHelperParam*)params)->kprime_over_k*v2/2);
    return pow(v2,1.5)*pow(r,6) * (coef_I0(params)*(8*Q2+2) + coef_In(params, 0.5)*3 + coef_In(params, 1.5)*(8*Q2-1))*8.0/45; 
}
double DifLLlocEfix(double r, void* params)
{
    double Phi=-((CHelperParam*)params)->model->pot(r);
    double E=((CHelperParam*)params)->E;
    double v2=2*(Phi-E);
    if(v2<=0) return 0;
    return pow(r,4.0)*((CHelperParam*)params)->model->dens(r);
}
// orbit-averaged dif.coefs at given E, R
void CMassModel::DifCoefAvER(double E, double R, double *DEE, double *DRR, double *DLL)
{
    gsl_root_fsolver* solv_vr = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_integration_workspace* ws_difcoefs = gsl_integration_workspace_alloc(1000);
    // find orbit extent (r_min, r_max)
    double rc=rcirc(-E), rmin, rmax;
    gsl_function F;
    CHelperParam params;
    F.function=&find_vr2;
    F.params=&params;
    params.model=this;
    params.E=E;
    params.L2=rc*rc* 2*(-pot(rc)-E) * R;
    double rc1=rcirc(-E*0.9999);
    params.kprime_over_k = (1 - rc1*rc1* 2*(-pot(rc1)-E*0.9999) / (rc*rc* 2*(-pot(rc)-E)))/(0.0001*E);
    int status=0, iter=0;
    double rlow=rc*R/100, rupp=rc;
    gsl_root_fsolver_set (solv_vr, &F, rlow, rupp);
    do{
        iter++;
        status = gsl_root_fsolver_iterate (solv_vr);
        rmin = gsl_root_fsolver_root (solv_vr);
        rlow = gsl_root_fsolver_x_lower (solv_vr);
        rupp = gsl_root_fsolver_x_upper (solv_vr);
        status = gsl_root_test_interval (rlow, rupp, 0, 1e-8);
    }
    while (status == GSL_CONTINUE && iter < 100);
    status=0, iter=0;
    rlow=rc, rupp=rc*4;
    gsl_root_fsolver_set (solv_vr, &F, rlow, rupp);
    do{
        iter++;
        status = gsl_root_fsolver_iterate (solv_vr);
        rmax = gsl_root_fsolver_root (solv_vr);
        rlow = gsl_root_fsolver_x_lower (solv_vr);
        rupp = gsl_root_fsolver_x_upper (solv_vr);
        status = gsl_root_test_interval (rlow, rupp, 0, 1e-8);
    }
    while (status == GSL_CONTINUE && iter < 100);
    // now integrate local dif.coefs from r_min to r_max
    F.function=&Dif1;
    double result, error;
    gsl_integration_qags (&F, rmin, rmax, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    double trad = result;
    if(DEE){
    F.function=&DifEElocal;
    gsl_integration_qags (&F, rmin, rmax, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    *DEE = result / trad;
    }
    if(DRR){
    F.function=&DifRRlocal;
    gsl_integration_qags (&F, rmin, rmax, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    *DRR = result*R*R / trad;
    }
    if(DLL){
    F.function=&DifLLlocal;
    gsl_integration_qags (&F, rmin, rmax, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    *DLL = result*params.L2 / trad;
    }
    gsl_root_fsolver_free(solv_vr);
    gsl_integration_workspace_free(ws_difcoefs);
}
// fully-averaged (ergodic) dif.coefs (depend on E only)
void CMassModel::DifCoefAvE(double E, double *DEE, double *DRR, double *DLL, double *DEL, double* B)
{
    gsl_root_fsolver* solv_vr = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_integration_workspace* ws_difcoefs = gsl_integration_workspace_alloc(1000);
    // find orbit extent (r_min, r_max)
    gsl_function F;
    CHelperParam params;
    F.params=&params;
    params.model=this;
    params.E=E;
    double rmx=rmax(-E), rc=rcirc(-E), rc1=rcirc(-E*0.9999), Lc2=rc*rc* 2*(-pot(rc)-E);
    params.kprime_over_k = (1 - rc1*rc1* 2*(-pot(rc1)-E*0.9999) / Lc2)/(0.0001*E);
    // now integrate local dif.coefs from r_min to r_max
    F.function=&Dif1E;
    double result, error;
    gsl_integration_qags (&F, 0, rmx, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    double phasevol = result;
    if(DEE){
    F.function=&DifEElocE;
    gsl_integration_qags (&F, 0, rmx, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    *DEE = result / phasevol;
    }
    if(DRR){
    F.function=&DifRRlocE;
    gsl_integration_qags (&F, 0, rmx, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    *DRR = result / phasevol / pow_2(Lc2);
    }
    if(DLL){
    F.function=&DifLLlocE;
    gsl_integration_qags (&F, 0, rmx, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    *DLL = result / phasevol;
    }
    if(DEL){
    F.function=&DifLLlocEfix;//&DifELlocE;
    gsl_integration_qags (&F, 0, rmx, 0, 1e-4, 1000, ws_difcoefs, &result, &error); 
    *DEL = result / phasevol;
    }
    if(B)
        *B = trad(-E)/(4/pow_2(lcirc(-E))*phasevol);
    gsl_root_fsolver_free(solv_vr);
    gsl_integration_workspace_free(ws_difcoefs);
}
// lim(DRR/R), R->0
double CMassModel::DRRrad(double E)
{
    double DRR;
    DifCoefAvER(E, 0.001, NULL, &DRR, NULL);
    return DRR/0.001 / 2;  // D_RR = 1/2 * <Delta R^2>
}

// create model from array of { r, m(r) } pairs
CMassModel::CMassModel(std::vector<double> &initradii, std::vector<double> &initmass)
{
    spline_loginmass=NULL;
    spline_logrmass=NULL;
    spline_logpotential=NULL;
    spline_logdenspot=NULL;
    spline_logdistrf=NULL;
    spline_logded=NULL;
    spline_logtrad=NULL;
    spline_logrcirc=NULL;
    spline_logrmax=NULL;
    acc_loginmass=NULL;
    acc_logrmass=NULL;
    acc_logpotential=NULL;
    acc_logdenspot=NULL;
    acc_logdistrf=NULL;
    acc_logded=NULL;
    acc_logtrad=NULL;
    acc_logrcirc=NULL;
    acc_logrmax=NULL;
    gsl_function F;
    CHelperParam params;
    params.model=this;
    F.params=&params;
    errorcode=0;
    if(initradii.size()!=initmass.size())
    {
        errorcode=-1;  // input arrays not equal
        return;
    }
    int Npoints=static_cast<int>(initradii.size());
    if(Npoints<5) {
        errorcode=-2;  // Too few points
        return;
    }
    totalmass=initmass.back();
    masscenter=initmass.front();

    logradii.resize(Npoints,0);
    loginmass.resize(Npoints,0);
    potential.resize(Npoints,0);
    logdistrf.resize(Npoints,0);
    logded.resize(Npoints,0);
    logdens.resize(Npoints,0);
    logpotential.resize(Npoints,0);
    logtrad.resize(Npoints,0);
    logrcirc.resize(Npoints,0);

    for(int i=1; i<Npoints-1; i++)
    {
        if(initradii[i]<=initradii[i-1]){
            errorcode=-3; // "Radii not in ascending order!\n";
            return;
        }
        if(initmass[i]<=initmass[i-1]){
            errorcode=-4; // "Mass not increasing with radius!\n";
            return;
        }
        logradii[i]=log(initradii[i]);
        loginmass[i]=log((initmass[i]-initmass[0])/(totalmass-initmass[i]));
    }
    // regularize inner and outer segments and check whether potential diverges at r=0
    logradii[0]=2*logradii[1]-logradii[2];
    loginmass[0]=2*loginmass[1]-loginmass[2];
    potfinite = masscenter==0 && (loginmass[1]-loginmass[0])/(logradii[1]-logradii[0]) > 1;
    logradii[Npoints-1]=2*logradii[Npoints-2]-logradii[Npoints-3];
    loginmass[Npoints-1]=2*loginmass[Npoints-2]-loginmass[Npoints-3];
    
    // check if density tends to constant at r=0
    alpha=log((initmass[2]-initmass[0])/(initmass[1]-initmass[0]))/log(initradii[2]/initradii[1]);
    if(potfinite && fabs(alpha-3) < 0.1)
    {
        densconst=true;
        params.p2=initmass[2]/initmass[1]*pow(initradii[1]/initradii[2],3.0);
        params.p3=initmass[3]/initmass[1]*pow(initradii[1]/initradii[3],3.0);
        params.r1=initradii[1];
        params.r2=initradii[2];
        params.r3=initradii[3];
        F.function = &find_alpha;
        int status=0, iter=0;
        double alow=0.1, aupp=5;
        if(F.function(alow, F.params) * F.function(aupp, F.params) < 0) {
            gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
            gsl_root_fsolver_set (s, &F, alow, aupp);
            do{
                iter++;
                status = gsl_root_fsolver_iterate (s);
                alpha = gsl_root_fsolver_root (s);
                alow = gsl_root_fsolver_x_lower (s);
                aupp = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (alow, aupp, 0, 1e-8);
            }
            while (status == GSL_CONTINUE && iter < 100);
            gsl_root_fsolver_free (s);
        } else alpha=2; // default for Plummer and similar types of cored potentials
        rho0 = initmass[1]*3/12.56637/pow(initradii[1],3.0);
        coefdrhodr = (1+alpha/3)*(1-params.p2)/(pow(initradii[2],alpha)-params.p2*pow(initradii[1],alpha));
        if(coefdrhodr<0) coefdrhodr=0;  // density should be monotonic with radius
#ifdef DEBUGPRINT
        std::cerr << "Constant-density core determined: rho(r) = " << rho0 << "*(1-" << coefdrhodphi << "*r^" << alpha << ")\n";
#endif
        coefdrhodphi = coefdrhodr*3/12.56637*alpha* pow(3/6.2832/rho0, alpha/2-1);
    } else densconst=false;

    // init spline to compute mass and density within radius
    acc_loginmass = gsl_interp_accel_alloc ();
    spline_loginmass = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    gsl_spline_init (spline_loginmass, &(logradii.front()), &(loginmass.front()), Npoints);

    acc_logrmass = gsl_interp_accel_alloc ();
    spline_logrmass = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    gsl_spline_init (spline_logrmass, &(loginmass.front()), &(logradii.front()), Npoints);

    // compute potential by integration of M(r)/r^2
    F.function=&int_pot;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    potential[0]=0;
    for(int i=1; i<Npoints-1; i++)
    {
        double result, error;
        gsl_integration_qags (&F, initradii[i-1], initradii[i], 0, 1e-7, 1000, w, &result, &error); 
        potential[i]=potential[i-1]+result;
    }
    // subtle correction for outermost potential: first determine the power-law index of density at r->inf
    double gammaout = 3 + log((totalmass-initmass[Npoints-3])/(totalmass-initmass[Npoints-2]))/log(initradii[Npoints-2]/initradii[Npoints-3]);
    potential[Npoints-1] = potential[Npoints-2] + (totalmass - (totalmass-initmass[Npoints-2])/(gammaout-2))/initradii[Npoints-2];
    potentialcenter=-potential[Npoints-1];
    for(int i=0; i<Npoints-1; i++)
    {
        potential[i]-=potential[Npoints-1] + (masscenter?masscenter/initradii[i]:0);
        logpotential[i]= potfinite ? log(-1/(1/potential[i]-1/potential[0])) : log(-potential[i]);
    }
    logpotential[0]= /*inmass[0] ? logpotential[1]+logradii[1]-logradii[0] :*/ 2*logpotential[1]-logpotential[2]; ///!!!
    logpotential[Npoints-1]=2*logpotential[Npoints-2]-logpotential[Npoints-3];
    
    // init splines: potential(radius), density(potential), radius(potential)
    for(int i=1; i<Npoints; i++)
        logdens[i]=log(dens(initradii[i]));
    logdens[0]=2*logdens[1]-logdens[2];
    logdens[Npoints-1]=logdens[Npoints-2]+(logdens[Npoints-2]-logdens[Npoints-3])/(logradii[Npoints-2]-logradii[Npoints-3])*(logradii[Npoints-1]-logradii[Npoints-2]);

    acc_logpotential = gsl_interp_accel_alloc ();
    spline_logpotential = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    gsl_spline_init (spline_logpotential, &(logradii.front()), &(logpotential.front()), Npoints);
 
    acc_logdenspot = gsl_interp_accel_alloc ();
    spline_logdenspot = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    std::vector<double> logd1(Npoints), logp1(Npoints);
    for(int i=0; i<Npoints; i++) { logd1[i]=logdens[Npoints-1-i]; logp1[i]=logpotential[Npoints-1-i]; }
    gsl_spline_init (spline_logdenspot, &(logp1.front()), &(logd1.front()), Npoints);

    acc_logrmax = gsl_interp_accel_alloc ();
    spline_logrmax = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    for(int i=0; i<Npoints; i++) { logd1[i]=logradii[Npoints-1-i]; }
    gsl_spline_init (spline_logrmax, &(logp1.front()), &(logd1.front()), Npoints);

    // compute integrated DF using Eddington inversion formula
    minEdfpositive=0;
    F.function=&int_distr;
    for(int i=1; i<Npoints-1; i++)
    {
        params.E=-potential[i];
        double result, error;
        gsl_integration_qags (&F, 0, params.E, 0, 1e-7, 1000, w, &result, &error); 
        logdistrf[i]=result<0 ? log(-result) : -42.0;
        if(result>=0) minEdfpositive=-params.E;
    }
    logdistrf[0]=2*logdistrf[1]-logdistrf[2];
    logdistrf[Npoints-1]=2*logdistrf[Npoints-2]-logdistrf[Npoints-3];
    if(minEdfpositive) errorcode=-5;  // "Warning, f(E) not positive!\n";

    // init spline for f(E)
    acc_logdistrf = gsl_interp_accel_alloc ();
    spline_logdistrf = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    for(int i=0; i<Npoints; i++) { logd1[i]=logdistrf[Npoints-1-i]; }
    gsl_spline_init (spline_logdistrf, &(logp1.front()), &(logd1.front()), Npoints);

    // compute differential energy distribution  g(E)
    F.function=&int_ded;
    for(int i=1; i<Npoints-1; i++)
    {
        params.E=potential[i];
        double result, error;
        gsl_integration_qags (&F, 0, initradii[i], 0, 1e-7, 1000, w, &result, &error); 
        logded[i]=log(result*157.91367);  // 16 pi^2
    }
    logded[0]=2*logded[1]-logded[2];
    logded[Npoints-1]=2*logded[Npoints-2]-logded[Npoints-3];

    // init spline for g(E)
    acc_logded = gsl_interp_accel_alloc ();
    spline_logded = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    for(int i=0; i<Npoints; i++) { logd1[i]=logded[Npoints-1-i]; }
    gsl_spline_init (spline_logded, &(logp1.front()), &(logd1.front()), Npoints);

    // compute radial period as a function of energy
    F.function=&int_trad;
    for(int i=1; i<Npoints-1; i++)
    {
        params.E=potential[i];
        double result, error;
        gsl_integration_qags (&F, 0, initradii[i], 0, 1e-7, 1000, w, &result, &error); 
        logtrad[i]=log(2*result);
    }
    logtrad[0]=2*logtrad[1]-logtrad[2];
    logtrad[Npoints-1]=2*logtrad[Npoints-2]-logtrad[Npoints-3];

    // init spline for T_rad(E)
    acc_logtrad = gsl_interp_accel_alloc ();
    spline_logtrad = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    for(int i=0; i<Npoints; i++) { logd1[i]=logtrad[Npoints-1-i]; }
    gsl_spline_init (spline_logtrad, &(logp1.front()), &(logd1.front()), Npoints);

    // compute radius of circular orbit as a function of energy
    F.function=&find_rcirc;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    for(int i=1; i<Npoints-1; i++)
    {
        params.E=potential[i];
        int status=0, iter=0;
        double rlow=initradii[i]/4, rupp=initradii[i], result;
        gsl_root_fsolver_set (s, &F, rlow, rupp);
        do{
            iter++;
            status = gsl_root_fsolver_iterate (s);
            result = gsl_root_fsolver_root (s);
            rlow = gsl_root_fsolver_x_lower (s);
            rupp = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (rlow, rupp, 0, 1e-8);
        }
        while (status == GSL_CONTINUE && iter < 100);
        logrcirc[i]=log(result);
    }
    gsl_root_fsolver_free (s);
    logrcirc[0]=2*logrcirc[1]-logrcirc[2];
    logrcirc[Npoints-1]=2*logrcirc[Npoints-2]-logrcirc[Npoints-3];

    // init spline for r_circ(E)
    acc_logrcirc = gsl_interp_accel_alloc ();
    spline_logrcirc = gsl_spline_alloc (gsl_interp_cspline, Npoints);
    for(int i=0; i<Npoints; i++) { logd1[i]=logrcirc[Npoints-1-i]; }
    gsl_spline_init (spline_logrcirc, &(logp1.front()), &(logd1.front()), Npoints);

    gsl_integration_workspace_free (w);
}
CMassModel::~CMassModel()
{
    gsl_spline_free(spline_loginmass);
    gsl_spline_free(spline_logrmass);
    gsl_spline_free(spline_logpotential);
    gsl_spline_free(spline_logdenspot);
    gsl_spline_free(spline_logdistrf);
    gsl_spline_free(spline_logded);
    gsl_spline_free(spline_logtrad);
    gsl_spline_free(spline_logrcirc);
    gsl_spline_free(spline_logrmax);
    gsl_interp_accel_free(acc_loginmass);
    gsl_interp_accel_free(acc_logrmass);
    gsl_interp_accel_free(acc_logpotential);
    gsl_interp_accel_free(acc_logdenspot);
    gsl_interp_accel_free(acc_logdistrf);
    gsl_interp_accel_free(acc_logded);
    gsl_interp_accel_free(acc_logtrad);
    gsl_interp_accel_free(acc_logrcirc);
    gsl_interp_accel_free(acc_logrmax);
}

}