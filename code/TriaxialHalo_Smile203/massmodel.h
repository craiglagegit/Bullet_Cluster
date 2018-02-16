// this file defines class for dealing with spherical isotropic mass models (built on a given mass(r) profile), computing all necessary dynamical quantities such as T(r), f(E) etc.
#pragma once
#include <vector>
#include <gsl/gsl_spline.h>

namespace smile{

// class defining spherical mass model
class CMassModel{
public:
    int errorcode;  // signals error from initialization
    CMassModel(std::vector<double> &initradii, std::vector<double> &initmas);  // create model from array of { r, m(r) } pairs
    ~CMassModel();
    double mass(double rad);   // return enclosed mass at given radius
    double rad(double mass);   // find radius for enclosed mass
    double dens(double rad);   // density at given radius
    double pot(double rad);    // potential at given radius
    double drhodphi(double E); // d(rho)/d(phi) at given energy E=phi
    double distrf(double E);   // distribution function of energy (negative)
    double ded(double E);      // differential energy distribution  g(E)
    double masse(double E);    // mass of particles having energies <E
    double trad(double E);     // period of radial orbit at energy E
    double rmax(double E);     // max radius at given energy
    double rcirc(double E);    // radius of circular orbit at energy E
    double lcirc(double E);    // ang.momentum at circular orbit
    double veldisp(double rad);// velocity dispersion at given radius
    double densproj(double rad);//projected density at given projected radius
    void DifCoefAvER(double E, double R, double *DEE, double *DRR, double *DLL);   // diffusion coefs depending on E and R
    void DifCoefAvE (double E, double *DEE, double *DRR, double *DLL, double *DEL=NULL, double* B=NULL);             // R-averaged diffusion coefs depending on E
    double DRRrad(double E);   // lim(DRR/R), R->0
    double samplev(double r);  // return (v,r) from distribution function
    bool hasDensityCore(double &rho0_, double &drhodr_, double &alpha_);  // informs whether model has no BH and constant-density core
    double getTotalMass() { return totalmass; }
private:
    std::vector<double> potential, logradii, loginmass, logpotential, logdens, logdistrf, logded, logtrad, logrcirc;
    gsl_interp_accel *acc_loginmass, *acc_logrmass, *acc_logpotential, *acc_logdenspot, *acc_logdistrf, *acc_logded, *acc_logtrad, *acc_logrcirc, *acc_logrmax;
    gsl_spline *spline_loginmass, *spline_logrmass, *spline_logpotential, *spline_logdenspot, *spline_logdistrf, *spline_logded, *spline_logtrad, *spline_logrcirc, *spline_logrmax;
    bool potfinite;    // whether potential is finite at origin
    double potentialcenter;  // central potential not accounting for SMBH
    double totalmass;  // total mass at infinity
    double masscenter;  // BH mass if it exists
    bool densconst;    // precautions for constant-density cores
    double alpha, rho0, coefdrhodr, coefdrhodphi;
    double minEdfpositive;
};

}