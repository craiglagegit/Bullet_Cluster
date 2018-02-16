#include "orbit.h"
#include <cassert>
#include <algorithm>
#include <cmath>
#include <set>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include "potential.h"
#include "stringconv.h"

namespace smile{
 
CConfigOrbit configOrbit = {
    1e-10, // accuracyRel
    1e-15, // accuracyAbs
    0.3,   // accTreecode
    false, // treecodeSymmetrizeTimestep
    0.05,  // adaptiveTimeThreshold
    true,  // freqStrictOrder
    0.98,  // minFreqX
    0.99   // minFreqYZ
};

template<typename NumT> std::string CPericenterInformation<NumT>::toString() const
{
    return "Lmin\xB2="+convertToString(fitIntercept,4)+" + "+convertToString(fitSlope,3)+"*w [sct="+convertToString(fitScatter,3)+",sig="+convertToString(fitSignificance,5)+"]\nLcirc\xB2="+convertToString(Lcirc2,3)+"\n";
}
template std::string CPericenterInformation<float >::toString() const;
template std::string CPericenterInformation<double>::toString() const;

template<typename NumT> std::string COrbitInformation<NumT>::toString() const
{
    return Description 
        + "\nlfx=" + convertToString(lf[0],8)
        + "\nlfy=" + convertToString(lf[1],8)
        + "\nlfz=" + convertToString(lf[2],8)
        + "\nFreqDiff=" + convertToString(lfccdiff,3) 
        + (lambda>=0 ? "\nLyapunovExp=" + convertToString(lambda,3) : "")
        + "\nInertia=" + convertToString(inertia[0],3) + "," +  convertToString(inertia[1],3) + "," +  convertToString(inertia[2],3)
        + "\nLz=" + convertToString(Lavg[2],3)+"\xB1"+convertToString(Lvar[2],3) 
        + (Lavg[0]!=0 ? ", Lx=" + convertToString(Lavg[0],3)+"\xB1"+convertToString(Lvar[0],3) : "") 
        + "\nE=" + convertToString(Einit) + ", deltaE=" + convertToString(Ediff,3) + "\n";
}

template<typename NumT> std::string CTrajSampleInformation<NumT>::toString() const
{
    if(trajSample.empty()) return std::string(); 
    else return convertToString(trajSample.size())+" sampling points\n";
}

template<typename NumT> std::string CPoincareInformation<NumT>::toString() const
{
    if(ps.empty()) return std::string(); 
    else return convertToString(ps.size())+" points in Poincare section\n";
}
template std::string CPoincareInformation<float >::toString() const;
template std::string CPoincareInformation<double>::toString() const;

////-----------------   Below follows COrbit class  -------------////

template<typename NumT> COrbit::COrbit(const COrbitInitData<NumT> &_InitData, const vectorRuntimeFncCreators* creators) :
  InitData(completeInitData(_InitData)),
  potential(InitData.potential),
  N_dim(potential->N_dim)
{
    hout=0;  // timestep size - autodetect
    endCond=InitData.initCond;
    intTime=0;
    state=OS_INITIALIZED;
    if(InitData.calcLyapunov)
    {
        needrenorm=allowrenorm=false;
        // initialize deviation vector
        double w2=0;
        for(unsigned int k=0; k<N_dim; k++) 
        { 
            devVector.Pos[k]=rand()*1.0/RAND_MAX; 
            devVector.Vel[k]=rand()*1.0/RAND_MAX; 
            w2+=pow_2(devVector.Pos[k]) + pow_2(devVector.Vel[k]); 
        }
        double wmag=sqrt(w2);
#ifdef LYAPUNOVVAREQ
        lnwadd=0;
        for(unsigned int k=0; k<N_dim; k++) 
        {
            devVector.Pos[k]/=wmag;
            devVector.Vel[k]/=wmag;
        }
#else
        lnwadd = -log(LYAP_DEV_INIT);
        for(unsigned int k=0; k<N_dim; k++) 
        {
            devVector.Pos[k]=InitData.initCond.Pos[k] + devVector.Pos[k]/wmag * LYAP_DEV_INIT;
            devVector.Vel[k]=InitData.initCond.Vel[k] + devVector.Vel[k]/wmag * LYAP_DEV_INIT;
        }
#endif
    }
    if(creators!=NULL)  //< construct timestep function objects using creator classes
        for(vectorRuntimeFncCreators::const_iterator iter=creators->begin(); iter!=creators->end(); iter++)
        {
            CBasicOrbitRuntimeFnc* fnc=(*iter)->createRuntimeFnc(this);
            if(fnc!=NULL) RuntimeFnc.push_back(fnc);
        }
}
// explicit instantiations of constructor for two variants of template argument
template COrbit::COrbit(const COrbitInitData<float > &_InitData, const vectorRuntimeFncCreators* creators);
template COrbit::COrbit(const COrbitInitData<double> &_InitData, const vectorRuntimeFncCreators* creators);

template<typename NumT> COrbit::COrbit(const CPotential* _potential, const double _timeStep, const double _timeUnit, const std::vector< CPosVelPoint<NumT> > &_traj) :
    InitData(completeInitData(COrbitInitData<double>(_potential, _timeStep, _timeUnit, /*initCond*/ _traj[0], /*calcLyapunov*/ false))),
    potential(_potential),
    N_dim(potential->N_dim)
{
    hout=0;
    state=OS_DONE;
    endCond=_traj[_traj.size()-1];
    intTime=InitData.timeStep*_traj.size();
    RuntimeFnc.push_back( new COrbitRuntimeTrajectory(this, _traj));
    RuntimeFnc.push_back( new COrbitRuntimePericenter(this, _traj));
    for(vectorRuntimeFncs::iterator iter=RuntimeFnc.begin(); iter!=RuntimeFnc.end(); iter++)
        (*iter)->Finish();
};
// explicit instantiations of constructor for two variants of template argument
template COrbit::COrbit(const CPotential* _potential, const double _timeStep, const double _timeUnit, const std::vector< CPosVelPoint<float > > &_traj);
template COrbit::COrbit(const CPotential* _potential, const double _timeStep, const double _timeUnit, const std::vector< CPosVelPoint<double> > &_traj);

COrbit::~COrbit()
{
    // if there were any timestep functions attached, destroy 'em all!!
    for(vectorRuntimeFncs::iterator iter=RuntimeFnc.begin(); iter!=RuntimeFnc.end(); iter++)
        delete (*iter);
}

std::string COrbit::toString(bool runtimeFncInfo) const 
{
    std::string strInfo;
    if(runtimeFncInfo)  
        for(size_t rf=0; rf<getInfoNum(); rf++)
        {
            const CBasicInformation* info=getNewInfo(rf);
            strInfo += info->toString();
            delete info;
        }
    strInfo+="fnc.eval="+convertToString(nfcn)+", naccp="+convertToString(naccpt)+", nrej="+convertToString(nrejct)+"\n";
    return strInfo;
}

double COrbit::getInterpolatedTrajectory(unsigned int c, double t) const // returns interpolated trajectory (c-th coordinate);
{
    if(t<xold || t>xout) return 0;
    if(InitData.calcLyapunov && c==N_dim*2)  // return log(w)
    {
        double wx =contd8(N_dim*2,t);
        double wy =contd8(N_dim*2+1,t);
        double wz =(N_dim==3)?contd8(N_dim*2+2,t):0;
        double wvx=contd8(N_dim*3,t);
        double wvy=contd8(N_dim*3+1,t);
        double wvz=(N_dim==3)?contd8(N_dim*3+2,t):0;
#ifndef LYAPUNOVVAREQ
        wx -= contd8(0, t);
        wy -= contd8(1, t);
        wz -= (N_dim==3)?contd8(2, t):0;
        wvx-= contd8(N_dim, t);
        wvy-= contd8(N_dim+1, t);
        wvz-= (N_dim==3)?contd8(5, t):0;
#endif
        return (log(pow_2(wx) + pow_2(wy) + pow_2(wz) + pow_2(wvx) + pow_2(wvy) + pow_2(wvz))/2 + lnwadd);
    }
    assert(c<N_dim*2);
    if(leapfrog)
        return contlf(c, t);
    else
        return contd8(c, t);
}

// initialize variables, then integrate using dop853 routine
void COrbit::integrateToTime(double time)
{
    double t,vars[N_DIM*4];
    double tend;           /* final t-value */
    double rtoler=configOrbit.accuracyRel; /* relative error tolerance */
    double atoler=configOrbit.accuracyAbs; /* absolute error tolerance */
    if(potential->PotentialType()==CDensity::PT_NB) 
    {
        leapfrog=true;     /* treecode has discontinuous potential, so low-order leapfrog integrator performs better */
        if(N_dim!=N_DIM || InitData.calcLyapunov) return;  // allow run only in 3d without Lyapunov exp
    }
    else leapfrog=false;
    int itoler=0;          /* switch for rtoler and atoler */
    int iout=2;            /* switch for calling solout */
    double uround=0.0;     /* rounding unit */
    double safe=0.0;       /* safety factor */
    double fac1=0.0;       /* parameters for step size selection */
    double fac2=0.0;
    double beta=0.0;       /* for stabilized step size control */
    double hmax=0.0;       /* maximal step size */
    double h=hout;         /* initial step size */
    long nmax=NUMSTEP_MAX; /* maximal number of allowed steps */
    int meth=0;            /* switch for the choice of the coefficients */
    long nstiff=0;         /* test for stiffness */
    unsigned neq;          /* number of equations */
    unsigned nrdens;       /* number of components for which dense output is required */
    unsigned licont=0;     /* declared length of icont */
    unsigned* icont=NULL;  /* indexes for which dense output is required  NULL=all */
    //FILE* fileout=stderr;    /* messages stream */

    t = intTime;  // current value 
    tend = time;  // value at the end of integration
    if(time<=0 || time<=intTime || InitData.timeUnit<=0) return;  // apparent error
    neq=nrdens = N_dim*2*(InitData.calcLyapunov?2:1);
    for(unsigned int d=0; d<N_dim; d++)
    {
        vars[d]      = endCond.Pos[d];  // current final values (they are equal to initial values at the beginning of integration)
        vars[d+N_dim]= endCond.Vel[d];
        //traj.PosVelData.reserve( std::max<int>(0,(int)ceil(tend/InitData.timeStep-1e-8)));
    }
    drdtsign= vars[0]*vars[N_dim]+vars[1]*vars[N_dim+1]+(N_dim==3?vars[2]*vars[5]:0) >0 ? 1 : -1;

    if(InitData.calcLyapunov)
    {
        double w2=0;
        for(unsigned int d=0; d<N_dim; d++)
        {
            vars[N_dim*2+d]=devVector.Pos[d];  // stored value
            vars[N_dim*3+d]=devVector.Vel[d];  // stored value
#ifdef LYAPUNOVVAREQ
            w2+=pow_2(devVector.Pos[d]);
            w2+=pow_2(devVector.Vel[d]);
#else
            w2+=pow_2(devVector.Pos[d]-vars[d]) + pow_2(devVector.Vel[d]-vars[d+N_dim]);
#endif
        }
        //lnw.push_back(log(w2)/2 + lnwadd);
    }
    /*if(grid && t>0)   ///!!! regression
    {
        latestcelltimes.clear();
        latestcelltimes.resize(grid->cellsInModel+1, 0);
    }*/

    state=OS_RUNNING;

    dop853(neq, /*fcn,*/t, vars, tend, &rtoler, &atoler, itoler, /*solout,*/iout,/*fileout,*/
        uround, safe, fac1,fac2, beta, hmax,h , nmax, meth, nstiff, nrdens, icont, licont);
    //latestIntTime = xout-t;  // duration of integration interval
    intTime = xout;  // latest time value reached in the integrator
    for(unsigned int c=0; c<N_dim; c++)
    {
        endCond.Pos[c]=vars[c];
        endCond.Vel[c]=vars[c+N_dim];
    }
}

void COrbit::integrateAdaptiveTime(double minTime, double /*maxTime*/)
{  ///!!! regression 
    //if(grid==NULL)    // used only in the context of creating Schwarzschild model orbit library
    {
        integrateToTime(minTime);
        return;
    }
/*    if(maxTime<minTime) maxTime=minTime;
    double interval=minTime/2;
    integrateToTime(interval);
    if(time==0) die;
    double delta=1;
    int numint=1;
    bool goon=true;
    while(goon)
    {
        numint++;
        interval = std::min<double>((minTime/2)*pow(1.1, numint-2), maxTime-intTime);
        integrateToTime(intTime+interval);
        delta=diffCellTimes();
        goon = intTime < maxTime-1e-5 && delta*3.0/(2+numint)>configOrbit.adaptiveTimeThreshold;
    }*/
}

void COrbit::finishIntegration()
{
    state=OS_DONE;
    for(vectorRuntimeFncs::iterator iter=RuntimeFnc.begin(); iter!=RuntimeFnc.end(); iter++)
        (*iter)->Finish();
}

// performs regular-timestep output (store trajectory, populate Poincare section, etc). Called from dop853 integration routine
void COrbit::solout(long /*nr*/, double told, double t, double y[], unsigned /*n*/)
{
    // watch for peri/apocenter passage (only for Lyapunov deviation vector renormalization which should be made at apocenter to avoid unacceptably large error growth)
    int prevsign=drdtsign;
    drdtsign= y[0]*y[N_dim]+y[1]*y[N_dim+1]+(N_dim==3?y[2]*y[5]:0) >0 ? 1 : -1;
    if(InitData.calcLyapunov && prevsign>0 && drdtsign<=0 && needrenorm) allowrenorm=true; // passed apocenter

    // perform specific tasks at each timestep if they are given
    for(vectorRuntimeFncs::iterator iter=RuntimeFnc.begin(); iter!=RuntimeFnc.end(); iter++)
        (*iter)->Timestep(told, t, y);
}

void COrbit::fcn(unsigned int n, double t, double *v, double *f)
{
    potential->DiffEq(n, t, v, f);
}

// leap-frog integrator for frozen-N treecode potential
int COrbit::intlf(unsigned int n, double x, double *y, double xend)
{
    double f[4*N_DIM+2];
    double hprev=0, h=hprev;
    fcn(n, x, y, f);
    naccpt=nrejct=nstep=0; nfcn=1;
    while(x<xend && state!=OS_NEEDTOTERMINATE)
    {
        for(unsigned int d=0; d<N_DIM; d++)
        {
            yprev[d] = y[d];
            yprev[d+N_DIM] = y[d+N_DIM];
            aprev[d] = f[d+N_DIM];
        }
        xold=x;
        // first estimate timestep based on two criteria: minimum distance to a particle and maximum acceleration; these were computed during tree-walk
        double vel = sqrt(pow_2(y[N_DIM])+pow_2(y[N_DIM+1])+pow_2(y[N_DIM+2]));
        double acc = sqrt(pow_2(f[N_DIM])+pow_2(f[N_DIM+1])+pow_2(f[N_DIM+2]));
        double ts1 = vel>0 ? f[N_DIM*2]/vel : xend;
        double ts2 = acc>0 ? sqrt(f[N_DIM*2]/acc) : xend;
        if(potential->Mbh>0) ts1=std::min<double>(ts1, sqrt(pow_2(y[0])+pow_2(y[1])+pow_2(y[2]))/vel*TREECODE_TIMESTEP_BH_FACTOR);
        h = std::min<double> (ts1, ts2) * configOrbit.accTreecode;
        if(configOrbit.treecodeSymmetrizeTimestep){
            // iteratively estimate timestep according to (Hut, Makino & McMillan 1995)
            for(unsigned int d=0; d<N_DIM; d++)  // first step for position
                y[d] = yprev[d] + yprev[d+N_DIM]*h + aprev[d]*h*h/2;
            fcn(n, x+h, y, f);
            nfcn++;
            for(unsigned int d=0; d<N_DIM; d++)
                y[d+N_DIM] = yprev[d+N_DIM] + (aprev[d]+f[d+N_DIM])*h/2;
            vel = sqrt(pow_2(y[N_DIM])+pow_2(y[N_DIM+1])+pow_2(y[N_DIM+2]));
            acc = sqrt(pow_2(f[N_DIM])+pow_2(f[N_DIM+1])+pow_2(f[N_DIM+2]));
            ts1 = vel>0 ? f[N_DIM*2]/vel : xend;
            ts2 = acc>0 ? sqrt(f[N_DIM*2]/acc) : xend;
            if(potential->Mbh>0) ts1=std::min<double>(ts1, sqrt(pow_2(y[0])+pow_2(y[1])+pow_2(y[2]))/vel*TREECODE_TIMESTEP_BH_FACTOR);
            double h1 = std::min<double> (ts1, ts2) * configOrbit.accTreecode;
            h = (h+h1)/2;
        }
        if(x+h>xend) h=xend-x;
        x+=h;
        for(unsigned int d=0; d<N_DIM; d++)  // first step for position
            y[d] = yprev[d] + yprev[d+N_DIM]*h + aprev[d]*h*h/2;
        fcn(n, x, y, f);
        nfcn++;
        for(unsigned int d=0; d<N_DIM; d++)
        {
            anew[d] = f[d+N_DIM];
            y[d+N_DIM] = yprev[d+N_DIM] + (aprev[d]+anew[d])*h/2;
        }
        hprev=h;
        hout=h;
        xout=x;
        solout (naccpt+1, xold, x, y, n);
        naccpt++; nstep++;
        if(!gsl_finite(h) || fabs(h) <= fabs(x) * 1E-14) state=OS_NEEDTOTERMINATE;   // timestep too small
        if(nstep>NUMSTEP_MAX) state=OS_NEEDTOTERMINATE;  // too many steps
    }
    return 1;
}

double COrbit::contlf(unsigned int d, double x) const
{
    double dh=x-xold;
    if(d<N_DIM)  // interpolate position
        return yprev[d] + yprev[d+N_DIM]*dh + aprev[d]*dh*dh/2/* + (anew[d]-aprev[d])/hout*dh*dh*dh/6*/;
    else   // interpolate velocity
        return yprev[d] + aprev[d-N_DIM]*dh + (anew[d-N_DIM]-aprev[d-N_DIM])/hout*dh*dh/2;
}
    
///////////////////////////////////////////////////////
////  below follow dop853.c encapsulated routines  ////

double COrbit::hinit (unsigned int n, /*FcnEqDiff fcn,*/ double x, double* y,
              double posneg, double* f0, double* f1, double* yy1, int iord,
              double hmax, double* atoler, double* rtoler, int itoler)
{
  double   dnf, dny, atoli, rtoli, sk, h, h1, der2, der12, sqre;
  unsigned i;

  dnf = 0.0;
  dny = 0.0;
  atoli = atoler[0];
  rtoli = rtoler[0];

  if (!itoler)
    for (i = 0; i < n; i++)
    {
      sk = atoli + rtoli * fabs(y[i]);
      sqre = f0[i] / sk;
      dnf += sqre*sqre;
      sqre = y[i] / sk;
      dny += sqre*sqre;
    }
  else
    for (i = 0; i < n; i++)
    {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqre = f0[i] / sk;
      dnf += sqre*sqre;
      sqre = y[i] / sk;
      dny += sqre*sqre;
    }

  if ((dnf <= 1.0E-10) || (dny <= 1.0E-10))
    h = 1.0E-6;
  else
    h = sqrt (dny/dnf) * 0.01;

  h = std::min<double>(h, hmax);
  h *= posneg>=0?1:-1;

  /* perform an explicit Euler step */
  for (i = 0; i < n; i++)
    yy1[i] = y[i] + h * f0[i];
  fcn (n, x+h, yy1, f1);

  /* estimate the second derivative of the solution */
  der2 = 0.0;
  if (!itoler)
    for (i = 0; i < n; i++)
    {
      sk = atoli + rtoli * fabs(y[i]);
      sqre = (f1[i] - f0[i]) / sk;
      der2 += sqre*sqre;
    }
  else
    for (i = 0; i < n; i++)
    {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqre = (f1[i] - f0[i]) / sk;
      der2 += sqre*sqre;
    }
  der2 = sqrt (der2) / h;

  /* step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01 */
  der12 = std::max<double> (fabs(der2), sqrt(dnf));
  if (der12 <= 1.0E-15)
    h1 = std::max<double> (1.0E-6, fabs(h)*1.0E-3);
  else
    h1 = pow (0.01/der12, 1.0/(double)iord);
  h = std::min<double> (100.0 * h, std::min<double> (h1, hmax));

  return h*(posneg>=0?1:-1);

} /* hinit */


/* core integrator */
int COrbit::dopcor (unsigned n, /*FcnEqDiff fcn,*/ double x, double* y, double xend,
                   double hmax, double h, double* rtoler, double* atoler,
                   int itoler, /*FILE* fileout, SolTrait solout,*/ int iout,
                   unsigned int nmax, double uround, int /*meth*/, long nstiff, double safe,
                   double beta, double fac1, double fac2, unsigned* icont)
{
  double   facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
  double   atoli, rtoli, hlamb, err, sk, hnew, ydiff, bspl;
  double   stnum, stden, sqre, err2, erri, deno;
  int      iasti, iord, irtrn, reject, last, nonsti;
  unsigned i, j;
  double   c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c14, c15, c16;
  double   b1, b6, b7, b8, b9, b10, b11, b12, bhh1, bhh2, bhh3;
  double   er1, er6, er7, er8, er9, er10, er11, er12;
  double   a21, a31, a32, a41, a43, a51, a53, a54, a61, a64, a65, a71, a74, a75, a76;
  double   a81, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98;
  double   a101, a104, a105, a106, a107, a108, a109;
  double   a111, a114, a115, a116, a117, a118, a119, a1110;
  double   a121, a124, a125, a126, a127, a128, a129, a1210, a1211;
  double   a141, a147, a148, a149, a1410, a1411, a1412, a1413;
  double   a151, a156, a157, a158, a1511, a1512, a1513, a1514;
  double   a161, a166, a167, a168, a169, a1613, a1614, a1615;
  double   d41, d46, d47, d48, d49, d410, d411, d412, d413, d414, d415, d416;
  double   d51, d56, d57, d58, d59, d510, d511, d512, d513, d514, d515, d516;
  double   d61, d66, d67, d68, d69, d610, d611, d612, d613, d614, d615, d616;
  double   d71, d76, d77, d78, d79, d710, d711, d712, d713, d714, d715, d716;

  /* memory allocation: moved here from dop853 (change from the original code) */
  /* the vectors below will be automatically destroyed upon exit; double* pointers to correspoinding vectors are used in the original code */
  std::vector<double> vec_yy1(n), vec_k1(n), vec_k2(n), vec_k3(n), vec_k4(n), 
      vec_k5(n), vec_k6(n), vec_k7(n), vec_k8(n), vec_k9(n), vec_k10(n);
  double *yy1=&(vec_yy1.front());
  double *k1=&(vec_k1.front()), *k2=&(vec_k2.front()), *k3=&(vec_k3.front()), *k4=&(vec_k4.front()), *k5=&(vec_k5.front()),
         *k6=&(vec_k6.front()), *k7=&(vec_k7.front()), *k8=&(vec_k8.front()), *k9=&(vec_k9.front()), *k10=&(vec_k10.front());

  /* initialisations */
//  switch (meth)
  {
//    case 1:

      c2  = 0.526001519587677318785587544488E-01;
      c3  = 0.789002279381515978178381316732E-01;
      c4  = 0.118350341907227396726757197510E+00;
      c5  = 0.281649658092772603273242802490E+00;
      c6  = 0.333333333333333333333333333333E+00;
      c7  = 0.25E+00;
      c8  = 0.307692307692307692307692307692E+00;
      c9  = 0.651282051282051282051282051282E+00;
      c10 = 0.6E+00;
      c11 = 0.857142857142857142857142857142E+00;
      c14 = 0.1E+00;
      c15 = 0.2E+00;
      c16 = 0.777777777777777777777777777778E+00;

      b1 =   5.42937341165687622380535766363E-2;
      b6 =   4.45031289275240888144113950566E0;
      b7 =   1.89151789931450038304281599044E0;
      b8 =  -5.8012039600105847814672114227E0;
      b9 =   3.1116436695781989440891606237E-1;
      b10 = -1.52160949662516078556178806805E-1;
      b11 =  2.01365400804030348374776537501E-1;
      b12 =  4.47106157277725905176885569043E-2;

      bhh1 = 0.244094488188976377952755905512E+00;
      bhh2 = 0.733846688281611857341361741547E+00;
      bhh3 = 0.220588235294117647058823529412E-01;

      er1  =  0.1312004499419488073250102996E-01;
      er6  = -0.1225156446376204440720569753E+01;
      er7  = -0.4957589496572501915214079952E+00;
      er8  =  0.1664377182454986536961530415E+01;
      er9  = -0.3503288487499736816886487290E+00;
      er10 =  0.3341791187130174790297318841E+00;
      er11 =  0.8192320648511571246570742613E-01;
      er12 = -0.2235530786388629525884427845E-01;

      a21 =    5.26001519587677318785587544488E-2;
      a31 =    1.97250569845378994544595329183E-2;
      a32 =    5.91751709536136983633785987549E-2;
      a41 =    2.95875854768068491816892993775E-2;
      a43 =    8.87627564304205475450678981324E-2;
      a51 =    2.41365134159266685502369798665E-1;
      a53 =   -8.84549479328286085344864962717E-1;
      a54 =    9.24834003261792003115737966543E-1;
      a61 =    3.7037037037037037037037037037E-2;
      a64 =    1.70828608729473871279604482173E-1;
      a65 =    1.25467687566822425016691814123E-1;
      a71 =    3.7109375E-2;
      a74 =    1.70252211019544039314978060272E-1;
      a75 =    6.02165389804559606850219397283E-2;
      a76 =   -1.7578125E-2;

      a81 =    3.70920001185047927108779319836E-2;
      a84 =    1.70383925712239993810214054705E-1;
      a85 =    1.07262030446373284651809199168E-1;
      a86 =   -1.53194377486244017527936158236E-2;
      a87 =    8.27378916381402288758473766002E-3;
      a91 =    6.24110958716075717114429577812E-1;
      a94 =   -3.36089262944694129406857109825E0;
      a95 =   -8.68219346841726006818189891453E-1;
      a96 =    2.75920996994467083049415600797E1;
      a97 =    2.01540675504778934086186788979E1;
      a98 =   -4.34898841810699588477366255144E1;
      a101 =   4.77662536438264365890433908527E-1;
      a104 =  -2.48811461997166764192642586468E0;
      a105 =  -5.90290826836842996371446475743E-1;
      a106 =   2.12300514481811942347288949897E1;
      a107 =   1.52792336328824235832596922938E1;
      a108 =  -3.32882109689848629194453265587E1;
      a109 =  -2.03312017085086261358222928593E-2;

      a111 =  -9.3714243008598732571704021658E-1;
      a114 =   5.18637242884406370830023853209E0;
      a115 =   1.09143734899672957818500254654E0;
      a116 =  -8.14978701074692612513997267357E0;
      a117 =  -1.85200656599969598641566180701E1;
      a118 =   2.27394870993505042818970056734E1;
      a119 =   2.49360555267965238987089396762E0;
      a1110 = -3.0467644718982195003823669022E0;
      a121 =   2.27331014751653820792359768449E0;
      a124 =  -1.05344954667372501984066689879E1;
      a125 =  -2.00087205822486249909675718444E0;
      a126 =  -1.79589318631187989172765950534E1;
      a127 =   2.79488845294199600508499808837E1;
      a128 =  -2.85899827713502369474065508674E0;
      a129 =  -8.87285693353062954433549289258E0;
      a1210 =  1.23605671757943030647266201528E1;
      a1211 =  6.43392746015763530355970484046E-1;

      a141 =  5.61675022830479523392909219681E-2;
      a147 =  2.53500210216624811088794765333E-1;
      a148 = -2.46239037470802489917441475441E-1;
      a149 = -1.24191423263816360469010140626E-1;
      a1410 =  1.5329179827876569731206322685E-1;
      a1411 =  8.20105229563468988491666602057E-3;
      a1412 =  7.56789766054569976138603589584E-3;
      a1413 = -8.298E-3;

      a151 =  3.18346481635021405060768473261E-2;
      a156 =  2.83009096723667755288322961402E-2;
      a157 =  5.35419883074385676223797384372E-2;
      a158 = -5.49237485713909884646569340306E-2;
      a1511 = -1.08347328697249322858509316994E-4;
      a1512 =  3.82571090835658412954920192323E-4;
      a1513 = -3.40465008687404560802977114492E-4;
      a1514 =  1.41312443674632500278074618366E-1;
      a161 = -4.28896301583791923408573538692E-1;
      a166 = -4.69762141536116384314449447206E0;
      a167 =  7.68342119606259904184240953878E0;
      a168 =  4.06898981839711007970213554331E0;
      a169 =  3.56727187455281109270669543021E-1;
      a1613 = -1.39902416515901462129418009734E-3;
      a1614 =  2.9475147891527723389556272149E0;
      a1615 = -9.15095847217987001081870187138E0;

      d41  = -0.84289382761090128651353491142E+01;
      d46  =  0.56671495351937776962531783590E+00;
      d47  = -0.30689499459498916912797304727E+01;
      d48  =  0.23846676565120698287728149680E+01;
      d49  =  0.21170345824450282767155149946E+01;
      d410 = -0.87139158377797299206789907490E+00;
      d411 =  0.22404374302607882758541771650E+01;
      d412 =  0.63157877876946881815570249290E+00;
      d413 = -0.88990336451333310820698117400E-01;
      d414 =  0.18148505520854727256656404962E+02;
      d415 = -0.91946323924783554000451984436E+01;
      d416 = -0.44360363875948939664310572000E+01;

      d51  =  0.10427508642579134603413151009E+02;
      d56  =  0.24228349177525818288430175319E+03;
      d57  =  0.16520045171727028198505394887E+03;
      d58  = -0.37454675472269020279518312152E+03;
      d59  = -0.22113666853125306036270938578E+02;
      d510 =  0.77334326684722638389603898808E+01;
      d511 = -0.30674084731089398182061213626E+02;
      d512 = -0.93321305264302278729567221706E+01;
      d513 =  0.15697238121770843886131091075E+02;
      d514 = -0.31139403219565177677282850411E+02;
      d515 = -0.93529243588444783865713862664E+01;
      d516 =  0.35816841486394083752465898540E+02;

      d61 =  0.19985053242002433820987653617E+02;
      d66 = -0.38703730874935176555105901742E+03;
      d67 = -0.18917813819516756882830838328E+03;
      d68 =  0.52780815920542364900561016686E+03;
      d69 = -0.11573902539959630126141871134E+02;
      d610 =  0.68812326946963000169666922661E+01;
      d611 = -0.10006050966910838403183860980E+01;
      d612 =  0.77771377980534432092869265740E+00;
      d613 = -0.27782057523535084065932004339E+01;
      d614 = -0.60196695231264120758267380846E+02;
      d615 =  0.84320405506677161018159903784E+02;
      d616 =  0.11992291136182789328035130030E+02;

      d71  = -0.25693933462703749003312586129E+02;
      d76  = -0.15418974869023643374053993627E+03;
      d77  = -0.23152937917604549567536039109E+03;
      d78  =  0.35763911791061412378285349910E+03;
      d79  =  0.93405324183624310003907691704E+02;
      d710 = -0.37458323136451633156875139351E+02;
      d711 =  0.10409964950896230045147246184E+03;
      d712 =  0.29840293426660503123344363579E+02;
      d713 = -0.43533456590011143754432175058E+02;
      d714 =  0.96324553959188282948394950600E+02;
      d715 = -0.39177261675615439165231486172E+02;
      d716 = -0.14972683625798562581422125276E+03;

  }

  facold = 1.0E-4;
  expo1 = 1.0/8.0 - beta * 0.2;
  facc1 = 1.0 / fac1;
  facc2 = 1.0 / fac2;
  posneg = xend>x ? 1.0 : -1.0;

  /* initial preparations */
  atoli = atoler[0];
  rtoli = rtoler[0];
  last  = 0;
  hlamb = 0.0;
  iasti = 0;
  nonsti= 0;
  fcn (n, x, y, k1);
  hmax = fabs (hmax);
  iord = 8;
  if (h == 0.0)
    h = hinit (n, /*fcn,*/ x, y, posneg, k1, k2, k3, iord, hmax, atoler, rtoler, itoler);
  nfcn += 2;
  reject = 0;
  xold = x;
  
  if (iout)
  {
    irtrn = 1;
    hout = 1.0;
    xout = x;
    if(x>xold) solout (naccpt+1, xold, x, y, n); 
  }

  /* basic integration step */
  while (1)
  {
    if (state==OS_NEEDTOTERMINATE)   /* !!! alteration of original code */
    {
      xout = x;
      hout = h;
      return -42;
    }

    if (nstep > nmax)
    {
      /*if (fileout)
        fprintf (fileout, "Exit of dop853 at x = %.16e, more than nmax = %li are needed\r\n", x, nmax);*/
      xout = x;
      hout = h;
      return -2;
    }

    if (!gsl_finite(h) || 0.1 * fabs(h) <= fabs(x) * uround)
    {
      /*if (fileout)
        fprintf (fileout, "Exit of dop853 at x = %.16e, step size too small h = %.16e\r\n", x, h);*/
      xout = x;
      hout = h;
      return -3;
    }

    if ((x + 1.01*h - xend) * posneg > 0.0)
    {
      h = xend - x;
      last = 1;
    }

    nstep++;

    /* the twelve stages */
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * a21 * k1[i];
    fcn (n, x+c2*h, yy1, k2);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a31*k1[i] + a32*k2[i]);
    fcn (n, x+c3*h, yy1, k3);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a41*k1[i] + a43*k3[i]);
    fcn (n, x+c4*h, yy1, k4);
    for (i = 0; i <n; i++)
      yy1[i] = y[i] + h * (a51*k1[i] + a53*k3[i] + a54*k4[i]);
    fcn (n, x+c5*h, yy1, k5);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a61*k1[i] + a64*k4[i] + a65*k5[i]);
    fcn (n, x+c6*h, yy1, k6);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a71*k1[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
    fcn (n, x+c7*h, yy1, k7);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a81*k1[i] + a84*k4[i] + a85*k5[i] + a86*k6[i] +
                          a87*k7[i]);
    fcn (n, x+c8*h, yy1, k8);
    for (i = 0; i <n; i++)
      yy1[i] = y[i] + h * (a91*k1[i] + a94*k4[i] + a95*k5[i] + a96*k6[i] +
                          a97*k7[i] + a98*k8[i]);
    fcn (n, x+c9*h, yy1, k9);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a101*k1[i] + a104*k4[i] + a105*k5[i] + a106*k6[i] +
                          a107*k7[i] + a108*k8[i] + a109*k9[i]);
    fcn (n, x+c10*h, yy1, k10);
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a111*k1[i] + a114*k4[i] + a115*k5[i] + a116*k6[i] +
                          a117*k7[i] + a118*k8[i] + a119*k9[i] + a1110*k10[i]);
    fcn (n, x+c11*h, yy1, k2);
    xph = x + h;
    for (i = 0; i < n; i++)
      yy1[i] = y[i] + h * (a121*k1[i] + a124*k4[i] + a125*k5[i] + a126*k6[i] +
                          a127*k7[i] + a128*k8[i] + a129*k9[i] +
                          a1210*k10[i] + a1211*k2[i]);
    fcn (n, xph, yy1, k3);
    nfcn += 11;
    for (i = 0; i < n; i++)
    {
      k4[i] = b1*k1[i] + b6*k6[i] + b7*k7[i] + b8*k8[i] + b9*k9[i] +
              b10*k10[i] + b11*k2[i] + b12*k3[i];
      k5[i] = y[i] + h * k4[i];
    }
     
    /* error estimation */
    err = 0.0;
    err2 = 0.0;
    if (!itoler)
    {
      /*   !!!! alteration of original dop853 code !!!!   */
      rtoli=rtoler[0]; atoli=atoler[0];
      if(potential->Mbh!=0)
      {
          double r=sqrt(pow_2(y[0])+pow_2(y[1])+(N_dim==3?pow_2(y[2]):0.0));
          if(r < potential->Mbh*BH_RADIUS_INCREASE_ACCURACY)   // take measures to increase accuracy
          {
              rtoli=rtoler[0]*(r/(potential->Mbh*BH_RADIUS_INCREASE_ACCURACY));
              atoli=atoler[0]*(r/(potential->Mbh*BH_RADIUS_INCREASE_ACCURACY));
          }
      }
      for (i = 0; i < n; i++)
      {
        sk = atoli + rtoli * std::max<double> (fabs(y[i]), fabs(k5[i]));
        erri = k4[i] - bhh1*k1[i] - bhh2*k9[i] - bhh3*k3[i];
        sqre = erri / sk;
        err2 += sqre*sqre;
        erri = er1*k1[i] + er6*k6[i] + er7*k7[i] + er8*k8[i] + er9*k9[i] +
               er10 * k10[i] + er11*k2[i] + er12*k3[i];
        sqre = erri / sk;
        err += sqre*sqre;
      }
    }
    else
      for (i = 0; i < n; i++)
      {
        sk = atoler[i] + rtoler[i] * std::max<double> (fabs(y[i]), fabs(k5[i]));
        erri = k4[i] - bhh1*k1[i] - bhh2*k9[i] - bhh3*k3[i];
        sqre = erri / sk;
        err2 += sqre*sqre;
        erri = er1*k1[i] + er6*k6[i] + er7*k7[i] + er8*k8[i] + er9*k9[i] +
               er10 * k10[i] + er11*k2[i] + er12*k3[i];
        sqre = erri / sk;
        err += sqre*sqre;
      }
    deno = err + 0.01 * err2;
    if (deno <= 0.0)
      deno = 1.0;
    err = fabs(h) * err * sqrt (1.0 / (deno*(double)n));

    /* computation of hnew */
    fac11 = pow (err, expo1);
    /* Lund-stabilization */
    fac = fac11 / pow(facold,beta);
    /* we require fac1 <= hnew/h <= fac2 */
    fac = std::max<double> (facc2, std::min<double> (facc1, fac/safe));
    hnew = h / fac;

    if (err <= 1.0)
    {
      /* step accepted */

      facold = std::max<double> (err, 1.0E-4);
      naccpt++;
      fcn (n, xph, k5, k4);
      nfcn++;
      
      /* stiffness detection */
      if (!(naccpt % nstiff) || (iasti > 0))
      {
        stnum = 0.0;
        stden = 0.0;
        for (i = 0; i < n; i++)
        {
          sqre = k4[i] - k3[i];
          stnum += sqre*sqre;
          sqre = k5[i] - yy1[i];
          stden += sqre*sqre;
        }
        if (stden > 0.0)
          hlamb = h * sqrt (stnum / stden);
        if (hlamb > 6.1)
        {
          nonsti = 0;
          iasti++;
          if (iasti == 15)
            /*if (fileout)
              fprintf (fileout, "The problem seems to become stiff at x = %.16e\r\n", x);
            else*/
            {
              xout = x;
              hout = h;
              return -4;
            }
        }
        else
        {
          nonsti++;
          if (nonsti == 6)
            iasti = 0;
        }
      }
       
      /* final preparation for dense output */
      if (iout == 2)
      {
        /* save the first function evaluations */
        if (nrds == n)
          for (i = 0; i < n; i++)
          {
            rcont1[i] = y[i];
            ydiff = k5[i] - y[i];
            rcont2[i] = ydiff;
            bspl = h * k1[i] - ydiff;
            rcont3[i] = bspl;
            rcont4[i] = ydiff - h*k4[i] - bspl;
            rcont5[i] = d41*k1[i] + d46*k6[i] + d47*k7[i] + d48*k8[i] +
                        d49*k9[i] + d410*k10[i] + d411*k2[i] + d412*k3[i];
            rcont6[i] = d51*k1[i] + d56*k6[i] + d57*k7[i] + d58*k8[i] +
                        d59*k9[i] + d510*k10[i] + d511*k2[i] + d512*k3[i];
            rcont7[i] = d61*k1[i] + d66*k6[i] + d67*k7[i] + d68*k8[i] +
                        d69*k9[i] + d610*k10[i] + d611*k2[i] + d612*k3[i];
            rcont8[i] = d71*k1[i] + d76*k6[i] + d77*k7[i] + d78*k8[i] +
                        d79*k9[i] + d710*k10[i] + d711*k2[i] + d712*k3[i];
          }
        else
          for (j = 0; j < nrds; j++)
          {
            i = icont[j];
            rcont1[j] = y[i];
            ydiff = k5[i] - y[i];
            rcont2[j] = ydiff;
            bspl = h * k1[i] - ydiff;
            rcont3[j] = bspl;
            rcont4[j] = ydiff - h*k4[i] - bspl;
            rcont5[j] = d41*k1[i] + d46*k6[i] + d47*k7[i] + d48*k8[i] +
                        d49*k9[i] + d410*k10[i] + d411*k2[i] + d412*k3[i];
            rcont6[j] = d51*k1[i] + d56*k6[i] + d57*k7[i] + d58*k8[i] +
                        d59*k9[i] + d510*k10[i] + d511*k2[i] + d512*k3[i];
            rcont7[j] = d61*k1[i] + d66*k6[i] + d67*k7[i] + d68*k8[i] +
                        d69*k9[i] + d610*k10[i] + d611*k2[i] + d612*k3[i];
            rcont8[j] = d71*k1[i] + d76*k6[i] + d77*k7[i] + d78*k8[i] +
                        d79*k9[i] + d710*k10[i] + d711*k2[i] + d712*k3[i];
          }

        /* the next three function evaluations */
        for (i = 0; i < n; i++)
          yy1[i] = y[i] + h * (a141*k1[i] + a147*k7[i] + a148*k8[i] +
                              a149*k9[i] + a1410*k10[i] + a1411*k2[i] +
                              a1412*k3[i] + a1413*k4[i]);
        fcn (n, x+c14*h, yy1, k10);
        for (i = 0; i < n; i++)
          yy1[i] = y[i] + h * (a151*k1[i] + a156*k6[i] + a157*k7[i] + a158*k8[i] +
                              a1511*k2[i] + a1512*k3[i] + a1513*k4[i] +
                              a1514*k10[i]);
        fcn (n, x+c15*h, yy1, k2);
        for (i = 0; i < n; i++)
          yy1[i] = y[i] + h * (a161*k1[i] + a166*k6[i] + a167*k7[i] + a168*k8[i] +
                              a169*k9[i] + a1613*k4[i] + a1614*k10[i] +
                              a1615*k2[i]);
        fcn (n, x+c16*h, yy1, k3);
        nfcn += 3;

        /* final preparation */
        if (nrds == n)
          for (i = 0; i < n; i++)
          {
            rcont5[i] = h * (rcont5[i] + d413*k4[i] + d414*k10[i] +
                             d415*k2[i] + d416*k3[i]);
            rcont6[i] = h * (rcont6[i] + d513*k4[i] + d514*k10[i] +
                             d515*k2[i] + d516*k3[i]);
            rcont7[i] = h * (rcont7[i] + d613*k4[i] + d614*k10[i] +
                             d615*k2[i] + d616*k3[i]);
            rcont8[i] = h * (rcont8[i] + d713*k4[i] + d714*k10[i] +
                             d715*k2[i] + d716*k3[i]);
          }
        else
          for (j = 0; j < nrds; j++)
          {
            i = icont[j];
            rcont5[j] = h * (rcont5[j] + d413*k4[i] + d414*k10[i] +
                             d415*k2[i] + d416*k3[i]);
            rcont6[j] = h * (rcont6[j] + d513*k4[i] + d514*k10[i] +
                             d515*k2[i] + d516*k3[i]);
            rcont7[j] = h * (rcont7[j] + d613*k4[i] + d614*k10[i] +
                             d615*k2[i] + d616*k3[i]);
            rcont8[j] = h * (rcont8[j] + d713*k4[i] + d714*k10[i] +
                             d715*k2[i] + d716*k3[i]);
          }
      }

      //memcpy (k1, k4, n * sizeof(double)); 
      //memcpy (y, k5, n * sizeof(double));
      for(unsigned int indx=0; indx<n; indx++) k1[indx]=k4[indx]; 
      for(unsigned int indx=0; indx<n; indx++) y[indx]=k5[indx]; 
      xold = x;
      x = xph;

      if (iout)
      {
        hout = h;
        xout = x;
        solout (naccpt+1, xold, x, y, n);
      }

      /* normal exit */
      if (last)
      {
        hout=hnew;
        xout = x;
        return 1;
      }

      if (fabs(hnew) > hmax)
        hnew = posneg * hmax;
      if (reject)
        hnew = posneg * std::min<double> (fabs(hnew), fabs(h));

      reject = 0;

      /*   !!!! alteration of original dop853 code !!!!   */
      if(InitData.calcLyapunov)
      {
          // check if deviation vector requires renormalization
          double wnorm=0;
#ifdef LYAPUNOVVAREQ
          for(int indx=N_dim*2; indx<N_dim*4; indx++) wnorm+=pow_2(y[indx]);
          wnorm=sqrt(wnorm);
          if(wnorm>LYAP_VAREQ_DEV_MAX)
          {
            lnwadd+=log(wnorm);
            for(int indx=N_dim*2; indx<N_dim*4; indx++)
            {
                y[indx]/=wnorm;
                k1[indx]/=wnorm;
                k2[indx]/=wnorm;
                k3[indx]/=wnorm;
                k4[indx]/=wnorm;
                k5[indx]/=wnorm;
                k6[indx]/=wnorm;
                k7[indx]/=wnorm;
                k8[indx]/=wnorm;
                k9[indx]/=wnorm;
                k10[indx]/=wnorm;
            }
          }
#else
          for(unsigned int indx=N_dim*2; indx<N_dim*4; indx++)
              wnorm+=pow_2(y[indx-N_dim*2]-y[indx]); 
          wnorm=sqrt(wnorm)/LYAP_DEV_INIT;
          if(wnorm>LYAP_DEV_MAX/LYAP_DEV_INIT)
          {
            needrenorm=true;   // it will be performed later, when allowed
          }
          if(needrenorm && allowrenorm)
          {
            // ensure that we do not decrease normalization
            if(wnorm<2) wnorm=2;
            lnwadd+=log(wnorm);
            for(unsigned int indx=N_dim*2; indx<N_dim*4; indx++)
            {
                y[indx] = y[indx-N_dim*2] + (y[indx]-y[indx-N_dim*2])/wnorm;
                k1[indx] = k1[indx-N_dim*2] + (k1[indx]-k1[indx-N_dim*2])/wnorm;
                k2[indx] = k2[indx-N_dim*2] + (k2[indx]-k2[indx-N_dim*2])/wnorm;
                k3[indx] = k3[indx-N_dim*2] + (k3[indx]-k3[indx-N_dim*2])/wnorm;
                k4[indx] = k4[indx-N_dim*2] + (k4[indx]-k4[indx-N_dim*2])/wnorm;
                k5[indx] = k5[indx-N_dim*2] + (k5[indx]-k5[indx-N_dim*2])/wnorm;
                k6[indx] = k6[indx-N_dim*2] + (k6[indx]-k6[indx-N_dim*2])/wnorm;
                k7[indx] = k7[indx-N_dim*2] + (k7[indx]-k7[indx-N_dim*2])/wnorm;
                k8[indx] = k8[indx-N_dim*2] + (k8[indx]-k8[indx-N_dim*2])/wnorm;
                k9[indx] = k9[indx-N_dim*2] + (k9[indx]-k9[indx-N_dim*2])/wnorm;
                k10[indx] = k10[indx-N_dim*2] + (k10[indx]-k10[indx-N_dim*2])/wnorm;
            }
            needrenorm = allowrenorm = false;
          }
#endif
      }
    }
    else
    {
      /* step rejected */
      hnew = h / std::min<double> (facc1, fac11/safe);
      reject = 1;
      if (naccpt >= 1)
        nrejct=nrejct + 1;
      last = 0;
    }

    h = hnew;
  }

} /* dopcor */


/* front-end */
int COrbit::dop853
 (unsigned n, /*FcnEqDiff fcn,*/ double x, double* y, double xend, double* rtoler,
  double* atoler, int itoler, /*SolTrait solout,*/ int iout, /*FILE* fileout,*/ double uround,
  double safe, double fac1, double fac2, double beta, double hmax, double h,
  long nmax, int meth, long nstiff, unsigned nrdens, unsigned* icont, unsigned licont)
{
  if(leapfrog)
    return intlf(n, x, y, xend);
  int       arret, idid;
  unsigned  i;

  /* initialisations */
  nfcn = nstep = naccpt = nrejct = arret = 0;

  /* n, the dimension of the system */
  if (n == UINT_MAX)
  {
    /*if (fileout)
      fprintf (fileout, "System too big, max. n = %u\r\n", UINT_MAX-1);*/
    arret = 1;
  }

  /* nmax, the maximal number of steps */
  if (!nmax)
    nmax = NUMSTEP_MAX;
  else if (nmax <= 0)
  {
    /*if (fileout)
      fprintf (fileout, "Wrong input, nmax = %li\r\n", nmax);*/
    arret = 1;
  }

  /* meth, coefficients of the method */
  if (!meth)
    meth = 1;
  else if ((meth <= 0) || (meth >= 2))
  {
    /*if (fileout)
      fprintf (fileout, "Curious input, meth = %i\r\n", meth);*/
    arret = 1;
  }

  /* nstiff, parameter for stiffness detection */
  if (!nstiff)
    nstiff = 1000;
  else if (nstiff < 0)
    nstiff = nmax + 10;

  /* iout, switch for calling solout */
  if ((iout < 0) || (iout > 2))
  {
    /*if (fileout)
      fprintf (fileout, "Wrong input, iout = %i\r\n", iout);*/
    arret = 1;
  }

  std::vector<unsigned int> vec_indir;
  std::vector<double> vec_rcont[8];
  indir=NULL;
  rcont1=rcont2=rcont3=rcont4=rcont5=rcont6=rcont7=rcont8=NULL;
  /* nrdens, number of dense output components */
  if (nrdens > n)
  {
    /*if (fileout)
      fprintf (fileout, "Curious input, nrdens = %u\r\n", nrdens);*/
    arret = 1;
  }
  else if (nrdens)
  {
    vec_rcont[0].resize(nrdens);  rcont1=&(vec_rcont[0].front());
    vec_rcont[1].resize(nrdens);  rcont2=&(vec_rcont[1].front());
    vec_rcont[2].resize(nrdens);  rcont3=&(vec_rcont[2].front());
    vec_rcont[3].resize(nrdens);  rcont4=&(vec_rcont[3].front());
    vec_rcont[4].resize(nrdens);  rcont5=&(vec_rcont[4].front());
    vec_rcont[5].resize(nrdens);  rcont6=&(vec_rcont[5].front());
    vec_rcont[6].resize(nrdens);  rcont7=&(vec_rcont[6].front());
    vec_rcont[7].resize(nrdens);  rcont8=&(vec_rcont[7].front());
    if (nrdens < n)
    {
      vec_indir.resize(nrdens);
      indir = &(vec_indir.front());
    }
    /* control of length of icont */
    if (nrdens == n)
    {
      /*if (icont && fileout)
        fprintf (fileout, "Warning : when nrdens = n there is no need allocating memory for icont\r\n");*/
      nrds = n;
    }
    else if (licont < nrdens)
    {
      /*if (fileout)
        fprintf (fileout, "Insufficient storage for icont, min. licont = %u\r\n", nrdens);*/
      arret = 1;
    }
    else
    {
      /*if ((iout < 2) && fileout)
        fprintf (fileout, "Warning : put iout = 2 for dense output\r\n");*/
      nrds = nrdens;
      for (i = 0; i < n; i++)
        indir[i] = UINT_MAX;
      for (i = 0; i < nrdens; i++)
        indir[icont[i]] = i;
    }
  }
  else nrds=0;

  /* uround, smallest number satisfying 1.0+uround > 1.0 */
  if (uround == 0.0)
    uround = 2.3E-16;
  else if ((uround <= 1.0E-35) || (uround >= 1.0))
  {
    /*if (fileout)
      fprintf (fileout, "Which machine do you have ? Your uround was : %.16e\r\n", uround);*/
    arret = 1;
  }

  /* safety factor */
  if (safe == 0.0)
    safe = 0.9;
  else if ((safe >= 1.0) || (safe <= 1.0E-4))
  {
    /*if (fileout)
      fprintf (fileout, "Curious input for safety factor, safe = %.16e\r\n", safe);*/
    arret = 1;
  }

  /* fac1, fac2, parameters for step size selection */
  if (fac1 == 0.0)
    fac1 = 0.333;
  if (fac2 == 0.0)
    fac2 = 6.0;

  /* beta for step control stabilization */
  if (beta == 0.0)
    beta = 0.0;
  else if (beta < 0.0)
    beta = 0.0;
  else if (beta > 0.2)
  {
    /*if (fileout)
      fprintf (fileout, "Curious input for beta : beta = %.16e\r\n", beta);*/
    arret = 1;
  }

  /* maximal step size */
  if (hmax == 0.0)
    hmax = xend - x;


  /* when a failure has occured, we return -1 */
  if (arret)
  {
    return -1;
  }
  else
  {
    idid = dopcor (n, /*fcn,*/ x, y, xend, hmax, h, rtoler, atoler, itoler, /*fileout,
                   solout,*/ iout, nmax, uround, meth, nstiff, safe, beta, fac1, fac2, icont);
    return idid;
  }

} /* dop853 */


/* dense output function */
double COrbit::contd8 (unsigned ii, double x) const
{
  unsigned i;
  double   s, s1;

  i = UINT_MAX;

  if (!indir)
    i = ii;
  else
    i = indir[ii];

  if (i == UINT_MAX)
  {
    //printf ("No dense output available for %uth component", ii);
    return 0.0;
  }

  s = (x - xold) / hout;
  s1 = 1.0 - s;

  return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*(rcont5[i]+
         s*(rcont6[i]+s1*(rcont7[i]+s*rcont8[i]))))));
} /* contd8 */

//----------------------------------------------------------------------------//
// Runtime function for recording and analyzing trajectory by means of frequency analysis and Lyapunov exponent 

COrbitRuntimeTrajectory::COrbitRuntimeTrajectory(const smile::COrbit *_orbit) :
    CBasicOrbitRuntimeFnc(_orbit), InitData(orbit->getInitData()), potential(InitData.potential), N_dim(potential->N_dim)
{
    Einit = potential->Phi(InitData.initCond.Pos[0], InitData.initCond.Pos[1], InitData.initCond.Pos[2]) + 
        (pow_2(InitData.initCond.Vel[0])+pow_2(InitData.initCond.Vel[1])+pow_2(InitData.initCond.Vel[2]))/2.0;
    Ediff=0;
}

template<typename NumT> COrbitRuntimeTrajectory::COrbitRuntimeTrajectory(const smile::COrbit *_orbit, const std::vector< CPosVelPoint<NumT> > &_traj) :
    CBasicOrbitRuntimeFnc(_orbit), InitData(orbit->getInitData()), potential(InitData.potential), N_dim(potential->N_dim), traj(_traj.begin(), _traj.end())
{
    Einit = potential->Phi(InitData.initCond.Pos[0], InitData.initCond.Pos[1], InitData.initCond.Pos[2]) + 
        (pow_2(InitData.initCond.Vel[0])+pow_2(InitData.initCond.Vel[1])+pow_2(InitData.initCond.Vel[2]))/2.0;
    Ediff=0;  // not quite true, but we won't bother calculating it..
}

void COrbitRuntimeTrajectory::Timestep(const double told, const double t, const double y[])
{
    int i1=(int)(told/InitData.timeStep); 
    int i2=(int)(t/InitData.timeStep);

    // check energy conservation
    if(i2>i1 && i2%10==0)  // do not check too frequently (potential evaluation may be costly)
    {
        double Z=(N_dim==3)?y[2]:0;
        double vZ=(N_dim==3)?y[N_dim+2]:0;
        double Ecur = potential->Phi(y[0], y[1], Z) + (pow_2(y[N_dim]) + pow_2(y[N_dim+1]) + pow_2(vZ))/2;
        double Edifcur = fabs(Einit - Ecur);
        Ediff=std::max<double>(Edifcur, Ediff);
    }
    // store trajectory at regular intervals of time    
    for (int i=i1; i<=i2; i++)
    {
        double tt=InitData.timeStep*i;
        if ((tt>=told) && (tt<t))
        {
            CPosVelPoint<double> pt;
            for(unsigned int c=0; c<N_dim; c++)
            {
                pt.Pos[c]=orbit->getInterpolatedTrajectory(c, tt);
                pt.Vel[c]=orbit->getInterpolatedTrajectory(c+N_dim, tt);
            }
            traj.push_back(pt);
            if(InitData.calcLyapunov)
                lnw.push_back(orbit->getInterpolatedTrajectory(2*N_dim, tt));
        }
    }
}

template<typename NumT> COrbitInformation<NumT>* COrbitRuntimeTrajectory::analyzeOrbit(size_t start, size_t count, CSpectrum *spectrum, vectord *fa) const
{
    if(start>=traj.size())  // invalid call
        return new COrbitInformation<NumT>("Invalid call: start>=size");
    if(count<=0 || count>traj.size()-start) count = traj.size() - start;   // passing zero for count means analysing entire orbit
    if(count<=1)  // empty orbit
        return new COrbitInformation<NumT>("?");
    CSpectrum tempSpectrum[N_DIM];
    CSpectrum* sl = spectrum!=NULL ? spectrum : tempSpectrum;   // if outputting spectrum is not necessary, use temp. array for storing the information
    double Lavg_d[N_DIM]={0,0,0}, Lvar_d[N_DIM]={0,0,0};
    std::string Description = performFrequencyAnalysis(start, count, sl, fa, Lavg_d, Lvar_d);
    double lfccdiff = calcLFCCdiffusion(start, count);
    double lambda = InitData.calcLyapunov ? calcLambda() : -1;
    double inertia_d[3];
    calcInertiaOverInterval(start, count, inertia_d);
    NumT leadfreq[N_DIM]={static_cast<NumT>(sl[0][0]), static_cast<NumT>(sl[1][0]), static_cast<NumT>((N_dim==3)?sl[2][0].freq:0)};
    NumT inertia[N_DIM] ={static_cast<NumT>(inertia_d[0]), static_cast<NumT>(inertia_d[1]), static_cast<NumT>(inertia_d[2])};
    NumT Lavg[N_DIM] ={static_cast<NumT>(Lavg_d[0]), static_cast<NumT>(Lavg_d[1]), static_cast<NumT>(Lavg_d[2])};
    NumT Lvar[N_DIM] ={static_cast<NumT>(Lvar_d[0]), static_cast<NumT>(Lvar_d[1]), static_cast<NumT>(Lvar_d[2])};
    // note that Ediff is returned for entire orbit integration time, not for the given interval
    return new COrbitInformation<NumT>(Description, static_cast<NumT>(Einit), static_cast<NumT>(Ediff), leadfreq, static_cast<NumT>(lfccdiff), static_cast<NumT>(lambda), inertia, Lavg, Lvar);
}
// explicit instantiations of a template function
template COrbitInformation<double>* COrbitRuntimeTrajectory::analyzeOrbit<double>(size_t start, size_t count, CSpectrum *spectrum, vectord *fa) const;

CBasicInformation* COrbitRuntimeTrajectory::getData() const
{
    return analyzeOrbit<float>(0, 0);
}

// simple routine that tests if an integer is prime (in the FFT primes are great evil)
void tuneNonPrime(size_t &n)
{
    size_t m=n, f=2;
    while(f*f<=m)
        if(m%f==0) m/=f; else f++;
    if(m>=MAX_FFT_PRIME && n>=MAX_FFT_PRIME)
    {
        n--;
        tuneNonPrime(n);
    }
}

// perform FFT of orbit (output amplitudes to 'fa' array) and find 'maxnlines' lines, which will be placed into 'sl' array. 
// If fa=NULL, do not output frequency amplitudes
void COrbitRuntimeTrajectory::findLines(size_t start, size_t npoints, CSpectrum *sl, vectord *fa) const
{
    if(npoints==0) return;
    tuneNonPrime(npoints);
    vectord in(npoints);
    gsl_fft_real_wavetable * wt = gsl_fft_real_wavetable_alloc (npoints);
    gsl_fft_real_workspace * ws = gsl_fft_real_workspace_alloc (npoints);
    std::vector<complexd> fc(npoints/2+1);  // temp.array of complex amplitudes

    for(unsigned int c=0; c<N_dim; c++)
    {
        for(size_t i=0; i<npoints; i++)
            in[i]=traj[start+i].Pos[c];
        gsl_fft_real_transform (&(in.front()), 1, npoints, wt, ws);
        fc[0]=complexd(in[0],0)  / (0.5*npoints);
        for(size_t i = 1; i < npoints/2; i++)
            fc[i]=complexd(in[2*i-1], in[2 * i])/(0.5*npoints);
        if(npoints%2==0)
            fc[npoints/2]=complexd(in[npoints-1], 0)/(0.5*npoints);
        if(fa!=NULL)
        {  // reserve space in vector
            fa[c].resize(npoints/2+1, 0);
            for(size_t i=0; i<npoints/2+1; i++)
                fa[c][i]=abs(fc[i]);    // store amplitude of spectrum if needed (to display only) - and then forget it; only the true complex spectrum is used hereafter
        }

        // perform frequency analysis
        sl[c].clear();
        bool finish=false;
        double maxamp=0; double ma=0;
        while(!finish)
        {   // find next greatest line
            int m=0; 
            ma=0;
            for(int i=0; i<static_cast<int>(npoints/2)-2; i++)
                if(abs(fc[i])>ma)
                {
                    ma=abs(fc[i]);
                    m=i;
                }
            if(maxamp==0 && m>0) maxamp=ma;  // record max.amplitude if the frequency is non-zero
            if(ma>maxamp*MIN_SPECTRAL_LINE_AMPLITUDE)
            {   // add line
                // peak is somewhere near m.  
                if(m>1)
                {
                    int m1=m;
                    if(abs(fc[m-1])>abs(fc[m+1]))  m1--;
                    //int m2=((abs(fc[m-1])>abs(fc[m+1]))?(m-1):(m+1));
                    //if(m2<m1) { int mm=m2; m2=m1; m1=mm; }   // not necessary. but convenient to have m2=m1+1
                    double s;   // fractional index of peak frequency 
                    complexd Ac;// amplitude
                    complexd f1 = (fc[m1] - (fc[m1-1]+fc[m1+1])/2.0);
                    complexd f2 = (fc[m1+1] - (fc[m1]+fc[m1+2])/2.0);
                    // C=f2/f1; check out that C<>0 and C<>infinity
                    if(abs(f1) < abs(f2)*1e-4)
                    {   
                        m1++;
                        f1=f2;
                        f2 = (fc[m1+1] - (fc[m1]+fc[m1+2])/2.0);
                    }
                    if((abs(f2) < abs(f1)*1e-4) && (m1>0))
                    {   
                        m1--;
                        f2=f1;
                        f1 = (fc[m1] - (fc[m1-1]+fc[m1+1])/2.0);
                    }
                    complexd C = f2 / f1;
                    double Ca=C.real();  // in ideal case, Im C=0
                    double omega=2*M_PI/npoints;
                    bool Precise=true;
                    if((Ca<-1.0/FREQ_PRECISE_THRESHOLD) && (Ca>-FREQ_PRECISE_THRESHOLD)/*abs(C.imag())*0 < abs(Ca)*/)
                    {   // try Hunter's precise method
                        double sgn = (Ca<1)?+1:-1;
                        double t = (Ca*cos(omega)-1 + sgn * sqrt( pow_2(Ca*cos(omega)-1) - (1+2*Ca)*pow_2(sin(omega))))/sin(omega);
                        s = m1 + npoints/M_PI * atan(t);
                        complexd z = std::polar(1.0, 2*M_PI*(s-m1)/npoints); //complexd((1-pow_2(t))/(1+pow_2(t)), 2*t/(1+pow_2(t)));
                        complexd beta = std::polar(1.0, omega);
                        complexd P = (std::polar(1.0, 2*M_PI*s)-1.0) / (std::polar(1.0, 2*M_PI*(s-m1)/npoints)-1.0) / (1.0*npoints);
                        complexd g1 = P * (1.0 - 0.5*(z-1.0)*(1.0/(z*beta-1.0) + 1.0/(z*std::conj(beta)-1.0)));
                        Ac = f1/g1;
                    }
                    else
                    {   // Follow the method by Carpintero&Aguilar(1998)
                        Precise=false;
                        double alpha = atan(sin(omega/2)/(cos(omega/2)+abs(fc[m1]/fc[m1+1])));
                        s = m1 + npoints*alpha/M_PI;
                        complexd P=(std::polar(1.0, 2*alpha*npoints)-1.0) / (std::polar(1.0, 2*alpha)-1.0) / (1.0*npoints);
                        complexd Q=(std::polar(1.0, -2*M_PI*(s+m1))-1.0) / (std::polar(1.0, -2*M_PI*(s+m1)/npoints)-1.0) / (1.0*npoints);
                        Ac=fc[m1] / P;
                    }
                    sl[c].push_back(CSpectralLine(s*InitData.timeUnit/(npoints*InitData.timeStep)/*intTime*/, abs(Ac), std::arg(Ac), Precise));

                    if(sl[c].size()>=MAX_SPECTRAL_LINES) finish=true;    // not to iterate forever

                    // subtract line from spectrum
                    if(!finish)
                    {
                        if(s-floor(s) < 1e-5)    // frequency was very close to integer value, so subtract only this point
                        {
                            fc[(int)floor(s)]-=Ac;
                        } else
                        if(s-floor(s) > 1-1e-5)
                        {
                            fc[(int)floor(s)+1]-=Ac;
                        }
                        else   // subtract this line with its wings from the whole spectrum
                        {
                            complexd Amul = Ac / (1.0*npoints) * (std::polar(1.0, 2*M_PI*(s))-1.0);
                            complexd Acmul = std::conj(Ac) / (1.0*npoints) * (std::polar(1.0,-2*M_PI*(s))-1.0);
                            int jmin=std::max<int>(0, m-1000);
                            int jmax=std::min<int>(static_cast<int>(npoints/2)+1, m+1000);
                            for(int j=jmin; j<jmax; j++)
                            {
                                complexd linea = Amul / (std::polar(1.0, 2*M_PI*(s-j)/npoints)-1.0) +
                                    Acmul / (std::polar(1.0,-2*M_PI*(s+j)/npoints)-1.0);
                                fc[j] -= linea;
                            }
                        }
                    }
                }
                else
                {   // zero frequency -- just nullify one array element (do not record zero frequencies)
                    fc[m]=0;
                }
            }
            else finish=true;
        }
        std::sort(sl[c].begin(), sl[c].end());
        if(sl[c].size()==0)
            sl[c].push_back(CSpectralLine());  // zero
        int sw=0; CSpectralLine ftmp;
        switch(c)   // ensure that the first line in each coordinate has frequency >1 (for X), >X (for Y), and >Y (for Z); tolerate small deviations given in MIN_FREQ_*
        {
        case 0:
            if(configOrbit.freqStrictOrder) 
                while((sw<(int)sl[0].size()-1) && (sl[0][0].freq<configOrbit.minFreqX)) 
                { 
                    ftmp=sl[0][sw+1]; 
                    sl[0][sw+1]=sl[0][0]; 
                    sl[0][0]=ftmp; sw++; 
                }
            break;
        case 1:
            if(configOrbit.freqStrictOrder) 
                while((sw<(int)sl[1].size()-1) && (sl[1][0].freq<configOrbit.minFreqYZ*sl[0][0].freq)) 
                { 
                    ftmp=sl[1][sw+1]; 
                    sl[1][sw+1]=sl[1][0]; 
                    sl[1][0]=ftmp; sw++; 
                }
            break;
        case 2:
            if(configOrbit.freqStrictOrder && N_dim==3) 
                while((sw<(int)sl[2].size()-1) && (sl[2][0].freq<configOrbit.minFreqYZ*sl[1][0].freq)) 
                { 
                    ftmp=sl[2][sw+1]; 
                    sl[2][sw+1]=sl[2][0]; 
                    sl[2][0]=ftmp; sw++; 
                }
            break;
        }
    }
    gsl_fft_real_wavetable_free (wt);
    gsl_fft_real_workspace_free (ws);
}

void simplify(int &n1, int &n2)    // remove common factors
{
    int nmin=std::min<int>(n1, n2);
    int f=1;
    for(int k=1; k<nmin/2; k++)
        if((n1%k==0) && (n2%k==0)) f=k;
    n1/=f;
    n2/=f;
}

bool commensurable(double a, double b, double acc, int &na, int &nb)   // check that a:b ~=~ na:nb
{
    if(a*b==0) return(false);
    na=1; nb=0;
    while(na<=FREQ_RES_2)
    {
        int nb1=(int)floor(b/a*na+.5);
        if(nb1 && (fabs(nb1*a-na*b) < acc))
        {
            nb=nb1;
            return true;
        }
        else na++;
    }
    return false;
}

bool lineardepend(double f, double a, double b, double acc, int &ma, int &mb)   // check linear dependency: f ~=~ ma*a + mb*b
{
    ma=mb=0; 
    while(abs(ma)<=FREQ_RES_2)
    {
        mb=(int)floor((f-ma*a)/b+0.5);
        if(fabs(ma*a+mb*b-f) < acc) 
            return true;
        else 
            if(ma>0) ma=-ma; else ma=-ma+1;
    }
    return false;
}

bool lineardepend3(double f, double a, double b, double c, double acc, int &ma, int &mb, int &mc)   // check 3d linear dependency: f ~=~ ma*a + mb*b + mc*c
{
    ma=mb=mc=0; 
    while(abs(ma)<=FREQ_RES_2)
    {
      mb=0;
      while(abs(mb)<=FREQ_RES_2)
      {
        mc=(int)floor((f-ma*a-mb*b)/c+0.5);
        if(fabs(ma*a+mb*b+mc*c-f) < acc) 
            return true;
        else 
            if(mb>0) mb=-mb; else mb=-mb+1;
      }
      if(ma>0) ma=-ma; else ma=-ma+1;
    }
    return false;
}

bool commensurable3(double a, double b, double c, double acc, int &ma, int &mb, int &mc)
{
    for(int na=1; na<=FREQ_RES_3; na++)
        for(int nb=-FREQ_RES_3; nb<=FREQ_RES_3; nb++)
            for(int nc=-FREQ_RES_3; nc<=FREQ_RES_3; nc++)
                if(fabs(a*na+b*nb+c*nc) < acc)
                {
                    ma=na; mb=nb; mc=nc;
                    return true;
                }
    return false;
}

std::string COrbitRuntimeTrajectory::performFrequencyAnalysis(size_t start, size_t count, CSpectrum sl[N_DIM], vectord *fa, double Lavg[N_DIM], double Lvar[N_DIM]) const
{
    findLines(start, count, sl, fa);
    std::string Description;
    double Labs[N_DIM];  // abs.value of ang.mom.components
    double Cavg[N_DIM];  // avg.value of C-th coordinate
    double Cmax[N_DIM];  // max.value of this coordinate
    for(unsigned int d=0; d<N_DIM; d++)
    {
        Lavg[d]=0;  // mean value of d-axis component of angular momentum
        Lvar[d]=0;  // variation in this component (r.m.s.)
        Labs[d]=0;  // mean absolute value of d-axis component of ang.mom.
        Cavg[d]=0;
        Cmax[d]=0;
    }
    for(size_t i=start; i<start+count; i++)
    {
        for(unsigned int d=0; d<N_DIM; d++)
        {
            double L = traj[i].Pos[(d+1)%N_DIM]*traj[i].Vel[(d+2)%N_DIM] - traj[i].Pos[(d+2)%N_DIM]*traj[i].Vel[(d+1)%N_DIM];
            Lavg[d] += L;
            Labs[d] += fabs(L);
            Lvar[d] += pow_2(L);
            Cavg[d] += traj[i].Pos[d];
            Cmax[d] = std::max<double>(Cmax[d], traj[i].Pos[d]);
        }
    }
    for(unsigned int d=0; d<N_DIM; d++)
    {
        Lavg[d]/=count;
        Labs[d]/=count;
        Cavg[d]/=count;
        Lvar[d]=sqrt(std::max<double>(0, Lvar[d]/count-pow_2(Lavg[d])));
    }

    // determine orbit class
    bool tubez=(fabs(Lavg[2])>Labs[2]*OC_TUBE_AVGL);  // has definite sense of rotation about z axis (more robust indicator than 1:1 correspondence of leading freqencies)
    bool tubex=(N_dim==3) && (fabs(Lavg[0])>Labs[0]*OC_TUBE_AVGL);  // same for x axis
    bool tubey=(N_dim==3) && (fabs(Lavg[1])>Labs[1]*OC_TUBE_AVGL);  // same for y axis (although we don't expect to find y-axis tubes if axes are sorted properly)
    double lfX=sl[0][0];
    double lfY=sl[1][0];
    double acc = InitData.timeUnit/(InitData.timeStep*count) * FREQ_ACCURACY;  // accuracy of comparison
    if(acc<FREQ_DIFF_REG_MIN) acc=FREQ_DIFF_REG_MIN;  // don't be too demanding.. (actually the accuracy limit is determined by energy conservation error, for simplicity ignored here)
    if(N_dim==2)
    {
        // try to find nx:ny correspondence of leading frequencies
        int nx, ny;
        if(commensurable(lfX, lfY, acc, nx, ny))    // boxlet family or tube (in case lfX==lfY)
        {
            double lf = (lfX/nx + lfY/ny)/2;  // mean value of common denominator 
            // find additional (libration) frequency
            double lfa = 0;
            bool regular=true;
            int rank=1;
            int ma, mb;
            for(unsigned int c=0; c<N_dim; c++)
            {
                for(size_t i=1; regular&&(i<sl[c].size()); i++)
                {
                    if(commensurable(sl[c][i], lf, acc, ma, mb)) 
                    {
                        if(((ma==1) || (abs(ma-mb)==1)) && (mb>1))
                            rank=mb;
                    }
                    else
                    {
                        double lfa1=std::min<double>(sl[c][i].freq, fabs(sl[c][0].freq-sl[c][i].freq));
                        if(lfa)
                        {
                            if(!lineardepend(sl[c][i], lf, lfa, acc, ma, mb)) 
                            {
                                regular=false;
                            }
                            else 
                                if(lfa1<lfa) 
                                {
                                    lfa=lfa1;
                                }
                        }
                        else
                        {
                            lfa=lfa1;
                        }
                    }
                }
            }
            if(tubez)
                Description = "Tube";
            else if((nx==1)&&(ny==2))
                Description = "banana";
            else if((nx==2)&&(ny==3))
                Description = "fish";
            else if((nx==3)&&(ny==4))
                Description = "pretzel";
            else 
                Description = convertToString(nx)+":"+convertToString(ny)+" resonant boxlet";
            Description += " rank "+convertToString(rank);
            if(!regular)
                Description = "chaotic "+Description;
        }
        else 
        {
            // for regular box orbit the remaining lines must be linear-dependent on two leading frequencies
            bool regular=true;
            int ma, mb;
            for(int c=0; c<2; c++)
                for(size_t i=1; regular&&(i<sl[c].size()); i++)
                    if(!lineardepend(sl[c][i], lfX, lfY, acc, ma, mb)) regular=false;
            if(regular)
                Description = "box";
            else
                Description = "chaotic";
        }
    }
    else    // N_dim==3
    {
        // determine if there are resonances between each pair of coordinates
        bool res[3];
        int nres=0;          // number of resonances
        int na[3], nb[3];    // resonance quotients
        double lf[3]={0,0,0};// leading frequencies (1 for nres=3, 2 for nres=1, 3 for nres=0)
        int nlf=0;           // number of leading frequencies
        for(int d=0; d<3; d++)
        {
            res[d] = commensurable(sl[d][0], sl[(d+1)%3][0], acc, na[d], nb[d]);
            if(res[d]) nres++;
        }
        
        if(nres==2)  // check for consistency (there cannot be two resonances)
        {   // reconstruct missing resonance number dn
            int dn=!res[0] ? 0 : (!res[1] ? 1 : 2);
            na[dn]=nb[(dn+2)%3]*nb[(dn+1)%3];
            nb[dn]=na[(dn+1)%3]*na[(dn+2)%3];
            simplify(na[dn], nb[dn]);
            nres=3;
        }
        if(nres==0)
        {
            nlf=3;
            lf[0]=sl[0][0];
            lf[1]=sl[1][0];
            lf[2]=sl[2][0];
            int n1, n2, n3;
            if(commensurable3(lf[0],lf[1],lf[2],acc,n1,n2,n3))
            {
                Description = "("+convertToString(n1)+","+convertToString(n2)+","+convertToString(n3)+") thin orbit";
                nlf = 2;
            }
            else
            {
                // either box or pyramid (near black hole)
                if(fabs(Cavg[2])>fabs(Cmax[2]*OC_PYRAMID_AVGZ))
                    Description = "pyramid";
                else if(fabs(Cavg[0])>fabs(Cmax[0]*OC_PYRAMID_AVGZ))
                    Description = "x-pyramid";   // shouldn't exist normally
                else if(fabs(Cavg[1])>fabs(Cmax[1]*OC_PYRAMID_AVGZ))
                    Description = "y-pyramid";   // shouldn't exist normally
                else
                    Description = "box";
            }
        }
        if(nres==1)
        {   
            int de=res[0] ? 0 : (res[1] ? 1 : 2);
            nlf=2;
            lf[0] = (sl[de][0]/na[de] + sl[(de+1)%3][0]/nb[de])/2;
            lf[1] = sl[(de+2)%3][0];
            if(lf[1]==0)
                nlf=1;  // avoid degenerate cases
            switch(de)
            {
            case 0: Description = convertToString(na[0])+":"+convertToString(nb[0])+":*"; break;
            case 1: Description = "*:"+convertToString(na[1])+":"+convertToString(nb[1]); break;
            case 2: Description = convertToString(nb[2])+":*:"+convertToString(na[2]); break;
            }
            Description +=" resonance";
        }
        if(nres==3)
        {
            nlf=1;
            lf[0]=std::min<double>((sl[0][0].freq/na[0] + sl[1][0].freq/nb[0])/2, std::min<double>((sl[1][0].freq/na[1] + sl[2][0].freq/nb[1])/2, (sl[2][0].freq/na[2] + sl[0][0].freq/nb[2])/2));
            int nx=(int)floor(sl[0][0]/lf[0]+0.5); 
            int ny=(int)floor(sl[1][0]/lf[0]+0.5); 
            int nz=(int)floor(sl[2][0]/lf[0]+0.5); 
            Description = convertToString(nx)+":"+convertToString(ny)+":"+convertToString(nz)+" resonance";
        }
        if(tubex && !tubez && !tubey) 
        {
            Description = "LAT "+Description;
            // try to determine whether it's inner or outer LAT
            // ..remains to be implemented
        }
        if(tubez && !tubey && !tubex) 
        {
            if(fabs(Cavg[2])>fabs(Cmax[2]*OC_PYRAMID_AVGZ))
                Description = "SAU "+Description;  // saucer orbit (a variant of short-axis tube)
            else
                Description = "SAT "+Description;
        }
        if(tubey && !tubex && !tubez)
        {
            Description = "IAT "+Description;  // intermediate-axis tube, they simply don't exist.
        }
        // check that every other freq is commensurable with base
        int ma, mb, mc;
        bool regular=true;
        for(int c=0; c<3; c++)
            for(size_t i=1; regular&&(i<sl[c].size()); i++)
            {
                switch(nlf)
                {
                case 1:
                    if(!commensurable(sl[c][i], lf[0], acc, ma, mb)) 
                    {
                        nlf++;
                        lf[1] = std::min<double>(sl[c][i].freq, fabs(lf[0]-sl[c][i].freq));
                    }
                    break;
                case 2:
                    if(!lineardepend(sl[c][i], lf[0], lf[1], acc, ma, mb))
                    {
                        nlf++;
                        lf[2] = std::min<double>(sl[c][i].freq, std::min<double>(fabs(lf[0]-sl[c][i].freq), fabs(lf[1]-sl[c][i].freq)));
                    }
                    break;
                case 3:
                    if(!lineardepend3(sl[c][i], lf[0], lf[1], lf[2], acc, ma, mb, mc))
                        regular=false;
                }
            }
        if(!regular)
            Description = "chaotic "+Description;
    }
    return Description;
}

void COrbitRuntimeTrajectory::calcLFCCoverInterval(size_t start, size_t count, double leadfreq[]) const
{
    CSpectrum tempSpectrum[N_DIM];
    findLines(start, count, tempSpectrum, NULL);
    for(unsigned int c=0; c<N_dim; c++)
        leadfreq[c]=(tempSpectrum[c].size()>0)? tempSpectrum[c][0].freq : -1;
}

// calculate rate of diffusion of leading frequencies in cartesian coordinates (LFCCs), which are related to the degree of chaos in the orbit
// determines the LFCCs for first and second halves of orbit, computes difference, returns delta f/f averaged over all coordinates
double COrbitRuntimeTrajectory::calcLFCCdiffusion(size_t start, size_t count) const
{
    if(count<=0) count=traj.size()-start;  // whole orbit by default
    size_t count2=count/2;
    double lf1[N_DIM], lf2[N_DIM];  // leading frequencies of 1st and 2nd half for each coordinate
    calcLFCCoverInterval(start, count2, lf1);
    calcLFCCoverInterval(start+count2, count2, lf2);
    double FreqDiff = 0;
#if FREQ_DIFF_METHOD == 0
    for(int c=N_dim-1; c>=0 && FreqDiff==0; c++)
      if(lf1[c] + lf2[c] > 0 )
        FreqDiff = fabs(lf1[c] - lf2[c]) * 2 / (lf1[c] + lf2[c]);
#else
    int nf=0;
    for(unsigned int c=0; c<N_dim; c++)
        if(lf1[c] + lf2[c] > 0 )
#if FREQ_DIFF_METHOD == 2
          if(fabs(lf1[c] - lf2[c]) * 2 / (lf1[c] + lf2[c]) < 1)
#endif
        {
            FreqDiff += fabs(lf1[c] - lf2[c]) * 2 / (lf1[c] + lf2[c]);
            nf++;
        }
    if(nf>0) FreqDiff /= nf;
#endif
    /*
    double maxFreqDiff = pow_2(timeUnit/intTime);  // rough estimate of regular/chaotic divide
    if(maxFreqDiff<FREQ_DIFF_REG_MIN) maxFreqDiff=FREQ_DIFF_REG_MIN;

    if(FreqDiff<1*maxFreqDiff)
        DescriptionLFCCdiff = "regular";
    else 
    if(FreqDiff<2*maxFreqDiff)
        DescriptionLFCCdiff = "probably regular";
    else 
    if(FreqDiff<4*maxFreqDiff)
        DescriptionLFCCdiff = "probably chaotic";
    else 
        DescriptionLFCCdiff = "chaotic";*/
    return FreqDiff;
}

void COrbitRuntimeTrajectory::calcInertiaOverInterval(size_t start, size_t count, double inertia[]) const
{
    for(unsigned int c=0; c<N_DIM; c++)
        inertia[c] = 0;
    for(size_t i=0; i<count; i++)
    {
        for(unsigned int c=0; c<N_dim; c++)
            inertia[c] += pow_2(traj[i+start].Pos[c]);
    }
    for(unsigned int c=0; c<N_dim; c++)
        inertia[c] = sqrt(inertia[c]/count);
}

// calculates maximal Lyapunov characteristic exponent (if the deviation equation was integrated along with orbit)
double COrbitRuntimeTrajectory::calcLambda() const
{
    if(!InitData.calcLyapunov)
        return 0;
    double sumavg=0;     // to calculate average value of ln(|w|/t) over the linear-growth interval
    int cntavg=0;
    double sumlambda=0;  // to calculate average value of ln(|w|)/t over the exponential growth interval
    int cntlambda=0;
    int timeblock=static_cast<int>(InitData.timeUnit/InitData.timeStep*2);    // averaging interval (roughly orbital period * 2)
    std::vector<double> interval(timeblock);
    bool chaotic=false;
    for(int t=timeblock; t<static_cast<int>(traj.size()); t+=timeblock)
    {
        int cnt = std::min<int>(timeblock, static_cast<int>(traj.size())-t);
        if(!chaotic)
        {
            // compute median value over 'timeblock' points  to avoid outliers
            interval.resize(cnt);
            for(int i=0; i<cnt; i++)
                interval[i] = (lnw[i+t]-log((i+t)*InitData.timeStep));
            std::sort(interval.begin(), interval.end());
            double mid = interval[cnt/2];
            sumavg+=mid;
            cntavg++;
            if((cntavg>4) && (mid > sumavg/cntavg + 2))
                chaotic=true;     // from here on, assume that we got into t->infinity regime and start collecting data for estimating lambda
        }
        if(chaotic)
        {
            for(int i=0; i<cnt; i++)
            {
                sumlambda += lnw[i+t]/(i+t);
                cntlambda++;
            }
        }
    }
    if(chaotic)
        return (sumlambda / cntlambda * InitData.timeUnit / InitData.timeStep);
    else
        return 0;
}

// return maximal distance between adjacent points
double COrbitRuntimeTrajectory::getMaxDist(size_t start, size_t count) const
{
    if(start>=traj.size())  return -1; // invalid call
    if(count<=0 || count>traj.size()-start) count = traj.size() - start;   // passing zero for count means analysing entire orbit
    if(count<=1) return -1; // empty orbit
    double maxd=0;
    CPosPoint<double> prev=traj[start];
    for(size_t i=start+1; i<start+count; i++)
    {
        CPosPoint<double> cur=traj[i];
        maxd = std::max<double>(maxd, sqrt(pow_2(cur.Pos[0]-prev.Pos[0])+pow_2(cur.Pos[1]-prev.Pos[1])+pow_2(cur.Pos[2]-prev.Pos[2])));
        prev=cur;
    }
    return maxd;
}

//----------------------------------------------------------------------------//
// Timestep function for storing a number of sampling points from trajectory
COrbitRuntimeTrajSample::COrbitRuntimeTrajSample(const COrbit* _orbit, size_t _numSamplingPoints) :
  CBasicOrbitRuntimeFnc(_orbit), 
  numSamplingPoints(_numSamplingPoints), 
  timeStep(orbit->getInitData().timeStep),
  N_dim(orbit->getPotential()->N_dim)
{}; 

void COrbitRuntimeTrajSample::Timestep(const double told, const double t, const double /*y*/[])
{
    int i1=(int)(told/timeStep); 
    int i2=(int)(t/timeStep);

    // store trajectory at regular intervals of time    
    for (int i=i1; i<=i2; i++)
    {
        double tt=timeStep*i;
        if ((tt>=told) && (tt<t))
        {
            CPosVelPoint<float> pt;
            for(unsigned int c=0; c<N_dim; c++)
            {
                pt.Pos[c]=(float)orbit->getInterpolatedTrajectory(c, tt);
                pt.Vel[c]=(float)orbit->getInterpolatedTrajectory(c+N_dim, tt);
            }
            trajEntire.push_back(pt);
        }
    }
}

void COrbitRuntimeTrajSample::Finish()
{
    size_t numActual=trajEntire.size();
    if(numActual<=1) return;
    numSamplingPoints = std::min<unsigned int>(static_cast<unsigned int>(numSamplingPoints), static_cast<unsigned int>(numActual*0.99));
    std::set<size_t> samplingIndexes;
    for(size_t i=0; i<numSamplingPoints; i++)
    {
        size_t np;
        do{
            np=static_cast<size_t>(floor(rand()*1.0/RAND_MAX*(numActual-1)));
        }
        while(samplingIndexes.find(np)!=samplingIndexes.end());
        samplingIndexes.insert(np);   // try to avoid duplication
        trajSample.push_back(trajEntire[np]);
    }
}

CBasicOrbitRuntimeFnc* COrbitRuntimeTrajSampleCreator::createRuntimeFnc(const COrbit* orbit) const
{ 
    // find if this orbit is listed in the 'special' list
    if(!numSamplingPointsSpecial.empty())
    {
        CPosVelPoint<double> initCond(orbit->getInitData().initCond);
        for(size_t i=0; i<numSamplingPointsSpecial.size(); i++)
        {
            CPosVelPoint<double> specialCond(numSamplingPointsSpecial[i].first);
            const double eps=1e-6;  // accuracy of comparison (float rather than double)
            if( gsl_fcmp(initCond.Pos[0], specialCond.Pos[0], eps)==0 &&
                gsl_fcmp(initCond.Pos[1], specialCond.Pos[1], eps)==0 &&
                gsl_fcmp(initCond.Pos[2], specialCond.Pos[2], eps)==0 &&
                gsl_fcmp(initCond.Vel[0], specialCond.Vel[0], eps)==0 &&
                gsl_fcmp(initCond.Vel[1], specialCond.Vel[1], eps)==0 &&
                gsl_fcmp(initCond.Vel[2], specialCond.Vel[2], eps)==0)
                return new COrbitRuntimeTrajSample(orbit, numSamplingPointsSpecial[i].second);
        }
    }
    // otherwise return default 
    return new COrbitRuntimeTrajSample(orbit, numSamplingPointsDefault); 
}

//----------------------------------------------------------------------------//
// simple trajectory timestep function recording points in Poincare section
struct CPoincareParam {
    const COrbit* orbit;
    unsigned int psindex;
};

double findPS(double t, void* params)
{
    return static_cast<CPoincareParam*>(params)->orbit->getInterpolatedTrajectory(static_cast<CPoincareParam*>(params)->psindex,t);
}

void COrbitRuntimePoincare::Timestep(const double told, const double t, const double y[])
{
    // Poincare section
    unsigned int N_dim=orbit->getPotential()->N_dim;
    double yprev=(told>0)?orbit->getInterpolatedTrajectory(PS2, told):y[PS2];  // previous y-value
    if((ps.size()<MAX_POINCARE_POINTS) && (yprev<=0) && (y[PS2]>=0))  // put a point
    {
        if(y[PS2]>yprev)
        {
            // interpolate coordinates at Poincare section crossing
            gsl_function F;
            CPoincareParam params;
            F.function=&findPS;
            F.params=&params;
            params.orbit=orbit;
            params.psindex=PS2;
            gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
            gsl_root_fsolver_set (s, &F, told, t);
            int status=0, iter=0;
            double tcross, tlow, tupp;
            do{
                iter++;
                status = gsl_root_fsolver_iterate (s);
                tcross = gsl_root_fsolver_root (s);
                tlow = gsl_root_fsolver_x_lower (s);
                tupp = gsl_root_fsolver_x_upper (s);
                status = gsl_root_test_interval (tlow, tupp, 0, 1e-8*(t-told));
            }
            while (status == GSL_CONTINUE && iter < 100);
            gsl_root_fsolver_free (s);
            ps.push_back(pairf((float)orbit->getInterpolatedTrajectory(PS1,tcross), (float)orbit->getInterpolatedTrajectory(PS1+N_dim, tcross)));
        }
        else
        {
            ps.push_back(pairf((float)y[PS1], (float)y[PS1+N_dim]));
        }
    }
}

//----------------------------------------------------------------------------//
COrbitRuntimePericenter::COrbitRuntimePericenter(const COrbit* _orbit):
CBasicOrbitRuntimeFnc(_orbit)
{
    fitSlope=fitIntercept=fitScatter=fitSignificance=0;
    Lavg[0]=Lavg[1]=Lavg[2]=L2avg[0]=L2avg[1]=L2avg[2]=0;
    // init sign
    COrbitInitData<double> data = orbit->getInitData();
    drdtsign= data.initCond.Pos[0]*data.initCond.Vel[0]+data.initCond.Pos[1]*data.initCond.Vel[1]+data.initCond.Pos[2]*data.initCond.Vel[2] >0 ? 1 : -1;
    // find approximate L2 for a circular orbit
    Lcirc2=findLcirc2(orbit->getPotential(), orbit->getPotential()->totalEnergy(data.initCond));
}

template<typename NumT> COrbitRuntimePericenter::COrbitRuntimePericenter(const COrbit* _orbit, const std::vector< CPosVelPoint<NumT> > &traj):
CBasicOrbitRuntimeFnc(_orbit)
{
    fitSlope=fitIntercept=fitScatter=fitSignificance=0;
    Lavg[0]=Lavg[1]=Lavg[2]=L2avg[0]=L2avg[1]=L2avg[2]=0;
    drdtsign= traj[0].Pos[0]*traj[0].Vel[0]+traj[0].Pos[1]*traj[0].Vel[1]+traj[0].Pos[2]*traj[0].Vel[2] >0 ? 1 : -1;
    int N_dim=orbit->getPotential()->N_dim;  // shortcut
    double timeStep=orbit->getInitData().timeStep;
    for(size_t i=1; i<traj.size(); i++)
    {
        int prevsign=drdtsign;
        drdtsign= traj[i].Pos[0]*traj[i].Vel[0]+traj[i].Pos[1]*traj[i].Vel[1]+traj[i].Pos[2]*traj[i].Vel[2] >0 ? 1 : -1;
        if(drdtsign*prevsign<0 && drdtsign>0)
        {   // direction reversed - passed pericenter; do not try to refine exact time, just record current squared ang.mom.
            double Lx=(N_dim==3)? traj[i].Pos[1]*traj[i].Vel[2]-traj[i].Pos[2]*traj[i].Vel[1] : 0;
            double Ly=(N_dim==3)? traj[i].Pos[2]*traj[i].Vel[0]-traj[i].Pos[0]*traj[i].Vel[2] : 0;
            double Lz=traj[i].Pos[0]*traj[i].Vel[1]-traj[i].Pos[1]*traj[i].Vel[0];
            double L2=pow_2(Lx)+pow_2(Ly)+pow_2(Lz);
            pericentert.push_back(i*timeStep);
            pericenterr.push_back(sqrt(pow_2(traj[i].Pos[0])+pow_2(traj[i].Pos[1])+pow_2(traj[i].Pos[2])));
            pericenterl2.push_back(L2);
            Lavg[0]+=Lx; L2avg[0]+=pow_2(Lx);
            Lavg[1]+=Ly; L2avg[1]+=pow_2(Ly);
            Lavg[2]+=Lz; L2avg[2]+=pow_2(Lz);
        }
    }
    // find approximate L2 for a circular orbit
    Lcirc2=findLcirc2(orbit->getPotential(), orbit->getPotential()->totalEnergy(traj[0]));
}

struct CPericenterParam {
    const COrbit* orbit;
};

double finddrdt(double t, void* params)
{
    const COrbit* o=static_cast<CPericenterParam*>(params)->orbit;
    unsigned int N_dim=o->getPotential()->N_dim;
    return o->getInterpolatedTrajectory(0,t)*o->getInterpolatedTrajectory(N_dim,t) + o->getInterpolatedTrajectory(1,t)*o->getInterpolatedTrajectory(N_dim+1,t) + (N_dim==3?o->getInterpolatedTrajectory(2,t)*o->getInterpolatedTrajectory(5,t):0);
}

void COrbitRuntimePericenter::Timestep(const double told, const double t, const double y[])
{
    unsigned int N_dim=orbit->getPotential()->N_dim;  // shortcut
    int prevsign=drdtsign;
    drdtsign= y[0]*y[N_dim]+y[1]*y[N_dim+1]+(N_dim==3?y[2]*y[5]:0) >0 ? 1 : -1;
    if(t>told && (drdtsign*prevsign<0) && drdtsign>0)
    {
        // direction reversed -- passed pericenter; find exact passagee time
        gsl_function F;
        CPericenterParam params;
        params.orbit=orbit;
        F.function=&finddrdt;
        F.params=&params;
        gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        gsl_root_fsolver_set (s, &F, told, t);
        int status=0, iter=0;
        double tperi, tlow, tupp;
        do{
            iter++;
            status = gsl_root_fsolver_iterate (s);
            tperi = gsl_root_fsolver_root (s);
            tlow = gsl_root_fsolver_x_lower (s);
            tupp = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (tlow, tupp, 0, 1e-3*(t-told));
        }
        while (status == GSL_CONTINUE && iter < 100);
        gsl_root_fsolver_free (s);
        pericentert.push_back(tperi);
        double yperi[N_DIM*2];
        for(unsigned int i=0; i<N_dim*2; i++)
            yperi[i]=orbit->getInterpolatedTrajectory(i, tperi);
        pericenterr.push_back(sqrt(pow_2(yperi[0])+pow_2(yperi[1])+pow_2(yperi[2])));
        double Lx=(N_dim==3)? yperi[1]*yperi[5]-yperi[2]*yperi[4] : 0;
        double Ly=(N_dim==3)? yperi[2]*yperi[3]-yperi[0]*yperi[5] : 0;
        double Lz=(N_dim==3)? yperi[0]*yperi[4]-yperi[1]*yperi[3] : yperi[0]*yperi[3]-yperi[1]*yperi[2];
        double L2=pow_2(Lx)+pow_2(Ly)+pow_2(Lz);
        pericenterl2.push_back(L2);
        Lavg[0]+=Lx; L2avg[0]+=pow_2(Lx);
        Lavg[1]+=Ly; L2avg[1]+=pow_2(Ly);
        Lavg[2]+=Lz; L2avg[2]+=pow_2(Lz);
    }
}

void COrbitRuntimePericenter::Finish()
{
    vectord lp(pericenterl2);  // !! here one may substitute either pericenter radius or pericenter ang.mom.sq., but the latter is more general
    std::sort(lp.begin(), lp.end());
    int maxn=static_cast<int>(lp.size());
    int nump=static_cast<int>(maxn*PERICENTER_FIT_FRACTION);
    if(nump<PERICENTER_FIT_MIN_POINTS) nump=PERICENTER_FIT_MIN_POINTS;
    if(nump>PERICENTER_FIT_MAX_POINTS) nump=PERICENTER_FIT_MAX_POINTS;
    if(nump>=maxn) nump=maxn-1;
    if(nump<3) { fitScatter=fitSignificance=fitSlope=fitIntercept=0; return; }
    // find out whether there are any conserved components of ang.mom.
    double L2min=0;  // minimum possible L^2 due to conserved quantities
    for(unsigned int d=0; d<N_DIM; d++)
    {
        Lavg[d]/=lp.size();
        L2avg[d]/=lp.size();
        const double eps=1e-5;  // relative accuracy of comparison (std.dev. in L[d] < eps*|L[d]| )
        if(gsl_fcmp(pow_2(Lavg[d]), L2avg[d], eps*eps)==0)
            L2min+=L2avg[d];
    }
    if(L2min>0)  // subtract it from the values for fitting
    {
        if(L2min>lp[0]) L2min=lp[0];  // to be sure that first point is not below that...
        for(int i=0; i<nump; i++)
            lp[i]-=L2min;
    }
    vectord np(nump)/*, weight(nump)*/;
    for(int i=0; i<nump; i++)
    {
        np[i]=(i+1.0)/maxn;
        /*weight[i] = (i<PERICENTER_FIT_MIN_POINTS ? 1.0 : (PERICENTER_FIT_MIN_POINTS+0.0)/i);*/  // not used -- may be improved in the future to give more weight to innermost points
    }
    double cov00, cov01, cov11, sumsq;
    // fit linear regression with and without constant term
    gsl_fit_linear(&(np.front()), 1, &(lp.front()), 1, nump, &fitIntercept, &fitSlope, &cov00, &cov01, &cov11, &sumsq);
    // compute expected dispersion assuming that chi^2/(N-2)=1
    double sigma_expected = fitSlope/maxn * sqrt(nump*1.0);
    double sigma_actual = sqrt(sumsq/(nump-2));
#if 0
    double sigma2_used = pow_2(std::max(sigma_expected, sigma_actual));
#else
    double sigma2_used = pow_2(sigma_actual);
#endif
    double chi2 = sumsq/sigma2_used;    // equals to nump-2 or less

    double fitSlopeM, cov11M, sumsqM;
    gsl_fit_mul(&(np.front()), 1, &(lp.front()), 1, nump, &fitSlopeM, &cov11M, &sumsqM);
    double deltachi2 = //(sumsqM/sumsq - 1)*(nump-2);
        sumsqM/sigma2_used - chi2;
    if(deltachi2 < 11.8 || fitIntercept<0)  // <3sigma
    {
        fitSlope=fitSlopeM;
        fitIntercept=0;
        sumsq=sumsqM;
    }
    // compute ratio of actual over expected dispersion
    fitScatter = sigma_actual/sigma_expected;
    // compute how significant (how many sigmas) is intercept!=0 with respect to zero intercept
    double stddevcount=0;
    if(deltachi2<50)  // less than ~7sigma - directly compute equivalent number of 'gaussian' sigmas from chi^2 tail
    {
        stddevcount = gsl_cdf_gaussian_Qinv(gsl_cdf_chisq_Q(deltachi2, 2)/2, 1);
    }
    else
    {
        double intrcept, stddev;
        gsl_fit_linear_est(0, fitIntercept, fitSlope, cov00, cov01, cov11, &intrcept, &stddev);
        stddevcount = intrcept/stddev;  // distance in std.dev.
    }
    fitSignificance=stddevcount;
    fitIntercept += L2min;  // add back conserved value
}

/**** Utility function to find a 1:1 periodic orbit in a plane specified by coord1:coord2 ****/
bool findPeriodicOrbit(const CPotential* potential, double E, int coord1, int coord2, double& cross1)
{
    if(coord1>coord2) { int tmp=coord1; coord1=coord2; coord1=tmp; }
    double Xlo=0, Xhi=potential->longaxisradius(E);
    double tol=Xhi*1e-8;
    double Torb=potential->longaxisperiod(E, Xhi);
    COrbitInitData<double> IC(potential, Torb/10, Torb, CPosVelPoint<double>(), false);
    cross1=0;
    int niter=0;
    vectorRuntimeFncCreators crts;
    crts.push_back(new COrbitRuntimePoincareCreator(coord1, coord2));
    do{
        double v2;
        do{
            IC.initCond.Pos[coord1]=(Xlo+Xhi)/2; 
            v2=2*(E-potential->Phi(IC.initCond.Pos[0],IC.initCond.Pos[1],IC.initCond.Pos[2]));
            if(v2<0) Xhi=IC.initCond.Pos[coord1];
            else IC.initCond.Vel[coord2]=sqrt(v2);
        } 
        while(v2<0);
        COrbit* orb=new COrbit(IC, &crts);
        orb->integrateToTime(Torb*2);   // should perform more than one revolution
        orb->finishIntegration();
        const CPoincareInformation<float>* pi = static_cast<const CPoincareInformation<float>*>(orb->getNewInfo(0));
        if(pi->size()<2)
        {
            delete crts[0];
            return false;
        }
        if(pi->at(1).first > IC.initCond.Pos[coord1])
        {   Xlo=IC.initCond.Pos[coord1]; }
        else
        {   Xhi=IC.initCond.Pos[coord1]; }
        niter++;
        delete pi;
        delete orb;
    }
    while(niter<25 && Xhi-Xlo>tol && Xhi>tol);
    cross1=(Xhi>tol) ? IC.initCond.Pos[coord1] : 0;
    delete crts[0];
    return (Xhi>tol);
}

struct CRcircParam{
    const CPotential* P;
    double E;
    unsigned int direction;
};

double findRcirc(double r, void* params)
{
    double X=0, Y=0, Z=0;
    switch(static_cast<CRcircParam*>(params)->direction)
    {
    case 0: X=r; break;
    case 1: Y=r; break;
    case 2: Z=r; break;
    }
    const CPotential* P=static_cast<CRcircParam*>(params)->P;
    double V2 = 2*( static_cast<CRcircParam*>(params)->E - P->Phi(X, Y, Z) );
    double cv[N_DIM*2]={X,Y,Z,0,0,0};
    double force[N_DIM*2]={0,0,0,0,0,0};
    P->DiffEq(2*P->N_dim, 0, cv, force);
    return V2/r + force[static_cast<CRcircParam*>(params)->direction + P->N_dim];
}

/**** Utility function to find approximate L^2 of a circular orbit with given energy ****/
double findLcirc2(const CPotential* potential, double E)
{
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_function F;
    CRcircParam params;
    params.E=E;
    params.P=potential;
    F.function=&findRcirc;
    F.params=&params;
    // average over 2 or 3 coordinates
    double Lcirc2=0;
    unsigned int numAvg=0;
    for(unsigned int d=0; d<potential->N_dim; d++)
    {
        // find r such that v^2/r = F(r) and v^2/2+Phi(r)=E
        params.direction=d;
        double rupp=potential->findintersection( E, d==0?1.:0, d==1?1.:0, d==2?1.:0);
        double rlow=rupp/4;
        double rcirc=0;
        if(findRcirc(rlow, &params) * findRcirc(rupp, &params)>=0) continue;  // endpoints do not enclose root
        gsl_root_fsolver_set (s, &F, rlow, rupp);
        int status=0, iter=0;
        do{
            iter++;
            status = gsl_root_fsolver_iterate (s);
            rcirc= gsl_root_fsolver_root (s);
            rlow = gsl_root_fsolver_x_lower (s);
            rupp = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (rlow, rupp, 0, 1e-4);
        }
        while (status == GSL_CONTINUE && iter < 100);
        Lcirc2 += 2*(E - potential->Phi( d==0?rcirc:0, d==1?rcirc:0, d==2?rcirc:0))*pow_2(rcirc);
        numAvg++;
    }
    gsl_root_fsolver_free (s);
    if(numAvg) Lcirc2/=numAvg;
    return Lcirc2;
}

}