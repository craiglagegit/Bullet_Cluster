#include "orbitlib.h"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include "orbit.h"
#include "potential.h"
#include "massmodel.h"
#include "stringconv.h"

#ifdef DEBUGPRINT
#include <iostream>
#endif
namespace smile{

// helper function that creates an equivalent spherical mass model for given potential; returns NULL on fail
CMassModel* createMassModel(const CPotential* potential, int numNodes=50, double innerShellMass=0, double outerShellMass=0)
{
    double totalMass=potential->totalMass();
    if(innerShellMass==0) innerShellMass=totalMass*0.1/numNodes;
    if(outerShellMass==0) outerShellMass=totalMass*(1-0.1/numNodes);
    if(potential->Mbh>0) innerShellMass=std::min<double>(innerShellMass, potential->Mbh/10);
    double rout = totalMass<0 ? 10 : potential->getRadiusByMass(outerShellMass);   // totalMass<0 should not really happen! means that mass is infinite
    double rin = std::min<double>(0.1/numNodes, potential->getRadiusByMass(innerShellMass));  // somewhat arbitrary choice for min/max radii, but probably reasonable
    if(rin<0 && totalMass<0) rin=rout/numNodes;   // trying to create a grid for an infinite-mass model.. silly enough but we could try
    if(rout<0 || rin<0 || rout<rin) return NULL;  // some weird error
    if(potential->Mbh>0) rin=std::min<double>(rin, potential->Mbh/10);

    vectord radii(numNodes), inmass(numNodes);
    createNonuniformGrid(radii, numNodes-1, rin, rout, true); 
    radii.push_back(rout*10);
    for(int i=1; i<numNodes; i++)
        inmass[i]=potential->Mass(radii[i])+potential->Mbh;
    inmass[0]=potential->Mbh;
    // ensure that mass strictly increases with radius
    for(int i=1; i<numNodes; i++)
        if(inmass[i]<=inmass[i-1])
        {
            radii.erase(radii.begin()+i);
            inmass.erase(inmass.begin()+i);
            numNodes--;
            i--;
        }
    CMassModel* massModel=new CMassModel(radii, inmass);
    if(massModel->errorcode<0 && massModel->errorcode!=-5) { delete massModel; return NULL; }   // error code=-5 means that DF is negative somewhere, however that is non-critical for the present task
    return massModel;
}

//function that rounds double to NUM_DIGITS_ROUND significant digits. Used to initialize ICs
template<typename NumT> double xroundd(NumT x)
{
    if(x==0) return 0;
    return StringVariant(x, NUM_DIGITS_ROUND).toDouble();
}
template<typename NumT> float xroundf(NumT x)
{
    if(x==0) return 0;
    return StringVariant(x, NUM_DIGITS_ROUND).toFloat();
}

////--------------  COrbitDesc class - shortened orbit data -------------////
template<typename NumT> COrbitDesc::COrbitDesc(const COrbitInitData<NumT> &_InitData, const NumT _intTime, const NumT _maxTime, const vectorRuntimeFncCreators* _creators) :
InitData(completeInitData(COrbitInitData<float>(
    _InitData.potential, 
    xroundf(_InitData.timeStep), xroundf(_InitData.timeUnit),
    CPosVelPoint<float>(xroundf(_InitData.initCond.Pos[0]),xroundf(_InitData.initCond.Pos[1]),xroundf(_InitData.initCond.Pos[2]),
                        xroundf(_InitData.initCond.Vel[0]),xroundf(_InitData.initCond.Vel[1]),xroundf(_InitData.initCond.Vel[2])),
    _InitData.calcLyapunov))),
    creators(_creators)
{
    state=OS_INITIALIZED;
    orbit=NULL;
    intTime=xroundf(_intTime);
    maxTime=std::max<float>(intTime, xroundf(_maxTime));
    orbitWeight=0;
}

COrbitDesc::COrbitDesc(const COrbitInitData<float> &_InitData, const float _intTime, const float _maxTime, const float _weight, const vectorInformation& _info) :
    InitData(completeInitData(_InitData)), creators(NULL)
{
    state=OS_DONE;
    orbit=NULL;
    intTime=_intTime;
    maxTime=_maxTime;
    orbitWeight=_weight;
    info=_info;
}

COrbitDesc::~COrbitDesc()
{
    if(orbit!=NULL)
        delete orbit;
    for(vectorInformation::const_iterator iter=info.begin(); iter!=info.end(); iter++)
        delete (*iter);
}

void COrbitDesc::run()
{
    if(state!=OS_INITIALIZED && state!=OS_PREPARING) 
        return;  // shouldn't happen though...
    state=OS_RUNNING;
    for(vectorInformation::const_iterator iter=info.begin(); iter!=info.end(); iter++)  // delete all existing data
        delete (*iter);
    info.clear();
    try{
        orbit = new COrbit(InitData, creators);
        if(maxTime>intTime)
            orbit->integrateAdaptiveTime(intTime, maxTime);
        else
            orbit->integrateToTime(intTime);
        {
            orbit->finishIntegration();
            intTime=(float)orbit->getIntTime();  // may be longer than desired
            // get information from any additional timestep functions that were attached (including orbit analysis, trajectory sampler, etc)
            if(creators!=NULL) 
                for(size_t i=0; i<creators->size(); i++)
                    info.push_back(orbit->getNewInfo(i));
        }
        delete(orbit); 
        orbit=NULL;
    }
    catch(const std::bad_alloc&) {
        my_error("Not enough memory in orbit integration");
        info.push_back(new COrbitInformation<float>("Failed"));
        // do nothing..
    }
    state=OS_DONE;
}

void COrbitDesc::halt()
{
    if(state!=OS_PREPARING && state!=OS_RUNNING)
        return;
    if(orbit!=NULL)
    {
        state=OS_NEEDTOTERMINATE;
        orbit->halt();
    }
}

std::string COrbitDesc::toString() const 
{
    std::string strInfo;
    for(size_t i=0; i<getInfoNum(); i++)
        strInfo += getInfo(i)->toString();
    strInfo += "Int.time=" + convertToString(static_cast<int>(intTime/InitData.timeUnit))+" Torb" +
        "\nOrbit weight=" + convertToString(orbitWeight, 3) + "\n";
    return strInfo;
}

////--------------  COrbitLibrary class - array of COrbitDesc -------------////

COrbitLibrary::~COrbitLibrary() 
{
    for(std::vector<COrbitDesc*>::iterator iter=OrbitList.begin(); iter!=OrbitList.end(); ++iter)
        delete (*iter);
}

unsigned int COrbitLibrary::numComplete() const
{
    unsigned int nComplete=0;
    for(std::vector<COrbitDesc*>::const_iterator iter=OrbitList.begin(); iter!=OrbitList.end(); ++iter)
        if( (*iter)->getState()==COrbitDesc::OS_DONE) nComplete++;
    return nComplete;
}

int COrbitLibrary::findUnstartedOrbit() const
{
    for(size_t ind=0; ind<OrbitList.size(); ind++)
        if( OrbitList[ind]->getState()==COrbitDesc::OS_INITIALIZED)
        {
            OrbitList[ind]->setState(COrbitDesc::OS_PREPARING);   ///!!!??? how comes???const
            return static_cast<int>(ind);
        }
    return -1;
}

void COrbitLibrary::halt()
{
    for(std::vector<COrbitDesc*>::iterator iter=OrbitList.begin(); iter!=OrbitList.end(); ++iter)
        (*iter)->halt();
}

void COrbitLibrary::prepareInitialConditionsE(unsigned int NStationary, unsigned int NPrincipalPlane, unsigned int NYalpha, unsigned int NRandom, double E, double intTimeStepsPerPeriod, double intTimeInPeriods, double intTimeMaxAdaptive, bool calcLyapunov)
{
    double Xmax = potential->longaxisradius(E);
    if(Xmax<0) return;  // error...
    double Torbx= xroundd(potential->longaxisperiod(E, Xmax));
    double intTime = xroundd(intTimeInPeriods * Torbx);
    double maxTime = xroundd(intTimeMaxAdaptive * Torbx);
    double intTimeStep = xroundd(Torbx / intTimeStepsPerPeriod);
    //int rotation = (Omega>0)?2:1;    // !!! for rotating models, need to duplicate number of orbits since they are no more symmetric
    // Prepare initial conditions
    int nfms=NStationary;            // number of stationary start-space elements
    int nfmp=NPrincipalPlane/potential->N_dim;  // number of principal-plane start-space elements
    COrbitInitData<double> IC(potential, intTimeStep, Torbx, CPosVelPoint<double>(), calcLyapunov);  // used to initialize initial conditions
    if(potential->N_dim==3)
    {
        // stationary
#ifdef SCHWARZSCHILD_START_SPACE
        int ngrd=(int)sqrt(nfms/3.0);
#else
        int ngrd=sqrt(1.0*nfms*M_PI/2);
#endif
        for(int ith=0; ith<ngrd; ith++)
        {
#ifdef SCHWARZSCHILD_START_SPACE
            double theta = tan((ith+0.5)/ngrd * M_PI/4);
            for(int iph=0; iph<ngrd*3; iph++)
            {
                double phi = tan((iph/3+0.5)/ngrd * M_PI/4);
                double X=0, Y=0, Z=0;
                switch(iph%3)
                {
                case 0:
                    X=Xmax;
                    Y=Xmax*theta;
                    Z=Xmax*phi;
                    break;
                case 1:
                    Y=Xmax;
                    Z=Xmax*theta;
                    X=Xmax*phi;
                    break;
                case 2:
                    Z=Xmax;
                    X=Xmax*theta;
                    Y=Xmax*phi;
                    break;
                }
#else
            double theta = (1 - (ith+0.5)/ngrd) * M_PI/2;
            int nph = ngrd * sin(theta) + 1 ;
            for(int iph=0; iph<nph; iph++)
            {
                double phi = (iph+0.5)/nph * M_PI/2;
                double X=0, Y=0, Z=0;

                X=Xmax*sin(theta)*cos(phi);
                Y=Xmax*sin(theta)*sin(phi);
                Z=Xmax*cos(theta);
#endif
                double r=potential->findintersection(E, X, Y, Z);
                IC.initCond.Pos[0]=X*r; IC.initCond.Pos[1]=Y*r; IC.initCond.Pos[2]=Z*r;
                IC.initCond.Vel[0]=IC.initCond.Vel[1]=IC.initCond.Vel[2]=0;
                OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
            }
        }
        // principal-plane
        if(NPrincipalPlane>0)
        {
#ifdef ANNULUS_PP_START_SPACE
            double rmin=Xmax, rminc;
            findPeriodicOrbit(potential, E, 0, 1, rminc);
            if(rminc<rmin) rmin=rminc;
            if(rmin>0) findPeriodicOrbit(potential, E, 1, 2, rminc);
            if(rminc<rmin) rmin=rminc;
            if(rmin>0) findPeriodicOrbit(potential, E, 0, 2, rminc);
            if(rminc<rmin) rmin=rminc;
            rminc=rmin/Xmax;
#else
            double rminc=0;
#endif
            // use elliptic annulus with inner radius rmin and outer xmax
            ngrd = (int)(sqrt(4*nfmp/M_PI/(1+rminc))+0.5);
            for(int ir=0; ir<ngrd; ir++)
            {
                double r = (ir+0.5)/ngrd * (1-rminc) + rminc;
                int nph = (int)(r * ngrd * M_PI/2 + 0.22);
                for(int iph=0; iph<nph; iph++)
                {
                    double phi = (iph+0.5)/nph * M_PI/2;
                    double rmax=potential->findintersection(E, Xmax*cos(phi), Xmax*sin(phi), 0);
                    IC.initCond.Pos[0]=r*rmax*Xmax*cos(phi); IC.initCond.Pos[1]=r*rmax*Xmax*sin(phi); IC.initCond.Pos[2]=0;
                    IC.initCond.Vel[0]=IC.initCond.Vel[1]=0; IC.initCond.Vel[2]=sqrt(2*(E-potential->Phi(IC.initCond.Pos[0], IC.initCond.Pos[1], IC.initCond.Pos[2])));
                    if(IC.initCond.Vel[2]>0)
                        OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
                    rmax=potential->findintersection(E, 0, Xmax*cos(phi), Xmax*sin(phi));
                    IC.initCond.Pos[0]=0; IC.initCond.Pos[1]=r*rmax*Xmax*cos(phi); IC.initCond.Pos[2]=r*rmax*Xmax*sin(phi);
                    IC.initCond.Vel[1]=IC.initCond.Vel[2]=0; IC.initCond.Vel[0]=sqrt(2*(E-potential->Phi(IC.initCond.Pos[0], IC.initCond.Pos[1], IC.initCond.Pos[2])));
                    if(IC.initCond.Vel[0]>0)
                        OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
                    rmax=potential->findintersection(E, Xmax*sin(phi), 0, Xmax*cos(phi));
                    IC.initCond.Pos[0]=r*rmax*Xmax*sin(phi); IC.initCond.Pos[1]=0; IC.initCond.Pos[2]=r*rmax*Xmax*cos(phi);
                    IC.initCond.Vel[0]=IC.initCond.Vel[2]=0; IC.initCond.Vel[1]=sqrt(2*(E-potential->Phi(IC.initCond.Pos[0], IC.initCond.Pos[1], IC.initCond.Pos[2])));
                    if(IC.initCond.Vel[1]>0)
                        OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
                }
            }
        }
    }
    else   // N_dim==2
    {
        for(int iph=0; iph<nfms; iph++)
        {
            double phi = (iph+0.5)/nfms * M_PI/2;
            double X=Xmax*cos(phi);
            double Y=Xmax*sin(phi);
            double r=potential->findintersection(E, X, Y, 0);
            IC.initCond.Pos[0]=X*r; IC.initCond.Pos[1]=Y*r; IC.initCond.Pos[2]=0;
            IC.initCond.Vel[0]=IC.initCond.Vel[1]=IC.initCond.Vel[2]=0;
            OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
        }

        for(int ir=0; ir<nfmp; ir++)
        {
            double r = (ir+0.5)/nfmp;
            IC.initCond.Pos[0]=r*Xmax; IC.initCond.Pos[1]=0; IC.initCond.Pos[2]=0;
            IC.initCond.Vel[0]=IC.initCond.Vel[2]=0; IC.initCond.Vel[1]=sqrt(2*(E-potential->Phi(IC.initCond.Pos[0], IC.initCond.Pos[1], IC.initCond.Pos[2])));
            if(IC.initCond.Vel[1]>0)
                OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
            IC.initCond.Pos[0]=0; IC.initCond.Pos[1]=r*Xmax; IC.initCond.Pos[2]=0;
            IC.initCond.Vel[1]=IC.initCond.Vel[2]=0; IC.initCond.Vel[0]=sqrt(2*(E-potential->Phi(IC.initCond.Pos[0], IC.initCond.Pos[1], IC.initCond.Pos[2])));
            if(IC.initCond.Vel[0]>0)
                OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
        }

    }

    if(NYalpha>0)
    {
        // Y-alpha start space 
        int ny=(int)sqrt(NYalpha*1.0), na=NYalpha/ny;
        double maxx=potential->longaxisradius(E);
        double maxy=0;
        for(int niter=0; niter<50 && maxx-maxy>maxx*1e-8; niter++)
        {
            double Ey=potential->Phi(0,(maxy+maxx)/2,0);
            if(Ey>E) maxx=(maxx+maxy)/2; else maxy=(maxx+maxy)/2;
        }   // thus we have found approximately the intermediate axis radius at given energy;
        for(int iy=0; iy<ny; iy++)
        {
            double Y = (iy+0.5)/ny * maxy;
            double V=E-potential->Phi(0,Y,0);
            if(V<0) continue;   // should not happen though..
            V=sqrt(2*V);
            IC.initCond.Pos[0]=0; IC.initCond.Pos[1]=Y; IC.initCond.Pos[2]=0; IC.initCond.Vel[1]=0;
            for(int ia=0; ia<na; ia++)
            {
                double alpha = (ia+0.5)/na * M_PI;
                IC.initCond.Vel[0]=V*cos(alpha); IC.initCond.Vel[2]=V*sin(alpha);
                OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
            }
        }
    }

    // random initial conditions (if desired)
    if(NRandom<=0) return;

    const int npsampling = std::max<int>(NRandom*10, 1000);
    std::vector<double> Xs(npsampling), Ys(npsampling), Zs(npsampling), Ws(npsampling);
    double X, Y, Z, totalWeight=0;
    for(int s=0; s<npsampling; s++)
    {   // sample accessible configuration space with appropriate weight
        double v2;
        do{
            X=rand()*Xmax/RAND_MAX;
            Y=rand()*Xmax/RAND_MAX;
            Z=rand()*Xmax/RAND_MAX;
            v2=2*(E-potential->Phi(X,Y,Z));
        } while(v2<0);
        double W=sqrt(v2);
        Xs[s]=X; Ys[s]=Y; Zs[s]=Z; Ws[s]=W; totalWeight+=W;
      
    }
    // now select points from sampled
    for(unsigned int np=0; np<NRandom; np++)
    {
        double indw=rand() * totalWeight / RAND_MAX;
        int i=0;
        while(i<npsampling-1 && indw>=Ws[i]) 
        {
            indw-=Ws[i];
            i++;
        }
        IC.initCond.Pos[0]=Xs[i]; IC.initCond.Pos[1]=Ys[i]; IC.initCond.Pos[2]=Zs[i]; 
        double V=Ws[i];
        double costh=(potential->N_dim==3)? (rand()*2.0/RAND_MAX-1) : 0;
        double sinth=sqrt(1-pow_2(costh));
        double phi=rand()*2*M_PI/RAND_MAX;
        IC.initCond.Vel[0]=V*sinth*cos(phi); IC.initCond.Vel[1]=V*sinth*sin(phi); IC.initCond.Vel[2]=V*costh;
        OrbitList.push_back(new COrbitDesc(IC, intTime, maxTime, creators));
        // remove this point from sample list, to avoid duplication
        totalWeight-=Ws[i];
        Ws[i]=0;  // it won't occur again
    }
}

// sorting operator
bool compareOrbitWeight(COrbitDesc* elem1, COrbitDesc* elem2)
{
    return elem1->getWeight()<elem2->getWeight();
}

struct CDensityFinderParam {
    const CDensity* P;
    double dens;
    double X, Y, Z; 
};
double densityFinder(double r, void* params)
{
    return ((CDensityFinderParam*)params)->P->Rho(
        r*((CDensityFinderParam*)params)->X, r*((CDensityFinderParam*)params)->Y, r*((CDensityFinderParam*)params)->Z) 
        - ((CDensityFinderParam*)params)->dens;
}

// prepare IC for the whole model (using random initial conditions spanning whole range of energies)
void COrbitLibrary::prepareInitialConditions(unsigned int NRandom, double maxRadius, double intTimeStepsPerPeriod, double intTimeInPeriods, double intTimeMaxAdaptive, bool calcLyapunov)
{
    for(std::vector<COrbitDesc*>::iterator iter=OrbitList.begin(); iter!=OrbitList.end(); ++iter)
        delete (*iter);
    OrbitList.clear();
    if(maxRadius<=0 && !potential->finite()){
        my_error("Error, cannot create initial conditions for an infinite-mass model without cutoff in radius");
        return;   // cannot sample entire space for non-finite-mass models
    }
    try{
        const size_t NUM_RADIAL_POINTS_SPHERICAL_MODEL=50;  // parameter to pass to createMassModel controlling number of spline control points
        const size_t NUM_RADIAL_POINTS_AXIS_RATIO=10;       // evaluate shape at this many points in radius
        // create spherical mass model with its distribution function, which will be used to draw initial conditions
        CMassModel* massModel = createMassModel(potential, NUM_RADIAL_POINTS_SPHERICAL_MODEL, 0, maxRadius>0 ? potential->Mass(maxRadius):0);
        if(massModel==NULL) {
            my_error("Unknown error in creating equivalent spherical mass model");
            return;
        }
        double maxEnergy = maxRadius>0 ? potential->Phi(0,0,maxRadius) : 0;  // if maxRadius is given, use it; otherwise use r=infinity provided that potential is zero at infinity
        double totalMassWithoutBH = maxRadius>0 ? potential->Mass(maxRadius) : potential->totalMass();
        if(totalMassWithoutBH<0) {
            my_error("Error, total mass is not finite");
            return;
        }
        totalMassWithoutBH = std::min<double>(totalMassWithoutBH, massModel->getTotalMass()-potential->Mbh);   // to ensure correct initialization of outermost initial conditions
        // evaluate shape as function of radius (defined as axis ratios of equidensity ellipsoids)
        // NOTE: this works well only if density is decreasing with radius. A fail-safe approach is to assign axes ratio of unity if anything goes wrong.
        vectord radii_shape(NUM_RADIAL_POINTS_AXIS_RATIO), axes_a(NUM_RADIAL_POINTS_AXIS_RATIO), axes_b(NUM_RADIAL_POINTS_AXIS_RATIO), axes_c(NUM_RADIAL_POINTS_AXIS_RATIO);
        gsl_function F;
        CDensityFinderParam PP;
        PP.P=potential;
        F.function=&densityFinder;
        F.params=&PP;
        gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
#ifdef DEBUGPRINT
        std::cerr << "radius\taxis_x\taxis_y\taxis_z\n";
#endif
        for(size_t i=0; i<NUM_RADIAL_POINTS_AXIS_RATIO; i++)
        {
            radii_shape[i]=massModel->rad((i+0.5)/NUM_RADIAL_POINTS_AXIS_RATIO * totalMassWithoutBH + potential->Mbh);
            PP.dens=massModel->dens(radii_shape[i]);
            if(radii_shape[i]<0 || !gsl_finite(radii_shape[i]) || PP.dens<0 || !gsl_finite(PP.dens)) 
            {
                axes_a[i]=1.0; axes_b[i]=1.0; axes_c[i]=1.0;
                continue;  // weird things do happen, just skip it.
            }
            // find the radii along x,y,z axes where density in real model is equal to that of spherical model
            double Rlow, Rupp, Rroot;
            int status=0, iter=0;
            const double ARlow=0.1, ARupp=10;    // we won't be able to deal correctly with highly flattened models
            // X axis
            PP.X=radii_shape[i]; PP.Y=0; PP.Z=0;
            Rlow=ARlow; Rupp=ARupp;
            if(F.function(Rlow, F.params) * F.function(Rupp, F.params) < 0) {  // otherwise endpoints do not enclose root, assume axis ratio to be unity
                gsl_root_fsolver_set (s, &F, Rlow, Rupp);
                status=0; iter=0;
                do{
                    iter++;
                    status = gsl_root_fsolver_iterate (s);
                    Rroot= gsl_root_fsolver_root (s);
                    Rlow = gsl_root_fsolver_x_lower (s);
                    Rupp = gsl_root_fsolver_x_upper (s);
                    status = gsl_root_test_interval (Rlow, Rupp, 0, 1e-3);
                }
                while (status == GSL_CONTINUE && iter < 100);
            } else Rroot=1.0;
            if(i>0 && Rroot*radii_shape[i] < axes_a[i-1]*radii_shape[i-1])
                Rroot = axes_a[i-1]*radii_shape[i-1]/radii_shape[i] * 1.01;  // ensure monotonicity
            axes_a[i]=Rroot;
            // Y axis
            PP.X=0; PP.Y=radii_shape[i]; PP.Z=0;
            Rlow=ARlow; Rupp=ARupp;
            if(F.function(Rlow, F.params) * F.function(Rupp, F.params) < 0) {
                gsl_root_fsolver_set (s, &F, Rlow, Rupp);
                status=0; iter=0;
                do{
                    iter++;
                    status = gsl_root_fsolver_iterate (s);
                    Rroot= gsl_root_fsolver_root (s);
                    Rlow = gsl_root_fsolver_x_lower (s);
                    Rupp = gsl_root_fsolver_x_upper (s);
                    status = gsl_root_test_interval (Rlow, Rupp, 0, 1e-3);
                }
                while (status == GSL_CONTINUE && iter < 100);
            } else Rroot=1.0;
            if(i>0 && Rroot*radii_shape[i] < axes_b[i-1]*radii_shape[i-1])
                Rroot = axes_b[i-1]*radii_shape[i-1]/radii_shape[i] * 1.01;  // ensure monotonicity
            axes_b[i]=Rroot;
            // Z axis
            PP.X=0; PP.Y=0; PP.Z=radii_shape[i];
            Rlow=ARlow; Rupp=ARupp;
            if(F.function(Rlow, F.params) * F.function(Rupp, F.params) < 0) {
                gsl_root_fsolver_set (s, &F, Rlow, Rupp);
                status=0; iter=0;
                do{
                    iter++;
                    status = gsl_root_fsolver_iterate (s);
                    Rroot= gsl_root_fsolver_root (s);
                    Rlow = gsl_root_fsolver_x_lower (s);
                    Rupp = gsl_root_fsolver_x_upper (s);
                    status = gsl_root_test_interval (Rlow, Rupp, 0, 1e-3);
                }
                while (status == GSL_CONTINUE && iter < 100);
            } else Rroot=1.0;
            if(i>0 && Rroot*radii_shape[i] < axes_c[i-1]*radii_shape[i-1])
                Rroot = axes_c[i-1]*radii_shape[i-1]/radii_shape[i] * 1.01;  // ensure monotonicity
            axes_c[i]=Rroot;
#ifdef DEBUGPRINT
            std::cerr << radii_shape[i] << "\t" << axes_a[i] << "\t" << axes_b[i] << "\t" << axes_c[i] << "\n";
#endif
        }
        gsl_root_fsolver_free (s);
        // now use the mass model for sampling initial conditions
        int numOrbits=NRandom;
        gsl_rng* randn=gsl_rng_alloc(gsl_rng_default);
        CPosVelPoint<double> IC;
        for(int o=0; o<numOrbits; o++)
        {   // sample accessible configuration space with appropriate weight
            double radius=0, veloc, encmass, E, poten, costh, sinth, phi;
            encmass = (o+0.5)/numOrbits * totalMassWithoutBH + potential->Mbh;
            radius=massModel->rad(encmass);
            if(radius<0 || !gsl_finite(radius)) 
                continue;  // weird things do happen, just skip it.
            // find axis ratios at this radius
            size_t indbin=0;
            while(indbin<NUM_RADIAL_POINTS_AXIS_RATIO && radius > radii_shape[indbin]) indbin++;
            double axis_a=1, axis_b=1, axis_c=1;
            if(indbin==0) {
                axis_a = axes_a[0];
                axis_b = axes_b[0];
                axis_c = axes_c[0];
            } else if(indbin>=NUM_RADIAL_POINTS_AXIS_RATIO) {
                axis_a = axes_a.back();
                axis_b = axes_b.back();
                axis_c = axes_c.back();
            } else {  // linear interpolation
                axis_a = (radii_shape[indbin]*axes_a[indbin-1]-radii_shape[indbin-1]*axes_a[indbin] + radius*(axes_a[indbin]-axes_a[indbin-1]))/(radii_shape[indbin]-radii_shape[indbin-1]);
                axis_b = (radii_shape[indbin]*axes_b[indbin-1]-radii_shape[indbin-1]*axes_b[indbin] + radius*(axes_b[indbin]-axes_b[indbin-1]))/(radii_shape[indbin]-radii_shape[indbin-1]);
                axis_c = (radii_shape[indbin]*axes_c[indbin-1]-radii_shape[indbin-1]*axes_c[indbin] + radius*(axes_c[indbin]-axes_c[indbin-1]))/(radii_shape[indbin]-radii_shape[indbin-1]);
            }
            int ntrial=0;
            do{
                costh=gsl_rng_uniform(randn)*2.0-1;
                sinth=sqrt(1-costh*costh);
                phi=gsl_rng_uniform(randn)*2*M_PI;
                IC.Pos[0]=axis_a * radius*sinth*cos(phi); 
                IC.Pos[1]=axis_b * radius*sinth*sin(phi); 
                IC.Pos[2]=axis_c * radius*costh;
                poten=potential->Phi(IC.Pos[0], IC.Pos[1], IC.Pos[2]);
                ntrial++;
            } while(poten>=maxEnergy && ntrial<1e3);
            // correction factor for escape velocity
            double velcorr = sqrt(poten / massModel->pot(radius));
            if(gsl_isnan(velcorr)) velcorr=1.0;
            ntrial=0;
            do{
                veloc=massModel->samplev(radius) * velcorr;
                ntrial++;
                E=poten+veloc*veloc/2;
            } while(E >= maxEnergy && ntrial<1e3);  // to avoid infinite loop in the case of weird potential 
            costh=gsl_rng_uniform(randn)*2.0-1;
            sinth=sqrt(1-costh*costh);
            phi=gsl_rng_uniform(randn)*2*M_PI;
            IC.Vel[0]=veloc*sinth*cos(phi); 
            IC.Vel[1]=veloc*sinth*sin(phi); 
            IC.Vel[2]=veloc*costh;

            // faster and more approximate evaluation of Torb for PT_NB, although it underestimates actual period for a smoothed Nbody potential. But anyway it's not a good idea to run S.M. on frozen Nbody potentials..
            double Torbx= ( potential->PotentialType()==CDensity::PT_NB ? massModel->trad(E)*2 : potential->longaxisperiod(E));
            if(Torbx<0)
                continue;  // shouldn't happen but if it does then skip this orbit
            double intTime = (intTimeInPeriods * Torbx);
            double maxTime = (intTimeMaxAdaptive * Torbx);
            double intTimeStep = (Torbx / intTimeStepsPerPeriod);
            OrbitList.push_back(new COrbitDesc(COrbitInitData<double>(potential, intTimeStep, Torbx, IC, calcLyapunov), intTime, maxTime, creators));
            OrbitList.back()->setWeight((float)E);   // temporarily store total energy in the "weight" field, for sorting
        }
        gsl_rng_free(randn); 
        delete massModel;
        // sort orbits by energy
        std::sort(OrbitList.begin(), OrbitList.end(), compareOrbitWeight);
        for(size_t o=0; o<OrbitList.size(); o++) 
            OrbitList[o]->setWeight(0);   // clear this field
    }
    catch(...) { my_error("Unknown exception in generating initial conditions"); }
}

void COrbitLibrary::modifyInitialConditions(const CPotential* _potential, const vectorRuntimeFncCreators* _creators, double intTimeStepsPerPeriod, double intTimeInPeriods, double intTimeMaxAdaptive, bool calcLyapunov, bool modifyAll)
{
    if(!_potential) return;
    if(OrbitList.size()==0) return;
    potential=_potential;
    creators=_creators;
    COrbitDesc* orbnew;
    for(std::vector<COrbitDesc*>::iterator iter=OrbitList.begin(); iter!=OrbitList.end(); ++iter)
        if(modifyAll || (*iter)->getState()==COrbitDesc::OS_INITIALIZED) 
    {
        COrbitInitData<float> InitData((*iter)->getInitData());
        InitData.potential=potential;
        InitData.calcLyapunov=calcLyapunov;
        float intTime=(*iter)->getIntTime();
        float maxTime=(*iter)->getMaxTime();
        if(intTimeInPeriods>0 && modifyAll)
        {
            InitData.timeStep=(float)(InitData.timeUnit/intTimeStepsPerPeriod);
            orbnew=new COrbitDesc(InitData, (float)(intTimeInPeriods*InitData.timeUnit), (float)(intTimeMaxAdaptive*InitData.timeUnit), creators);
        }
        else  // modify only creators and potential
        {
            orbnew=new COrbitDesc(InitData, intTime, maxTime, creators);
            orbnew->setWeight( (*iter)->getWeight());
        }
        delete(*iter);
        *iter=orbnew;
    }
}

void COrbitLibrary::removeUnused()
{
    if(OrbitList.size()==0) return;
    for(std::vector<COrbitDesc*>::iterator iter=OrbitList.end()-1; iter!=OrbitList.begin(); )
        if((*iter)->getState()!=COrbitDesc::OS_DONE)
        {
            delete (*iter);
            OrbitList.erase(iter--);
        }
        else
            --iter;
}

std::string COrbitLibrary::getOrbitPopulation(const CBasicOrbitFilteringFnc* filter, const CBasicOrbitFilteringFnc* chaos, bool useWeight) const
{
    if(OrbitList.size()==0) return "Error: no orbits";
    std::map<std::string, double> families, familiesr, familiesc;
    double totalweight=0, weightreg=0, weightch=0;
    for(size_t o=0; o<OrbitList.size(); o++)
    {
        double w = useWeight ? OrbitList[o]->getWeight() : 1.0+rand()*0.00001/RAND_MAX;   // add random to avoid duplicates
        const COrbitInformation<float>* infoOrbit=static_cast<const COrbitInformation<float>*>(OrbitList[o]->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
        double acceptFactor=filter!=NULL ? (filter->eval(OrbitList[o])>=0.5 ? 1 : 0) : 1;
        double chaosFactor =chaos!=NULL ? chaos->eval(OrbitList[o]) : 1;
        if(chaosFactor<0) chaosFactor=0;
        if(chaosFactor>1) chaosFactor=1;
        w*=acceptFactor;
        if(w==0 || infoOrbit==NULL) continue;
        std::string family=infoOrbit->getDescription();
        if(family.substr(0, 8)=="chaotic ") family=family.erase(0, 8);
        if(family.find("SAT")!=std::string::npos && family.find("resonance")!=std::string::npos) family="SAT";  else // merge all subgroups like 'SAT 6:6:7' into one
        if(family.find("LAT")!=std::string::npos && family.find("resonance")!=std::string::npos) family="LAT";  
        families[family] +=w;
        familiesc[family]+=w*chaosFactor;
        familiesr[family]+=w*(1-chaosFactor);
        totalweight+=w;
        weightreg  +=w*(1-chaosFactor);
        weightch   +=w*chaosFactor;
    }
    std::vector<double> weights;
    for(std::map<std::string, double>::const_iterator iterm=families.begin(); iterm!=families.end(); iterm++) weights.push_back(iterm->second);
    std::sort(weights.begin(), weights.end(), std::greater<double>());
    int nf=0; 
    if(totalweight==0) return "No weights assigned to orbits";
    std::string info = StringVariant(weightch/totalweight*100,3).toString()+"% chaotic\n";
    std::vector<double>::const_iterator iterw=weights.begin();
    while(nf<8 && iterw!=weights.end())
    {
        // find key in map
        std::string f;
        for(std::map<std::string, double>::const_iterator iterm=families.begin(); f.empty() && iterm!=families.end(); iterm++) 
            if(iterm->second == *iterw) f=iterm->first;
        info +=  f + ": " + StringVariant((familiesr[f]+familiesc[f])/totalweight*100,3).toString() + "% (" + StringVariant(familiesc[f]/totalweight*100,3).toString() + "%)\n";
        nf++;
        ++iterw;
    }
    return info;
}

std::string COrbitLibrary::exportNbody(unsigned int numPoints, const CBasicOrbitFilteringFnc* massRefineFactor, CPointMassSetFloat* &result, CPointCountSetFloat* &numSamplingPointsMin) const
{
    unsigned int maxSamplingPoints=0;
    unsigned int numOrbitsUsed=0;
    unsigned int numOrbitsReintegrate=0;
    if(numSamplingPointsMin!=NULL)
        numSamplingPointsMin->clear();
    double sumWeight=0;   // total weight of all orbits (adjusted by refinement factors)
    double totalMass=0;   // total mass of model (sum of orbit weights, unadjusted)
    // mass refinement if necessary
    vectord orbitMassFactor(OrbitList.size(), 1.0);
    if(massRefineFactor!=NULL)
    {
        for(size_t o=0; o<OrbitList.size(); o++)
        {
            orbitMassFactor[o] = massRefineFactor->eval(OrbitList[o]);
            if(orbitMassFactor[o]<=0)
            {
                return("Error, refinement factor is non-positive");
            }
            sumWeight += OrbitList[o]->getWeight() / orbitMassFactor[o];
            totalMass += OrbitList[o]->getWeight();
        }
    }
    else 
    {
        for(size_t o=0; o<OrbitList.size(); o++)
            sumWeight += OrbitList[o]->getWeight();
        totalMass=sumWeight;
    }
    if(sumWeight<=0 || totalMass<=0) 
    {
        return("Error, total orbit weight is zero, need to assign weights first (run Schwarzschild model optimization)");
        // in principle, could use equal weights in this case, but it makes little sense
    }
    if(result==NULL) result=new CPointMassSetFloat;
    if(potential->Mbh>0)
    {
        result->push_back(std::pair<CPosVelPoint<float>,float>(CPosVelPoint<float>(), static_cast<float>(potential->Mbh)));
        numPoints--;
    }
    double outSumWeight=0;  // accumulator for sum of orbit weights up to current orbit (adjusted by refinement factors)
    int outPointIndex=0;    // accumulator for sum of points emitted up to current orbit
    double outTotalMass=0;  // accumulator for total mass in points up to current orbit
    // controlling center mass position, net momentum and rotation
    double CX=0, CY=0, CZ=0, VX=0, VY=0, VZ=0, LX=0, LY=0, LZ=0, Lhm=0, Chm=0, Vhm=0;
    for(size_t o=0; o<OrbitList.size(); o++)
    {
        outSumWeight += OrbitList[o]->getWeight() / orbitMassFactor[o];
        int newPointIndex = static_cast<int>(outSumWeight/sumWeight * numPoints);
        int pointsPerOrbit = newPointIndex - outPointIndex;   // the number of sampling points necessary for an orbit
        outPointIndex=newPointIndex;
        if(pointsPerOrbit==0) continue;
        double pointMass = orbitMassFactor[o] * sumWeight / numPoints;  // individual mass of points from this orbit
        const CTrajSampleInformation<float>* trajSampleInfo = static_cast<const CTrajSampleInformation<float>*>(OrbitList[o]->getInfoByType(CBasicInformation::IT_TRAJ_SAMPLE));
        CPosVelDataFloat onePoint(1, OrbitList[o]->getInitData().initCond);  // if there is no trajectory sample stored, at least have one point corresponding to initial conditions
        const CPosVelDataFloat* trajSample = (trajSampleInfo!=NULL && !trajSampleInfo->getTraj().empty() ? &trajSampleInfo->getTraj() : &onePoint); 
        int pointsInSample=static_cast<int>(trajSample->size());

        maxSamplingPoints = std::max<unsigned int>(maxSamplingPoints, pointsPerOrbit);
        numOrbitsUsed++;
        if(pointsPerOrbit>pointsInSample)
        {
            numOrbitsReintegrate++;
            OrbitList[o]->setState(COrbitDesc::OS_INITIALIZED);
            // communicate desired number of points
            if(numSamplingPointsMin==NULL)
                numSamplingPointsMin=new std::vector<std::pair<CPosVelPoint<float>, unsigned int> >;
            numSamplingPointsMin->push_back(std::pair<CPosVelPoint<float>, unsigned int>(OrbitList[o]->getInitData().initCond, pointsPerOrbit+1));
        }
        for(int p=0; p<pointsPerOrbit; p++)
        {
            outTotalMass+=pointMass;
            CPosVelPoint<float> xv = trajSample->at(p % pointsInSample);  // loop through entire sample more than once if number of points necessary is greater than sample size
            /*{   //!!! do not emit points with very extreme kinetic to total energy ratio
            double W=potential->Phi(xv.Pos[0], xv.Pos[1], xv.Pos[2]);
            double K=(pow_2(xv.Vel[0])+pow_2(xv.Vel[1])+pow_2(xv.Vel[2]))/2;
            if(W<0 && (K+W>=0 || -K/(K+W)>sqrt(numPoints*1.0)*0.2)) continue;  // otherwise it will distort virial equilibrium too much
            }*/
            double lx=pointMass*(xv.Pos[1]*xv.Vel[2]-xv.Pos[2]*xv.Vel[1]);
            double ly=pointMass*(xv.Pos[2]*xv.Vel[0]-xv.Pos[0]*xv.Vel[2]); 
            double lz=pointMass*(xv.Pos[0]*xv.Vel[1]-xv.Pos[1]*xv.Vel[0]);
            // decide whether to flip the orbit, and about which plane.
            int nflip=0;
            double absL = fabs(LX+lx)+fabs(LY+ly)+fabs(LZ+lz) + 
                (fabs(CX+pointMass*xv.Pos[0])+fabs(CY+pointMass*xv.Pos[1])+fabs(CZ+pointMass*xv.Pos[2])) *
                (fabs(VX+pointMass*xv.Vel[0])+fabs(VY+pointMass*xv.Vel[1])+fabs(VZ+pointMass*xv.Vel[2])) / NBEXPORT_XV_FACTOR;
            double absLx= fabs(LX+lx)+fabs(LY-ly)+fabs(LZ-lz) + 
                (fabs(CX-pointMass*xv.Pos[0])+fabs(CY+pointMass*xv.Pos[1])+fabs(CZ+pointMass*xv.Pos[2])) *
                (fabs(VX-pointMass*xv.Vel[0])+fabs(VY+pointMass*xv.Vel[1])+fabs(VZ+pointMass*xv.Vel[2])) / NBEXPORT_XV_FACTOR;
            if(absLx<absL) { absL=absLx; nflip=1; }
            double absLy= fabs(LX-lx)+fabs(LY+ly)+fabs(LZ-lz) + 
                (fabs(CX+pointMass*xv.Pos[0])+fabs(CY-pointMass*xv.Pos[1])+fabs(CZ+pointMass*xv.Pos[2])) *
                (fabs(VX+pointMass*xv.Vel[0])+fabs(VY-pointMass*xv.Vel[1])+fabs(VZ+pointMass*xv.Vel[2])) / NBEXPORT_XV_FACTOR;
            if(absLy<absL) { absL=absLy; nflip=2; }
            double absLz= fabs(LX-lx)+fabs(LY-ly)+fabs(LZ+lz) + 
                (fabs(CX+pointMass*xv.Pos[0])+fabs(CY+pointMass*xv.Pos[1])+fabs(CZ-pointMass*xv.Pos[2])) *
                (fabs(VX+pointMass*xv.Vel[0])+fabs(VY+pointMass*xv.Vel[1])+fabs(VZ-pointMass*xv.Vel[2])) / NBEXPORT_XV_FACTOR;
            if(absLz<absL) { absL=absLz; nflip=3; }
            if(nflip>0)
            {   
                xv.Pos[nflip-1] *= -1;
                xv.Vel[nflip-1] *= -1; 
            }
            result->push_back(std::pair<CPosVelPoint<float>,float>(xv, static_cast<float>(pointMass)));
            // accumulate total momentum and angular momentum
            LX+=pointMass*(xv.Pos[1]*xv.Vel[2]-xv.Pos[2]*xv.Vel[1]);
            LY+=pointMass*(xv.Pos[2]*xv.Vel[0]-xv.Pos[0]*xv.Vel[2]);
            LZ+=pointMass*(xv.Pos[0]*xv.Vel[1]-xv.Pos[1]*xv.Vel[0]);
            CX+=pointMass*xv.Pos[0];
            CY+=pointMass*xv.Pos[1];
            CZ+=pointMass*xv.Pos[2];
            VX+=pointMass*xv.Vel[0];
            VY+=pointMass*xv.Vel[1];
            VZ+=pointMass*xv.Vel[2];
        }
        if(outTotalMass>0.5*totalMass && Lhm==0 && Vhm==0 && Chm==0)
        {   // store half-mass values
            Chm = sqrt(CX*CX+CY*CY+CZ*CZ)/outTotalMass;
            Vhm = sqrt(VX*VX+VY*VY+VZ*VZ)/outTotalMass;
            Lhm = sqrt(LX*LX+LY*LY+LZ*LZ)/outTotalMass;
        }
    }
    CX/=outTotalMass;
    CY/=outTotalMass;
    CZ/=outTotalMass;
    VX/=outTotalMass;
    VY/=outTotalMass;
    VZ/=outTotalMass;
    LX/=outTotalMass;
    LY/=outTotalMass;
    LZ/=outTotalMass;
    for(size_t i=0; i<result->size(); i++)
    {   // cancel total momentum
        result->at(i).first.Vel[0] -= static_cast<float>(VX);
        result->at(i).first.Vel[1] -= static_cast<float>(VY);
        result->at(i).first.Vel[2] -= static_cast<float>(VZ);
    }
    // return textual representation of some statistics
    std::string info = convertToString(result->size())+" points created\n" +
        convertToString(outTotalMass) + " total mass\n" + 
        convertToString(numOrbitsUsed)+" orbits used\n" +
        "max "+convertToString(maxSamplingPoints)+" sampling points per orbit\n" + 
        "Avg.X="+convertToString(sqrt(CX*CX+CY*CY+CZ*CZ),3) + " (" + convertToString(Chm,3) + 
        ")\nAvg.V="+convertToString(sqrt(VX*VX+VY*VY+VZ*VZ),3) + " (" + convertToString(Vhm,3) + ") [cancelled]" +
        "\nAvg.L="+convertToString(sqrt(LX*LX+LY*LY+LZ*LZ),3) + " (" + convertToString(Lhm,3) + ")";
    return info;
}

void CPointMassSetHandler::computeVirialRatio(double* _T, double* _W)
{
    state=PS_RUNNING;
    CPotentialNB treePot(3, 0 /*Mbh*/, 0 /*eps*/, 0.5 /*theta*/, *points);
    double T=0, W=0;
    for(size_t ind=0; ind<points->size() && state!=PS_NEEDTOTERMINATE; ind++)
    {
        const CPosVelPoint<float>& pt=points->at(ind).first;
        double v2 = pow_2(pt.Vel[0])+pow_2(pt.Vel[1])+pow_2(pt.Vel[2]);
        T += points->at(ind).second * v2/2.;
        W -= points->at(ind).second * treePot.Phi(pt.Pos[0], pt.Pos[1], pt.Pos[2])/2.;
    }
    state=PS_DONE;
    *_T=T; *_W=W;
}

/**** orbit evaluating function for mass refinement in exportNbody ****/
CMassRefinementFnc::CMassRefinementFnc(const smile::COrbitLibrary *orbitlib, unsigned int numBins, const CBasicOrbitFilteringFnc* _rankFnc) :
  rankFnc(_rankFnc)
{
    // divide all orbits in orbitlib into several bins, based on the value returned by rank evaluating function; mass factors of adjacent bins differ by a factor of two
    unsigned int numOrbits=orbitlib->size();
    if(numOrbits==0 || numBins==0) return;   // not an error, but no ranking will be performed
    vectord orbitRank(numOrbits);
    for(unsigned int o=0; o<numOrbits; o++)
        orbitRank[o] = rankFnc->eval(orbitlib->getOrbitDesc(o));
    std::sort(orbitRank.begin(), orbitRank.end());
    binBoundary.resize(numBins+1);
    for(unsigned int bin=0; bin<=numBins; bin++)
        binBoundary[bin] = orbitRank[ static_cast<unsigned int>((numOrbits-1) / pow(2.0, 1.0*bin)) ];
}

double CMassRefinementFnc::eval(const smile::COrbitDesc *orbit) const
{
    if(binBoundary.empty() || orbit==NULL) return 1;
    double orbitRank=rankFnc->eval(orbit);
    unsigned int numBins=static_cast<unsigned int>(binBoundary.size());
    double factor=0;  // shouldn't remain zero, signals an error
    for(unsigned int bin=0; bin<numBins; bin++)
        if(orbitRank <= binBoundary[bin]) factor=pow(2.0, numBins-1.0-bin);
    return factor;
}

/**** an example of orbit rank evaluation based on total energy ****/
double CEnergyOrbitFilteringFnc::eval(const smile::COrbitDesc *orbit) const
{
    const COrbitInformation<float>* info=static_cast<const COrbitInformation<float>*>(orbit->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
    if(info!=NULL) 
        return info->getEinit();
    else return 0;   // shouldn't happen because this function is supposed to be called when exporting Nbody model, which obviously should be done after orbit integration
}

}  // namespace