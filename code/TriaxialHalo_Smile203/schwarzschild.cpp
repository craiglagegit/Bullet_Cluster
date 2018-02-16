#include "schwarzschild.h"
#include "potential.h"
#include "orbitlib.h"
#include "orbit.h"
#include "stringconv.h"

#include <cmath>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <cassert>

namespace smile{

/// class for recording time spent in each shell during orbit integration
COrbitRuntimeBasicShellSchw::COrbitRuntimeBasicShellSchw(const COrbit* _orbit, const CBasicShellSchwModel* _model): 
CBasicOrbitRuntimeFnc(_orbit), model(_model)
{
    assert(orbit->getPotential()->N_dim==N_DIM);
    shellTime.assign(model->getNumShellsKinem()+1, 0); 
    shellVr.assign(model->getNumShellsKinem()+1, 0); 
    shellVt.assign(model->getNumShellsKinem()+1, 0); 
    COrbitInitData<double> ID=orbit->getInitData();
    shellprev = model->whichShellR(sqrt(pow_2(ID.initCond.Pos[0])+pow_2(ID.initCond.Pos[1])+pow_2(ID.initCond.Pos[2])));
    timeAccuracy = 1e-6*ID.timeStep;
}

void COrbitRuntimeBasicShellSchw::Timestep(const double told, const double t, const double y[])
{
    double r2=pow_2(y[0]) + pow_2(y[1]) + pow_2(y[2]);   // assume N_dim==3
    double Lx=y[1]*y[5]-y[2]*y[4];
    double Ly=y[2]*y[3]-y[0]*y[5];
    double Lz=y[0]*y[4]-y[1]*y[3];
    double L=sqrt(pow_2(Lx)+pow_2(Ly)+pow_2(Lz));
    double vr2= pow_2(y[0]*y[3]+y[1]*y[4]+y[2]*y[5])/r2;
    double vt2= pow_2(L)/r2;
    if(r2==0) { vr2=pow_2(y[0])+pow_2(y[1])+pow_2(y[2]); vt2=0; }
    unsigned int shell=model->whichShellR(sqrt(r2));
    if(shell==shellprev)
    {
        shellTime[shell] += (t-told);
        shellVr[shell] += (t-told)*vr2;
        shellVt[shell] += (t-told)*vt2;
    }
    else  // calculate fractional time spent in shell
        calcShellTime(told, t, shellprev, shell, vr2, vt2);
    shellprev=shell;
}

void COrbitRuntimeBasicShellSchw::calcShellTime(double tlow, double tupp, unsigned int slow, unsigned int supp, double vr2, double vt2)
{
    double tl=tlow, tu=tupp;
    int niter=0;
    double tmid;
    do{
        tmid=(tu+tl)/2;
        unsigned int smid=model->whichShellR(sqrt(pow_2(orbit->getInterpolatedTrajectory(0,tmid)) + pow_2(orbit->getInterpolatedTrajectory(1,tmid)) + pow_2(orbit->getInterpolatedTrajectory(2,tmid))));
        if(smid==slow)
            tl=tmid;
        else if(smid==supp)
            tu=tmid;
        else
        {   // neither of the two shells - go deeper into recursion
            calcShellTime(tlow, tmid, slow, smid, vr2, vt2);
            calcShellTime(tmid, tupp, smid, supp, vr2, vt2);
            return;
        }
        niter++;
    }
    while((tu-tl>timeAccuracy) && (niter<20));
    shellTime[slow] += tmid-tlow;
    shellTime[supp] += tupp-tmid;
    shellVr[slow]   += (tmid-tlow)*vr2;
    shellVr[supp]   += (tupp-tmid)*vr2;
    shellVt[slow]   += (tmid-tlow)*vt2;
    shellVt[supp]   += (tupp-tmid)*vt2;
}

void COrbitRuntimeBasicShellSchw::Finish()
{
    double intTime=orbit->getIntTime();
    if(intTime>0)
    {
        for(size_t s=0; s<shellTime.size(); s++)
        {
            if(shellTime[s]>0)
            {
                shellVr[s] /= shellTime[s];
                shellVt[s] /= shellTime[s];
            }
            shellTime[s] /= intTime;    // normalize to unit time
        }
    }
    shellVr.pop_back();  // eliminate last shell which goes to infty
    shellVt.pop_back();  
    shellTime.pop_back();
}

/// class for recording time spent in each cell during orbit integration for classic Schwarzschild modelling
COrbitRuntimeSchwClassic::COrbitRuntimeSchwClassic(const COrbit* _orbit, const CSchwModelClassic* _model): 
COrbitRuntimeBasicShellSchw(_orbit, _model), model(_model)
{
    cellTime.assign(model->getNumCoefsDens()+1, 0);
    COrbitInitData<double> ID=orbit->getInitData();
    cellprev = model->whichCell(ID.initCond.Pos[0],ID.initCond.Pos[1],ID.initCond.Pos[2]);
}

void COrbitRuntimeSchwClassic::Timestep(const double told, const double t, const double y[])
{
    COrbitRuntimeBasicShellSchw::Timestep(told, t, y);
    unsigned int cell=model->whichCell(y[0], y[1], y[2]);
    if(cell==cellprev)
    {
        cellTime[cell] += (t-told);
    }
    else  // calculate fractional time spent in cell
        calcCellTime(told, t, cellprev, cell);
    cellprev=cell;
}

void COrbitRuntimeSchwClassic::calcCellTime(double tlow, double tupp, unsigned int clow, unsigned int cupp)
{
    double tl=tlow, tu=tupp;
    int niter=0;
    double tmid;
    do{
        tmid=(tu+tl)/2;
        unsigned int cmid=model->whichCell(orbit->getInterpolatedTrajectory(0,tmid), orbit->getInterpolatedTrajectory(1,tmid), orbit->getInterpolatedTrajectory(2,tmid));  // assume that N_dim=3
        if(cmid==clow)
            tl=tmid;
        else if(cmid==cupp)
            tu=tmid;
        else
        {   // neither of the two cells - go deeper into recursion
            calcCellTime(tlow, tmid, clow, cmid);
            calcCellTime(tmid, tupp, cmid, cupp);
            return;
        }
        niter++;
    }
    while((tu-tl>timeAccuracy) && (niter<20));
    cellTime[clow]  += tmid-tlow;
    cellTime[cupp]  += tupp-tmid;
}

void COrbitRuntimeSchwClassic::Finish()
{
    COrbitRuntimeBasicShellSchw::Finish();
    double intTime=orbit->getIntTime();
    if(intTime>0)
        for(size_t c=0; c<model->getNumCoefsDens(); c++)
            cellTime[c] /= intTime;    // normalize to unit time
    cellTime.pop_back();  // eliminate last cell which goes to infty
}

/// class for recording density data in terms of spherical-harmonic expansion
COrbitRuntimeSchwSH::COrbitRuntimeSchwSH(const COrbit* _orbit, const CSchwModelSH* _model) : 
  COrbitRuntimeBasicShellSchw(_orbit, _model), model(_model) 
{
    timeStep=orbit->getInitData().timeStep;
}

void COrbitRuntimeSchwSH::Timestep(const double told, const double t, const double y[])
{
    COrbitRuntimeBasicShellSchw::Timestep(told, t, y);
    if(t-told<timeStep)
        traj.push_back(std::pair<CPosVelPoint<double>,double>(CPosVelPoint<double>(y[0], y[1], y[2], y[3], y[4], y[5]), t-told));
    else
    {
        int numIntervals=static_cast<int>(ceil((t-told)/timeStep));
        for(int i=numIntervals-1; i>=0; i--)
        {
            double t1=t-i*(t-told)/numIntervals;
            traj.push_back(std::pair<CPosVelPoint<double>,double>(CPosVelPoint<double>(
                orbit->getInterpolatedTrajectory(0, t1),
                orbit->getInterpolatedTrajectory(1, t1),
                orbit->getInterpolatedTrajectory(2, t1),
                0, 0, 0), (t-told)/numIntervals));
        }
    }
}

void COrbitRuntimeSchwSH::Finish()
{
    COrbitRuntimeBasicShellSchw::Finish();
    double intTime=orbit->getIntTime();
    if(intTime>0) 
        for(size_t i=0; i<traj.size(); i++)
            traj[i].second /= intTime;
    SHcoefs = model->computeDensityData(traj);
    traj.clear();
}

/**--------------  CBasicShellSchwModel class - parent class for variants of Schw.modelling approaches with shell partitioning of configuration space -------------*/

CBasicShellSchwModel::CBasicShellSchwModel(const CPotential* potential, unsigned int numShells, double innerShellMass, double outerShellMass) :
  CBasicSchwModel()
{
    shellRadius.assign(numShells, 0);
    shellEnergy.assign(numShells, 0);
    if(outerShellMass==0) outerShellMass=1.0-1.0/(numShells+1);
    if(innerShellMass==0) innerShellMass=outerShellMass/numShells;
    createNonuniformGrid(massInteriorToShell, numShells, innerShellMass, outerShellMass, false); 
    totalMass = potential->totalMass();
    if(totalMass<=0) { my_error("Error, density model does not seem to have finite mass!"); totalMass=1; }
    for(unsigned int s=0; s<numShells; s++)
    {
        massInteriorToShell[s] *= totalMass;
        shellRadius[s] = potential->getRadiusByMass(massInteriorToShell[s]);
        shellEnergy[s] = potential->Phi(shellRadius[s], 0, 0);
    }
}

int CBasicShellSchwModel::solveOptimizationProblem(COrbitLibrary* orbits, const CBasicOptimizationSolver* solver, const CBasicOrbitFilteringFnc* orbitPenaltyFnc) const
{
    if(orbits==NULL || orbits->size()==0) return -1;
    // internal arrays used in optimization procedure
    std::vector<vectorn> linearMatrix;    // matrix M of linear equations for the optimization problem  "M w = rhs"
    vectorn rhs;                          // the rhs of equations
    vectorn constraintWeight;             // 'importance' of linear constraints (in the form of penalty for slack variables L, M)
    vectorn orbitWeight;                  // the solution "w"
    vectorn weightFactor;                 // penalties for individual orbits (=chaoticWeightFactor/numOrbits for chaotic orbits, 0 for regular ones)
    // 1st stage: preparation
    my_message("Preparing model");
    try{
        orbitWeight.clear();
        weightFactor.clear();
        linearMatrix.clear();
        const size_t numCells=densityData.size();  // shortcuts
        const size_t numShells=shellRadius.size();
        vectorn row(1 + numCells + (configSchw.sm_constrainBeta ? numShells : 0));
        for(unsigned int o=0; o<orbits->size(); o++)
        {
            // load grid information data from orbit desc
            const CSchwInformation* infoGrid=static_cast<const CSchwInformation*> (orbits->getOrbitDesc(o)->getInfoByType(CBasicInformation::IT_SCHW));
            const COrbitInformation<float>* infoOrb=static_cast<const COrbitInformation<float>*>(orbits->getOrbitDesc(o)->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
            if(infoGrid==NULL || infoOrb==NULL)
            {
                my_error("Error, missing necessary orbit-cell-weight data");
                return false;
            }
            if(infoGrid->getDensityData().size() != numCells || infoGrid->getShellVr().size() != numShells || infoGrid->getShellVt().size() != numShells)
            {
                my_error("Error, grid dimensions are incompatible with orbit-cell-weight data");
                return false;
            }
            for(size_t c=0; c<numCells; c++)
                row[c]=static_cast<smnumtype>(infoGrid->getDensityData(c) / getNormFactor(c));
            if(configSchw.sm_constrainBeta)
            {
                for(size_t s=0; s<numShells; s++)
                {
                    double beta = configSchw.sm_betaIn + (configSchw.sm_betaOut-configSchw.sm_betaIn)*massInteriorToShell[s]/totalMass;
                    row[numCells+s] = static_cast<smnumtype>((infoGrid->getShellVt(s)-2*(1-beta)*infoGrid->getShellVr(s)) * infoGrid->getShellTime(s));
                }
            }
            row.back()=1.0;  // sum of all orbit weights should be equal to total model mass
            linearMatrix.push_back(row);
            if(orbitPenaltyFnc!=NULL)
                weightFactor.push_back(static_cast<smnumtype>(orbitPenaltyFnc->eval(orbits->getOrbitDesc(o))));
            orbitWeight.push_back(static_cast<smnumtype>(orbits->getOrbitDesc(o)->getWeight()));
        }
        // initialize rhs
        rhs.clear();
        for(size_t c=0; c<numCells; c++)
            rhs.push_back(static_cast<smnumtype>(densityData[c] / getNormFactor(c)));
        constraintWeight.assign(rhs.size(), 1.0);
        if(configSchw.sm_constrainBeta)
            for(size_t sh=0; sh<numShells; sh++)
            {
                rhs.push_back(0.0);
                constraintWeight.push_back(static_cast<smnumtype>(SCHW_CONSTRAINT_BETA_WEIGHT));
            }
        rhs.push_back(static_cast<smnumtype>(orbits->getPotential()->totalMass()));  // sum of all orbit weights should be equal to total model mass
        constraintWeight.push_back(1.0);
    }
    catch(const std::bad_alloc&) {
        my_error("Not enough memory in preparing model");
        return -1;
    }
    // 2nd stage: optimization itself
    my_message("Calling solver");
    int result = solver->callSolver(linearMatrix, rhs, constraintWeight, weightFactor, configSchw.sm_maxWeight, &orbitWeight);
    // 3rd stage: post-processing
    my_message("Finishing model");
    for(unsigned int o=0; o<orbits->size(); o++)
    {
        if(orbitWeight[o] < SCHW_MIN_ORBIT_WEIGHT_USED/orbits->size()) 
            orbitWeight[o]=0;
        orbits->setOrbitWeight(o, orbitWeight[o]);
    }
    return result;
}

unsigned int CBasicShellSchwModel::whichShellR(double radius) const
{
    unsigned int shnum=0;
    while(shnum<shellRadius.size() && radius>=shellRadius[shnum]) shnum++;
    return shnum;
}

unsigned int CBasicShellSchwModel::whichShellE(double energy) const
{
    unsigned int shnum=0;
    while(shnum<shellRadius.size() && energy>=shellEnergy[shnum]) shnum++;
    return shnum;
}

std::string CBasicShellSchwModel::getStatistics(const COrbitLibrary* orbits) const
{
    if(orbits==NULL || orbits->size()==0) return "Error in getStatistics: no orbits!";
    unsigned int numOrbitsNonzero=0;   // number of orbits with non-zero weight
    unsigned int numCoefsInfeasible=0; // number of coefficients that did not meet requirement values
    double totalDensityPenalty=0;      // sum of densityPenalty for all infeasible coefs
    double entropy=0;
    double totalWeight=0;
    std::vector<double> weights(orbits->size(), 0);
    for(unsigned int o=0; o<orbits->size(); o++)
    {
        double orbitWeight=orbits->getOrbitDesc(o)->getWeight();
        totalWeight+=orbitWeight;
        if(orbitWeight>0) numOrbitsNonzero++;
        weights[o]=orbitWeight;
    }
    std::sort(weights.begin(), weights.end(), std::greater<double>());
    double avgWeight=totalWeight/orbits->size();
    double sumHighWeights=0;
    unsigned int numOrbitsHalfMass=0;
    if(avgWeight>0) {
        for(unsigned int o=0; o<orbits->size(); o++)  // compute entropy assuming flat (uniform) priors for weights
        {
            double orbitWeight=orbits->getOrbitDesc(o)->getWeight();
            if(orbitWeight>0) entropy -= orbitWeight*log(orbitWeight/avgWeight);
            sumHighWeights+=weights[o];
            if(sumHighWeights>=totalWeight*0.5 && numOrbitsHalfMass==0) numOrbitsHalfMass=o;
        }
        entropy /= totalWeight;
    }
    double maxToAvgWeight = weights.front()/avgWeight;
    double fracOrbitsHalfMass = numOrbitsHalfMass*1.0/orbits->size();

    vectord densityPenalty(calcPenalty(orbits));
    for(unsigned int c=0; c<static_cast<unsigned int>(densityPenalty.size()); c++)
        if(densityPenalty[c]!=0) 
        {
            numCoefsInfeasible++;
            totalDensityPenalty+=densityPenalty[c]/getNormFactor(c);
        }
    return
        convertToString(numOrbitsNonzero)+" orbits out of "+convertToString(orbits->size())+" with nonzero weights\n" +
        (totalWeight>0 ?   "Total weight: " + convertToString(totalWeight) + "\n" +
        convertToString(fracOrbitsHalfMass*100, 3) + "% orbits make up 50% total mass" +
        (fracOrbitsHalfMass>0.2 ? " (good)" : fracOrbitsHalfMass>0.1 ? " (fair)" : " (poor)") + "\n" +
        "entropy: " + convertToString(entropy, 3) + 
        (entropy>=-0.5 ? " (good)" : entropy>=-1.0 ? " (fair)" : " (poor)") + "\n" +
        "Max/avg weight: " + convertToString(maxToAvgWeight, 3) + 
        (maxToAvgWeight<=10 ? " (good)" : maxToAvgWeight<=50 ? " (fair)" : " (poor)") + "\n" : "" ) +
        (numCoefsInfeasible>0 ? "WARNING: " +
        convertToString(numCoefsInfeasible)+" coefs out of "+convertToString(densityData.size())+" infeasible\n" +
        convertToString(totalDensityPenalty, 3)+" total difference" : "Model is feasible");
}

vectord CBasicShellSchwModel::calcPenalty(const COrbitLibrary* orbits) const
{
    vectord densityPenalty(densityData.size(), 0);
    for(unsigned int o=0; o<orbits->size(); o++)
    {
        // load grid information data from orbit desc
        const CSchwInformation* infoGrid=static_cast<const CSchwInformation*> (orbits->getOrbitDesc(o)->getInfoByType(CBasicInformation::IT_SCHW));
        if(infoGrid!=NULL && infoGrid->getDensityData().size()==densityData.size())
        {
            double weight=orbits->getOrbitDesc(o)->getWeight();
            for(size_t c=0; c<densityData.size(); c++)
                densityPenalty[c] += weight * infoGrid->getDensityData(c);
        }
    }
    for(unsigned int c=0; c<densityData.size(); c++)
    {
        densityPenalty[c] -= densityData[c];
        if(fabs(densityPenalty[c]/densityData[c])<=SCHW_MAX_CELLMASS_REL_DIFFERENCE) 
            densityPenalty[c]=0;
    }
    return densityPenalty;
}

/// Class for performing Schwarzschild modelling

CSchwModelClassic::CSchwModelClassic(const CPotential* potential, unsigned int _numShells, double innerShellMass, double outerShellMass, unsigned int _linesPerSegment): 
  CBasicShellSchwModel(potential, _numShells, innerShellMass, outerShellMass),
  linesPerSegment(_linesPerSegment), 
  cellsPerSegment(linesPerSegment*linesPerSegment),
  cellsPerShell(3*cellsPerSegment),
  numCells(_numShells*cellsPerShell)
{
    cellBoundary.resize(linesPerSegment+1);
    for(size_t i=0; i<linesPerSegment; i++)
        cellBoundary[i] = tan(M_PI/4/linesPerSegment*i);
    cellBoundary[linesPerSegment] = 1.0;  // exact value
    initCellMass(potential);
}

void CSchwModelClassic::initCellMass(const CPotential* potential)
{
    // calculate mass in cells
    densityData.assign(numCells, 0);
    if(potential->PotentialType()!=CDensity::PT_NB)
    for(unsigned int c=0; c<numCells; c++)
    {
        unsigned int nshell = c / cellsPerShell;
        unsigned int indseg = c % cellsPerShell;
        unsigned int nseg = indseg / cellsPerSegment;
        unsigned int ind1 =(indseg % cellsPerSegment) % linesPerSegment;
        unsigned int ind2 =(indseg % cellsPerSegment) / linesPerSegment;
        double dr1 = cellBoundary[ind1+1] - cellBoundary[ind1];
        double dr2 = cellBoundary[ind2+1] - cellBoundary[ind2];
        double mass=0;
        for(unsigned int i1=0; i1<NUM_SUBDIV_ANGULAR; i1++)
            for(unsigned int i2=0; i2<NUM_SUBDIV_ANGULAR; i2++)
            {
                double r1l = cellBoundary[ind1] + dr1*i1/NUM_SUBDIV_ANGULAR;
                double r1h = cellBoundary[ind1] + dr1*(i1+1)/NUM_SUBDIV_ANGULAR;
                double r2l = cellBoundary[ind2] + dr2*i2/NUM_SUBDIV_ANGULAR;
                double r2h = cellBoundary[ind2] + dr2*(i2+1)/NUM_SUBDIV_ANGULAR;
                double phill = asin(r1l*r2l/sqrt((1+pow_2(r1l))*(1+pow_2(r2l))));
                double philh = asin(r1l*r2h/sqrt((1+pow_2(r1l))*(1+pow_2(r2h))));
                double phihl = asin(r1h*r2l/sqrt((1+pow_2(r1h))*(1+pow_2(r2l))));
                double phihh = asin(r1h*r2h/sqrt((1+pow_2(r1h))*(1+pow_2(r2h))));
                double area = (phill+phihh-philh-phihl);
                double ratio1 = (r1l+r1h)/2;
                double ratio2 = (r2l+r2h)/2;
                double r=0, rprev, rho=0, rhoprev;
                for(unsigned int ir=0; ir<=NUM_SUBDIV_RADIAL; ir++)
                {
                    rprev=r; rhoprev=rho;
                    double shellRprev=nshell==0 ? 0 : shellRadius[nshell-1];
                    r = shellRprev + (shellRadius[nshell] - shellRprev) * ir / NUM_SUBDIV_RADIAL;
                    if(r==0) r=shellRadius[0]*0.01/NUM_SUBDIV_RADIAL;
                    double coord0 = r/sqrt(1 + pow_2(ratio1) + pow_2(ratio2));
                    double coord1 = ratio1*coord0;
                    double coord2 = ratio2*coord0;
                    double X=0, Y=0, Z=0;
                    switch(nseg)
                    {
                    case 0:
                        X=coord0;
                        Y=coord1;
                        Z=coord2;
                        break;
                    case 1:
                        Y=coord0;
                        Z=coord1;
                        X=coord2;
                        break;
                    case 2:
                        Z=coord0;
                        X=coord1;
                        Y=coord2;
                        break;
                    }
                    rho = potential->Rho(X, Y, Z);
                    if(ir>0)
                    {
                        double gamma = - log(rho/rhoprev)/log(r/rprev);
                        if(!gsl_finite(gamma)) gamma=0; // avoid NaN
                        if(gamma<0)
                            gamma=0; 
                        double intrad = rho*pow(r,gamma) * (gamma!=3 ? (pow(r, 3-gamma)-pow(rprev, 3-gamma))/(3-gamma) : log(r/rprev));
                        if(!gsl_finite(intrad)) intrad=0;
                        mass +=  area * intrad * 8;
                    }
                }
            }
        densityData[c] = mass;
    }
    else  // frozen Nbody potential - simply count bodies in each cell
    {
        const CPotentialNB* potentialNB=static_cast<const CPotentialNB*>(potential);
        for(size_t ind=0; ind<potentialNB->bodyNum(); ind++)
        {
            CPosPoint<double> pt=potentialNB->bodyPos<double>(ind);
            unsigned int cell=whichCell(pt.Pos[0], pt.Pos[1], pt.Pos[2]);
            double ptmass=potentialNB->bodyMass(ind);
            if(cell<numCells) densityData[cell] += ptmass;
        }
    }
}

void CSchwModelClassic::cellCenter(const unsigned int C, double& X, double& Y, double& Z) const
{
    unsigned int nshell = C / cellsPerShell;
    if(nshell>=shellRadius.size()) 
    {
        X=0; Y=0; Z=0;
        return;
    }
    unsigned int indseg = C % cellsPerShell;
    unsigned int nseg = indseg / cellsPerSegment;
    unsigned int ind1 =(indseg % cellsPerSegment) % linesPerSegment;
    unsigned int ind2 =(indseg % cellsPerSegment) / linesPerSegment;
    double ratio1 = (cellBoundary[ind1] + cellBoundary[ind1+1])/2;
    double ratio2 = (cellBoundary[ind2] + cellBoundary[ind2+1])/2;
    double r = ((nshell>0 ? shellRadius[nshell-1] : 0) + shellRadius[nshell])/2;
    double coord0 = r/sqrt(1 + pow_2(ratio1) + pow_2(ratio2));
    double coord1 = ratio1*coord0;
    double coord2 = ratio2*coord0;
    switch(nseg)
    {
    case 0:
        X=coord0;
        Y=coord1;
        Z=coord2;
        break;
    case 1:
        Y=coord0;
        Z=coord1;
        X=coord2;
        break;
    case 2:
        Z=coord0;
        X=coord1;
        Y=coord2;
        break;
    }
    if(whichCell(X,Y,Z) != C)
        my_error("Error in determining cell center, shouldn't happen!");
}

unsigned int CSchwModelClassic::whichCell(const double _X, const double _Y, const double _Z) const
{
    double X=fabs(_X), Y=fabs(_Y), Z=fabs(_Z);
    int nseg=-1;
    unsigned int nshell=0;
    double coord0=0, coord1=0, coord2=0;
    if((X>=Y)&&(X>Z))
    {
        nseg=0;
        coord0=X;
        coord1=Y;
        coord2=Z;
    } else
    if((Y>=Z)&&(Y>X))
    {
        nseg=1;
        coord0=Y;
        coord1=Z;
        coord2=X;
    } else
    if((Z>=X)&&(Z>=Y))
    {
        nseg=2;
        coord0=Z;
        coord1=X;
        coord2=Y;
    }
    if(nseg==-1)
    {
        my_error("Error in determining cell index, shouldn't happen!");
        return 0;  // should not happen, but who knows?...
    }
    double ratio1 = coord1/coord0;
    double ratio2 = coord2/coord0;
    unsigned int ind1=0, ind2=0;
    for(unsigned int i=1; i<linesPerSegment; i++)
    {
        if(ratio1>=cellBoundary[i]) ind1=i;
        if(ratio2>=cellBoundary[i]) ind2=i;
    }
    unsigned int indseg = ind1 + ind2*linesPerSegment + nseg*cellsPerSegment;
    nshell=whichShellR(sqrt(pow_2(X)+pow_2(Y)+pow_2(Z)));
    unsigned int ncell = nshell*cellsPerShell + indseg;
    return (ncell<numCells)?ncell:numCells;
}

double CSchwModelClassic::getNormFactor(unsigned int index) const
{ 
    unsigned int nshell=index/cellsPerShell; 
    return (massInteriorToShell[nshell]-(nshell>0?massInteriorToShell[nshell-1]:0))/cellsPerShell; 
}

//-------------------------------------------------------------//
CSchwModelSHGrid::CSchwModelSHGrid(const smile::CPotential *potential, unsigned int _numShells, double innerShellMass, double outerShellMass, unsigned int _Nang) :
  CSchwModelSH(potential, _numShells, innerShellMass, outerShellMass, _Nang),
  numGrid(_numShells)   // for simplicity, assign equal number of grid nodes for velocity dispersion bins and potential expansion grid, however, this restriction may be lifted in future
{
    std::vector< vectord > tmpCoefs;
    GridRadii.push_back(0);  // additional grid node inserted at origin
    GridRadii.insert(GridRadii.end(), shellRadius.begin(), shellRadius.end());  // copy from grid of velocity bins
    if(potential->PotentialType()==CPotential::PT_SPLINE)
        static_cast<const CPotentialSpline*>(potential)->getCoefs(&GridRadii, &tmpCoefs, false);   // HACK!! to speed up initialization
    else
    {
        CPotentialSpline pot(N_DIM, 0, numGrid, Nang, potential, &GridRadii);
        pot.getCoefs(&GridRadii, &tmpCoefs, true);
    }
    assert(GridRadii.size()==numGrid+1 && tmpCoefs.size()==numGrid+1);
    // take only even l,m>=0 coefficients
    unsigned int nAngCoefs = getNumCoefsAtRadius();
    densityData.assign(1 + numGrid * nAngCoefs, 0);  // one coef for 0th node, nAngCoefs coefs for each of other nodes
    densityData[0] = tmpCoefs[0][0];
    for(unsigned int radInd=0; radInd<numGrid; radInd++)
        for(unsigned int l=0; l<=Nang/2; l++)
            for(unsigned int m=0; m<=l; m++)
            {
                unsigned int tmpInd=2*l*(2*l+1)+2*m;
                densityData[1 + radInd*nAngCoefs + l*(l+1)/2 + m] = 
                    tmpInd<tmpCoefs[radInd+1].size() ? tmpCoefs[radInd+1][tmpInd] : 0;  // may happen because of hack above, i.e. tmpcoefs array may have less items that required by the model
            }
}

vectord CSchwModelSHGrid::computeDensityData(const smile::CPointMassSetDouble &traj) const
{
    std::vector< vectord > tmpCoefs;
    unsigned int nAngCoefs = getNumCoefsAtRadius();
    vectord data(1 + numGrid * nAngCoefs, 0);
    if(traj.size()<=1) return data;  // zerozero
    CPotentialSpline pot(Nang, traj, &GridRadii, &tmpCoefs);
    assert(tmpCoefs.size()==numGrid+1);
    // take only even l,m>=0 coefficients
    data[0] = tmpCoefs[0][0];
    for(unsigned int radInd=0; radInd<numGrid; radInd++)
        for(unsigned int l=0; l<=Nang/2; l++)
            for(unsigned int m=0; m<=l; m++)
            {
                data[1 + radInd*nAngCoefs + l*(l+1)/2 + m] = tmpCoefs[radInd+1][2*l*(2*l+1)+2*m];
            }
    return data;
}

double CSchwModelSHGrid::getNormFactor(unsigned int index) const
{
    if(index==0) return densityData[0];
    unsigned int nAngCoefs = getNumCoefsAtRadius();
    unsigned int n=(index-1)/nAngCoefs;
    unsigned int lm=(index-1)%nAngCoefs;
    int l=static_cast<int>(sqrt(2.*lm-floor(sqrt(2.*lm))));
    double norm=0;
    for(int m=0; m<=l; m++)
        norm+=pow_2(densityData[1 + n*nAngCoefs + l*(l+1)/2 + m]);
    if(norm>1e-20) return sqrt(norm) * pow_2(l+1.); 
    else if(densityData[1+n*nAngCoefs]>1e-10) return densityData[1+n*nAngCoefs];
    else return 1;
}

//-------------------------------------------------------------//
CSchwModelSHBSE::CSchwModelSHBSE(const smile::CPotential *potential, unsigned int _numShells, double innerShellMass, double outerShellMass, unsigned int _Nang, unsigned int _Nrad, double _Alpha) :
  CSchwModelSH(potential, _numShells, innerShellMass, outerShellMass, _Nang), Nrad(_Nrad)
{
    std::vector< vectord > tmpCoefs;
    if(potential->PotentialType()==CPotential::PT_BSE && 
        (_Alpha==0 || static_cast<const CPotentialBSE*>(potential)->getAlpha()==_Alpha))
    {  // HACK: to speed up initialization, assign coefficients from existing density model instead of computing costly integrals
        static_cast<const CPotentialBSE*>(potential)->getCoefs(&tmpCoefs);
        Alpha=static_cast<const CPotentialBSE*>(potential)->getAlpha();
        // the number of coefficients in potential may be greater or less than needed in S.M.
    }
    else
    {
        CPotentialBSE pot(N_DIM, 0, _Alpha, Nrad, Nang, potential);
        Alpha=pot.getAlpha();   // if _Alpha==0, it means autodetect; here we get this autodetected value (or the original supplied value if it was >=0.5)
        pot.getCoefs(&tmpCoefs);
        assert(tmpCoefs.size()==Nrad+1);
    }
    // take only even l,m>=0 coefficients
    unsigned int nAngCoefs = getNumCoefsAtRadius();
    densityData.assign((Nrad+1) * nAngCoefs, 0);  // Nang*(Nang+1)/2 coefs for each radial basis function
    for(unsigned int radInd=0; radInd<=Nrad; radInd++)
        for(unsigned int l=0; l<=Nang/2; l++)
            for(unsigned int m=0; m<=l; m++)
            {
                unsigned int tmpInd=2*l*(2*l+1)+2*m;
                densityData[radInd*nAngCoefs + l*(l+1)/2 + m] = 
                    radInd<tmpCoefs.size() && tmpInd<tmpCoefs[radInd].size() ? tmpCoefs[radInd][tmpInd] : 0;  // may happen because of hack above, i.e. tmpcoefs array may have less items that required by the model
            }
}

vectord CSchwModelSHBSE::computeDensityData(const smile::CPointMassSetDouble &traj) const
{
    std::vector< vectord > tmpCoefs;
    unsigned int nAngCoefs = getNumCoefsAtRadius();
    vectord data((Nrad+1) * nAngCoefs, 0);
    if(traj.size()<=1) return data;  // zerozero
    CPotentialBSE pot(N_DIM, 0, Alpha, Nrad, Nang, traj);
    pot.getCoefs(&tmpCoefs);
    assert(tmpCoefs.size()==Nrad+1);
    // take only even l,m>=0 coefficients
    for(unsigned int radInd=0; radInd<=Nrad; radInd++)
        for(unsigned int l=0; l<=Nang/2; l++)
            for(unsigned int m=0; m<=l; m++)
            {
                data[radInd*nAngCoefs + l*(l+1)/2 + m] = tmpCoefs[radInd][2*l*(2*l+1)+2*m];
            }
    return data;
}

double CSchwModelSHBSE::getNormFactor(unsigned int index) const
{
    unsigned int n=index/getNumCoefsAtRadius();
    unsigned int lm=index%getNumCoefsAtRadius();
    int l=static_cast<int>(sqrt(2.*lm-floor(sqrt(2.*lm))));
    /*double norm=0;
    for(int m=0; m<=l; m++)
        norm+=pow_2(densityData[n*Nang*(Nang+1)/2 + l*(l+1)/2 + m]);
    if(norm>1e-20) return sqrt(norm);*/
    // some arbitrary, not very good heuristical approximation
    double coef_n0=pow(n/(Alpha+1)+1, -4-2*log(Alpha)/log(2.)) * pow_2(n+1.);
    if(l==0) return coef_n0;
    double A=2.8*Alpha-0.6;
    double B=-0.1*Alpha-0.024;
    double C=-0.01*Alpha-0.05;
    double D=-0.11*Alpha*Alpha+0.07*Alpha-1.9;
    return coef_n0*exp(A*l + B*n + C*l*n + D) * pow_2(l+1.);
}

//-------------------------------------------------------------//
double CChaosOrbitFilteringFnc::eval(const smile::COrbitDesc *orbit) const
{
    if(orbit==NULL) return 0;
    const COrbitInformation<float>* info=static_cast<const COrbitInformation<float>*>(orbit->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
    if(info!=NULL) 
    {
        double wf=log10(info->getlfccdiff())-log10(chaoticMinFreqDiff) + 0.5;   // `soft' transition from reg to chaos, based on freq.diff.rate
        if(wf<0) wf=0; if(wf>1) wf=1;
        if(info->getlambda()>chaoticMinLambda) wf=1;
        return chaoticWeightFactor*wf;
    }
    else return 0;
}

double CShellOrbitFilteringFnc::eval(const smile::COrbitDesc *orbit) const
{
    if(orbit==NULL) return 0;
    if(model==NULL || numShell==0) return 1;
    double E=0;
    const COrbitInformation<float>* info=static_cast<const COrbitInformation<float>*>(orbit->getInfoByType(CBasicInformation::IT_TRAJ_ANALYSIS));
    if(info!=NULL) E=info->getEinit();
    // need to find which energy shell does this orbit belong to
    if(E==0)  // need to figure out orbit energy
    {
        const COrbitInitData<float>& InitData = orbit->getInitData();
        E = InitData.potential->totalEnergy(InitData.initCond);
    }
    if((numShell>model->getNumShellsKinem() && E>model->getShellEnergy(model->getNumShellsKinem()-1)) ||
       (E<=model->getShellEnergy(numShell-1) && (numShell==1 || E>model->getShellEnergy(numShell-2)))) 
       return 1;
    else return 0;
}

CConfigSchw configSchw = 
{
    CBasicSchwModel::MT_CLASSIC,
    0,     // double sm_maxWeight
    false, // bool sm_constrainBeta
    0, 0,  // double sm_betaIn, sm_betaOut
};

/** Correspondence between Schw.model names and types **/
std::string getSchwModelNameByType(CBasicSchwModel::MODELTYPE ModelType)
{
    switch(ModelType) {
    case CBasicSchwModel::MT_CLASSIC: return CSchwModelClassic::myName();
    case CBasicSchwModel::MT_SHGRID:  return CSchwModelSHGrid::myName();
    case CBasicSchwModel::MT_SHBSE:   return CSchwModelSHBSE::myName();
    default: return "Error";
    }
};

CBasicSchwModel::MODELTYPE getSchwModelTypeByName(std::string ModelName)
{
    if(ModelName==CSchwModelClassic::myName()) return CBasicSchwModel::MT_CLASSIC;
    if(ModelName==CSchwModelSHGrid ::myName()) return CBasicSchwModel::MT_SHGRID;
    if(ModelName==CSchwModelSHBSE  ::myName()) return CBasicSchwModel::MT_SHBSE;
    return CBasicSchwModel::MT_CLASSIC;  // default
};

}