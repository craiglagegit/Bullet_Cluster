// this file defines COrbitDesc - a class with 'shortened' orbit data (initial conditions and results of integration), which is used as a basic item in orbit library
// and COrbitLibrary which keeps a collection of COrbitDesc and has functions to generate initial conditions for a set of orbits, save and load library
#pragma once
#include "common.h"

namespace smile{

class COrbitDesc
{
public:
enum ORBITSTATE{
    OS_INITIALIZED,
    OS_PREPARING,
    OS_RUNNING,
    OS_NEEDTOTERMINATE,
    OS_DONE
};
    template<typename NumT> COrbitDesc(const COrbitInitData<NumT> &_InitData, const NumT _intTime, const NumT _maxTime, const vectorRuntimeFncCreators* _creators);
    COrbitDesc(const COrbitInitData<float> &_InitData, const float _intTime, const float _maxTime, const float _weight, const vectorInformation& _info);
    ~COrbitDesc();
    void run();                   // create internal COrbit instance and run orbit integration and analysis
    void halt();                  // if running, immediately stops computation
    // get/set functions
    const COrbitInitData<float>& getInitData() const { return InitData; };
    float getIntTime() const { return intTime; };
    float getMaxTime() const { return maxTime; };
    unsigned int getInfoNum() const { return static_cast<unsigned int>(info.size()); };
    const CBasicInformation* getInfo(unsigned int index) const { return info.at(index); };  // returning a pointer, though const, to a private data member!
    const CBasicInformation* getInfoByType(CBasicInformation::INFOTYPE InfoType) const {    // convenience function to get a pointer to information of a given type, or NULL if it doesn't exist
        for(unsigned int i=0; i<info.size(); i++)
            if(info[i]->infoType()==InfoType) return info[i];
        return NULL; };
    std::string toString() const;    // returns string representation of all information associated with the orbit, plus the weights and integration time
    double getWeight() const { return orbitWeight; };
    void setWeight(double _weight) { orbitWeight=_weight; }
    ORBITSTATE getState() const { return state; };
    void setState(ORBITSTATE _state) { if((state==OS_INITIALIZED && _state==OS_PREPARING)||(state==OS_DONE && _state==OS_INITIALIZED)) state=_state; };
private:
    const COrbitInitData<float> InitData;
    COrbit* orbit;                // created in run(), may be accessed via halt() to stop computation
    float intTime, maxTime;         // integration time
    const vectorRuntimeFncCreators* creators;
    vectorInformation info;
    double orbitWeight;           // weight in Schwarzschild model
    volatile ORBITSTATE state;
};

class COrbitLibrary
{
public:
    explicit COrbitLibrary(const CPotential* _potential, const vectorRuntimeFncCreators* _creators) : potential(_potential), creators(_creators) {};
    ~COrbitLibrary();
    void prepareInitialConditionsE(unsigned int NStationary, unsigned int NPrincipalPlane, unsigned int NYalpha, unsigned int NRandom, double E, double intTimeStepsPerPeriod, double intTimeInPeriods, double intTimeMaxAdaptive, bool calcLyapunov);   // create initial conditions for given energy
    void prepareInitialConditions(unsigned int NRandom, double maxRadiusMass, double intTimeStepsPerPeriod, double intTimeInPeriods, double intTimeMaxAdaptive, bool calcLyapunov);  // create initial conditions for the whole Schwarzschild model, second parameter is fraction of total mass which gives the max radius and energy of created particles (typically ~0.9(..) )
    void modifyInitialConditions(const CPotential* _potential, const vectorRuntimeFncCreators* _creators, double intTimeStepsPerPeriod, double intTimeInPeriods, double intTimeMaxAdaptive, bool calcLyapunov, bool modifyAll);        // replace potential and integration time/timestep for each orbit
    void removeUnused();               // if integration of the whole library was terminated prematurely, delete unfinished orbits
    std::string getOrbitPopulation(const CBasicOrbitFilteringFnc* filter, const CBasicOrbitFilteringFnc* chaos, bool useWeight) const; // return a string containing orbit population (5 most abundant orbit families) within given energy shell or for entire orbit library, weighted, if necessary, by orbit weights, and for regular and chaotic orbits separately
    double totalMass() const;          // returns sum of all orbit weights (if they were assigned)
    /** tries to create Nbody model from orbit library (with trajectory sample data); 
        if number of sampling points was not enough for some orbits, put duplicate points in the output for orbits with insufficient sampling, 
        mark those orbits as unstarted and return list of orbits for reintegration 
        (an associative array numSamplingPointsMin, which relates orbit initial data to the number of sampling points).
        massRefineFnc may specify unequal mass refinement (the mass of points sampled from an orbit is proportional to 
        the number returned by this function), if no function given then just assign equal mass to all points.
        Return some statistics about number of orbits and points.
    **/
    std::string exportNbody(unsigned int numPoints, 
        const CBasicOrbitFilteringFnc* massRefineFactor, 
        CPointMassSetFloat* &result,
        CPointCountSetFloat* &numSamplingPointsMin
        ) const;  
    // get functions
    const CPotential* getPotential() const { return potential; };
    unsigned int size() const { return static_cast<unsigned int>(OrbitList.size()); };
    const COrbitDesc* getOrbitDesc(unsigned int index) const { return OrbitList.at(index); }; ///!!! not very good to return even a const pointer to private member...
    unsigned int numComplete() const;  // finds how many orbits are already finished
    int findUnstartedOrbit() const;    // finds first orbit with state=OS_INITIALIZED and set it to OS_PREPARING; if none exist return -1. Should be called with external lock in multi-threaded environment.
    // set functions
    void addOrbitDesc(COrbitDesc* const &od) { OrbitList.push_back(od); };
    void setOrbitWeight(unsigned int index, double weight) { OrbitList.at(index)->setWeight(weight); };
    void runOrbit(unsigned int index) { OrbitList.at(index)->run(); }; // integrates orbit with given index
    void halt();
private:
    const CPotential* potential;       // potential for all orbits
    const vectorRuntimeFncCreators* creators;
    std::vector<COrbitDesc*> OrbitList; // array of shortened orbit info
};
/*
template<typename NumT> class CPointMassSet
{
public:
    CPointMassSet() {};
    explicit CPointMassSet(CPotential* _potential) {};
    CPointMassSet( const CPointMassSet<float> &src) {
        PosVelData.clear();
        MassData.clear();
        for(std::vector< CPosVelPoint<float> >::const_iterator iter=src.PosVelData.begin(); iter!=src.PosVelData.end(); ++iter) 
            PosVelData.push_back(*iter);
        for(std::vector< float >::const_iterator iter=src.MassData.begin(); iter!=src.MassData.end(); ++iter) 
            MassData.push_back(*iter);
    }
    CPointMassSet( const CPointMassSet<double> &src) {
        PosVelData.clear();
        for(std::vector< CPosVelPoint<double> >::const_iterator iter=src.PosVelData.begin(); iter!=src.PosVelData.end(); ++iter) 
            PosVelData.push_back(*iter);
        for(std::vector< double >::const_iterator iter=src.MassData.begin(); iter!=src.MassData.end(); ++iter) 
            MassData.push_back(*iter);
    }
    size_t size() const { return PosVelData.size(); };
    void resize(const size_t newSize) { PosVelData.resize(newSize); MassData.resize(newSize,0); };
    CPosVelPoint<NumT>& point(const size_t index) { return PosVelData[index]; };
    const CPosVelPoint<NumT>& point(const size_t index) const { return PosVelData[index]; };
    CPosVelPoint<NumT>& mass(const size_t index) { return MassData[index]; };
    const CPosVelPoint<NumT>& mass(const size_t index) const { return MassData[index]; };
    void push_back(const CPosVelPoint<NumT>& point, const NumT mass) { PosVelData.push_back(point); MassData.push_back(mass); };
private:
    std::vector< CPosVelPoint<NumT> > PosVelData;
    std::vector< NumT > MassData;
};
*/
class CPointMassSetHandler
{
public:
    CPointMassSetHandler(const CPointMassSetFloat* _points) : points(_points), state(PS_IDLE) {};
    void computeVirialRatio(double* T, double* W);
    void halt() { if(state==PS_RUNNING) state=PS_NEEDTOTERMINATE; };
private:
enum PSTATE{
    PS_IDLE,
    PS_RUNNING,
    PS_NEEDTOTERMINATE,
    PS_DONE
};
    const CPointMassSetFloat* points;
    volatile PSTATE state;
};

class CMassRefinementFnc: public CBasicOrbitFilteringFnc
{
public:
    CMassRefinementFnc(const COrbitLibrary* orbitlib, unsigned int numBins, const CBasicOrbitFilteringFnc* _evalFnc);
    virtual double eval(const COrbitDesc* orbit) const;
private:
    const CBasicOrbitFilteringFnc* rankFnc;
    vectord binBoundary;
};

class CEnergyOrbitFilteringFnc: public CBasicOrbitFilteringFnc
{
public:
    virtual double eval(const COrbitDesc* orbit) const;
};

}