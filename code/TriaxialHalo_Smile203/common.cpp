#include "common.h"
#include "potential.h"
#include <cmath>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

namespace smile{

//----------------------------------------------------------------------------//
// Creation of grid with exponentially increasing cells
double findA(double A, void* params)
{
    double dynrange = ((double*)params)[0];
    double nnodes   = ((double*)params)[1];
    return (exp(A*nnodes)-1)/(exp(A)-1) - dynrange;
}
void createNonuniformGrid(vectord &grid, int nnodes, double xmin, double xmax, bool zeroelem)
{   // create grid so that x_k = B*(exp(A*k)-1)
    double A, B, dynrange=xmax/xmin;
    grid.resize(nnodes);
    int indexstart=zeroelem?1:0;
    if(zeroelem) { grid[0]=0; nnodes--; }
    if(nnodes>=dynrange)  // no need for non-uniform grid
    {
        for(int i=0; i<nnodes; i++)
            grid[i+indexstart] = xmin+(xmax-xmin)*i/(nnodes-1);
        return;
    }
    // solve for A:  dynrange = (exp(A*nnodes)-1)/(exp(A)-1)
    double params[2] = {dynrange, nnodes};
    gsl_function F;
    F.function=&findA;
    F.params=&params;
    double Alow=0.0001, Aupp=100./nnodes;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set (s, &F, Alow, Aupp);
    int status=0, iter=0;
    do{
        iter++;
        status = gsl_root_fsolver_iterate (s);
        A    = gsl_root_fsolver_root (s);
        Alow = gsl_root_fsolver_x_lower (s);
        Aupp = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (Alow, Aupp, 0, 1e-9);
    }
    while (status == GSL_CONTINUE && iter < 100);
    gsl_root_fsolver_free (s);
    B = xmin / (exp(A)-1);
    for(int i=0; i<nnodes; i++)
        grid[i+indexstart] = B*(exp(A*(i+1))-1);
    grid[nnodes-1+indexstart]=xmax;
}

static void my_stderr_show_message(const std::string &message)
{
    std::cerr << message <<"\n";
}

show_message_type* my_error_ptr   = &my_stderr_show_message;
show_message_type* my_message_ptr = &my_stderr_show_message;

/** helper function for initializing missing timestep/timeunit in incomplete initial conditions data **/
template<typename NumT> COrbitInitData<NumT> completeInitData(const COrbitInitData<NumT>& InitData)
{
    NumT timeStep=InitData.timeStep;
    NumT timeUnit=InitData.timeUnit;
    bool changed=false;
    if(timeUnit==0 && InitData.potential!=NULL)
    {
        changed=true;
        timeUnit = static_cast<NumT>(InitData.potential->longaxisperiod( InitData.potential->totalEnergy(InitData.initCond)));
    }
    if(timeStep==0)
    {
        changed=true;
        timeStep=timeUnit/10;  //< somewhat arbitrary number, but should work reasonably
    }
    if(changed)
        return COrbitInitData<NumT>(InitData.potential, timeStep, timeUnit, InitData.initCond, InitData.calcLyapunov);
    else
        return InitData;
}
template COrbitInitData<double> completeInitData(const COrbitInitData<double>& InitData);
template COrbitInitData<float > completeInitData(const COrbitInitData<float >& InitData);

}  // namespace