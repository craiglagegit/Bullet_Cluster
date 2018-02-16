#pragma once
#include "common.h"

namespace smile{

/* solves the optimization problem defined as  
   M sol = rhs, where M is N_v * N_c matrix, sol is N_v vector of variables to be found, rhs is N_c vector of constraints to be satisfied;
   tries to achieve exact solution but allows for deviation from it, which is penalized in the cost function;
   denote deviation dev[c] = | \Sum_o M[v][c] sol[v] - rhs[c] |
   with cost function = (quad.coef.)*\Sum_v sol[v]^2 + \Sum_v solWeight[v]*sol[v] + \Sum_c rhsWeight[c]*dev[c]
   if solved by linear programming, the first terms is omitted.
   Solution is returned in vector sol, which also may contain some initial guess for the solution as an input parameter (might be used in Lucy iterations, for example).
*/
class COptimizationSolverBPMPD: public CBasicOptimizationSolver
{
private:
    std::string execName, tempDir;
    bool QP;
public:
    COptimizationSolverBPMPD(const std::string &_execName, const std::string &_tempDir, const bool _QP): 
        execName(_execName), tempDir(_tempDir), QP(_QP) {};
    virtual int callSolver(
        const std::vector< vectorn > &linearMatrix, // matrix M of linear equations for the optimization problem  "M w = rhs"
        const vectorn &rhs,         // constraints to be satisfied
        const vectorn &rhsWeight,   // constraint importance
        const vectorn &solWeight,   // penalties for elements in sol vector
        const double maxSolValue,   // if nonzero, specifies upper limit on sol[v] in addition to lower limit of zero
        vectorn *sol) const;        // returns solution in this vector
    virtual std::string errorDescription(int result) const;
};

class COptimizationSolverGLPK: public CBasicOptimizationSolver
{
public:
    COptimizationSolverGLPK() {};
    virtual int callSolver(
        const std::vector< vectorn > &linearMatrix,
        const vectorn &rhs,
        const vectorn &rhsWeight,
        const vectorn &solWeight,
        const double maxSolValue,
        vectorn *sol) const;
    virtual std::string errorDescription(int result) const;
};

}