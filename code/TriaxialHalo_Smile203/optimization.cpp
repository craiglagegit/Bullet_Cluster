#include "optimization.h"

#include <QFile>
#include <QTextStream>
#include <QProcess>
#include <glpk.h>    // GNU linear programming kit - optimization problem solver
namespace smile{

/*------- LP solver based on GNU linear programming kit -------*/

static int glpk_output_hook(void * tmp, const char *data)
{
    std::string* infostr = static_cast<std::string*>(tmp);
    if(*data == '\r' || *data == '\n')
    {
        my_message(*infostr);
        infostr->clear();
    }
    else
        infostr->append(data);
    return 1;
}

int COptimizationSolverGLPK::callSolver(const std::vector<vectorn> &linearMatrix, const vectorn &rhs, const vectorn &rhsWeight, const vectorn &solWeight, const double maxSolValue, vectorn *sol) const
{
    try{
    unsigned int numVariables = static_cast<unsigned int>(linearMatrix.size());
    unsigned int numConstraints = static_cast<unsigned int>(rhs.size());
    if(rhsWeight.size()!=numConstraints) { my_error("Invalid call to solver"); return -999; }
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, numConstraints);
    for(unsigned int c=0; c<numConstraints; c++)
        glp_set_row_bnds(lp, c+1, GLP_FX, rhs[c], rhs[c]);
    glp_add_cols(lp, numVariables+2*numConstraints);
    if(maxSolValue>0)
        for(unsigned int c=0; c<numVariables; c++)
            glp_set_col_bnds(lp, c+1, GLP_DB, 0.0, maxSolValue);
    else
        for(unsigned int c=0; c<numVariables; c++)
            glp_set_col_bnds(lp, c+1, GLP_LO, 0.0, 0.0);
    for(unsigned int c=0; c<2*numConstraints; c++)  // lambda, mu - slack variables (should be zero after optimization)
    {
        glp_set_col_bnds(lp, c+1+numVariables, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, c+1+numVariables, rhsWeight[c/2]);
    }
    std::vector<int> ic(1);
    std::vector<int> io(1);
    std::vector<double> matrix(1);
    for(unsigned int o=0; o<numVariables; o++)
    {
        if(linearMatrix[o].size()!=numConstraints) { my_error("Invalid call to solver"); return -999; }
        for(unsigned int c=0; c<numConstraints; c++)
        {
            if(linearMatrix[o][c]!=0)
            {
                ic.push_back(c+1);
                io.push_back(o+1);
                matrix.push_back(linearMatrix[o][c]);
            }
        }
        if(solWeight.size() && solWeight[o]!=0)
            glp_set_obj_coef(lp, o+1, solWeight[o]/numVariables);
    }
    for(unsigned int c=0; c<numConstraints; c++)
    {

        ic.push_back(c+1);
        io.push_back(numVariables+2*c+1);
        matrix.push_back(1.0);
        ic.push_back(c+1);
        io.push_back(numVariables+2*c+2);
        matrix.push_back(-1.0);
    }
    glp_load_matrix(lp, static_cast<int>(ic.size())-1, &(ic.front()), &(io.front()), &(matrix.front()));

    my_message("Starting optimization");
    std::string infostr;
    glp_term_hook(glpk_output_hook, &infostr);

    //int status = glp_simplex(lp, NULL);
    int status = glp_interior(lp, NULL);
    
    sol->assign(numVariables, 0);
    for(unsigned int o=0; o<numVariables; o++)
    {
        double v = glp_ipt_col_prim(lp, o+1);
        //double v = glp_get_col_prim(lp, o+1);
        (*sol)[o] = v;
    }
    glp_delete_prob(lp); 
    return -status;  // 0 for success, <0 for failure
    }
    catch(const std::bad_alloc&) {
        return -1;
    }
}

std::string COptimizationSolverGLPK::errorDescription(int result) const
{
    switch(result)
    {
    case 0: return "Success"; break;
    case -1: return "Not enough memory in GLPK solver"; break;
    case GLP_EFAIL: return "Empty problem in GLPK solver"; break;
    case GLP_ENOCVG: return "Bad convergence in GLPK solver"; break;
    case GLP_EITLIM: return "Number of iterations exceeded limit in GLPK solver"; break;
    case GLP_EINSTAB: return "Numerical instability in GLPK solver"; break;
    default: return "Unknown error in GLPK solver";
    }
}

/*------- B P M P D  - based LP/QP solver -------*/
int COptimizationSolverBPMPD::callSolver(const std::vector<vectorn> &linearMatrix, const vectorn &rhs, const vectorn &rhsWeight, const vectorn &solWeight, const double maxSolValue, vectorn *sol) const
{
    unsigned int numVariables = static_cast<unsigned int>(linearMatrix.size());
    unsigned int numVarActive = 0;
    unsigned int numConstraints = static_cast<unsigned int>(rhs.size());
    if(rhsWeight.size()!=numConstraints) { my_error("Invalid call to solver"); return -999; }
    // first write out MPS file
    my_message("Writing file");

    QString dataFileName = QString::fromStdString(tempDir)+"data.mps";
    QFile file(dataFileName);
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    QString Str;
    out << "NAME  testmp" << endl;
    out << "ROWS" << endl << " N  OBJFUN" << endl;
    for(unsigned int c=0; c<numConstraints; c++)
        out << Str.sprintf(" E  R%07d", c+1) << endl;
    if(maxSolValue>0)
    {
        for(unsigned int o=0; o<numVariables; o++)
            out << Str.sprintf(" L  W%07d", o+1) << endl;
    }
    out << "COLUMNS";
    int parity=0;
    for(unsigned int o=0; o<numVariables; o++)
    {
        if(linearMatrix[o].size()!=numConstraints) { my_error("Invalid call to solver"); return -999; }
        bool active=false;
        if(solWeight.size()>0 && solWeight[o]!=0)
        {
            out << endl << Str.sprintf("    C%07d  OBJFUN    %15.13G ", o+1, solWeight[o]/numVariables);
            parity=1;
            active=true;
        }
        if(maxSolValue>0)
        {
            if(parity==0) 
                out << endl << Str.sprintf("    C%07d", o+1);
            out << Str.sprintf("  W%07d             1 ", o+1);
            parity=1-parity;
            active=true;
        }
        for(unsigned int c=0; c<numConstraints; c++)
            if(linearMatrix[o][c]!=0)
            {
                if(parity==0) 
                    out << endl << Str.sprintf("    C%07d", o+1);
                out << Str.sprintf("  R%07d  %15.13G ", c+1, linearMatrix[o][c]);
                parity = 1-parity;
                active=true;
            }
        parity=0;
        if(active) numVarActive++;
    }
    out << endl;
    for(unsigned int c=0; c<numConstraints; c++)
    {
        out << Str.sprintf("    L%07d  OBJFUN    %15.13G   R%07d  %15.13G", c+1, rhsWeight[c], c+1,  1.0) << endl;
        out << Str.sprintf("    M%07d  OBJFUN    %15.13G   R%07d  %15.13G", c+1, rhsWeight[c], c+1, -1.0) << endl;
    }
    out << "RHS";
    parity=0;
    for(unsigned int c=0; c<numConstraints; c++)
    {
        if(parity==0)
            out << endl << "    RHS1    ";
        out << Str.sprintf("  R%07d  %15.13G ", c+1, rhs[c]);
        parity=1-parity;
    }
    if(maxSolValue>0)
    {
        for(unsigned int o=0; o<numVariables; o++)
        {
            if(parity==0)
                out << endl << "    RHS1    ";
            out << Str.sprintf("  W%07d  %12G ", o+1, maxSolValue);
            parity=1-parity;
        }
    }
    if(QP)
    {   // add Quadratic optimization block
        out << endl << "QUADOBJ";
        for(unsigned int o=0; o<numVariables; o++)
            out << endl << Str.sprintf("    C%07d  C%07d  %12G", o+1, o+1, SCHW_CONSTRAINT_QUAD_WEIGHT);
    }
    out << endl << "ENDATA" << endl;
    file.close();

    // then run BPMPD
    my_message("Starting optimization");
    QProcess process;
    process.setWorkingDirectory(QString::fromStdString(tempDir));  // to read .par file correctly
    process.setProcessChannelMode(QProcess::ForwardedChannels);  // reading both stdout and stderr
    process.start(QString::fromStdString(execName), QStringList("data") );
    if(!process.waitForStarted()) return -2;  // report error
    QProcess::ProcessState st=process.state();
    QProcess::ProcessError er=process.error();
    if(st!=QProcess::Running) return -2-er;  // report error
    process.waitForFinished(-1);
    file.remove();
    // finally, read results
    my_message("Reading results");
    QFile fileo(QString::fromStdString(tempDir)+"data.out");
    int status=-4;
    if(fileo.open(QIODevice::ReadOnly))
    {
        sol->assign(numVariables, 0);
        QTextStream in(&fileo);
        QString S;
        unsigned int numColumnsResult=0;
        do{
            S = in.readLine();
            if(S.contains("Limit of the evaluation version is reached")) return -5;  // BPMPD: sadly we won't know the answer...
            bool ok=false;
            QStringList fields=S.trimmed().split(QRegExp("\\s+"));
            unsigned int c=fields[0].mid(1,7).toInt(&ok);
            if(fields[0].mid(0,1)=="C" && ok && c>0 && c<=numVariables)
            {
                double v=fields[1].toDouble(&ok);
                if(ok)
                {
                    (*sol)[c-1] = v;
                    numColumnsResult++;
                }
            }
        }
        while(!in.atEnd());
        status=(numColumnsResult==numVarActive ? 0 : -6);
        fileo.close();
        fileo.remove();
    }
    return status;
}

std::string COptimizationSolverBPMPD::errorDescription(int result) const
{
    switch(result)
    {
    case -1: return "Not enough memory in BPMPD solver"; break;
    case -2: return "BPMPD: Process failed to start"; break;
    case -3: return "BPMPD: Process crashed unexpectedly"; break;
    case -4: return "BPMPD: Result data file not found"; break;
    case -5: return "BPMPD: Problem size exceeds limit of evaluation version"; break;
    case -6: return "BPMPD: Unrecognized output format"; break;
    default: return "BPMPD: Unknown error";
    }
}

}
