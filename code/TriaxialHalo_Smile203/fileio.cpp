#include "fileio.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cassert>
#include "potential.h"
#include "orbit.h"
#include "orbitlib.h"
#include "schwarzschild.h"
#include "stringconv.h"

namespace smile{

void splitString(const std::string& src, const std::string& delim, std::vector<std::string> *result)
{
    result->clear();
    std::string str(src);
    std::string::size_type indx=str.find_first_not_of(delim);
    if(indx==std::string::npos) 
    {
        result->push_back("");   // ensure that result contains at least one element
        return;
    }
    if(indx>0)  // remove unnecessary delimiters at the beginning
        str=str.erase(0, indx);
    while(!str.empty())
    {
        indx=str.find_first_of(delim);
        if(indx==std::string::npos) indx=str.size();
        result->push_back(str.substr(0, indx));
        str=str.erase(0, indx);
        indx=str.find_first_not_of(delim);
        if(indx==std::string::npos) return;
        str=str.erase(0, indx);
    }
}

void CFileInputText::open(const std::string &_baseFileName)
{
    baseFileName=_baseFileName;
    if(isopen) 
    {
        strmText.close();
        strmSchwData.close();
        strmTrajSample.close();
        isopen=false;
    }
    status=true;
    strmText.open(baseFileName.c_str(), std::ios::in);
    if(!strmText) { status=false; return; }
    isopen=status;
}

CFileInputText::~CFileInputText()
{
    if(strmText.is_open()) strmText.close();
    if(strmSchwData.is_open()) strmSchwData.close();
    if(strmTrajSample.is_open()) strmTrajSample.close();
}

COrbit* CFileInputText::readOrbit(const CPotential* _potential)
{
    if(!isopen || !strmText) return NULL;
    CPosVelDataDouble trj;
    double time=0;
    int countlines=0;
    std::string buffer;
    std::vector<std::string> fields;
    while(std::getline(strmText, buffer) && !strmText.eof())
    {
        if(buffer.size()>0 && ((buffer[0]>='0' && buffer[0]<='9') || buffer[0]=='-' || buffer[0]=='+'))
        {
            splitString(buffer, " \t", &fields);
            if(fields.size()<7)
                return NULL;
            time=StringVariant(fields[0]).toDouble();
            trj.push_back(CPosVelPoint<double>(   // don't check for errors in conversion, assume value=0 instead
                StringVariant(fields[1]).toDouble(), StringVariant(fields[2]).toDouble(), StringVariant(fields[3]).toDouble(),
                StringVariant(fields[4]).toDouble(), StringVariant(fields[5]).toDouble(), StringVariant(fields[6]).toDouble() ));
            countlines++;
        }
    }
    if(countlines>1)
        return new COrbit(_potential, time/(countlines-1) /*timeStep*/, 0 /*timeUnit=autodetect*/, trj);
    else return NULL;
}

COrbitDesc* CFileInputText::readOrbitDesc(const CPotential* _potential, bool withModelData)
{
    if(!isopen || !strmText) return NULL;
    std::string buffer;
    std::vector<std::string> fields;
    if(std::getline(strmText, buffer) && buffer.size()>0 && ((buffer[0]>='0' && buffer[0]<='9') || buffer[0]=='-' || buffer[0]=='+'))
    {
        try{
        std::vector<std::string> fields;
        splitString(buffer, " \t", &fields);
        size_t nf=fields.size();
        if(nf<3)
            return NULL;
        vectorInformation info;
        COrbitInitData<float> InitData(_potential, 
            nf>8 ? StringVariant(fields[8]).toFloat() : 0, 
            nf>8 ? StringVariant(fields[7]).toFloat() : 0,
            CPosVelPoint<float>( StringVariant(fields[0]).toFloat(), StringVariant(fields[1]).toFloat(), StringVariant(fields[2]).toFloat(), 
            nf>5 ? StringVariant(fields[3]).toFloat() : 0, nf>5 ? StringVariant(fields[4]).toFloat() : 0, nf>5 ? StringVariant(fields[5]).toFloat() : 0), 
            nf>17 ? StringVariant(fields[17]).toFloat()>=0 /* calcLyapunov: if Lyapunov exp.!=-1, then it was at least calculated */ : 0);
        float weight= nf>6 ? StringVariant(fields[6]).toFloat() : 0;
        float intTime=nf>9 ? StringVariant(fields[9]).toFloat() : 0;
        float maxTime=nf>10? StringVariant(fields[10]).toFloat() : 0;
        float Einit  =nf>11? StringVariant(fields[11]).toFloat() : 0;
        float Ediff  =nf>12? StringVariant(fields[12]).toFloat() : 0;
        if(InitData.timeUnit==0)
        {
            if(Einit==0) Einit=static_cast<float>(_potential->totalEnergy(InitData.initCond));
            InitData.timeUnit =static_cast<float>(_potential->longaxisperiod(Einit));
        }
        float lf[N_DIM]={0,0,0};
        if(nf>15) {
            lf[0]=StringVariant(fields[13]).toFloat();
            lf[1]=StringVariant(fields[14]).toFloat();
            lf[2]=StringVariant(fields[15]).toFloat();
        }
        float lfdiff =nf>16? StringVariant(fields[16]).toFloat() : 0;
        float lambda =nf>17? StringVariant(fields[17]).toFloat() : -1;
        std::string Description=nf>18 ? fields[18] : "empty";
        std::string::size_type ind;
        while((ind = Description.find("_")) != std::string::npos)
            Description=Description.replace(ind, 1, " ");
        float inert[N_DIM]={0,0,0};
        if(nf>21) {
            inert[0]=StringVariant(fields[19]).toFloat();
            inert[1]=StringVariant(fields[20]).toFloat();
            inert[2]=StringVariant(fields[21]).toFloat();
        }
        float Lavg[N_DIM]={0,0,0};
        float Lvar[N_DIM]={0,0,0};
        if(nf>27) {
            Lavg[0]=StringVariant(fields[22]).toFloat();
            Lvar[0]=StringVariant(fields[23]).toFloat();
            Lavg[1]=StringVariant(fields[24]).toFloat();
            Lvar[1]=StringVariant(fields[25]).toFloat();
            Lavg[2]=StringVariant(fields[26]).toFloat();
            Lvar[2]=StringVariant(fields[27]).toFloat();
        }
        info.push_back(new COrbitInformation<float>(Description, Einit, Ediff, lf, lfdiff, lambda, inert, Lavg, Lvar));
        if(nf>32) {
            float fitIntercept=StringVariant(fields[28]).toFloat();
            float fitSlope=StringVariant(fields[29]).toFloat();
            float fitScatter=StringVariant(fields[30]).toFloat();
            float fitSignificance=StringVariant(fields[31]).toFloat();
            float Lcirc2=StringVariant(fields[32]).toFloat();
            info.push_back(new CPericenterInformation<float>(fitIntercept, fitSlope, fitScatter, fitSignificance, Lcirc2));
        }
        if(withModelData && strmSchwData.is_open() && strmSchwData.good())
        {
            unsigned int datasize=0;
            strmSchwData.read( reinterpret_cast<char*>(&datasize), sizeof(unsigned int) );
            vectorn densityData(datasize);
            if(datasize>0) strmSchwData.read( reinterpret_cast<char*>(&(densityData.front()) ), datasize*sizeof(smnumtype));
            strmSchwData.read( reinterpret_cast<char*>(&datasize), sizeof(unsigned int) );
            vectorn shellTime(datasize);
            if(datasize>0) strmSchwData.read( reinterpret_cast<char*>(&(shellTime.front()) ), datasize*sizeof(smnumtype));
            strmSchwData.read( reinterpret_cast<char*>(&datasize), sizeof(unsigned int) );
            vectorn shellVr(datasize);
            if(datasize>0) strmSchwData.read( reinterpret_cast<char*>(&(shellVr.front()) ), datasize*sizeof(smnumtype));
            strmSchwData.read( reinterpret_cast<char*>(&datasize), sizeof(unsigned int) );
            vectorn shellVt(datasize);
            if(datasize>0) strmSchwData.read( reinterpret_cast<char*>(&(shellVt.front()) ), datasize*sizeof(smnumtype));
            info.push_back(new CSchwInformation(densityData, shellTime, shellVr, shellVt));
            if(!strmSchwData) status=false;
        }
        if(withModelData && strmTrajSample.is_open() && strmTrajSample.good())
        {
            unsigned int nump=0;
            strmTrajSample.read( reinterpret_cast<char*>(&nump), sizeof(unsigned int) );
            if(nump>0)
            {
                std::vector< CPosVelPoint<float> > data(nump);
                if(nump>0) strmTrajSample.read( reinterpret_cast<char*>( &(data.front()) ), nump * sizeof(CPosVelPoint<float>));
                info.push_back( new CTrajSampleInformation<float>(CPosVelDataFloat(data)) );
            }
            if(!strmTrajSample) status=false;
        }
        if(!strmText) {
            status=false; 
            for(size_t ind=0; ind<info.size(); ind++) delete info[ind];
            return NULL;
        }
        return new COrbitDesc(InitData, intTime, maxTime, weight, info);
        }
        catch(const std::bad_alloc&) {
            return NULL;
        }
    }
    else return NULL;
}

COrbitLibrary* CFileInputText::readOrbitLibrary(const CPotential* _potential, bool withModelData)
{
    if(!status || !isopen) return NULL;
    COrbitLibrary* orbitlib=NULL;
    if(withModelData)
    {
        strmSchwData.open((baseFileName+".ocw").c_str(), std::ios::in | std::ios::binary);
        //if(!strmSchwData) status=false;
        strmTrajSample.open((baseFileName+".sam").c_str(), std::ios::in | std::ios::binary);
        //if(!strmSchwData) status=false;
    }
    while(strmText && !strmText.eof())
    {
        COrbitDesc* od=readOrbitDesc(_potential, withModelData);
        if(od!=NULL)
        {
            if(orbitlib==NULL) orbitlib=new COrbitLibrary(_potential, NULL);
            orbitlib->addOrbitDesc(od);
        }
    }
    return orbitlib;
}

const CPotential* CFileInputText::readPotential(CConfigPotential *configPotential)
{
    if(!status || !isopen) return NULL;
    if( configPotential==NULL || configPotential->N_dim!=N_DIM || 
        ! (configPotential->PotentialType == CPotential::PT_NB || 
        configPotential->PotentialType == CPotential::PT_SPLINE ||
        configPotential->PotentialType == CPotential::PT_BSE ||
        configPotential->PotentialType == CPotential::PT_SCALEFREESH ||
        configPotential->PotentialType == CPotential::PT_UNKNOWN ) ) 
        return NULL;   // only may load these types of potential
    CPointMassSetDouble points;
    std::string buffer;
    std::vector<std::string> fields;
    bool noMasses=false;
    while(std::getline(strmText, buffer) && !strmText.eof())
    {
        splitString(buffer, "# \t", &fields);
        if(fields[0] == "BSEcoefs" || fields[0] == "SFcoefs" || fields[0] == "SHEcoefs")
        {   // encountered header of potential coefs file, try to load coefficients
            CDensity::POTENTIALTYPE ptype = 
                fields[0] == "BSEcoefs" ? CDensity::PT_BSE :
                fields[0] == "SFcoefs" ? CDensity::PT_SCALEFREESH :
                fields[0] == "SHEcoefs" ? CDensity::PT_SPLINE : CDensity::PT_UNKNOWN;  // the last option shouldn't happen:)
            bool valid=true;
            if(!std::getline(strmText, buffer)) valid=false;
            splitString(buffer, "# \t", &fields);
            unsigned int ncoefsRadial=StringVariant(fields[0]).toInt();
            if(!std::getline(strmText, buffer)) valid=false;
            splitString(buffer, "# \t", &fields);
            unsigned int ncoefsAngular=StringVariant(fields[0]).toInt();
            if(!std::getline(strmText, buffer)) valid=false;
            splitString(buffer, "# \t", &fields);
            double param=StringVariant(fields[0]).toDouble();   // meaning of this parameter depends on potential type
            if(/*ncoefsRadial==0 || ncoefsAngular==0 || */  // zero values are possible, means just a single term in expansion
                (ptype == CDensity::PT_BSE && param<0.5) || 
                (ptype == CDensity::PT_SCALEFREESH && (param<0 || param>2 || ncoefsRadial!=0)) ||
                (ptype == CDensity::PT_SPLINE && ncoefsRadial<4) ) valid=false;
            std::vector< vectord > coefs;
            std::vector< double > radii;
            if(ncoefsRadial>MAX_NCOEFS_RADIAL) ncoefsRadial=MAX_NCOEFS_RADIAL;
            if(ncoefsAngular>MAX_NCOEFS_ANGULAR) ncoefsAngular=MAX_NCOEFS_ANGULAR;
            while(valid && std::getline(strmText, buffer))  // time, ignored
            {
                std::getline(strmText, buffer);  // comments, ignored
                radii.clear();
                coefs.clear();
                for(unsigned int n=0; valid && n<=ncoefsRadial; n++)
                {
                    std::getline(strmText, buffer);
                    splitString(buffer, "# \t", &fields);
                    radii.push_back(StringVariant(fields[0]).toDouble());
                    if((ptype == CDensity::PT_BSE && radii.back()!=n) || (ptype == CDensity::PT_SPLINE && n>0 && radii.back()<=radii[n-1])) 
                        valid=false;  // for BSE this field is basis function index, for spline the radii should be in increasing order
                    coefs.push_back( vectord() );
                    for(int l=0; l<=static_cast<int>(ncoefsAngular); l++)
                        for(int m=-l; m<=l; m++)
                        {
                            unsigned int fi=1+l*(l+1)+m;
                            coefs.back().push_back( fi<fields.size() ? StringVariant(fields[fi]).toDouble() : 0);
                        }
                }
            }
            if(!valid || (configPotential->PotentialType != CPotential::PT_UNKNOWN && ptype!=configPotential->PotentialType))
            {   // load coefs only if potential type is the same as requested, or if request is not specific
                my_message("Error loading potential "+PotentialNames[configPotential->PotentialType]+" coefs from file "+baseFileName);
                return NULL;
            }
            CPotential* pot=NULL;
            switch(ptype)
            {
            case CDensity::PT_BSE: 
                configPotential->Alpha=param;
                pot=new CPotentialBSE(N_DIM, configPotential->Mbh, /*Alpha*/param, coefs); 
                break;
            case CDensity::PT_SCALEFREESH: 
                configPotential->Gamma=param;
                pot=new CPotentialScaleFreeSH(N_DIM, configPotential->Mbh, /*Gamma*/param, coefs[0]); 
                break;
            case CDensity::PT_SPLINE: 
                pot=new CPotentialSpline(N_DIM, configPotential->Mbh, radii, coefs); 
                break;
            default:  return NULL;  // shouldn't happen
            }
            configPotential->Ncoefs_radial=ncoefsRadial;    
            configPotential->Ncoefs_angular=ncoefsAngular;
            configPotential->SymmetryType=pot->symmetry();
            return pot;
        }
        else  // load potential from a set of point masses
        if(fields.size()>=3 && ((fields[0][0]>='0' && fields[0][0]<='9') || fields[0][0]=='-' || fields[0][0]=='+'))
        {
            points.push_back(std::pair< CPosVelPoint<double>, double>(CPosVelPoint<double>(
                StringVariant(fields[0]).toDouble(), 
                StringVariant(fields[1]).toDouble(), 
                StringVariant(fields[2]).toDouble(), 
                0, 0, 0), 
                fields.size()>6 ? StringVariant(fields[6]).toDouble():(noMasses=true,1.0)));
        }
    }
    if(points.empty())
    {
        my_message("Error loading potential "+PotentialNames[configPotential->PotentialType]+" from file "+baseFileName);
        if(configPotential->PotentialType==CDensity::PT_NB) 
            return new CPotentialNB(N_DIM, configPotential->Mbh); 
        else if(configPotential->PotentialType==CDensity::PT_BSE) 
            return new CPotentialBSE(N_DIM, configPotential->Mbh, configPotential->Alpha, 1, 1);
        else if(configPotential->PotentialType==CDensity::PT_SPLINE) 
            return new CPotentialSpline(N_DIM, configPotential->Mbh, 1, 1);
        else return NULL;
    }
    if(noMasses)
    {
        double mass=1.0/points.size();
        for(size_t i=0; i<points.size(); i++)
            points[i].second=mass;
    }
    const CPotential* pot=NULL;
    if(configPotential->PotentialType==CDensity::PT_NB)
        pot = new CPotentialNB(N_DIM, configPotential->Mbh, configPotential->nbpotEps, configPotential->nbpotTol, points);
    else if(configPotential->PotentialType==CDensity::PT_SPLINE)
        pot = new CPotentialSpline(N_DIM, configPotential->Mbh, configPotential->Ncoefs_radial, configPotential->Ncoefs_angular, points, configPotential->SymmetryType);
    else if(configPotential->PotentialType==CDensity::PT_BSE)
        pot = new CPotentialBSE(N_DIM, configPotential->Mbh, configPotential->Alpha, configPotential->Ncoefs_radial, configPotential->Ncoefs_angular, points, configPotential->SymmetryType);
    else return NULL;  // unsupported type (unknown or scale-free)
    if(configPotential->PotentialType!=CDensity::PT_NB)
    {   // store coefficients in a text file, later may load this file instead for faster initialization
        CFileOutputText coefout(baseFileName+".coef");
        coefout.writePotential(pot);
    }
    return pot;
}


/// Output to text file (and possibly additional binary files with bulky orbit data ///

void CFileOutputText::open(const std::string &_baseFileName)
{
    baseFileName=_baseFileName;
    if(isopen) 
    {
        strmText.close();
        strmSchwData.close();
        strmTrajSample.close();
        isopen=false;
    }
    started=false;
    status=true;
    strmText.open(baseFileName.c_str(), std::ios::out);
    if(!strmText) { status=false; return; }
    isopen=status;
}

CFileOutputText::~CFileOutputText()
{
    if(strmText.is_open()) strmText.close();
    if(strmSchwData.is_open()) strmSchwData.close();
    if(strmTrajSample.is_open()) strmTrajSample.close();
}

bool CFileOutputText::writeOrbit(const COrbit* orbit)
{
    if(!status || !isopen) return false;
    const COrbitRuntimeTrajectory* trj = NULL;
    const COrbitInitData<double> InitData=orbit->getInitData();
    for(size_t rf=0; rf<orbit->getInfoNum(); rf++)
        if(orbit->getRuntimeFnc(rf)->FncType() == CBasicOrbitRuntimeFnc::FT_TRAJ_ANALYSIS)
            trj = static_cast<const COrbitRuntimeTrajectory*>(orbit->getRuntimeFnc(rf));
    if(trj==NULL || trj->getTrajSize()==0) 
        return false;
    if(!started)
    {
        strmText << "t\tx\ty\tz\t\vx\tvy\tvz" << std::endl;
        started=true;
    }
    for(size_t indx=0; indx<trj->getTrajSize(); indx++)
    {
        CPosVelPoint<double> pt(trj->getTraj(indx));
        strmText << (indx*InitData.timeStep) << "\t" << 
            pt.Pos[0] << "\t" << pt.Pos[1] << "\t" << pt.Pos[2] << "\t" << 
            pt.Vel[0] << "\t" << pt.Vel[1] << "\t" << pt.Vel[2] << std::endl;
        if(!strmText) status=false;
    }
    return status;
}

bool CFileOutputText::writeOrbitDesc(const COrbitDesc* value, bool withModelData)
{
    if(!status || !isopen) return false;
    const COrbitInformation<float>* infoOrbit=NULL;
    const CPericenterInformation<float>* infoPeri=NULL;
    const CSchwInformation* infoSchw=NULL;
    const CTrajSampleInformation<float>* infoTrajSample=NULL;
    for(size_t inf=0; inf<value->getInfoNum(); inf++)
    {
        const CBasicInformation* info=value->getInfo(inf);
        if(info->infoType()==CBasicInformation::IT_TRAJ_ANALYSIS) 
            infoOrbit=static_cast<const COrbitInformation<float>* >(info);
        if(info->infoType()==CBasicInformation::IT_TRAJ_SAMPLE) 
            infoTrajSample=static_cast<const CTrajSampleInformation<float>* >(info);
        if(info->infoType()==CBasicInformation::IT_SCHW)
            infoSchw=static_cast<const CSchwInformation* >(info);
        if(info->infoType()==CBasicInformation::IT_PERICENTER) 
            infoPeri=static_cast<const CPericenterInformation<float>* >(info);
    }
    if(!started)
    {
        // write header line
        strmText << "x\ty\tz\tvx\tvy\tvz\tweight\ttimeunit\ttimestep\tinttime\tmaxtime";
        if(infoOrbit!=NULL) strmText << "\tEtot\tEdif\tlfx\tlfy\tlfz\tlfdiff\tlambda\tdescription\tinertx\tinerty\tinertz\tLXavg\tLXvar\tLYavg\tLYvar\tLZavg\tLZvar";
        if(infoPeri!=NULL) strmText << "\tL2min\tL2slope\tscatter\tsignificance\tL2circ";
        strmText << std::endl;
        started=true;
        if(withModelData)   // try to open additional binary streams
        {
            if(infoSchw!=NULL)   // open only if 1st OrbitDesc contains schwarzschild data, otherwise don't write it for any subsequent orbits
                strmSchwData.open((baseFileName+".ocw").c_str(), std::ios::out | std::ios::binary);
            if(infoTrajSample!=NULL)   // open only if 1st OrbitDesc contains trajectory sample data, otherwise don't write it for any subsequent orbits
                strmTrajSample.open((baseFileName+".sam").c_str(), std::ios::out | std::ios::binary);
            if(!strmSchwData || !strmTrajSample) status=false;
        }
    }
    COrbitInitData<float> InitData=value->getInitData();
    float intTime=value->getIntTime();
    float maxTime=value->getMaxTime();
    if(maxTime==intTime) maxTime=0;  // default
    strmText <<  std::setprecision(8) << 
        InitData.initCond.Pos[0] << "\t" << InitData.initCond.Pos[1] << "\t" << InitData.initCond.Pos[2] << "\t" << 
        InitData.initCond.Vel[0] << "\t" << InitData.initCond.Vel[1] << "\t" << InitData.initCond.Vel[2] << "\t" << 
        std::setprecision(6) << value->getWeight() << "\t" <<
        InitData.timeUnit << "\t" << InitData.timeStep << "\t" << intTime << "\t" << maxTime;
    if(infoOrbit!=NULL)
    {
        std::string Description=infoOrbit->getDescription();
        std::string::size_type ind;
        while((ind = Description.find(" ")) != std::string::npos)
            Description=Description.replace(ind, 1, "_");
        strmText << "\t" << infoOrbit->getEinit() << "\t" << infoOrbit->getEdiff() << "\t" << 
            infoOrbit->getlf(0) << "\t" << infoOrbit->getlf(1) << "\t" << infoOrbit->getlf(2) << "\t" << 
            infoOrbit->getlfccdiff() << "\t" << (InitData.calcLyapunov ? infoOrbit->getlambda() : -1) << "\t" << Description << "\t" << 
            infoOrbit->getinertia(0) << "\t" << infoOrbit->getinertia(1) << "\t" << infoOrbit->getinertia(2) << "\t" <<
            infoOrbit->getLavg(0) << "\t" << infoOrbit->getLvar(0) << "\t" <<
            infoOrbit->getLavg(1) << "\t" << infoOrbit->getLvar(1) << "\t" <<
            infoOrbit->getLavg(2) << "\t" << infoOrbit->getLvar(2);
    }
    if(infoPeri!=NULL)
    {
        float fitIntercept=0, fitSlope=0, fitScatter=0, fitSignificance=0, Lcirc2=0;
        infoPeri->getParams(&fitIntercept, &fitSlope, &fitScatter, &fitSignificance, &Lcirc2);
        strmText << "\t" << fitIntercept << "\t" << fitSlope << "\t" << fitScatter << "\t" << fitSignificance << "\t" << Lcirc2;
    }
    strmText << std::endl;
    if(!strmText) status=false;
    if(withModelData && strmSchwData.is_open())
    {
        unsigned int datasize=0;
        if(infoSchw!=NULL)
        {
            datasize = static_cast<unsigned int>(infoSchw->getDensityData().size());
            strmSchwData.write( reinterpret_cast<const char*>(&datasize), sizeof(unsigned int) );
            if(datasize>0) strmSchwData.write( reinterpret_cast<const char*>(&(infoSchw->getDensityData().front()) ), datasize*sizeof(smnumtype));
            datasize = static_cast<unsigned int>(infoSchw->getShellTime().size());
            strmSchwData.write( reinterpret_cast<const char*>(&datasize), sizeof(unsigned int) );
            if(datasize>0) strmSchwData.write( reinterpret_cast<const char*>(&(infoSchw->getShellTime().front()) ), datasize*sizeof(smnumtype));
            datasize = static_cast<unsigned int>(infoSchw->getShellVr().size());
            strmSchwData.write( reinterpret_cast<const char*>(&datasize), sizeof(unsigned int) );
            if(datasize>0) strmSchwData.write( reinterpret_cast<const char*>(&(infoSchw->getShellVr().front()) ), datasize*sizeof(smnumtype));
            datasize = static_cast<unsigned int>(infoSchw->getShellVt().size());
            strmSchwData.write( reinterpret_cast<const char*>(&datasize), sizeof(unsigned int) );
            if(datasize>0) strmSchwData.write( reinterpret_cast<const char*>(&(infoSchw->getShellVt().front()) ), datasize*sizeof(smnumtype));
        }
        else
            strmSchwData.write( reinterpret_cast<const char*>(&datasize), sizeof(unsigned int) );
        if(!strmSchwData) status=false;
    }
    if(withModelData && strmTrajSample.is_open())
    {
        unsigned int nump=0;
        if(infoTrajSample!=NULL && infoTrajSample->getTraj().size()>0)
        {
            nump=static_cast<unsigned int>(infoTrajSample->getTraj().size());
            std::vector< CPosVelPoint<float> > data;
            data.reserve(nump);
            for(unsigned int i=0; i<nump; i++)
                data.push_back(infoTrajSample->getTraj()[i]);
            strmTrajSample.write( reinterpret_cast<const char*>(&nump), sizeof(unsigned int) );
            strmTrajSample.write( reinterpret_cast<const char*>( &(data.front()) ), nump * sizeof(CPosVelPoint<float>));
        }
        else
            strmTrajSample.write( reinterpret_cast<const char*>(&nump), sizeof(unsigned int) );
        if(!strmTrajSample) status=false;
    }
    return status;
}

bool CFileOutputText::writeOrbitLibrary(const COrbitLibrary* value, bool withModelData)
{
    if(!status || !isopen) return false;
    for(unsigned int o=0; status && o<value->size(); o++)
        writeOrbitDesc(value->getOrbitDesc(o), withModelData);
    return status;
}

bool CFileOutputText::writePotential(const CPotential *potential)
{
    if(!status || !isopen) return false;
    if(! (potential->PotentialType()==CDensity::PT_SCALEFREESH || potential->PotentialType()==CDensity::PT_BSE || potential->PotentialType()==CDensity::PT_SPLINE))
        return false;  // export potential coefs for BSE/Spline/Scale-FreeSH potentials only
    vectord indices;
    std::vector< vectord > coefs;
    int ncoefsAngular=0;
    switch(potential->PotentialType())
    {
    case CDensity::PT_SCALEFREESH: {
        const CPotentialScaleFreeSH* potSH=static_cast<const CPotentialScaleFreeSH*>(potential);
        indices.assign(1, 0);
        coefs.resize(1);
        potSH->getCoefs(&(coefs[0]));
        ncoefsAngular=potSH->getNcoefs_angular();
        strmText << "SFcoefs\t#header\n" << 0 << "\t#n_radial\n" << ncoefsAngular << "\t#n_angular\n" << 
            potential->getGamma() << "\t#gamma\n0\t#time\n";
        strmText << "#index";
        break; }
    case CDensity::PT_BSE: {
        const CPotentialBSE* potBSE=static_cast<const CPotentialBSE*>(potential);
        indices.resize(potBSE->getNcoefs_radial()+1);
        for(size_t i=0; i<indices.size(); i++) indices[i]=i*1.0;
        potBSE->getCoefs(&coefs);
        assert(coefs.size()==indices.size());
        ncoefsAngular=potBSE->getNcoefs_angular();
        strmText << "BSEcoefs\t#header\n" << potBSE->getNcoefs_radial() << "\t#n_radial\n" << ncoefsAngular << "\t#n_angular\n" << 
            potBSE->getAlpha() <<"\t#alpha\n0\t#time\n";
        strmText << "#index";
        break; }
    case CDensity::PT_SPLINE: {
        const CPotentialSpline* potSpline=static_cast<const CPotentialSpline*>(potential);
        potSpline->getCoefs(&indices, &coefs);
        assert(coefs.size()==indices.size());
        assert(indices[0]==0);  // leftmost radius is 0
        coefs[0].resize(1);     // retain only l=0 term for r=0, the rest is supposed to be zero
        ncoefsAngular=potSpline->getNcoefs_angular();
        strmText << "SHEcoefs\t#header\n" << potSpline->getNcoefs_radial() << "\t#n_radial\n" << ncoefsAngular << "\t#n_angular\n" << 
            0 <<"\t#unused\n0\t#time\n";
        strmText << "#radius";
        break; }
    default:
        assert("wrong potential type in writePotential{coefs}");
        return false;   // shouldn't occur, we checked potential type in the beginning
    }
    for(int l=0; l<=ncoefsAngular; l++)
        for(int m=-l; m<=l; m++)
            strmText << "\tl="<<l<<",m="<<m;  // header line
    strmText << "\n";
    for(size_t n=0; n<indices.size(); n++)
    {
        strmText << indices[n];
        strmText << "\t" << std::setprecision(14) << coefs[n][0] << std::setprecision(7);   // leading coeft should be high-accuracy at least for spline potential
        for(size_t i=1; i<coefs[n].size(); i++)
            strmText << "\t" << coefs[n][i];
        strmText << "\n";
    }
    if(!strmText) status=false;
    return status;
}

bool CFileOutputText::writePointMassSet(const CPointMassSetFloat* points)
{
    if(!status || !isopen || points==NULL) return false;
    if(!started)
    {
        strmText << "x\ty\tz\tvx\tvy\tvz\tm" << std::endl;
        started=true;
    }
    for(size_t indx=0; indx<points->size(); indx++)
    {
        const CPosVelPoint<float>& pt=points->at(indx).first;
        strmText << 
            pt.Pos[0] << "\t" << pt.Pos[1] << "\t" << pt.Pos[2] << "\t" << 
            pt.Vel[0] << "\t" << pt.Vel[1] << "\t" << pt.Vel[2] << "\t" <<
            points->at(indx).second << std::endl;
        if(!strmText) status=false;
    }
    return status;
}


void CFileOutputNemo::open(const std::string &fileName, const std::string &header)
{
    if(SnapshotWriter!=NULL) 
    {
        delete SnapshotWriter;  // file closed automatically
    }
    SnapshotWriter=new CNemoSnapshotWriter(fileName, false);
    if(SnapshotWriter->ok())
        SnapshotWriter->WriteHistory(header);
}

/// helper class that writes NEMO-compatible snapshot file
CNemoSnapshotWriter::CNemoSnapshotWriter(const std::string &filename, bool append)
{
    snap.open(filename.c_str(), std::ios::binary | (append?std::ios_base::app : std::ios_base::trunc));
    level=0;
}

CNemoSnapshotWriter::~CNemoSnapshotWriter()
{
   snap.close();
}

void CNemoSnapshotWriter::PutInt(const std::string &name,int var)
{
   snap.put(-110);
   snap.put(9);
   snap.put('i');
   PutZString(name);
   snap.write((const char*)&var,4);
}

void CNemoSnapshotWriter::PutDouble(const std::string &name,double d)
{
   snap.put(-110);
   snap.put(9);
   snap.put('d');
   PutZString(name);
   snap.write((const char*)&d,8);
}

void CNemoSnapshotWriter::PutDoubleArray(const std::string &name,int ndim,const int dim[],const double*ar)
{
   snap.put(-110);
   snap.put(11);
   snap.put('d');
   PutZString(name);
   snap.write((const char*)dim,ndim*4);
   int buf=0;
   snap.write((const char*)&buf,4);
   buf=1;
   for (int i=0;i<ndim;i++)
     buf*=dim[i];
   snap.write((const char*)ar,sizeof(double)*buf);
}

void CNemoSnapshotWriter::PutFloatArray(const std::string &name,int ndim,const int dim[],const float*ar)
{
   snap.put(-110);
   snap.put(11);
   snap.put('f');
   PutZString(name);
   snap.write((const char*)dim,ndim*4);
   int buf=0;
   snap.write((const char*)&buf,4);
   buf=1;
   for (int i=0;i<ndim;i++)
     buf*=dim[i];
   snap.write((const char*)ar,sizeof(float)*buf);
}

void CNemoSnapshotWriter::PutIntArray(const std::string &name,int ndim,const int dim[],const int* ar)
{
   snap.put(-110);
   snap.put(11);
   snap.put('i');
   PutZString(name);
   snap.write((const char*)dim,ndim*4);
   int buf=0;
   snap.write((const char*)&buf,4);
   buf=1;
   for (int i=0;i<ndim;i++)
     buf*=dim[i];
   snap.write((const char*)ar,sizeof(int)*buf);
}

void CNemoSnapshotWriter::PutCharArray(const std::string &name,int ndim,const int dim[],const char* ar)
{
   snap.put(-110);
   snap.put(11);
   snap.put('c');
   PutZString(name);
   snap.write((const char*)dim,ndim*4);
   int buf=0;
   snap.write((const char*)&buf,4);
   buf=1;
   for (int i=0;i<ndim;i++)
     buf*=dim[i];
   snap.write((const char*)ar,buf);
}

void CNemoSnapshotWriter::StartLevel(const std::string &name)
{
   level++;
   snap.put(-110);
   snap.put(9);
   snap.put('(');
   PutZString(name);
}

void CNemoSnapshotWriter::EndLevel()
{
   level--;
   if (level<0)
     { // should raise error
     }
   snap.put(-110);
   snap.put(9);
   snap.put(')');
   snap.put(0);
}

void CNemoSnapshotWriter::PutString(const std::string &name,const std::string &str)
{
   int len=static_cast<int>(str.size());
   PutCharArray(name,1,&len,str.c_str());
}

void CNemoSnapshotWriter::PutZString(const std::string &name)
{
   snap.put(0);
   snap.write(name.c_str(),static_cast<std::streamsize>(name.size()));
   snap.put(0);
}

void CNemoSnapshotWriter::WriteHistory(const std::string &new_h)
{
   if (new_h.size() > 0)
       PutString("History",new_h.c_str());
}
void CNemoSnapshotWriter::WritePhase(const CPointMassSetFloat* points, double time)
{       
   int nbody    = static_cast<int>(points->size());
   float* phase = new float[nbody * 6];
   float* mass  = new float[nbody];
   
   for(int i = 0 ; i < nbody ; i++)
   {
       for(unsigned int j = 0 ; j < N_DIM ; j++)
       {
           phase[i*2*N_DIM + j] = points->at(i).first.Pos[j];
           phase[(i*2+1)*N_DIM + j] = points->at(i).first.Vel[j];
       }
       mass[i] = points->at(i).second;
   }
   
   StartLevel("SnapShot");
   StartLevel("Parameters");
   PutInt("Nobj", nbody);
   PutDouble("Time", time);
   EndLevel();
   StartLevel("Particles");
   PutInt("CoordSystem",0201402);
   int tmp_dim[3];
   tmp_dim[0] = nbody;
   PutFloatArray("Mass",1,tmp_dim,mass);
   tmp_dim[0] = nbody;
   tmp_dim[1] = 2;
   tmp_dim[2] = 3;
   PutFloatArray("PhaseSpace",3,tmp_dim,phase);
   EndLevel();
   EndLevel();
   
   delete [] phase;
   delete [] mass;
}

}