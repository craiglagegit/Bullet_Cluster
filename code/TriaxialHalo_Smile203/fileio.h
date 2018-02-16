#pragma once
#include "common.h"
#include <fstream>

namespace smile{

/// routine that splits one line from a text file into several items, separated by any set of chars from "delim"; empty items do not count; 
/// resulting items are stored in vector given by "result" (must be a pointer to an existing vector)
void splitString(const std::string& src, const std::string& delim, std::vector<std::string> *result);

struct CConfigPotential;

class CFileInputText
{
public:
    CFileInputText() { status=isopen=false; };
    explicit CFileInputText(const std::string &_baseFileName) { status=isopen=false; open(_baseFileName); };
    ~CFileInputText();
    void open(const std::string &_baseFileName);
    bool ok() const { return status; }
    COrbit* readOrbit(const CPotential* _potential);   // loads trajectory from file and initializes orbit, returns NULL on failure
    COrbitDesc* readOrbitDesc(const CPotential* _potential, bool withModelData);   // loads orbit description from single line in a text file, and possibly also reads binary data for trajectory sample and Schwarzschild model
    COrbitLibrary* readOrbitLibrary(const CPotential* _potential, bool withModelData);   // loads entire orbit library
    const CPotential* readPotential(CConfigPotential *configPotential);   // loads either potential coefficients from a text file, or initializes potential from N point masses stored in a text file. Updates parameters in configPotential
private:
    std::string baseFileName;
    bool status;
    bool isopen;
    std::ifstream strmText;
    std::ifstream strmSchwData;
    std::ifstream strmTrajSample;
};

class CFileOutputText
{
public:
    CFileOutputText() { status=isopen=started=false; };
    explicit CFileOutputText(const std::string &_baseFileName) { status=isopen=started=false; open(_baseFileName); };
    void open(const std::string &_baseFileName);
    ~CFileOutputText();
    bool ok() const { return status; }
    bool writeOrbit(const COrbit* value);
    bool writeOrbitDesc(const COrbitDesc* value, bool withModelData);
    bool writeOrbitLibrary(const COrbitLibrary* value, bool withModelData);
    bool writePotential(const CPotential* potential);
    bool writePointMassSet(const CPointMassSetFloat* points);
private:
    std::string baseFileName;
    bool status;
    bool isopen;
    bool started;
    std::ofstream strmText;
    std::ofstream strmSchwData;
    std::ofstream strmTrajSample;
};

// helper class that writes NEMO-compatible snapshot file
class CNemoSnapshotWriter
{
public:
   CNemoSnapshotWriter(const std::string &filename, bool append=false);
   ~CNemoSnapshotWriter();
   void PutInt(const std::string &name,int i);
   void PutDouble(const std::string &name,double d);
   // write array of double; ndim - number of dimensions, dim - length of array for each dimension
   void PutDoubleArray(const std::string &name, int ndim, const int dim[], const double* data);
   void PutFloatArray(const std::string &name, int ndim, const int dim[], const float* data);
   void PutIntArray(const std::string &name,int ndim,const int dim[],const int* data);
   void PutCharArray(const std::string &name,int ndim,const int dim[],const char* data);
   void StartLevel(const std::string &name);
   void EndLevel();
   void PutString(const std::string &name,const std::string &str);
   void WriteHistory(const std::string &new_h);
   void WritePhase(const CPointMassSetFloat* points, double time);
   bool ok() const { return snap.good(); }
private:
   std::ofstream snap;
   int level;  //index of current level (root level is 0)
   void PutZString(const std::string &name); 
};

class CFileOutputNemo
{
public:
    CFileOutputNemo() { SnapshotWriter=NULL; };
    explicit CFileOutputNemo(const std::string &fileName, const std::string &header="") { SnapshotWriter=NULL; open(fileName, header); };
    void open(const std::string &fileName, const std::string &header="");
    ~CFileOutputNemo() { if(SnapshotWriter!=NULL) delete SnapshotWriter; };
    bool ok() const { if(SnapshotWriter!=NULL) return SnapshotWriter->ok(); else return false; };
    bool writePointMassSet(const CPointMassSetFloat* points) { 
        if(points!=NULL && SnapshotWriter!=NULL) {
            SnapshotWriter->WritePhase(points, 0); 
            return SnapshotWriter->ok();
        }
        return false;
    };
private:
    CNemoSnapshotWriter* SnapshotWriter;
};

}