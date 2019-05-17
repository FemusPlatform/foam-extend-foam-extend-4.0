#ifndef __ParCommunicator__
#define __ParCommunicator__

#include "fvMesh.H"
#include <vector>
using namespace Foam;

class ParCommunicator{
public:
    Foam::fvMesh & _OFmesh;
    const int _nProcs;
    
    std::map<std::string, std::vector<int>> _OtherProcMap;
    std::map<std::string, int> _ThisProcMap;
    std::map<std::string, int> _PatchNameToId;
    std::map<int, std::string> _PatchIdToName;
    
    ParCommunicator(Foam::fvMesh &FoamMesh);
    ~ParCommunicator(){};    
    void ListBoundaryPatches();
    void FillProcMap();   
    void PrintThisProcMap(int proc);
    void PrintOtherProcMapFromProc(int proc);
    void CommunicateMap();
};

#endif
