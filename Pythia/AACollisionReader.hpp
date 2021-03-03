#ifndef WAC_AACollisionReader
#define WAC_AACollisionReader
#include "TRandom1.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "AACollisionPythiaGenerator.hpp"
#include "PythiaEventReader.hpp"
#include "CollisionGeometry.hpp"

class AACollisionReader : public PythiaEventReader
{
public:


  //////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////
  AACollisionReader(const TString & name,
    TaskConfiguration * configuration,
    Event * event,
    EventFilter * ef,
    ParticleFilter * pf,
    LogLevel selectedLevel,
    CollisionGeometry *collisionGeo);
  virtual ~AACollisionReader();
  virtual void initialize(); //overloaded with both superclasses
  virtual void execute();//overloaded
  virtual void Init(int i, TTree * t);//overloaded
  virtual Long64_t LoadTree(int i, Long64_t entry);


  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////

  int nCollisions;
  bool removePhotons;
  bool ppOnly;
  //one tree for events from each type of collision
  //The branches can be the same because at any point we will only use one tree at a time
  //In the order of pp,pn,np,nn

  TTree  * ppfChain[4];   //!pointer to the analyzed TTree or TChain
  Int_t   ppfCurrent[4]; //!current Tree number in a TChain
  Long64_t  ppnentries[4];
  Long64_t  ppnbytes[4];
  Long64_t  ppnb[4];
  long  ppjentry[4];

  /*TTree  *pnfChain;   //!pointer to the analyzed TTree or TChain
  Int_t   pnfCurrent; //!current Tree number in a TChain
  TFile  *pninputDataFile;
  Long64_t pnnentries;
  Long64_t pnnbytes;
  Long64_t pnnb;
  long pnjentry;

  TTree  *npfChain;   //!pointer to the analyzed TTree or TChain
  Int_t   npfCurrent; //!current Tree number in a TChain
  TFile  *npinputDataFile;
  Long64_t npnentries;
  Long64_t npnbytes;
  Long64_t npnb;
  long npjentry;

  TTree  *nnfChain;   //!pointer to the analyzed TTree or TChain
  Int_t   nnfCurrent; //!current Tree number in a TChain
  TFile  *nninputDataFile;
  Long64_t nnnentries;
  Long64_t nnnbytes;
  Long64_t nnnb;
  long nnjentry;*/


  CollisionGeometry * collisionGeometry;

  ////////////////////////////////////////////////////////////
  //All the branches and leafs from superclass

  ClassDef(AACollisionReader,0)
};


#endif /*WAC_AACollisionReader*/