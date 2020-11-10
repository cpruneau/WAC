#ifndef WAC_AACollisionReader
#define WAC_AACollisionReader
#include "TRandom1.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaEventGenerator.hpp"
#include "PythiaEventReader.hpp"
#include "CollisionGeometry.hpp"
#include "PythiaEventReader.hpp" 

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
    CollisionGeometry * collisionGeo);
  virtual ~AACollisionReader();
  virtual void execute();
  virtual void initialize();
  virtual void Init(TTree *tree);


  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////

  int nMax;
  int nCollisions;

  CollisionGeometry * collisionGeometry;

  ////////////////////////////////////////////////////////////
  //Leaf types
  //Int_t           eventNo; // in superclass
  //Int_t           mult; // in superclass
  //Float_t         px[arraySize];   //[Mult] // in superclass
  //Float_t         py[arraySize];   //[Mult] // in superclass
  //Float_t         pz[arraySize];   //[Mult] // in superclass
  int           ist[arraySize]; //[Mult] 
  int           pdg[arraySize]; //[Mult] 
  double         pE[arraySize]; //[Mult] 
  //TClonesArray *particles;

  ////////////////////////////////////////////////////////
  //Branch types
  //TBranch        *b_eventNo;   //! // in superclass
  //TBranch        *b_mult;   //!    // in superclass
  //TBranch        *b_px;   //!      // in superclass
  //TBranch        *b_py;   //!      // in superclass
  //TBranch        *b_pz;   //!      // in superclass
  TBranch        *b_ist;   //!
  TBranch        *b_pdg;   //!
  TBranch        *b_pE;   //!

  TString branchName_eventNo;
  TString branchName_mult;
  TString branchName_px;
  TString branchName_py;
  TString branchName_pz;
  TString branchName_ist;
  TString branchName_pdg;
  TString branchName_pE;

  ////////////////////////////////////////////////////////////
  //Requirements: 
  //One tree per file
  // treess must have branches with the above data
  //All trees have the same name
  //All files have the same name except for a final numeric subscript to distinguish them 
  //No .root extension on the files name
  TString eventTreeName; 
  TString fileLocation; 
  int numFiles;  //total number of root files



  ClassDef(AACollisionReader,0)
};


#endif /*WAC_AACollisionReader*/