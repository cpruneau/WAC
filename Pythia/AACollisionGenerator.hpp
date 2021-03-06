#ifndef WAC_AACollisionGenerator
#define WAC_AACollisionGenerator
#include "TRandom1.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaEventGenerator.hpp"
#include "CollisionGeometry.hpp"

class AACollisionGenerator : public Task
{
public:


  //////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////
  AACollisionGenerator(const TString & name,
  TaskConfiguration * configuration,
  Event * event,
  EventFilter * ef,
  ParticleFilter * pf,
  CollisionGeometry * collisionGeo);
  virtual ~AACollisionGenerator();
  virtual void initialize();
  virtual void finalize();
  virtual void reset();
  void execute();


  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////
  EventFilter * eventFilter;
  ParticleFilter * particleFilter;

  int nMax;
  int nCollisions;

  CollisionGeometry * collisionGeometry;

  TPythia8* pythia8; 

  TRandom1* g;

  TClonesArray* particles;


  ClassDef(AACollisionGenerator,0)
};


#endif /*WAC_AACollisionGenerator*/