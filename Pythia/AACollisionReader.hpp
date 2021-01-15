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
  virtual void execute();//overloaded


  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////

  int nMax;
  int nCollisions;

  CollisionGeometry * collisionGeometry;

  ////////////////////////////////////////////////////////////
  //All the branches and leafs from superclass

  ClassDef(AACollisionReader,0)
};


#endif /*WAC_AACollisionReader*/