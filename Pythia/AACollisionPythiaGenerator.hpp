#ifndef WAC_AACollisionPythiaGenerator
#define WAC_AACollisionPythiaGenerator
#include "TRandom.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "CollisionGeometry.hpp"
#include "PythiaConfiguration.hpp"

class AACollisionPythiaGenerator : public Task
{
public:

  AACollisionPythiaGenerator(const TString & name,
                             PythiaConfiguration * configuration,
                             CollisionGeometry * collisionGeometry,
                             Event * event,
                             EventFilter * ef,
                             ParticleFilter * pf,
                             LogLevel selectedLevel);
  virtual ~AACollisionPythiaGenerator();
  virtual void initialize();
  virtual void finalize();
  void execute();

  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////
  CollisionGeometry * collisionGeometry;
  EventFilter       * eventFilter;
  ParticleFilter    * particleFilter;
  TPythia8          * ppPythia;
  TPythia8          * pnPythia;
  TPythia8          * npPythia;
  TPythia8          * nnPythia;
  TClonesArray      * particles;
  int    clonesArraySize;
  bool   removePhotons;
  bool   ppOnly; 
  ClassDef(AACollisionPythiaGenerator,0)
};


#endif /*WAC_AACollisionPythiaGenerator*/
