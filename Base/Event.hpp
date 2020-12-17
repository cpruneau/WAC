// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/

#ifndef Event_hpp
#define Event_hpp
#include "Particle.hpp"

class Event
{
protected:

  Event();
  static Event * event;

  public:

  virtual ~Event();
  virtual void clear();
  virtual void reset();
  virtual void printProperties(ostream & output);

  Particle * getParticleAt(int index);

  //////////////////////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////////////////////
  long   index;
  int    eventNumber;
  int    nProjectile;
  int    nTarget;
  int    nParticleTotal;
  int    nParticles; // number of particles to analyze
  double centrality;
  double multiplicity;
  double impactParameter;
  double other;

  int    nFilters;
  double *nFiltered;
  double *eFiltered;
  double *qFiltered;
  double *bFiltered;

  static Event * getEvent();

  ClassDef(Event,0)

};




#endif /* Event_hpp */
