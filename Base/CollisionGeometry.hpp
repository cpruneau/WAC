// Author: Claude Pruneau   05/08/2020

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.  
 **********************************************************************/
#ifndef WAC_CollisionGeometry
#define WAC_CollisionGeometry
#include "Nucleon.hpp"
#include "Nucleus.hpp"
#include "NucleusGenerator.hpp"

// ***************************************************************************
// CollisionGeometry
//
// type          : proton or neutron (can be used in PYTHIA)
//                 all length units in fm
// x,y,z         : cartesian coordinates relative to its container (nucleus)
// r, phi, theta : spherical coordinates relative to its container (nucleus)
//                 redundant but memory is cheap
// radius        : radius of the CollisionGeometry used in collisions
// ***************************************************************************
class CollisionGeometry : public TObject 
{
public:

  Nucleus * nucleusA;
  Nucleus * nucleusB;
  int nBinaryMax;
  Nucleon ** nucleonsA;
  Nucleon ** nucleonsB;
  double b;
  double * x;
  double * y;
  double * z;
  double * t;
  int nBinary;       // number of binary collisions
  int nParticipant;  // number of participants
  int nPP;           // number of proton-proton collisions
  int nPN;           // number of proton-neutron collisions
  int nNN;           // number of neutron-neutron collisions

  CollisionGeometry(int nProtonsA, int nNeutronsA, int nProtonsB, int nNeutronsB);

  virtual ~CollisionGeometry()
  {
  // no ops.
  }

  virtual void initialize();
  virtual void reset();

  Nucleus * getNucleusA();
  Nucleus * getNucleusB();

  virtual void addNNCollision(Nucleon* nucleonA,Nucleon* nucleonB, double xAvg, double yAvg, double zAvg, double tAvg);
  virtual void countParticipants();

  static CollisionGeometry * getCollisionGeometry();

  ClassDef(CollisionGeometry,0)
  
};

#endif /* CollisionGeometry_hpp */
