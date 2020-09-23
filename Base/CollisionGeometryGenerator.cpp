// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */
#include "CollisionGeometryGenerator.hpp"
ClassImp(CollisionGeometryGenerator);

CollisionGeometryGenerator::CollisionGeometryGenerator(const TString & name,
                                                       CollisionGeometryConfiguration * configuration,
                                                       CollisionGeometry * _collisionGeometry,
                                                       NucleusGenerator * _nucleusGeneratorA,
                                                       NucleusGenerator * _nucleusGeneratorB)
:
Task(name, configuration, nullptr),
collisionGeometry(_collisionGeometry),
nucleusGeneratorA(_nucleusGeneratorA),
nucleusGeneratorB(_nucleusGeneratorB)
{
  minB   = configuration->minB;
  minBSq = minB*minB;
  maxB   = configuration->maxB;
  maxBSq = maxB*maxB;
  nnCrossSection  = configuration->nnCrossSection;
  maxNNDistanceSq = nnCrossSection/3.1415927;
  if (reportDebug()) cout << "CollisionGeometryGenerator::CollisionGeometryGenerator(...) No ops" << endl;
}


void CollisionGeometryGenerator::reset()
{
 collisionGeometry->reset();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read an ampt event from file
// Copy the event into Event for convenience...
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CollisionGeometryGenerator::execute()
{
  if (reportDebug()) cout << "CollisionGeometryGenerator::execute() Started" << endl;

  Nucleus * nucleusA = collisionGeometry->getNucleusA();
  Nucleus * nucleusB = collisionGeometry->getNucleusB();
  nucleusGeneratorA->generate(nucleusA);
  nucleusGeneratorB->generate(nucleusB);

  double rr = gRandom->Rndm();
  double b = sqrt(minBSq + rr*(maxBSq-minBSq));
  collisionGeometry->b = b;
  double xA, yA, zA;
  double xB, yB, zB;
  double xAvg, yAvg, zAvg, tAvg, dx, dy, dz, dSq;

  for (int i1=0; i1<nucleusA->nNucleons; i1++)
    {
    Nucleon * nucleonA = nucleusA->getNucleon(i1);
    xA = nucleonA->x;
    yA = nucleonA->y;
    zA = nucleonA->z;
    for (int i2=0; i2<nucleusB->nNucleons; i2++)
      {
      Nucleon * nucleonB = nucleusB->getNucleon(i2);
      xB = nucleonB->x;
      yB = nucleonB->y;
      zB = nucleonB->z;
      dx = xB - xA + b;
      dy = yB - yA;
      dz = zB - zA;
      dSq = dx*dx + dy*dy;
      if (dSq<maxNNDistanceSq)
        {
        // this is a collision
        xAvg = (xA + xB)/2.0;
        yAvg = (yA + yB)/2.0;
        zAvg = (zA + zB)/2.0;
        tAvg = 0.0; // for now...
        nucleonA->setWounded(true);
        nucleonB->setWounded(true);
        collisionGeometry->addNNCollision(nucleonA,nucleonB,xAvg,yAvg,zAvg,tAvg);
        }
      }
    }
  collisionGeometry->countParticipants();

  if (reportDebug()) cout << "CollisionGeometryGenerator::execute() Completed" << endl;
}
