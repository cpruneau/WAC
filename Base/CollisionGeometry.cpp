// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Nucleon
 \ingroup WAC

 Class defining Nucleon
 */
#include <iostream>
#include "CollisionGeometry.hpp"

using namespace std;

ClassImp(CollisionGeometry);

CollisionGeometry::CollisionGeometry(int nProtonsA, int nNeutronsA, int nProtonsB, int nNeutronsB)
:
nucleusA(nullptr),
nucleusB(nullptr),
nBinaryMax(0),
nucleonsA(0),
nucleonsB(0),
b(0),
x(0),
y(0),
z(0),
t(0),
nBinary(0),
nParticipant(0),
nPP(0),
nPN(0),
nNN(0),
binaryMoments(),
participantMoments()
{
  nucleusA = new Nucleus(nProtonsA,nNeutronsA);
  nucleusB = new Nucleus(nProtonsB,nNeutronsB);
  initialize();
}

void CollisionGeometry::initialize()
{
  nBinaryMax = nucleusA->nNucleons * nucleusB->nNucleons;
  if (nBinaryMax<100)
    nBinaryMax = 100;
  else if (nBinaryMax>40000)
    nBinaryMax = 40000;
  nucleonsA = new Nucleon * [nBinaryMax];
  nucleonsB = new Nucleon * [nBinaryMax];
  x = new double[nBinaryMax];
  y = new double[nBinaryMax];
  z = new double[nBinaryMax];
  t = new double[nBinaryMax];
}

void CollisionGeometry::reset()
{
  b = 0;
  nucleusA->reset();
  nucleusB->reset();

  // this loop not necessary because nBinary is set to zero...
//  for (int iNN=0; iNN<nBinary; iNN++)
//    {
//    nucleonsA[nBinary] = nullptr;
//    nucleonsB[nBinary] = nullptr;
//    x[nBinary] = 0;
//    y[nBinary] = 0;
//    z[nBinary] = 0;
//    t[nBinary] = 0;
//    }
  nBinary = 0;
  nParticipant = 0;
  nPP = 0;
  nPN = 0;
  nNN = 0;
}

void CollisionGeometry::addNNCollision(Nucleon* nucleonA,Nucleon* nucleonB, double xAvg, double yAvg, double zAvg, double tAvg)
{
  if (nBinary>=nBinaryMax)
    {
    cout << "<E> CollisionGeometry::addNNCollision(...) Container size exceeded." << endl;
    return;
    }
  nucleonsA[nBinary] = nucleonA;
  nucleonsB[nBinary] = nucleonB;
  x[nBinary] = xAvg;
  y[nBinary] = yAvg;
  z[nBinary] = zAvg;
  t[nBinary] = tAvg;
  nBinary++;
  if (nucleonA->isProton())
    {
    if (nucleonB->isProton())
      nPP++;
    else
      nPN++;
    }
  else
    {
    if (nucleonB->isProton())
      nPN++;
    else
      nNN++;
    }

}


void CollisionGeometry::countParticipants()
{
  nParticipant = nucleusA->countWounded() + nucleusB->countWounded();
}

Nucleus * CollisionGeometry::getNucleusA()
{
  return nucleusA;
}

Nucleus * CollisionGeometry::getNucleusB()
{
  return nucleusB;
}

CollisionGeometry::NucleonNucleonType CollisionGeometry::getNNTpye(int index)
{
  Nucleon::NucleonType nucleonTypeA = nucleonsA[index]->nucleonType;
  Nucleon::NucleonType nucleonTypeB = nucleonsB[index]->nucleonType;
  if (nucleonTypeA == Nucleon::Proton)
    {
    if (nucleonTypeB == Nucleon::Proton)
      return ProtonProton;
    else
      return ProtonNeutron;
    }
  else
    {
    if (nucleonTypeB == Nucleon::Proton)
      return NeutronProton;
    else
      return NeutronNeutron;
    }
}
