// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Nucleus
 \ingroup WAC

 Class defining Nucleus
 */

#include "Nucleus.hpp"
ClassImp(Nucleus);

Nucleus::Nucleus(int _nProtons, int _nNeutrons)
:
TObject(),
nProtons(_nProtons),
nNeutrons(_nNeutrons),
nNucleons(_nProtons+_nNeutrons),
nucleons(0)
{
  initialize();
}

Nucleus::Nucleus()
:
nProtons(82),
nNeutrons(126),
nNucleons(208),
nucleons(0)
{
  initialize();
}

//virtual void initialize();
//virtual Nucleon * getNucleon(int index);
//virtual void reset();
//virtual int  countWounded();

// create the container but do not assign any positions or properties.
void Nucleus::initialize()
{
  nucleons = new Nucleon*[nNucleons];
  for (int iNucleon=0; iNucleon<nNucleons; iNucleon++)
    {
    nucleons[iNucleon] = new Nucleon();
    }
}

Nucleon * Nucleus::getNucleon(int index)
{
  if (index<0)
    index = 0;
  else if (index>=nNucleons)
    index = nNucleons-1;
  return nucleons[index];
}

void Nucleus::reset()
{
//  x = -99999.0;
//  y = -99999.0;
//  z = -99999.0;
  for (int iNucleon=0; iNucleon<nNucleons; iNucleon++)
    {
    nucleons[iNucleon]->reset();
    }
}

int Nucleus::countWounded()
{
  int wounded = 0;
  for (int iNucleon=0; iNucleon<nNucleons; iNucleon++)
    {
    wounded += nucleons[iNucleon]->isWounded();
    }
  return wounded;
}
