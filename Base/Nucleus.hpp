// Author: Claude Pruneau   05/08/2020

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.  
 **********************************************************************/
#ifndef WAC_Nucleus
#define WAC_Nucleus
#include "TObject.h"
#include "Nucleon.hpp"

// ***************************************************************************
// Nucleus
//
// nNucleons     : number of nucleons
// nProtons      : number of protons
// nNeutrons     : number of neutron
//
// x,y,z         : cartesian coordinates relative to its container (collision)
// radius        : nominal radius of the nucleus used for generation purposes
// ***************************************************************************
class Nucleus  : public  TObject 
{
public:

  int nProtons;
  int nNeutrons;
  int nNucleons;

protected:

  Nucleon ** nucleons;

  public:


  Nucleus(int nProtons, int nNeutrons);
  Nucleus();
  virtual  ~Nucleus()
  {
  // no ops.
  }

  virtual void initialize();
  virtual Nucleon * getNucleon(int index);
  virtual void reset();
  virtual int  countWounded();

  //void setPosition(double _x, double _y, double _z);


  ClassDef(Nucleus,0)
  
};

#endif /* Nucleus_hpp */
