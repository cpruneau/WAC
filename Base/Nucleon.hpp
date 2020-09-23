// Author: Claude Pruneau   05/08/2020

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.  
 **********************************************************************/
#ifndef WAC_Nucleon
#define WAC_Nucleon
#include "TObject.h"

// ***************************************************************************
// Nucleon
//
// type          : proton or neutron (can be used in PYTHIA)
//                 all length units in fm
// x,y,z         : cartesian coordinates relative to its container (nucleus)
// r, phi, theta : spherical coordinates relative to its container (nucleus)
//                 redundant but memory is cheap
// radius        : radius of the nucleon used in collisions
// ***************************************************************************
class Nucleon : public  TObject 
{
public:

  enum NucleonType   { Proton, Neutron, Unknown};

  NucleonType  nucleonType;
  double x,y,z, r, phi, theta; // memory is cheap...
  bool   wounded;

  Nucleon();
  //Nucleon(double r, double theta, double phi, NucleonType  nucleonType=Proton);
  virtual ~Nucleon()
  {
  // no ops.
  }

  virtual void reset();
  virtual void setNucleonType(NucleonType _nucleonType);
  virtual void setXYZ(double _x, double _y, double _z);
  virtual void setRThetaPhi(double _r, double _theta, double _phi);
  virtual void setRCosThetaPhi(double _r, double _cosTheta, double _phi);

  virtual void setWounded(bool _wounded)
  {
  wounded = _wounded;
  }
  virtual bool isWounded() const
  {
  return wounded;
  }

  ClassDef(Nucleon,0)
  
};

#endif /* Nucleon_hpp */
