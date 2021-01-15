// Author: Claude Pruneau   05/08/2020

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.  
 **********************************************************************/
#ifndef WAC_GeometryMoments
#define WAC_GeometryMoments
#include "TString.h"

// ***************************************************************************
// GeometryMoments
//
// Just a container of moments of the geometry of AA collisions
// ***************************************************************************
class GeometryMoments
{
public:

  GeometryMoments();
  GeometryMoments(const GeometryMoments & source);
  virtual  ~GeometryMoments()
  {
  // no ops.
  }
  GeometryMoments& operator=(const GeometryMoments & source);
  void reset();

  double meanX;
  double meanY;
  double varX;
  double varY;
  double varXY;
  double epsX;
  double epsY;
  double epsMod;
  double psi2;
  double area;
  double cphiN[10];
  double sphiN[10];
  double rN[10];
  double psiN[10];
  double eccN[10];

  ClassDef(GeometryMoments,0)
};

#endif /* GeometryMoments_hpp */
