// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class GeometryMoments
 \ingroup WAC

 Class defining GeometryMoments
 */

#include "GeometryMoments.hpp"
ClassImp(GeometryMoments);

GeometryMoments::GeometryMoments()
:
meanX( 0 ),
meanY( 0 ),
varX( 0 ),
varY( 0 ),
varXY( 0 ),
epsX( 0 ),
epsY( 0 ),
epsMod( 0 ),
psi2( 0 ),
area( 0 )
{
  reset();
}

GeometryMoments::GeometryMoments(const GeometryMoments & source)
:
meanX( source.meanX ),
meanY( source.meanX ),
varX ( source.varX ),
varY ( source.varY ),
varXY( source.varXY ),
epsX ( source.epsX ),
epsY ( source.epsY ),
epsMod( source.epsMod ),
psi2  ( source.psi2 ),
area  ( source.area )
{
  for (int iN=0;iN<10;iN++)
  {
  cphiN[iN] = source.cphiN[iN];
  sphiN[iN] = source.sphiN[iN];
  rN[iN]    = source.rN[iN];
  psiN[iN]  = source.psiN[iN];
  eccN[iN]  = source.eccN[iN];
  }
}

GeometryMoments& GeometryMoments::operator=(const GeometryMoments & source)
{
  if (&source != this)
    {
    meanX = source.meanX;
    meanY = source.meanX;
    varX  = source.varX;
    varY  = source.varY;
    varXY = source.varXY;
    epsX  = source.epsX;
    epsY  = source.epsY;
    epsMod = source.epsMod;
    psi2   = source.psi2;
    area   = source.area;
    for (int iN=0;iN<10;iN++)
      {
      cphiN[iN] = source.cphiN[iN];
      sphiN[iN] = source.sphiN[iN];
      rN[iN]    = source.rN[iN];
      psiN[iN]  = source.psiN[iN];
      eccN[iN]  = source.eccN[iN];
      }
    }
  return *this;
}

void GeometryMoments::reset()
{
  meanX = 0.0;
  meanY = 0.0;
  varX = 0.0;
  varY = 0.0;
  varXY = 0.0;
  epsX = 0.0;
  epsY = 0.0;
  epsMod = 0.0;
  psi2 = 0.0;
  area = 0.0;
  for (int iN=0;iN<10;iN++)
  {
  cphiN[iN] = 0.0;
  sphiN[iN] = 0.0;
  rN[iN]    = 0.0;
  psiN[iN]  = 0.0;
  eccN[iN]  = 0.0;
  }
}
