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

#include "Nucleon.hpp"
#include "TMath.h"

ClassImp(Nucleon);

//double Nucleon::radius = 1.2;

Nucleon::Nucleon()
:
 nucleonType(Unknown),x(0.0),y(0.0),z(0.0), r(0.0), phi(0.0), theta(0.0), wounded(false)
{
}

//Nucleon::Nucleon(double _r, double _theta, double _phi, NucleonType  _nucleonType)
//:
//nucleonType(_nucleonType),x(0.0),y(0.0),z(0.0), r(0.0), phi(0.0), theta(0.0)
//{
//  setRThetaPhi( _r,  _theta,  _phi);
//  wounded = false;
//}

void Nucleon::reset()
{
  nucleonType = Unknown;
  x = y = z = r = phi = theta = 0.0;
  wounded = false;
}

void Nucleon::setNucleonType(NucleonType _nucleonType)
{
  nucleonType = _nucleonType;
}

void Nucleon::setXYZ(double _x, double _y, double _z)
{
  x = _x;
  y = _y;
  z = _z;
  r     = sqrt(x*x+y*y+z*z);
  phi   = atan2(y,x);
  theta = atan2(sqrt(x*x+y*y),z);
  wounded = false;
}


void Nucleon::setRThetaPhi(double _r, double _theta, double _phi)
{
  r     = _r;
  theta = _theta;
  phi   = _phi;
  double cosTheta = cos(theta);
  double sinTheta = sin(theta);
  x = r*sinTheta*cos(phi);
  y = r*sinTheta*sin(phi);
  z = r*cosTheta;
  wounded = false;
}

void Nucleon::setRCosThetaPhi(double _r, double _cosTheta, double _phi)
{
  r     = _r;
  phi   = _phi;
  double sinTheta = sqrt(1-_cosTheta*_cosTheta);
  theta = acos(_cosTheta);
  x = r*sinTheta*cos(phi);
  y = r*sinTheta*sin(phi);
  z = r*_cosTheta;
  wounded = false;
}

void Nucleon::shift(double dx, double dy, double dz)
{
  x += dx;
  y += dy;
  z += dz;
  r     = sqrt(x*x+y*y+z*z);
  phi   = atan2(y,x);
  theta = atan2(sqrt(x*x+y*y),z);
}

double Nucleon::distanceXYZSq(Nucleon  * otherNucleon)
{
  double dx = x - otherNucleon->x;
  double dy = y - otherNucleon->y;
  double dz = z - otherNucleon->z;
  return dx*dx + dy*dy + dz*dz;
}

double Nucleon::distanceXYZ(Nucleon * otherNucleon)
{
  double dx = x - otherNucleon->x;
  double dy = y - otherNucleon->y;
  double dz = z - otherNucleon->z;
  return sqrt(dx*dx + dy*dy + dz*dz);
}


double Nucleon::distanceXYSq(Nucleon * otherNucleon)
{
  double dx = x - otherNucleon->x;
  double dy = y - otherNucleon->y;
  return dx*dx + dy*dy;
}

double Nucleon::distanceXY(Nucleon * otherNucleon)
{
  double dx = x - otherNucleon->x;
  double dy = y - otherNucleon->y;
  return sqrt(dx*dx + dy*dy);
}
