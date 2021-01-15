// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Particle
 \ingroup WAC

 Class defining Particle
 */
#include <ostream>
#include "TLorentzVector.h"
#include "Particle.hpp"

Particle::Particle()
:
sourceIndex(0),
pid(0),
charge(0),
px(0.0), py(0.0), pz(0.0), e(0.0), pt(0.0), y(0.0), eta(0.0), phi(0.0), ixEtaPhi(-1), ixYPhi(-1)
{
  // no op
}

/////////////////////////////////////
// DTOR
/////////////////////////////////////
Particle::~Particle(){}

/////////////////////////////////////
// Copy the passed particle in this object
/////////////////////////////////////
Particle::Particle(const Particle& other)
:
sourceIndex( other.sourceIndex ),
pid      ( other.pid ),
charge   ( other.charge ),
px       ( other.px ),
py       ( other.py ),
pz       ( other.pz ),
e        ( other.e ),
pt       ( other.pt ),
y        ( other.y ),
eta      ( other.eta ),
phi      ( other.phi ),
ixEtaPhi ( other.ixEtaPhi ),
ixYPhi   ( other.ixYPhi )
{
  // no ops
}


Particle & Particle::operator=(const Particle & other)
{
  if (&other != this)
    {
    sourceIndex = other.sourceIndex;
    pid       = other.pid;
    charge    = other.charge;
    px        = other.px;
    py        = other.py;
    pz        = other.pz;
    e         = other.e;
    pt        = other.pt;
    y         = other.y;
    eta       = other.eta;
    phi       = other.phi;
    ixEtaPhi  = other.ixEtaPhi;
    ixYPhi    = other.ixYPhi;
    }
  return *this;
}
//////////////////////////////////////////////////////////////////////////
// Print the properties of this particle at the given output
//////////////////////////////////////////////////////////////////////////
void Particle::printProperties(ostream & output)
{
  double mass;
  double massSq = e*e - pt*pt - pz*pz;
  if (massSq<0.0)
    {
    massSq = 0.0;
    mass   = 0.0;
    }
  else
    {
    mass = sqrt(massSq);
    }
  output << "  Source: " << sourceIndex << endl;
  output << "     PID: " << pid << endl;
  output << "  charge: " << charge << endl;
  output << "    mass: " << mass << endl;
  output << "      px: " << px << endl;
  output << "      py: " << py << endl;
  output << "      pt: " << pt << endl;
  output << "      pz: " << pz << endl;
  output << "       e: " << e << endl;
  output << "       y: " << y << endl;
  output << "     eta: " << eta << endl;
  output << "     phi: " << phi << endl;
}

void Particle::setSourceIndex(int _sourceIndex)
{
  sourceIndex = _sourceIndex;
}

void Particle::setPxPyPzE(double p_x, double p_y, double p_z, double p_e)
{
  px  = p_x;
  py  = p_y;
  pz  = p_z;
  e   = p_e;
  pt  = sqrt(px*px + py*py);
  phi = atan2(py,px);
  if (phi<0) phi += 2.0*TMath::Pi();
  double theta = atan2(pt,pz);
  eta = -log(tan(theta/2.0));
  double plus  = e+pz;
  double minus = e-pz;
  if (plus==0)
    y = -1.0E50;
  else if (minus==0)
    y = 1.0E50;
  else
    y = 0.5 * log( plus/minus );

}

void Particle::setPidPxPyPzE(long thePid, long theCharge, double p_x, double p_y, double p_z, double p_e)
{
  pid = thePid;
  charge = theCharge;
  px  = p_x;
  py  = p_y;
  pz  = p_z;
  e   = p_e;
  pt  = sqrt(px*px + py*py);
  phi = atan2(py,px);
  if (phi<0) phi += 2.0*TMath::Pi();
  double theta = atan2(pt,pz);
  eta = -log(tan(theta/2.0));
  double plus  = e+pz;
  double minus = e-pz;
  if (plus==0)
    y = -1.0E50;
  else if (minus==0)
    y = 1.0E50;
  else
    y = 0.5 * log( plus/minus );
}

void Particle::setPidPtPhiYEta(long _id,long _ch,double _pT,double _phi,double _y,double _eta) {
  pid = _id;
  charge = _ch;
  phi = TVector2::Phi_0_2pi(_phi);
  pt = TMath::Abs(_pT);
  eta = _eta;
  y = _y;
  px = pt*cos(phi);
  py = pt*sin(phi);
  pz = pt/TMath::Tan(2.0*TMath::ATan(TMath::Exp(-eta)));
  e = 1.0/TMath::TanH(_y);
  ixEtaPhi = -1;
  ixYPhi = -1;
}


void Particle::boost(double ax, double ay, double az)
{
  TLorentzVector vec;
  vec.SetPxPyPzE(px,py,pz,e);
  vec.Boost(ax,ay,az);
  px = vec.Px();
  py = vec.Py();
  pz = vec.Pz();
  e  = vec.E();
  pt  = sqrt(px*px + py*py);
  phi = atan2(py,px);
  if (phi<0) phi += TMath::TwoPi();
  double theta = atan2(pt,pz);
  eta = -log(tan(theta/2.0));
  double plus  = e+pz;
  double minus = e-pz;
  if (plus==0)
    y = -1.0E50;
  else if (minus==0)
    y = 1.0E50;
  else
    y = 0.5 * log( plus/minus );
}

void Particle::boostRapidity(double boost)
{
  y += boost;
  double mt = sqrt(e*e - pz*pz);
  pz = mt * sinh(y);
  e  = mt * cosh(y);
  //
  //pt  = sqrt(px*px + py*py);
  //phi = atan2(py,px);
  double theta = atan2(pt,pz);
  eta = -log(tan(theta/2.0));
  double plus  = e+pz;
  double minus = e-pz;
  if (plus==0)
    y = -1.0E50;
  else if (minus==0)
    y = 1.0E50;
  else
    y = 0.5 * log( plus/minus );
}

int Particle::factorySize = 5000;
Factory<Particle> * Particle::factory = 0;
Factory<Particle> * Particle::getFactory()
{
  if (!factory)
    {
    factory = new Factory<Particle>();
    factory->initialize(factorySize);
    }
  return factory;
}
