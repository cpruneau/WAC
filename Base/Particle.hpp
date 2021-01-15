// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_Particle
#define WAC_Particle
#include "Factory.hpp"
using namespace std;


/////////////////////////////////////
// Class Particle
/////////////////////////////////////
class Particle
{
public:
  
  Particle();
  ~Particle();
  Particle(const Particle& other);
  Particle & operator=(const Particle & other);
  void printProperties(ostream & output);

  void setSourceIndex(int sourceIndex);
  void setPxPyPzE(double p_x, double p_y, double p_z, double p_e);
  void setPidPxPyPzE(long pid, long charge, double p_x, double p_y, double p_z, double p_e);
  void setPidPtPhiYEta(long _id,long _ch,double _pT,double _phi,double _y,double _eta);

  void boost(double ax, double ay, double az);
  void boostRapidity(double boost);

  int getSourceIndex() const
  {
  return sourceIndex;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  int   sourceIndex;
  int   pid;
  float charge,baryon;
  float px, py, pz, e, pt, y, eta, phi;
  int   ixEtaPhi, ixYPhi;

  static int factorySize;
  static Factory<Particle> * factory;
  static Factory<Particle> * getFactory();
};


#endif /* WAC_Particle */
