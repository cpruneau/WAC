// Author: Claude Pruneau   05/08/2020

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_NucleusGenerator
#define WAC_NucleusGenerator
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"
#include "TH1.h"
#include "Nucleus.hpp"
#include "Nucleon.hpp"

// ***************************************************************************
// NucleusGenerator
//
// nNucleons     : number of nucleons
// nProtons      : number of protons
// nNeutrons     : number of neutron
//
// x,y,z         : cartesian coordinates relative to its container (collision)
// radius        : nominal radius of the NucleusGenerator used for generation purposes
// ***************************************************************************
class NucleusGenerator  : public  TObject
{
public:

  enum GeneratorType   { Uniform, WoodsSaxon, Exponential, Gaussian, DoubleGaussian };

  GeneratorType generatorType;
  TString generatorName;

  int    nR;
  double minR, maxR;
  double parA;
  double parB;
  double parC;
  TH1 * rDensity;
  TH1 * rProfile;
  TH1 * rProfileGen;

  bool useRecentering;
  bool useNucleonExclusion;
  double exclusionRadius;
  double exclusionRadiusSq;


  NucleusGenerator(const TString & _generatorName,
                   GeneratorType _generatorType, double _parA, double _parB, double _parC,
                   int _nR, double _minR, double _maxR);
  virtual ~NucleusGenerator()
  {     }
  
  virtual void initialize();
  virtual void generate(Nucleus * nucleus);
  virtual void generate(double & r, double & cosTheta, double & phi);
  virtual void saveHistorgrams();

  ClassDef(NucleusGenerator,0)
  
};

#endif /* NucleusGenerator_hpp */
