// Author: Victor Gonzalez on 02/10/19 on an original idea from Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */

#include "GeneralizedGaussianConfiguration.hpp"

ClassImp(GeneralizedGaussianConfiguration);

GeneralizedGaussianConfiguration::GeneralizedGaussianConfiguration(const TString & name,
                                                                   const TString &  version)
:
 TaskConfiguration(name,"GeneralizedGaussianConfiguration", version),
 amplitude(0.0),
 gammaeta(0.0),
 gammaphi(0.0),
 omegaeta(0.0),
 omegaphi(0.0)
{
  // no ops
}

GeneralizedGaussianConfiguration::GeneralizedGaussianConfiguration(const GeneralizedGaussianConfiguration & source)
:
 TaskConfiguration(source),
 amplitude( source.amplitude ),
 gammaeta( source.gammaeta ),
 gammaphi( source.gammaphi ),
 omegaeta( source.omegaeta ),
 omegaphi( source.omegaphi )
{
  // no ops
}

GeneralizedGaussianConfiguration & GeneralizedGaussianConfiguration::operator=(const GeneralizedGaussianConfiguration & source)
{
  if (this != &source)
    {
    TaskConfiguration::operator=(source);
    amplitude = source.amplitude;
    gammaeta =  source.gammaeta;
    gammaphi =  source.gammaphi;
    omegaeta =  source.omegaeta;
    }
  return *this;
}

void GeneralizedGaussianConfiguration::printConfiguration(ostream & os)
{
  os
  << " ------------------------------------------------------------------------------------------" << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                                 name: " << getName()    << endl
  << "                              version: " << getVersion() << endl
  << "                            amplitude: " << amplitude     << endl
  << "                      gamma delta eta: " << gammaeta      << endl
  << "                      gamma delta phi: " << gammaphi      << endl
  << "                      omega delta eta: " << omegaeta      << endl
  << "                      omega delta phi: " << omegaphi      << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                      DELTA ETA WIDTH: " << TMath::Sqrt(omegaeta*omegaeta*TMath::Gamma(3./gammaeta)/TMath::Gamma(1./gammaeta)) << endl
  << "                      DELTA PHI WIDTH: " << TMath::Sqrt(omegaphi*omegaphi*TMath::Gamma(3./gammaphi)/TMath::Gamma(1./gammaphi)) << endl
  << " ------------------------------------------------------------------------------------------" << endl;
}
