// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/
/**
 \class Property
 \ingroup WAC

 Utility class used to defineRadialBoostConfiguration
 */

#include "RadialBoostConfiguration.hpp"
#include "TMath.h"

ClassImp(RadialBoostConfiguration);


RadialBoostConfiguration::RadialBoostConfiguration(const TString & name,
                                             const TString & type,
                                             const TString & version)
:
TaskConfiguration(name,type,version),
 param_a(0.9),    param_b(1.0),
 nBins_phi(72),   min_phi(0.0),  max_phi( TMath::TwoPi() ),
 nBins_r(100),    min_r(0.0),    max_r(1.0),
 nBins_beta(100), min_beta(0.0), max_beta(1.0)
{
  // no  ops
}

RadialBoostConfiguration::RadialBoostConfiguration(const RadialBoostConfiguration & source)
:
TaskConfiguration( source ),
param_a(source.param_a),       param_b(source.param_b),
nBins_phi(source.nBins_phi),   min_phi(source.min_phi),   max_phi(source.max_phi),
nBins_r(source.nBins_r),       min_r(source.min_r),       max_r(source.max_r),
nBins_beta(source.nBins_beta), min_beta(source.min_beta), max_beta(source.max_beta)
{
  // no  ops
}

RadialBoostConfiguration & RadialBoostConfiguration::operator=(const RadialBoostConfiguration & source)
{
  if (this!=&source)
    {
    TaskConfiguration::operator=( source );
    param_a    = source.param_a;
    param_b    = source.param_b;
    nBins_phi  = source.nBins_phi;
    min_phi    = source.min_phi;
    max_phi    = source.max_phi;
    nBins_r    = source.nBins_r;
    min_r      = source.min_r;
    max_r      = source.max_r;
    nBins_beta = source.nBins_beta;
    min_beta   = source.min_beta;
    max_beta   = source.max_beta;
    }
  return *this;
}



////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void RadialBoostConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  os
  << "               Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                  param_a: " << param_a << endl
  << "                  param_b: " << param_b << endl
  << "                nBins_phi: " << nBins_phi << endl
  << "                  min_phi: " << min_phi << endl
  << "                  max_phi: " << max_phi << endl
  << "                  nBins_r: " << nBins_r << endl
  << "                    min_r: " << min_r << endl
  << "                    max_r: " << max_r << endl
  << "               nBins_beta: " << nBins_beta << endl
  << "                 min_beta: " << min_beta << endl
  << "                 max_beta: " << max_beta << endl;
}
