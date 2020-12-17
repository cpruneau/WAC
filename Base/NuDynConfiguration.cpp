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

 Utility class used to defineNuDynConfiguration
 */

#include "NuDynConfiguration.hpp"

ClassImp(NuDynConfiguration);


NuDynConfiguration::NuDynConfiguration(const TString & name,
                                       const TString & type,
                                       const TString & version)
:
TaskConfiguration(name,type,version),
multiplicityType(NuDynConfiguration::AcceptedMultiplicity),
nBins_mult(20),
min_mult(0.0),
max_mult(200.0)
{
  // no ops
}

NuDynConfiguration::NuDynConfiguration(const NuDynConfiguration & source)
:
TaskConfiguration( source ),
multiplicityType(  source.multiplicityType ),
nBins_mult( source.nBins_mult ),
min_mult(   source.min_mult ),
max_mult(   source.max_mult )
{
}

NuDynConfiguration::~NuDynConfiguration()
{

}

NuDynConfiguration & NuDynConfiguration::operator=(const NuDynConfiguration & source)
{
  if (this!=&source)
    {
    TaskConfiguration::operator=( source );
    multiplicityType = source.multiplicityType;
    nBins_mult  =  source.nBins_mult;
    min_mult  =  source.min_mult;
    max_mult  =  source.max_mult;
    }
  return *this;
}



////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void NuDynConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  os
  << "   NuDyn Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "  multiplicityType : " << multiplicityType << endl
  << "        nBins_mult : " << nBins_mult << endl
  << "          min_mult : " << min_mult << endl
  << "          max_mult : " << max_mult << endl;
}
