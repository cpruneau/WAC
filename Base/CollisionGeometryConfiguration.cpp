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

 Utility class used to defineCollisionGeometryConfiguration
 */

#include "CollisionGeometryConfiguration.hpp"

ClassImp(CollisionGeometryConfiguration);


CollisionGeometryConfiguration::CollisionGeometryConfiguration(const TString & name,
                                                               const TString & type,
                                                               const TString & version)
:
TaskConfiguration(name,type,version),
histoBaseName(""),
minB(0.0), maxB(18.0),
nnCrossSection(5555),
nBins_b(100), min_b(0.0), max_b(18.0),
nBins_nPart(500), min_nPart(0), max_nPart(500),
nBins_nBinary(500), min_nBinary(0), max_nBinary(40000)
{
  
}

////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void CollisionGeometryConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  os
  << "    Analysis   Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "      Histogram Base Name: " << histoBaseName << endl
  << "                     minB: " << minB << endl
  << "                     maxB: " << minB << endl
  << "              nBins_nPart: " << nBins_nPart << endl
  << "                min_nPart: " << min_nPart << endl
  << "                max_nPart: " << max_nPart << endl
  << "            nBins_nBinary: " << nBins_nBinary << endl
  << "              min_nBinary: " << min_nBinary << endl
  << "              max_nBinary: " << max_nBinary << endl;
}
