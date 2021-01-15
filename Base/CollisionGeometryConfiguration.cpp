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
nProtonsA(0),
nNeutronsA(0),
aGeneratorType(NucleusGenerator::WoodsSaxon),
aParA(0.0),
aParB(0),
aParC(0),
aNR(0),
aMinR(0),
aMaxR(0),
nProtonsB(0),
nNeutronsB(0),
bGeneratorType(NucleusGenerator::WoodsSaxon),
bParA(0),
bParB(0),
bParC(0),
bNR(0),
bMinR(0),
bMaxR(0),
minB(0.0), maxB(0),
nnCrossSection(0),
nBins_b(0),
min_b(0),
max_b(0),
nBins_nPart(0),
min_nPart(0),
max_nPart(0),
nBins_nBinary(0),
min_nBinary(0),
max_nBinary(0),
nBins_bxSect(0),
bValues(nullptr)
{
  setDefaultConfiguration();
}

CollisionGeometryConfiguration::CollisionGeometryConfiguration(const CollisionGeometryConfiguration & source)
:
TaskConfiguration(source),
histoBaseName( source.histoBaseName  ),
nProtonsA( source.nProtonsA  ),
nNeutronsA( source.nNeutronsA  ),
aGeneratorType( source.aGeneratorType  ),
aParA( source.aParA  ),
aParB( source.aParB  ),
aParC( source.aParC  ),
aNR( source.aNR  ),
aMinR( source.aMinR  ),
aMaxR( source.aMaxR  ),
nProtonsB( source.nProtonsB  ),
nNeutronsB( source.nNeutronsB  ),
bGeneratorType( source.bGeneratorType  ),
bParA( source.bParA  ),
bParB( source.bParB  ),
bParC( source.bParC  ),
bNR( source.bNR  ),
bMinR( source.bMinR  ),
bMaxR( source.bMaxR  ),
minB( source.minB  ),
maxB( source.maxB  ),
nnCrossSection( source.nnCrossSection  ),
nBins_b( source.nBins_b  ),
min_b( source.min_b  ),
max_b( source.max_b  ),
nBins_nPart( source.nBins_nPart  ),
min_nPart( source.min_nPart  ),
max_nPart( source.max_nPart  ),
nBins_nBinary( source.nBins_nBinary  ),
min_nBinary( source.min_nBinary  ),
max_nBinary( source.max_nBinary  ),
nBins_bxSect( source.nBins_bxSect  ),
bValues(nullptr)
{
  setBXsectPartition(source.nBins_bxSect, source.bValues);
}

CollisionGeometryConfiguration & CollisionGeometryConfiguration::operator=(const CollisionGeometryConfiguration & source)
{
  if (this != &source)
    {
    TaskConfiguration::operator=(source);
    histoBaseName  =  source.histoBaseName;
    nProtonsA      =  source.nProtonsA;
    nNeutronsA     =  source.nNeutronsA;
    aGeneratorType =  source.aGeneratorType;
    aParA  =  source.aParA;
    aParB  =  source.aParB;
    aParC  =  source.aParC;
    aNR    =  source.aNR;
    aMinR  =  source.aMinR;
    aMaxR  =  source.aMaxR;
    nProtonsB       =  source.nProtonsB;
    nNeutronsB      =  source.nNeutronsB;
    bGeneratorType  =  source.bGeneratorType;
    bParA  =  source.bParA;
    bParB  =  source.bParB;
    bParC  =  source.bParC;
    bNR    =  source.bNR;
    bMinR  =  source.bMinR;
    bMaxR  =  source.bMaxR;
    minB   =  source.minB;
    maxB   =  source.maxB;
    nnCrossSection  =  source.nnCrossSection;
    nBins_b  =  source.nBins_b;
    min_b    =  source.min_b;
    max_b    =  source.max_b;
    nBins_nPart  =  source.nBins_nPart;
    min_nPart    =  source.min_nPart;
    max_nPart    =  source.max_nPart;
    nBins_nBinary  =  source.nBins_nBinary;
    min_nBinary    =  source.min_nBinary;
    max_nBinary    =  source.max_nBinary;
    setBXsectPartition(source.nBins_bxSect, source.bValues);
    }
  return *this;
}

CollisionGeometryConfiguration::~CollisionGeometryConfiguration()
{
  if (!bValues) delete bValues;
}

void CollisionGeometryConfiguration::setBXsectPartition(int nBins, double * sourceValues)
{
  // sanity check, do not change what we own already
  if (this->bValues == sourceValues) return;
  if (!bValues) delete bValues;
  nBins_bxSect = nBins;
  bValues = new double[nBins_bxSect+1];
  for (int iBin=0; iBin<nBins_bxSect; iBin++)
  {
  bValues[iBin] = sourceValues[iBin];
  }
}

void CollisionGeometryConfiguration::setDefaultBXsectPartition()
{
  // following is for PbPb at 2.74 TeV
  // but needs to updated for calculations
  // at other energies and for other nuclei.
  // You can overide this function in a subclass to
  // specificy values appropriate for another nucleus.
  setPbPb2740GeVBXsectPartition();
}

void CollisionGeometryConfiguration::setPbPb2740GeVBXsectPartition()
{
  // following is for PbPb at 2.74 TeV
  if (!bValues) delete bValues;
  nBins_bxSect = 11;
  bValues = new double[nBins_bxSect+1];
  bValues[0] = 0.0;
  bValues[1] = 3.39927;
  bValues[2] = 4.79929;
  bValues[3] = 6.79951;
  bValues[4] = 8.39898;
  bValues[5] = 9.79733;
  bValues[6] = 10.9965;
  bValues[7] = 11.9972;
  bValues[8] = 12.9963;
  bValues[9] = 13.7987;
  bValues[10] = 14.7986;
  bValues[11] = 20.0;
}

// Set the default configuration as PbPb at 2.74 TeV
void CollisionGeometryConfiguration::setDefaultConfiguration()
{
  setPbPb2740GeVConfiguration();
}

void CollisionGeometryConfiguration::setPbPb2740GeVConfiguration()
{
  nProtonsA         = 126;
  nNeutronsA        = 208-82;
  nProtonsB         = 82;
  nNeutronsB        = 126;
  outputPath        = getenv("WAC_INPUT_PATH");
  outputPath        = getenv("WAC_OUTPUT_PATH");
  rootOuputFileName =  "/CollisionGeometry";
  histoBaseName     =  "geom";
  minB = 0.0;
  maxB = 21.0;

  aParA = 6.62;
  aParB = 0.546;
  aParC = 0.0;
  aNR   = 11000;
  aMinR = 0.0;
  aMaxR = 11.0;
  bParA = 6.62;
  bParB = 0.546;
  bParC = 0.0;
  bNR   = 11000;
  bMinR = 0.0;
  bMaxR = 11.0;
  nnCrossSection = 6.4;
  nBins_b = 120;
  min_b   = 0.0;
  max_b   = 24.0;
  nBins_nPart = 100;
  min_nPart   = 0;
  max_nPart   = 500;
  nBins_nBinary = 600;
  min_nBinary   = 0;
  max_nBinary   = 3000;
  calculateDerivedHistograms = true;
  setPbPb2740GeVBXsectPartition();
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
  << "      Histogram Base Name : " << histoBaseName << endl
  << "                nProtonsA : " << nProtonsA << endl
  << "               nNeutronsA : " << nNeutronsA << endl
  << "           aGeneratorType : " << aGeneratorType << endl
  << "                    aParA : " << aParA << endl
  << "                    aParB : " << aParB << endl
  << "                    aParC : " << aParC << endl
  << "                      aNR : " << aNR << endl
  << "                    aMinR : " << aMinR << endl
  << "                    aMaxR : " << aMaxR << endl
  << "                nProtonsB : " << nProtonsB << endl
  << "               nNeutronsB : " << nNeutronsB << endl
  << "           bGeneratorType : " << bGeneratorType << endl
  << "                    bParA : " << bParA << endl
  << "                    bParB : " << bParB << endl
  << "                    bParC : " << bParC << endl
  << "                      bNR : " << bNR << endl
  << "                    bMinR : " << bMinR << endl
  << "                    bMaxR : " << bMaxR << endl
  << "                     minB : " << minB << endl
  << "                     maxB : " << minB << endl
  << "           nnCrossSection : " << nnCrossSection << endl
  << "                  nBins_b : " << nBins_b << endl
  << "                    min_b : " << min_b << endl
  << "                    max_b : " << max_b << endl
  << "              nBins_nPart : " << nBins_nPart << endl
  << "                min_nPart : " << min_nPart << endl
  << "                max_nPart : " << max_nPart << endl
  << "            nBins_nBinary : " << nBins_nBinary << endl
  << "              min_nBinary : " << min_nBinary << endl
  << "              max_nBinary : " << max_nBinary << endl
  << "             nBins_bxSect : " << nBins_bxSect << endl;
  for (int iBin=0; iBin<nBins_bxSect+1; iBin++)
  {
  os << "                     iBin : " << iBin << "   b: " <<  bValues[iBin] << " fm"<< endl;
  }
}
