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

 Utility class used to defineParticleAnalyzerConfiguration
 */

#include "ParticleAnalyzerConfiguration.hpp"
#include "TMath.h"

ClassImp(ParticleAnalyzerConfiguration);

ParticleAnalyzerConfiguration::ParticleAnalyzerConfiguration(const TString & name,
                                                             const TString & type,
                                                             const TString & version)
:
TaskConfiguration(name,type,version),
nBins_pt(40),           min_pt(0.0),    max_pt(2.0),   range_pt(2.0),
nBins_eta(20),          min_eta(-2.0),  max_eta(2.0),  range_eta(4.0),
nBins_y(20),            min_y(-2.0),    max_y(-2.0),   range_y(4.0),
nBins_phi(36),          min_phi(0.0),   max_phi(),     range_phi(0.0),
nBins_phiEta(720),
nBins_phiEtaPt(7200),
nBins_phiY(720),
nBins_phiYPt(7200),
fill3D(false),
fillY(true)
{
  min_phi   = 0.0;
  max_phi   = TMath::TwoPi();
  validate();
}

ParticleAnalyzerConfiguration::ParticleAnalyzerConfiguration(const ParticleAnalyzerConfiguration & source)
:
TaskConfiguration( source ),
nBins_pt( source.nBins_pt ),           min_pt( source.min_pt ),    max_pt( source.max_pt ),   range_pt( source.range_pt ),
nBins_eta( source.nBins_eta ),          min_eta( source.min_eta ),  max_eta( source.max_eta ),  range_eta( source.range_eta ),
nBins_y( source.nBins_y ),            min_y( source.min_y ),    max_y( source.max_y ),   range_y( source.range_y ),
nBins_phi( source.nBins_phi ),          min_phi( source.min_phi ),   max_phi( source.max_phi ),  range_phi( source.range_phi ),
nBins_phiEta( source.nBins_phiEta ),
nBins_phiEtaPt( source.nBins_phiEtaPt ),
nBins_phiY( source.nBins_phiY ),
nBins_phiYPt( source.nBins_phiYPt ),
fill3D( source.fill3D ),
fillY( source.fillY )
{
}

ParticleAnalyzerConfiguration & ParticleAnalyzerConfiguration::operator=(const ParticleAnalyzerConfiguration & source)
{
  if (this!=&source)
    {
    TaskConfiguration::operator=( source );
    nBins_pt  =  source.nBins_pt;
    min_pt  =  source.min_pt;
    max_pt  =  source.max_pt;
    range_pt  =  source.range_pt;
    nBins_eta  =  source.nBins_eta;
    min_eta  =  source.min_eta;
    max_eta  =  source.max_eta;
    range_eta  =  source.range_eta;
    nBins_y  =  source.nBins_y;
    min_y  =  source.min_y;
    max_y  =  source.max_y;
    range_y  =  source.range_y;
    nBins_phi  =  source.nBins_phi;
    min_phi  =  source.min_phi;
    max_phi  =  source.max_phi;
    range_phi  =  source.range_phi;
    nBins_phiEta  =  source.nBins_phiEta;
    nBins_phiEtaPt  =  source.nBins_phiEtaPt;
    nBins_phiY  =  source.nBins_phiY;
    nBins_phiYPt  =  source.nBins_phiYPt;
    fill3D  =  source.fill3D;
    fillY  =  source.fillY;
    }
  return *this;
}

// makes sure the difference variables are correct
void ParticleAnalyzerConfiguration::validate()
{
  range_pt   = max_pt  - min_pt;
  range_y    = max_y   - min_y;
  range_eta  = max_eta - min_eta;
  range_phi  = max_phi - min_phi;
}

////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void ParticleAnalyzerConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  os
  << "    Analysis   Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                  nBinsPt: " << nBins_pt    << endl
  << "                    minPt: " << min_pt      << endl
  << "                    maxPt: " << max_pt      << endl
  << "                 nBinsEta: " << nBins_eta   << endl
  << "                   minEta: " << min_eta     << endl
  << "                   maxEta: " << max_eta     << endl
  << "                     minY: " << min_y       << endl
  << "                     maxY: " << max_y       << endl
  << "                 nBinsPhi: " << nBins_phi   << endl
  << "                   minPhi: " << min_phi     << endl
  << "                   maxPhi: " << max_phi     << endl
  << "                   fill3D: " << fill3D      << endl
  << "                    fillY: " << fillY       << endl;
}

int ParticleAnalyzerConfiguration::getIxEtaPhi(double eta, double phi) {

  if(!(eta<min_eta || eta>max_eta)) {
    int iEta = int( double(nBins_eta)*(eta-min_eta)/range_eta );
    int iPhi = int( double(nBins_phi)*(phi-min_phi)/range_phi );
    if (iEta>=0 && iPhi>=0 && iEta<nBins_eta && iPhi<nBins_phi)
      {
      return nBins_phi*iEta + iPhi;
      }
    return -1;
  }
  return -1;
}

int ParticleAnalyzerConfiguration::getIxYPhi(double y, double phi) {

  if(!(y<min_y || y>max_y))
    {
    int iY = int( double(nBins_y)*(y-min_y)/range_y );
    int iPhi = int( double(nBins_phi)*(phi-min_phi)/range_phi );

    if (iY>=0 && iPhi>=0 && iY<nBins_y && iPhi<nBins_phi)
      {
      return nBins_phi*iY + iPhi;
      }
    return -1;
    }
  return -1;
}
