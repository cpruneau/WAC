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

 Utility class used to defineParticlePairAnalyzerConfiguration
 */

#include "ParticlePairAnalyzerConfiguration.hpp"
#include "TMath.h"

ClassImp(ParticlePairAnalyzerConfiguration);


ParticlePairAnalyzerConfiguration::ParticlePairAnalyzerConfiguration(const TString & name,
                                             const TString & type,
                                             const TString & version)
:
ParticleAnalyzerConfiguration(name,type,version),
binCorrPM(1.0),
binCorrMM(1.0),
binCorrPP(1.0),
nBins_DeltaPlong(10),   min_DeltaPlong(-1.0), max_DeltaPlong(1.0),  range_DeltaPlong(2.0),
nBins_DeltaPside(10),   min_DeltaPside(-1.0), max_DeltaPside(1.0),  range_DeltaPside(2.0),
nBins_DeltaPout(10),    min_DeltaPout(-1.0),  max_DeltaPout(1.0),   range_DeltaPout(2.0),
nBins_Dphi(36),         min_Dphi(0.0),        max_Dphi(0.0),        width_Dphi(0.0),
nBins_Dphi_shft(36),    min_Dphi_shft(0.0),   max_Dphi_shft(0.0),
nBins_Deta(39),         min_Deta(-2.0),       max_Deta(2.0),        width_Deta(4.0/39.0),
nBins_Dy(39),           min_Dy(-2.0),         max_Dy(2.0),          width_Dy(4.0/39.0),
fillPairs(true),
fill6D(false)
{
  validate();
}

// makes sure the difference variables are correct
void ParticlePairAnalyzerConfiguration::validate()
{
  ParticleAnalyzerConfiguration::validate();

  nBins_Deta = 2*nBins_eta-1;
  min_Deta   = min_eta - max_eta;
  max_Deta   = max_eta - min_eta;
  width_Deta = (max_Deta-min_Deta)/double(nBins_Deta);

  nBins_Dy = 2*nBins_y-1;
  min_Dy   = min_y - max_y;
  max_Dy   = max_y - min_y;
  width_Dy = (max_Dy-min_Dy)/double(nBins_Dy);


  max_Dphi   = TMath::TwoPi();
  width_Dphi = max_phi/double(nBins_Dphi);
  int nBin = 0.25 * double(nBins_Dphi);
  min_Dphi_shft = - TMath::TwoPi() * double(nBin)/double(nBins_Dphi);
  max_Dphi_shft =   TMath::TwoPi() * double(nBins_Dphi-nBin)/double(nBins_Dphi);
}

ParticlePairAnalyzerConfiguration::ParticlePairAnalyzerConfiguration(const ParticlePairAnalyzerConfiguration & source)
:
ParticleAnalyzerConfiguration( source ),
binCorrPM( source.binCorrPM ),
binCorrMM( source.binCorrMM ),
binCorrPP( source.binCorrPP ),
nBins_DeltaPlong( source.nBins_DeltaPlong ),   min_DeltaPlong( source.min_DeltaPlong ), max_DeltaPlong( source.max_DeltaPlong ),  range_DeltaPlong( source.range_DeltaPlong ),
nBins_DeltaPside( source.nBins_DeltaPside ),   min_DeltaPside( source.min_DeltaPside ), max_DeltaPside( source.max_DeltaPside ),  range_DeltaPside( source.range_DeltaPside ),
nBins_DeltaPout( source.nBins_DeltaPout ),    min_DeltaPout( source.min_DeltaPout ),  max_DeltaPout( source.max_DeltaPout ),   range_DeltaPout( source.range_DeltaPout ),
nBins_Dphi( source.nBins_Dphi ),         min_Dphi( source.min_Dphi ),        max_Dphi( source.max_Dphi ),        width_Dphi( source.width_Dphi ),
nBins_Dphi_shft( source.nBins_Dphi_shft ),    min_Dphi_shft( source.min_Dphi_shft ),   max_Dphi_shft( source.max_Dphi_shft ),
nBins_Deta( source.nBins_Deta ),         min_Deta( source.min_Deta ),       max_Deta( source.max_Deta ),        width_Deta( source.width_Deta ),
nBins_Dy( source.nBins_Dy ),           min_Dy( source.min_Dy ),         max_Dy( source.max_Dy ),          width_Dy( source.width_Dy ),
fillPairs( source.fillPairs ),
fill6D( source.fill6D )
{
  // no ops
}

ParticlePairAnalyzerConfiguration & ParticlePairAnalyzerConfiguration::operator=(const ParticlePairAnalyzerConfiguration & source)
{
  if (this!=&source)
    {
    ParticleAnalyzerConfiguration::operator=( source );
    binCorrPM  =  source.binCorrPM;
    binCorrMM  =  source.binCorrMM;
    binCorrPP  =  source.binCorrPP;
    nBins_DeltaPlong  =  source.nBins_DeltaPlong;
    min_DeltaPlong  =  source.min_DeltaPlong;
    max_DeltaPlong  =  source.max_DeltaPlong;
    range_DeltaPlong  =  source.range_DeltaPlong;
    nBins_DeltaPside  =  source.nBins_DeltaPside;
    min_DeltaPside  =  source.min_DeltaPside;
    max_DeltaPside  =  source.max_DeltaPside;
    range_DeltaPside  =  source.range_DeltaPside;
    nBins_DeltaPout  =  source.nBins_DeltaPout;
    min_DeltaPout  =  source.min_DeltaPout;
    max_DeltaPout  =  source.max_DeltaPout;
    range_DeltaPout  =  source.range_DeltaPout;
    nBins_Dphi  =  source.nBins_Dphi;
    min_Dphi  =  source.min_Dphi;
    max_Dphi  =  source.max_Dphi;
    width_Dphi  =  source.width_Dphi;
    nBins_Dphi_shft  =  source.nBins_Dphi_shft;
    min_Dphi_shft  =  source.min_Dphi_shft;
    max_Dphi_shft  =  source.max_Dphi_shft;
    nBins_Deta  =  source.nBins_Deta;
    min_Deta  =  source.min_Deta;
    max_Deta  =  source.max_Deta;
    width_Deta  =  source.width_Deta;
    nBins_Dy  =  source.nBins_Dy;
    min_Dy  =  source.min_Dy;
    max_Dy  =  source.max_Dy;
    width_Dy  =  source.width_Dy;
    fillPairs  =  source.fillPairs;
    fill6D  =  source.fill6D;
    fillQ3D  =  source.fillQ3D;
    }
  return *this;
}



////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void ParticlePairAnalyzerConfiguration::printConfiguration(ostream & os)
{
  ParticleAnalyzerConfiguration::printConfiguration(os);
  os
  << " ------------------------------------------------------------------------------------------" << endl
  << "                fillPairs: " << fillPairs        << endl
  << "                   fill6D: " << fill6D           << endl
  << "                  fillQ3D: " << fillQ3D          << endl
  << "        +- bin correction: " << binCorrPM        << endl
  << "        ++ bin correction: " << binCorrPP        << endl
  << "        -- bin correction: " << binCorrMM        << endl;
}
