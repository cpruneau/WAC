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

 Utility class used to defineAnalysisConfiguration
 */

#include "AnalysisConfiguration.hpp"

ClassImp(AnalysisConfiguration);


AnalysisConfiguration::AnalysisConfiguration(const TString & name,
                                             const TString & type,
                                             const TString & version)
:
TaskConfiguration(name,type,version),
histoBaseName(""),
nBins_nPartTot(10),     min_nPartTot(0.0),  max_nPartTot(100.0),
nBins_nCell(100),        min_nCell(0.0),     max_nCell(100.0),
nBins_nPartPerCell(10), min_nPartPerCell(0.0), max_nPartPerCell(10.0),
nBins_pt(40),           min_pt(0.0),    max_pt(2.0),   range_pt(2.0),
nBins_eta(20),          min_eta(-2.0),  max_eta(2.0),  range_eta(4.0),
nBins_y(20),            min_y(-2.0),    max_y(-2.0),   range_y(4.0),
nBins_phi(36),          min_phi(0.0),   max_phi(0.0),  range_phi(0.0),
nBins_phiEta(720),
nBins_phiEtaPt(7200),
nBins_phiY(720),
nBins_phiYPt(7200),
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
fill3D(false),
fill6D(false),
fillQ3D(false),
fillY(true),
nuDynVsMult(true),
nuDynVsCent(false),
ptCorrelatorVsMult(true),
ptCorrelatorVsCent(false),
nBins_mult(20),  min_mult(0.0), max_mult(200.0),
nBins_cent(20),  min_cent(0.0), max_cent(100.0)
{
  max_phi   = TMath::TwoPi();
  range_phi = TMath::TwoPi();
  max_Dphi  = max_phi;
  width_Dphi = max_phi/double(nBins_Dphi);
  min_Dphi_shft = -TMath::Pi()/4.0;
  max_Dphi_shft = 3.0*TMath::Pi()/4.0;
}

////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void AnalysisConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  os
  << "    Analysis   Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "      Histogram Base Name: " << histoBaseName << endl
  << "                  nBinsPt: " << nBins_pt    << endl
  << "                    minPt: " << min_pt      << endl
  << "                    maxPt: " << max_pt      << endl
  << "                 nBinsEta: " << nBins_eta   << endl
  << "                   minEta: " << min_eta     << endl
  << "                   maxEta: " << max_eta     << endl
  << "                     minY: " << min_y     << endl
  << "                     maxY: " << max_y     << endl
  << "                 nBinsPhi: " << nBins_phi   << endl
  << "                   minPhi: " << min_phi     << endl
  << "                   maxPhi: " << max_phi     << endl
  << "                fillPairs: " << fillPairs        << endl
  << "                   fill3D: " << fill3D           << endl
  << "                   fill6D: " << fill6D           << endl
  << "                  fillQ3D: " << fillQ3D          << endl
  << "                    fillY: " << fillY            << endl
  << "        +- bin correction: " << binCorrPM        << endl
  << "        ++ bin correction: " << binCorrPP        << endl
  << "        -- bin correction: " << binCorrMM        << endl;
}

int AnalysisConfiguration::getIxEtaPhi(double eta, double phi) {

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

int AnalysisConfiguration::getIxYPhi(double y, double phi) {

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
