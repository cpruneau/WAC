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

 Utility class used to defineGlobalAnalyzerConfiguration
 */

#include "GlobalAnalyzerConfiguration.hpp"
#include "TMath.h"

ClassImp(GlobalAnalyzerConfiguration);

GlobalAnalyzerConfiguration::GlobalAnalyzerConfiguration(const TString & name,
                                             const TString & type,
                                             const TString & version)
:
TaskConfiguration(name,type,version),
nBins_n(100), nBins_n2(20), min_n(0.0), max_n(1000.0), range_n(1000.0),
nBins_e(100), nBins_e2(20), min_e(0.0), max_e(1000.0), range_e(1000.0),
nBins_q(100), nBins_q2(20), min_q(0.0), max_q(1000.0), range_q(1000.0),
nBins_b(100), nBins_b2(20), min_b(0.0), max_b(1000.0), range_b(1000.0),
countParticles(true),
setEvent(true),
fillCorrelationHistos(true)
{
}

GlobalAnalyzerConfiguration::GlobalAnalyzerConfiguration(const GlobalAnalyzerConfiguration & source)
:
TaskConfiguration( source ),
nBins_n(source.nBins_n), nBins_n2(source.nBins_n2), min_n(source.min_n), max_n(source.max_n), range_n(source.range_n),
nBins_e(source.nBins_e), nBins_e2(source.nBins_e2), min_e(source.min_e), max_e(source.max_e), range_e(source.range_e),
nBins_q(source.nBins_q), nBins_q2(source.nBins_q2), min_q(source.min_q), max_q(source.max_q), range_q(source.range_q),
nBins_b(source.nBins_b), nBins_b2(source.nBins_b2), min_b(source.min_b), max_b(source.max_b), range_b(source.range_b),
countParticles(source.countParticles),
setEvent(source.setEvent),
fillCorrelationHistos(source.fillCorrelationHistos)
{
}

GlobalAnalyzerConfiguration & GlobalAnalyzerConfiguration::operator=(const GlobalAnalyzerConfiguration & source)
{
  if (this!=&source)
    {
    TaskConfiguration::operator=( source );
    nBins_n = source.nBins_n;  nBins_n2 = source.nBins_n2;  min_n = source.min_n;  max_n = source.max_n;  range_n = source.range_n;
    nBins_e = source.nBins_e;  nBins_e2 = source.nBins_e2;  min_e = source.min_e;  max_e = source.max_e;  range_e = source.range_e;
    nBins_q = source.nBins_q;  nBins_q2 = source.nBins_q2;  min_q = source.min_q;  max_q = source.max_q;  range_q = source.range_q;
    nBins_b = source.nBins_b;  nBins_b2 = source.nBins_b2;  min_b = source.min_b;  max_b = source.max_b;  range_b = source.range_b;
    fillCorrelationHistos = source.fillCorrelationHistos;
    countParticles = source.countParticles;
    setEvent       = source.setEvent;
    }
  return *this;
}

// makes sure the difference variables are correct
void GlobalAnalyzerConfiguration::validate()
{
  range_n = max_n - min_n;
  range_e = max_e - min_e;
  range_q = max_q - min_q;
  range_b = max_b - min_b;
}

////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void GlobalAnalyzerConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  os
  << "    Analysis   Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "    nBins_n: " << nBins_n << " nBins_n2: " << nBins_n2 << " min_n: " << min_n << " max_n: " << max_n << " range_n: " << range_n << endl
  << "    nBins_e: " << nBins_e << " nBins_e2: " << nBins_e2 << " min_e: " << min_e << " max_e: " << max_e << " range_e: " << range_e << endl
  << "    nBins_q: " << nBins_q << " nBins_q2: " << nBins_q2 << " min_q: " << min_q << " max_q: " << max_q << " range_q: " << range_q << endl
  << "    nBins_b: " << nBins_b << " nBins_b2: " << nBins_b2 << " min_b: " << min_b << " max_b: " << max_b << " range_b: " << range_b << endl
  << "    fillCorrelationHistos: " << fillCorrelationHistos << endl
  << "           countParticles: " << countParticles << endl
  << "                 setEvent: " << setEvent << endl;
}
