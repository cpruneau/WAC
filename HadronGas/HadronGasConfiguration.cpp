// Copyright 2016 Chun Shen

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "HadronGasConfiguration.hpp"

ClassImp(HadronGasConfiguration);

using namespace std;

HadronGasConfiguration::HadronGasConfiguration(const TString & name,
                                               const TString & type,
                                               const TString & version)
:
TaskConfiguration(name,type,version),
plotSystemProperties(true),
nT(50),     minT(100.0),  maxT(200.0),
nMu(20),    minMu(0.0),   maxMu(40),
nMuB(20),   minMuB(0.0),  maxMuB(20),
nMuS(20),   minMuS(0.0),  maxMuS(20),
nMuQ(20),   minMuQ(0.0),  maxMuQ(20),
plotParticleProperties(false),
nTypes(319),
nStableTypes(23)
{
  // no ops
}

HadronGasConfiguration::HadronGasConfiguration(const HadronGasConfiguration & source)
:
plotSystemProperties(source.plotSystemProperties),
nT(source.nT),       minT(source.minT),      maxT(source.maxT),
nMu(source.nMu),     minMu(source.minMu),    maxMu(source.maxMu),
nMuB(source.nMuB),   minMuB(source.minMuB),  maxMuB(source.maxMuB),
nMuS(source.nMuS),   minMuS(source.minMuS),  maxMuS(source.maxMuS),
nMuQ(source.nMuQ),   minMuQ(source.minMuQ),  maxMuQ(source.maxMuQ),
plotParticleProperties(source.plotParticleProperties),
nTypes(source.nTypes),
nStableTypes(source.nStableTypes)
{
  // no ops
}

HadronGasConfiguration::~HadronGasConfiguration()
{
  // no ops
}

HadronGasConfiguration & HadronGasConfiguration::operator=(const HadronGasConfiguration & source)
{
  if (this!=&source)
    {
    plotSystemProperties = source.plotSystemProperties;
    nT   = source.nT;     minT   = source.minT;   maxT   = source.maxT;
    nMu  = source.nMu;    minMu  = source.minMu;  maxMu  = source.maxMu;
    nMuB = source.nMuB;   minMuB = source.minMuB; maxMuB = source.maxMuB;
    nMuS = source.nMuS;   minMuS = source.minMuS; maxMuS = source.maxMuS;
    nMuQ = source.nMuQ;   minMuQ = source.minMuQ; maxMuQ = source.maxMuQ;
    plotParticleProperties = source.plotParticleProperties;
    nTypes       = source.nTypes;
    nStableTypes = source.nStableTypes;
    }
  return *this;
}

