// Adapted to root envinromemt 2021 C. Pruneau
//
// Copyright 2016 Chun Shen
#ifndef WAC_HadronGasConfiguration
#define WAC_HadronGasConfiguration
#include "TaskConfiguration.hpp"

using namespace std;

class HadronGasConfiguration : public TaskConfiguration
{
public:

  HadronGasConfiguration(const TString & name,
                         const TString & type,
                         const TString & version);
  HadronGasConfiguration(const HadronGasConfiguration & config);
  virtual ~HadronGasConfiguration();
  HadronGasConfiguration & operator=(const HadronGasConfiguration & config);

  ostream & printConfiguration(ostream & os);

  bool plotSystemProperties;
  int nT;    double minT,    maxT;
  int nMass; double minMass, maxMass;
  int nMu;   double minMu,   maxMu;
  int nMuB;  double minMuB,  maxMuB;
  int nMuS;  double minMuS,  maxMuS;
  int nMuQ;  double minMuQ,  maxMuQ;

  bool plotParticleProperties;
  int nTypes;
  int nStableTypes;

  ClassDef(HadronGasConfiguration,0)
};


#endif  // WAC_HadronGasConfiguration

