//
//  BWModelConfiguration.hpp
//  MyMC
//
//  Created by Victor Gonzalez on 20190206 based in the work of Claude Pruneau and the ALICE AliPhysics SW framework
//  Copyright © 2017-2019 Claude Pruneau. All rights reserved.
//

#ifndef BWModelConfiguration_hpp
#define BWModelConfiguration_hpp
#include "TaskConfiguration.hpp"
#include <TF1.h>

class BWModelConfiguration : public TaskConfiguration
{
public:

  BWModelConfiguration(const TString & name, const TString & type, const TString & version);
  ~BWModelConfiguration();
  virtual void printConfiguration(ostream & os);

  //////////////////////////////////////////////
  //  event generation parameters
  //////////////////////////////////////////////
  double  totalMultiplicityMean;
  double  totalMultiplicitySigma;
  double  netChargeMean;
  double  netChargeSigma;


  //////////////////////////////////////////////
  //  particle generation parameters
  //////////////////////////////////////////////
  
  double particleEta_Min;
  double particleEta_Max;
  double pT_Min;
  double pT_Max;

  bool   useAllKinds;
  bool   useBGBlastWave;
  double betaAvgBGBW;
  double betaProfileBGBW;
  double longitudinalFlow;
  bool   useFlow;
  
  // all kind of particles
  double temperatureAllKinds;
  double directedFlowAllKinds;
  double ellipticFlowAllKinds;
  double triangularFlowAllKinds;
  double quadrangularFlowAllKinds;
  double pentangularFlowAllKinds;

  // pions
  double pionPercentage;
  double temperaturePions;
  double directedFlowPions;
  double ellipticFlowPions;
  double triangularFlowPions;
  double quadrangularFlowPions;
  double pentangularFlowPions;
  
  // kaons
  double kaonPercentage;
  double temperatureKaons;
  double directedFlowKaons;
  double ellipticFlowKaons;
  double triangularFlowKaons;
  double quadrangularFlowKaons;
  double pentangularFlowKaons;

  // protons
  double protonPercentage;
  double temperatureProtons;
  double directedFlowProtons;
  double ellipticFlowProtons;
  double triangularFlowProtons;
  double quadrangularFlowProtons;
  double pentangularFlowProtons;

    ClassDef(BWModelConfiguration,0)
};


#endif /* BWModelConfiguration_hpp */
