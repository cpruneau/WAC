//
//  BWModelConfiguration.hpp
//  Created by Victor Gonzalez on 20190206 based in the work of Claude Pruneau and the ALICE AliPhysics SW framework
//  Copyright © 2017-2019 Claude Pruneau. All rights reserved.
//

#include <TF1.h>
#include "BWModelConfiguration.hpp"
ClassImp(BWModelConfiguration);

BWModelConfiguration::BWModelConfiguration(const TString & name, const TString & type, const TString & version)
:
TaskConfiguration(name,type, version),
totalMultiplicityMean(2000.0),
totalMultiplicitySigma(20.0),
netChargeMean(5.0),
netChargeSigma(4.0),
particleEta_Min(-3.0),
particleEta_Max(3.0),
pT_Min(0.2),
pT_Max(2.0),
useAllKinds(true),
useBGBlastWave(false),
betaAvgBGBW(0.0),
betaProfileBGBW(0.0),
longitudinalFlow(0.0),
useFlow(false),
temperatureAllKinds(),
directedFlowAllKinds(0.0),
ellipticFlowAllKinds(0.0),
triangularFlowAllKinds(0.0),
quadrangularFlowAllKinds(0.0),
pentangularFlowAllKinds(0.0),
temperaturePions(),
pionPercentage(),
directedFlowPions(0.0),
ellipticFlowPions(0.0),
triangularFlowPions(0.0),
quadrangularFlowPions(0.0),
pentangularFlowPions(0.0),
temperatureKaons(),
kaonPercentage(),
directedFlowKaons(0.0),
ellipticFlowKaons(0.0),
triangularFlowKaons(0.0),
quadrangularFlowKaons(0.0),
pentangularFlowKaons(0.0),
temperatureProtons(),
protonPercentage(),
directedFlowProtons(0.0),
ellipticFlowProtons(0.0),
triangularFlowProtons(0.0),
quadrangularFlowProtons(0.0),
pentangularFlowProtons(0.0)
{

}

BWModelConfiguration::~BWModelConfiguration()
{}

void BWModelConfiguration::printConfiguration(ostream & os)
{
  os
  << " ------------------------------------------------------------------------------------------" << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                           modelName: " << getName() << endl
  << "                        modelVersion: " << getVersion() << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "         Event Generation Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "               totalMultiplicityMean: " << totalMultiplicityMean   << endl
  << "              totalMultiplicitySigma: " << totalMultiplicitySigma  << endl
  << "                       netChargeMean: " << netChargeMean           << endl
  << "                      netChargeSigma: " << netChargeSigma          << endl
  << " ------------------------------------------------------------------------------------------" << endl;
  os
  << "      Particle Generation Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                     particleEta_Min: " << particleEta_Min  << endl
  << "                     particleEta_Max: " << particleEta_Max  << endl
  << "                              pT_Min: " << pT_Min           << endl
  << "                              pT_Max: " << pT_Max           << endl << endl

  << "                         useAllKinds: " << useAllKinds      << endl
  << "                      useBGBlastWave: " << useBGBlastWave   << endl
  << "                         betaAvgBGBW: " << betaAvgBGBW      << endl
  << "                     betaProfileBGBW: " << betaProfileBGBW  << endl
  << "                    longitudinalFlow: " << longitudinalFlow << endl
  << "                             useFlow: " << useFlow          << endl
  << " ------------------------------------------------------------------------------------------" << endl;
  os
  << "                           All kinds: " << endl << endl

  << "                 temperatureAllKinds: " << temperatureAllKinds        << endl
  << "                directedFlowAllKinds: " << directedFlowAllKinds       << endl
  << "                ellipticFlowAllKinds: " << ellipticFlowAllKinds       << endl
  << "              triangularFlowAllKinds: " << triangularFlowAllKinds     << endl
  << "            quadrangularFlowAllKinds: " << quadrangularFlowAllKinds   << endl
  << "             pentangularFlowAllKinds: " << pentangularFlowAllKinds    << endl << endl;
  os
  << "                               Pions: " << endl << endl

  << "                      pionPercentage: " << pionPercentage          << endl
  << "                    temperaturePions: " << temperaturePions        << endl
  << "                   directedFlowPions: " << directedFlowPions       << endl
  << "                   ellipticFlowPions: " << ellipticFlowPions       << endl
  << "                 triangularFlowPions: " << triangularFlowPions     << endl
  << "               quadrangularFlowPions: " << quadrangularFlowPions   << endl
  << "                pentangularFlowPions: " << pentangularFlowPions    << endl << endl;
  os
  << "                               Kaons: " << endl << endl

  << "                      kaonPercentage: " << kaonPercentage          << endl
  << "                    temperatureKaons: " << temperatureKaons        << endl
  << "                   directedFlowKaons: " << directedFlowKaons       << endl
  << "                   ellipticFlowKaons: " << ellipticFlowKaons       << endl
  << "                 triangularFlowKaons: " << triangularFlowKaons     << endl
  << "               quadrangularFlowKaons: " << quadrangularFlowKaons   << endl
  << "                pentangularFlowKaons: " << pentangularFlowKaons    << endl << endl;
  os
  << "                             Protons: " << endl << endl

  << "                    protonPercentage: " << protonPercentage          << endl
  << "                  temperatureProtons: " << temperatureProtons        << endl
  << "                 directedFlowProtons: " << directedFlowProtons       << endl
  << "                 ellipticFlowProtons: " << ellipticFlowProtons       << endl
  << "               triangularFlowProtons: " << triangularFlowProtons     << endl
  << "             quadrangularFlowProtons: " << quadrangularFlowProtons   << endl
  << "              pentangularFlowProtons: " << pentangularFlowProtons    << endl;
  os
  << " ------------------------------------------------------------------------------------------" << endl
  << " ------------------------------------------------------------------------------------------" << endl;
}


