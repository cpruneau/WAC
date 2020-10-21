#include "HeavyIonConfiguration.hpp"

ClassImp(HeavyIonConfiguration);


HeavyIonConfiguration::HeavyIonConfiguration(const TString & name,
                                             const TString & type,
                                             const TString & version)
:
AnalysisConfiguration(name,type,version),
hardBoost(true),
nCollisionsMax(1),
param_a(0.5), param_b(1),
maxOrder(0), totEvents(0)
{

}

////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void HeavyIonConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  AnalysisConfiguration::printConfiguration(os);
  os
  << "    Heavy Ion Collision Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                hardBoost: " << hardBoost   << endl
  << "                softBoost: " << !hardBoost  << endl
  << "           nCollisionsMax: " << nCollisionsMax      << endl
  << "                  param_a: " << param_a  << endl
  << "                  param_b: " << param_b  << endl
  << "                 maxOrder: " << maxOrder  << endl
  << "                totEvents: " << totEvents  << endl;
}
