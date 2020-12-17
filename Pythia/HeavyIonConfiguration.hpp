#ifndef WAC_HeavyIonConfiguration
#define WAC_HeavyIonConfiguration
#include "NuDynConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class HeavyIonConfiguration : public NuDynConfiguration
{
public:
  
  HeavyIonConfiguration(const TString & name,
                        const TString & type,
                        const TString & version);
  virtual ~HeavyIonConfiguration(){}
  void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  bool hardBoost; //corresponds to a hard boost scenario. If not hard boost the soft boost

  int nCollisionsMax; //holds the max number of binary collisions of an event

  double param_a; double param_b; //parameters in deciding the boost of particles (beta = a * r^b)

  int maxOrder; //order of correlation functions

  int totEvents; //total number of events


  // event wise parameters
  

ClassDef(HeavyIonConfiguration,0)
};

#endif /* WAC_AnalysisConfiguration */
