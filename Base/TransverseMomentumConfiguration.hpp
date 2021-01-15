#ifndef WAC_TransverseMomentumConfiguration
#define WAC_TransverseMomentumConfiguration
#include "TMath.h"
#include "TaskConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class TransverseMomentumConfiguration : public TaskConfiguration
{
public:
  TransverseMomentumConfiguration(const TString & name,
                        const TString & type,
                        const TString & version);
  TransverseMomentumConfiguration(const TransverseMomentumConfiguration & source);
  virtual ~TransverseMomentumConfiguration();
  TransverseMomentumConfiguration & operator=(const TransverseMomentumConfiguration & source);
  void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////

  int nBins_mult;  double min_mult; double max_mult;
  int nBins_cent;  double min_cent; double max_cent;
    bool ptCorrelatorVsMult; bool ptCorrelatorVsCent;
    int totEvents;
    int maxOrder; //order of correlation functions
    int numTypes; //number of particle types to consider >= maxOrder

ClassDef(TransverseMomentumConfiguration,0)
};

#endif /* WAC_TransverseMomentumConfiguration */