#ifndef WAC_PTCorrelator
#define WAC_PTCorrelator
#include <TParameter.h>
#include "TFile.h"
#include "TList.h"
#include "Task.hpp"
#include "Event.hpp"
#include "EventFilter.hpp"
#include "PTHistos.hpp"

class PTCorrelator : public Task
{
public:


  //////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////
  PTCorrelator(const TString &  name,
            TransverseMomentumConfiguration * configuration,
            Event * event,
            EventFilter * eventFilter,
            ParticleFilter ** particleFilter1);
  virtual ~PTCorrelator();
  virtual void execute();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void saveHistograms(TFile * outputFile);
  virtual void scaleHistograms(double factor);
  virtual void storeEventInfo();
  virtual void fillTransverseMomentumValues();
  virtual void fillYieldValues();
  virtual void calculateTransverseMomentumMoments();

  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////
  PTHistos * histos;
  EventFilter     * eventFilter;
  ParticleFilter  ** particleFilters;
  TString * partNames;
  
  //total # of events
  int maxEvents;

  //maximum order of correlation functions (2,3,4 ...) 
  //order is at most the number of particle filters
  int maxOrder;

  //temp variable 
  int correlatorIndex;

  //store the moments of the transverse momentum per event
  // first index is event, second is filter number, third is moment (1st, 2nd ...) maxorder^2 entries per event
  double *** transverseMomentumMoments;

  //store the yields of each combination of particles per event
  //first index is event, second index is combination(in recursive order)
  double ** yields; 

  //store the multiplicty of all the events
  double * multiplicity;

  //store the centrality of all the events
  double * centrality;

  //stores the number of particles in each event
  double * numParticles;


  ClassDef(PTCorrelator,0)
};


#endif /*WAC_PTCorrelator*/
