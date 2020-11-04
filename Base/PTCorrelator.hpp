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
            TaskConfiguration * configuration,
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
  virtual void calculateTransverseMomentumMoments();

  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////
  PTHistos * histos;
  EventFilter     * eventFilter;
  ParticleFilter  ** particleFilters;
  TString * partNames;
  
  int ** counts; //first index is event, 

  //total # of events
  int maxEvents;

  //maximum order of correlation functions (2,3,4 ...) 
  //order is at most the number of particle filters
  int maxOrder;

  //temp variable 
  int correlatorIndex;

  //store the moments of the transverse momentum per event
  // first index is event, second is moment (1st, 2nd ...)
  double ** transverseMomentumMoments;


  //stores the number of particles in each event
  int * numParticles;


  //store the transverse momentum of all the particles
  //first index is event number, second index is particle number
  double ** pT;

  //store the acceptances by the filters of all the particles
  //first index is event number, second number is filter number, third is particle number
  bool *** acceptances;


  //store the multiplicty of all the events
  double * multiplicity;

  //store the centrality of all the events
  double * centrality;


  ClassDef(PTCorrelator,0)
};


#endif /*WAC_PTCorrelator*/
