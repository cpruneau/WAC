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
  //virtual void loadBaseHistograms(TFile * inputFile);
  virtual void saveHistograms(TFile * outputFile);
  //  virtual void addHistogramsToExtList(TList *list, bool all=false);
  virtual void scaleHistograms(double factor);
  virtual void calculateAverage();
  virtual double calculateS(int * filters, int order, int curFilterIndex, int & count, int * particles);
  virtual void fillSValues(int depth, int filterIndex, int * filters, int & count, int *particles);
  virtual void storeEventInfo();


  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////
  PTHistos * histos;
  EventFilter     * eventFilter;
  ParticleFilter  ** particleFilters;
  TString * partNames;

  //Average pT of of particles in the event
  double eventAveragept;

  //overall event average counts of particles, pairs, triples ... NOT IN ORDER (They are in "recursive order")
  double* avgCounts;
  int ** counts; //first index is event, 

  //total # of events
  int maxEvents;

  //maximum order of correlation functions (2,3,4 ...) 
  //order is at most the number of particle filters
  int maxOrder;

  //temp variable used
  int correlatorIndex;


  //store the transverse momentum of all the particles
  //first index is event number, second index is particle number
  double ** pT;

  //store the acceptances by the filters of all the particles
  //first index is event number, second number is filter number, third is particle number
  bool *** acceptances;

  //store the overall event avgpT of each of the particle filters
  double * avgpT;

  //store the multiplicty of all the events
  double * multiplicity;

  //store the centrality of all the events
  double * centrality;

  // store the S values for all the events. NOT IN ORDER (They are in "recursive order")
  //first index is event number, second index is correlator number(in "recursive order")
  double ** S;



  ClassDef(PTCorrelator,0)
};


#endif /*WAC_PTCorrelator*/