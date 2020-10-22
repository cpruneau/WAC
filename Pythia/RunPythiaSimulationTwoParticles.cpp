//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "Event.hpp"
#include "AnalysisConfiguration.hpp"
#include "TwoPartCorrelationAnalyzer.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"
#include "TwoPartCorrelationAnalyzer.hpp"
//#include <time.h>

int main()
{
  time_t begin,end; // time_t is a datatype to store time values.
  time (&begin); // note time before execution


  //  long nEventsRequested = 100;
  long nEventsRequested = 100;
  int  nEventsReport    = 10000;

  // ==========================
  // Event Section
  // ==========================
  Event * event = Event::getEvent();

  // ==========================
  // Generator Section
  // ==========================
  int nOptions = 0;
  TString ** pythiaOptions  = new TString* [50];
  pythiaOptions[nOptions++] = new TString("Init:showChangedSettings = on");      // list changed settings
  pythiaOptions[nOptions++] = new TString("Init:showChangedParticleData = off"); // list changed particle data
  pythiaOptions[nOptions++] = new TString("Next:numberCount = 10000");            // print message every n events
  pythiaOptions[nOptions++] = new TString("Next:numberShowInfo = 1");            // print event information n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowProcess = 0");         // print process record n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowEvent = 0");
  pythiaOptions[nOptions++] = new TString("SoftQCD:all = on");                   // Allow total sigma = elastic/SD/DD/ND
                                                                                 //pythiaOptions[nOptions++] = new TString("HardQCD:all = on");
  PythiaConfiguration * pc = new PythiaConfiguration(2212 /* p */,
                                                     2212 /* p */,
                                                     14000.0, /* energy in GeV */
                                                     nOptions,
                                                     pythiaOptions);
  EventFilter     * eventFilterGen    = new EventFilter(EventFilter::MinBias,0.0,0.0);
  ParticleFilter  * particleFilterGen = new ParticleFilter(ParticleFilter::Hadron,
                                                           ParticleFilter::Charged,
                                                           0.2,100.0,
                                                           -6.0,6.0,
                                                           -10.0,10.0);
  Task * generator = new PythiaEventGenerator("PYTHIA",pc, event,eventFilterGen,particleFilterGen);


  // ==========================
  // Analysis Section
  // ==========================
  AnalysisConfiguration * ac = new AnalysisConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms         = false;
  ac->createHistograms       = true;
  ac->scaleHistograms        = true;
  ac->calculateDerivedHistograms  = false;
  ac->saveHistograms         = true;
  ac->resetHistograms        = false;
  ac->clearHistograms        = false;
  ac->forceHistogramsRewrite = true;
  ac->inputPath              = "/Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/";
  ac->rootInputFileName      = "";
  ac->outputPath             = "/Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/";
  ac->rootOuputFileName      =  "PYTHIA_softOnHardOff_Pairs_";
  ac->histoBaseName = "TEST";

  ac->nBins_pt    = 100;
  ac->min_pt      = 0.0;
  ac->max_pt      = 100.0;
  ac->nBins_eta   = 20;
  ac->min_eta     = -1;
  ac->max_eta     = 1;
  ac->nBins_y     = 20;
  ac->min_y       = -2;
  ac->max_y       = 2;
  ac->nBins_phi   = 36;
  ac->min_phi     = 0.0;
  ac->max_phi     = 2.0*3.1415927;

  AnalysisConfiguration * acWide = new AnalysisConfiguration(*ac);
  acWide->nBins_eta   = 120;
  acWide->min_eta     = -8;
  acWide->max_eta     =  8;
  acWide->nBins_y     = 120;
  acWide->min_y       = -8;
  acWide->max_y       =  8;


  TString taskName;
  int nAnalysisTasks    = 20;
  Task ** analysisTasks = new Task*[nAnalysisTasks];

  EventFilter     * eventFilter       = new EventFilter(EventFilter::MinBias,0.0,0.0);
  int nParticleFilters = 12;
  ParticleFilter  ** particleFilters = new ParticleFilter*[nParticleFilters];
  particleFilters[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[1]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[2]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[3]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Charged,  ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[4]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[5]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[6]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Charged,  ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[7]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[8]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[9]   = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Charged,  ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[10]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[11]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative, ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);

  int iTask = 0;
  analysisTasks[iTask++]  = new TwoPartCorrelationAnalyzer("NarrowHPHM",   ac, event, eventFilter,particleFilters[1],particleFilters[2] ); // h+ vs h-
  analysisTasks[iTask++]  = new TwoPartCorrelationAnalyzer("NarrowPiPPiM", ac, event, eventFilter,particleFilters[4],particleFilters[5] ); // Pi+ vs Pi-
  analysisTasks[iTask++]  = new TwoPartCorrelationAnalyzer("NarrowKPKM",   ac, event, eventFilter,particleFilters[7],particleFilters[8] ); // K+ vs K-
  analysisTasks[iTask++]  = new TwoPartCorrelationAnalyzer("NarrowPPPM",   ac, event, eventFilter,particleFilters[10],particleFilters[11] ); // P+ vs P-
  analysisTasks[iTask++]  = new TwoPartCorrelationAnalyzer("NarrowPiCKC",  ac, event, eventFilter,particleFilters[3],particleFilters[6] ); // PiC vs KC
  analysisTasks[iTask++]  = new TwoPartCorrelationAnalyzer("NarrowPiCPC",  ac, event, eventFilter,particleFilters[3],particleFilters[9] ); // PiC vs KC
  analysisTasks[iTask++]  = new TwoPartCorrelationAnalyzer("NarrowKCPC",   ac, event, eventFilter,particleFilters[6],particleFilters[9] ); // KC vs KC
  nAnalysisTasks = iTask;

//  ParticleFilter  ** particleFiltersWide = new ParticleFilter*[nParticleFilters];
//  particleFiltersWide[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[1]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[2]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[3]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Charged,  acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[4]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[5]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[6]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Charged,  acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[7]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[8]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[9]   = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Charged,  acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[10]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  particleFiltersWide[11]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative, acWide->min_pt+0.001,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
//  analysisTasks[1]         = new ParticleAnalyzer("Wide",   acWide, event, eventFilter, nParticleFilters, particleFiltersWide);

  // ==========================
  // Event Loop
  // ==========================

  EventLoop * eventLoop = new EventLoop();
  eventLoop->addTask( generator );
  for (int iAnalysisTask=0;iAnalysisTask<nAnalysisTasks;iAnalysisTask++)
  {
    eventLoop->addTask( analysisTasks[iAnalysisTask] );
  }
  eventLoop->run(nEventsRequested,nEventsReport);

  cout << "<INFO> PYTHIA Model Analysis - Pair Differential Correlations Histograms  - Completed" << endl;
  time (&end); // note time after execution
  double difference = difftime (end,begin);
  cout << "<INFO> in " <<  difference << " seconds";
}
