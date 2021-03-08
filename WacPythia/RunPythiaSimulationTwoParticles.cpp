//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "TString.h"
#include "Event.hpp"
#include "EventFilter.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "ParticleFilterAliceV0.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"
#include "PythiaEventReader.hpp"
#include "TwoPartCorrelationAnalyzer.hpp"
#include "GlobalAnalyzer.hpp"
#include "GlobalAnalyzerConfiguration.hpp"

int main()
{
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  EventLoop * eventLoop = new EventLoop("RunPythiaSimulationTwoParticles");
  eventLoop->setNEventRequested(10000000);
  eventLoop->setNEventReported(10000);
  eventLoop->setReportLevel(messageLevel);
  eventLoop->setNEventPartialSave(1000000);
  eventLoop->setPartialSave(true);
  eventLoop->setSubsampleAnalysis(true);
  eventLoop->setReportLevel(messageLevel);
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
                                                     7000.0, /* energy in GeV */
                                                     nOptions,
                                                     pythiaOptions);
  pc->dataOutputUsed      = false;
  pc->dataConversionToWac = true;
  pc->dataOutputFileName  = "Pythia_pp_7000.root";
  pc->dataOutputTreeName  = "PythiaTree";
  pc->dataOutputPath      = getenv("WAC_OUTPUT_DATA_PATH");
  pc->dataInputUsed       = true;
  pc->dataInputFileName   = "Pythia_pp_7000_10Mevents.root";
  pc->dataInputTreeName   = "PythiaTree";
  pc->dataInputPath       = getenv("WAC_INPUT_DATA_PATH");
  pc->saveHistograms      = false;



  EventFilter     * eventFilterGen    = new EventFilter(EventFilter::MinBias,0.0,0.0);
  ParticleFilter  * particleFilterGen = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  0.001,100.0, -6.0,6.0, -10.0,10.0);
  bool generateEvents = false;
  if (generateEvents)
    {
    eventLoop->addTask(  new PythiaEventGenerator("PYTHIA",pc, event,eventFilterGen,particleFilterGen, messageLevel) );
    }
  else
    {
    eventLoop->addTask( new PythiaEventReader("PYTHIA",pc, event,eventFilterGen,particleFilterGen, messageLevel) );
    }

  // ==========================
  // Analysis Section
  // ==========================
  GlobalAnalyzerConfiguration * gac = new GlobalAnalyzerConfiguration("PYTHIA","PYTHIA","1.0");
  gac->loadHistograms         = false;
  gac->createHistograms       = true;
  gac->scaleHistograms        = true;
  gac->saveHistograms         = true;
  gac->forceHistogramsRewrite = true;
  gac->inputPath              = getenv("WAC_INPUT_PATH");;
  gac->outputPath             = getenv("WAC_OUTPUT_PATH");;
  gac->rootOuputFileName      =  "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservablesForTwoPartCorr_";
  gac->nBins_n = 200; gac->nBins_n2 = 100; gac->min_n = 0.0; gac->max_n = 200.0; gac->range_n = 200.0;
  gac->nBins_e = 200; gac->nBins_e2 = 100; gac->min_e = 0.0; gac->max_e = 2000.0; gac->range_e = 2000.0;
  gac->nBins_q = 50;  gac->nBins_q2 = 50;  gac->min_q = -25.0; gac->max_q = 25.0; gac->range_q = 50.0;
  gac->nBins_b = 50;  gac->nBins_b2 = 50;  gac->min_b = -25.0; gac->max_b = 25.0; gac->range_b = 50.0;
  gac->countParticles = true;
  gac->setEvent       = true;
  gac->fillCorrelationHistos = true;

  EventFilter     * eventFilter = new EventFilter(EventFilter::MinBias,0.0,0.0);
  int nParticleFilters = 2;
  ParticleFilter  ** particleFilters = new ParticleFilter*[nParticleFilters];
  particleFilters[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  0.20, 100.0, -0.8, 0.8, 1.0, -1.0);
  particleFilters[1]   = new ParticleFilterAliceV0(ParticleFilterAliceV0::V0M);
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC", gac, event, eventFilter, nParticleFilters, particleFilters, messageLevel)      );

  ParticlePairAnalyzerConfiguration * ppac = new ParticlePairAnalyzerConfiguration("TwoPartWithPYTHIA","TwoPartWithPYTHIA","1.0");
  ppac->loadHistograms         = false;
  ppac->createHistograms       = true;
  ppac->scaleHistograms        = true;
  ppac->calculateDerivedHistograms  = true;
  ppac->saveHistograms         = true;
  ppac->resetHistograms        = false;
  ppac->clearHistograms        = false;
  ppac->forceHistogramsRewrite = true;
  ppac->inputPath              = getenv("WAC_INPUT_PATH");;
  ppac->rootInputFileName      = "";
  ppac->outputPath             = getenv("WAC_OUTPUT_PATH");
  ppac->rootOuputFileName      =  "PYTHIA_pp_7TeV_softOnHardOff_Pairs_";

  ppac->nBins_pt    = 18;
  ppac->min_pt      = 0.2;
  ppac->max_pt      = 2.0;
  ppac->nBins_eta   = 16;
  ppac->min_eta     = -0.8;
  ppac->max_eta     =  0.8;
  ppac->nBins_phi   = 72;
  ppac->min_phi     = 0.0;
  ppac->max_phi     = 2.0*3.1415927;
  ppac->validate();
  
  EventFilter * eventFilter_0_100  = new EventFilter(EventFilter::MinBias,0.0,0.0);
  EventFilter * eventFilter_0_5    = new EventFilter(EventFilter::V0M,66.0,1000.0);
  EventFilter * eventFilter_5_10   = new EventFilter(EventFilter::V0M,53.0,  66.0);
  EventFilter * eventFilter_10_20  = new EventFilter(EventFilter::V0M,38.0,  53.0);
  EventFilter * eventFilter_20_30  = new EventFilter(EventFilter::V0M,28.0,  38.0);
  EventFilter * eventFilter_30_40  = new EventFilter(EventFilter::V0M,20.0,  28.0);
  EventFilter * eventFilter_40_50  = new EventFilter(EventFilter::V0M,16.0,  20.0);
  EventFilter * eventFilter_50_60  = new EventFilter(EventFilter::V0M,12.0,  16.0);
  EventFilter * eventFilter_60_70  = new EventFilter(EventFilter::V0M,10.0,  12.0);
  EventFilter * eventFilter_70_80  = new EventFilter(EventFilter::V0M, 7.0,  10.0);
  EventFilter * eventFilter_80_90  = new EventFilter(EventFilter::V0M, 2.0,   7.0);
  EventFilter * eventFilter_90_100 = new EventFilter(EventFilter::V0M, 0.0,   2.0);

  
  nParticleFilters = 12;
  ParticleFilter  ** particleFilters2 = new ParticleFilter*[nParticleFilters];
  particleFilters2[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[1]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[2]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[3]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Charged,  ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[4]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[5]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[6]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Charged,  ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[7]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[8]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[9]   = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Charged,  ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[10]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);
  particleFilters2[11]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative, ppac->min_pt,ppac->max_pt,ppac->min_eta,ppac->max_eta, ppac->min_y,ppac->max_y);

  int option = 1;

  switch (option)
    {
      case 0: // inclusive for many types of pairs.
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM",   ppac, event, eventFilter_0_100,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_PiPPiM", ppac, event, eventFilter_0_100,particleFilters2[4],particleFilters2[5],messageLevel ) ); // Pi+ vs Pi-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_KPKM",   ppac, event, eventFilter_0_100,particleFilters2[7],particleFilters2[8],messageLevel ) ); // K+ vs K-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_PPPM",   ppac, event, eventFilter_0_100,particleFilters2[10],particleFilters2[11],messageLevel ) ); // P+ vs P-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_PiCKC",  ppac, event, eventFilter_0_100,particleFilters2[3],particleFilters2[6],messageLevel ) ); // PiC vs KC
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_PiCPC",  ppac, event, eventFilter_0_100,particleFilters2[3],particleFilters2[9],messageLevel ) ); // PiC vs KC
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_KCPC",   ppac, event, eventFilter_0_100,particleFilters2[6],particleFilters2[9],messageLevel ) ); // KC vs KC
      break;
      case 1: // inclusive and for various x-sect slices
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_0_100",   ppac, event, eventFilter_0_100,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_0_5",     ppac, event, eventFilter_0_5,  particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_5_10",    ppac, event, eventFilter_5_10, particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_10_20",   ppac, event, eventFilter_10_20,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_20_30",   ppac, event, eventFilter_20_30,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_30_40",   ppac, event, eventFilter_30_40,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_40_50",   ppac, event, eventFilter_40_50,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_50_60",   ppac, event, eventFilter_50_60,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_60_70",   ppac, event, eventFilter_60_70,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_70_80",   ppac, event, eventFilter_70_80,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_80_90",   ppac, event, eventFilter_80_90,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      eventLoop->addTask( new TwoPartCorrelationAnalyzer("Narrow_HPHM_90_100",  ppac, event, eventFilter_90_100,particleFilters2[1],particleFilters2[2],messageLevel ) ); // h+ vs h-
      break;
    }
   eventLoop->run();
}
