//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "TString.h"
#include "Event.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"
#include "PythiaEventReader.hpp"
#include "PythiaConfiguration.hpp"
#include "GlobalAnalyzer.hpp"
#include "GlobalAnalyzerConfiguration.hpp"
#include "GlobalHistos.hpp"
#include "ParticleFilterAliceV0.hpp"

int main()
{
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  cout << "<INFO> PYTHIA Model Analysis - GlobalObservables" << endl;
  EventLoop * eventLoop = new EventLoop("RunPythiaGlobalObservables");
  eventLoop->setNEventRequested(10000000);
  eventLoop->setNEventReported(100000);
  eventLoop->setReportLevel(messageLevel);
  //eventLoop->setNEventPartialSave(1000000);
  //eventLoop->setPartialSave(true);
  //eventLoop->setSubsampleAnalysis(true);
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
  GlobalAnalyzerConfiguration * ac = new GlobalAnalyzerConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms         = false;
  ac->createHistograms       = true;
  ac->scaleHistograms        = true;
  ac->saveHistograms         = true;
  ac->forceHistogramsRewrite = true;
  ac->inputPath              = getenv("WAC_INPUT_PATH");;
  ac->outputPath             = getenv("WAC_OUTPUT_PATH");;
  ac->rootOuputFileName      =  "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_";
  ac->nBins_n = 200; ac->nBins_n2 = 100; ac->min_n = 0.0; ac->max_n = 200.0; ac->range_n = 200.0;
  ac->nBins_e = 200; ac->nBins_e2 = 100; ac->min_e = 0.0; ac->max_e = 2000.0; ac->range_e = 2000.0;
  ac->nBins_q = 50;  ac->nBins_q2 = 50;  ac->min_q = -25.0; ac->max_q = 25.0; ac->range_q = 50.0;
  ac->nBins_b = 50;  ac->nBins_b2 = 50;  ac->min_b = -25.0; ac->max_b = 25.0; ac->range_b = 50.0;
  ac->countParticles = true;
  ac->setEvent       = true;
  ac->fillCorrelationHistos = true;
  GlobalAnalyzerConfiguration * acFiltered = new GlobalAnalyzerConfiguration(*ac);
  acFiltered->countParticles = false;
  acFiltered->setEvent       = false;
  acFiltered->fillCorrelationHistos = true;

  EventFilter     * eventFilter = new EventFilter(EventFilter::MinBias,0.0,0.0);
  int nParticleFilters = 2;
  ParticleFilter  ** particleFilters = new ParticleFilter*[nParticleFilters];
  particleFilters[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  0.20, 100.0, -0.8, 0.8, 1.0, -1.0);
  particleFilters[1]   = new ParticleFilterAliceV0(ParticleFilterAliceV0::V0M);
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC", ac,     event, eventFilter, nParticleFilters, particleFilters, messageLevel)      );

//  0.0507493   66   66.0007
//  0.103505   53   53.0035
//  0.201895   38   38.0019
//  0.300461   28   28.0005
//  0.414744   20   20.0147
//  0.503491   16   16.0035
//  0.631744   12   12.0317
//  0.70471   10   10.0047
//  0.805348   7   7.00535
//  0.913169   2   2.01317

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
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_0_5",   acFiltered,event, eventFilter_0_5,    nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_5_10",  acFiltered,event, eventFilter_5_10,   nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_10_20", acFiltered,event, eventFilter_10_20,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_20_30", acFiltered,event, eventFilter_20_30,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_30_40", acFiltered,event, eventFilter_30_40,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_40_50", acFiltered,event, eventFilter_40_50,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_50_60", acFiltered,event, eventFilter_50_60,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_60_70", acFiltered,event, eventFilter_60_70,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_70_80", acFiltered,event, eventFilter_70_80,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_80_90", acFiltered,event, eventFilter_80_90,  nParticleFilters,particleFilters,messageLevel ));
  eventLoop->addTask( new GlobalAnalyzer("V0vsTPC_90_100",acFiltered,event, eventFilter_90_100, nParticleFilters,particleFilters,messageLevel ));

  eventLoop->run();
}
