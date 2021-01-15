//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TString.h"
#include "TRandom.h"

//#include <iostream>
//#include <fstream>
//#include <TStyle.h>
//#include <TROOT.h>
//#include "Event.hpp"
//#include "EventLoop.hpp"
//#include "EventFilter.hpp"
//#include "ParticleFilter.hpp"
//#include "ParticleAnalyzerConfiguration.hpp"
//#include "ParticleAnalyzer.hpp"
//#include "PythiaConfiguration.hpp"
//#include "PythiaEventGenerator.hpp"
//#include "PythiaEventReader.hpp"

//int main()
int PythiaAnalysisSingleParticleAliceV0()
{
  TString includesPath = "/Users/claudeapruneau/Documents/GitHub/WAC/Base/";
  gSystem->Load(includesPath+"CanvasCollection.hpp");
  gSystem->Load(includesPath+"CanvasConfiguration.hpp");
  gSystem->Load(includesPath+"TaskConfiguration.hpp");
  gSystem->Load(includesPath+"EnhancedGraph.hpp");
  gSystem->Load(includesPath+"Factory.hpp");
  gSystem->Load(includesPath+"FunctionCollection.hpp");
  gSystem->Load(includesPath+"GraphConfiguration.hpp");
  gSystem->Load(includesPath+"HistogramCollection.hpp");
  gSystem->Load(includesPath+"Histograms.hpp");
  gSystem->Load(includesPath+"RandomGenerators.hpp");
  gSystem->Load(includesPath+"RapidityGenerator.hpp");
  gSystem->Load(includesPath+"Task.hpp");
  gSystem->Load(includesPath+"TaskCollection.hpp");
  gSystem->Load(includesPath+"Property.hpp");
  gSystem->Load(includesPath+"MessageLogger.hpp");
  gSystem->Load(includesPath+"AnalysisConfiguration.hpp");
  gSystem->Load(includesPath+"Event.hpp");
  gSystem->Load(includesPath+"EventFilter.hpp");
  gSystem->Load(includesPath+"EventHistos.hpp");
  gSystem->Load(includesPath+"EventLoop.hpp");
  gSystem->Load(includesPath+"GeneralizedGaussianConfiguration.hpp");
  gSystem->Load(includesPath+"GeneralizedGaussianProfile.hpp");
  gSystem->Load(includesPath+"GeneratorConfiguration.hpp");
  gSystem->Load(includesPath+"TwoPartCorrelationAnalyzer.hpp");
  gSystem->Load(includesPath+"Particle.hpp");
  gSystem->Load(includesPath+"ParticleFilter.hpp");
  gSystem->Load(includesPath+"ParticleHistos.hpp");
  gSystem->Load(includesPath+"ParticlePairCombinedHistos.hpp");
  gSystem->Load(includesPath+"ParticlePairDerivedHistos.hpp");
  gSystem->Load(includesPath+"ParticlePairFilter.hpp");
  gSystem->Load(includesPath+"ParticlePairHistos.hpp");
  gSystem->Load(includesPath+"SubSampleStatCalculator.hpp");
  gSystem->Load("libBase.dylib");


  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  cout << "<INFO> PythiaAnalysisSingleParticleAliceV0" << endl;
  EventLoop * eventLoop = new EventLoop("PythiaAnalysisSingleParticleAliceV0");
  eventLoop->setNEventRequested(1000);
  eventLoop->setNEventReported(100000);
  eventLoop->setReportLevel(messageLevel);
  eventLoop->setNEventPartialSave(1000000);
  eventLoop->setPartialSave(true);
  eventLoop->setSubsampleAnalysis(true);
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
  ParticleAnalyzerConfiguration * acV0A = new ParticleAnalyzerConfiguration("V0A","V0A","1.0");
  acV0A->loadHistograms         = false;
  acV0A->createHistograms       = true;
  acV0A->scaleHistograms        = true;
  acV0A->calculateDerivedHistograms  = true;
  acV0A->saveHistograms         = true;
  acV0A->resetHistograms        = false;
  acV0A->clearHistograms        = false;
  acV0A->forceHistogramsRewrite = true;
  acV0A->inputPath              = getenv("WAC_INPUT_PATH");;
  acV0A->outputPath             = getenv("WAC_OUTPUT_PATH");;
  acV0A->rootOuputFileName      =  "PYTHIA_pp_7TeV_softOnHardOff_SinglesV0A_";

  acV0A->nBins_pt    = 100;
  acV0A->min_pt      = 0.001;
  acV0A->max_pt      = 100.0;
  acV0A->nBins_eta   = 20;
  acV0A->min_eta     = -1;
  acV0A->max_eta     = 1;
  acV0A->nBins_y     = 20;
  acV0A->min_y       = -8;
  acV0A->max_y       = 8;
  acV0A->nBins_phi   = 36;
  acV0A->min_phi     = 0.0;
  acV0A->max_phi     = 2.0*3.1415927;

  ParticleAnalyzerConfiguration * acV0C = new ParticleAnalyzerConfiguration(*ac);
  acWide->nBins_eta   = 120;
  acWide->min_eta     = -8;
  acWide->max_eta     =  8;
  acWide->nBins_y     = 120;
  acWide->min_y       = -8;
  acWide->max_y       =  8;

  EventFilter     * eventFilter       = new EventFilter(EventFilter::MinBias,0.0,0.0);
  int nParticleFilters = 1;
  ParticleFilter  ** particleFiltersV0A = new ParticleFilter*[nParticleFilters];
  ParticleFilter  ** particleFiltersV0C = new ParticleFilter*[nParticleFilters];
  ParticleFilter  ** particleFiltersV0M = new ParticleFilter*[nParticleFilters];

  particleFiltersV0A[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  acV0A->min_pt,acV0A->max_pt,acV0A->min_eta,acV0A->max_eta, acV0A->min_y,acV0A->max_y);
  particleFiltersV0C[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  acV0C->min_pt,acV0C->max_pt,acV0C->min_eta,acV0C->max_eta, acV0C->min_y,acV0C->max_y);
  particleFiltersV0M[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  acV0M->min_pt,acV0M->max_pt,acV0M->min_eta,acV0M->max_eta, acV0C->min_y,acV0C->max_y);

  ParticleFilter  ** particleFiltersWide = new ParticleFilter*[nParticleFilters];
  eventLoop->addTask( new ParticleAnalyzer("V0A", acV0A, event, eventFilter, nParticleFilters, particleFiltersV0A, messageLevel)  );
  eventLoop->addTask( new ParticleAnalyzer("V0C", acV0C, event, eventFilter, nParticleFilters, particleFiltersV0C, messageLevel)  );
  eventLoop->addTask( new ParticleAnalyzer("V0M", acV0M, event, eventFilter, nParticleFilters, particleFiltersV0M, messageLevel)  );
  eventLoop->run();
}
