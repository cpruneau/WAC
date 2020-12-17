//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "Event.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "ParticleAnalyzerConfiguration.hpp"
#include "ParticleAnalyzer.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"
#include "PythiaEventReader.hpp"

int main()
{
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  cout << "<INFO> PYTHIA Model Analysis - Single Particle Histograms" << endl;
  EventLoop * eventLoop = new EventLoop("RunPythiaSimulationSingleParticle");
  eventLoop->setNEventRequested(10000000);
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
  ParticleAnalyzerConfiguration * ac = new ParticleAnalyzerConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms         = false;
  ac->createHistograms       = true;
  ac->scaleHistograms        = true;
  ac->calculateDerivedHistograms  = true;
  ac->saveHistograms         = true;
  ac->resetHistograms        = false;
  ac->clearHistograms        = false;
  ac->forceHistogramsRewrite = true;
  ac->inputPath              = getenv("WAC_INPUT_PATH");;
  ac->outputPath             = getenv("WAC_OUTPUT_PATH");;
  ac->rootOuputFileName      =  "PYTHIA_pp_7TeV_softOnHardOff_Singles_";

  ac->nBins_pt    = 100;
  ac->min_pt      = 0.001;
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

  ParticleAnalyzerConfiguration * acWide = new ParticleAnalyzerConfiguration(*ac);
  acWide->nBins_eta   = 120;
  acWide->min_eta     = -8;
  acWide->max_eta     =  8;
  acWide->nBins_y     = 120;
  acWide->min_y       = -8;
  acWide->max_y       =  8;

  EventFilter     * eventFilter       = new EventFilter(EventFilter::MinBias,0.0,0.0);
  int nParticleFilters = 12;
  ParticleFilter  ** particleFilters = new ParticleFilter*[nParticleFilters];
  particleFilters[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[1]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[2]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[3]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Charged,  ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[4]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[5]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[6]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Charged,  ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[7]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[8]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[9]   = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Charged,  ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[10]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);
  particleFilters[11]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative, ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y);

  ParticleFilter  ** particleFiltersWide = new ParticleFilter*[nParticleFilters];
  particleFiltersWide[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[1]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[2]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[3]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Charged,  acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[4]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[5]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[6]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Charged,  acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[7]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[8]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[9]   = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Charged,  acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[10]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  particleFiltersWide[11]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative, acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
  eventLoop->addTask( new ParticleAnalyzer("Narrow", ac,     event, eventFilter, nParticleFilters, particleFilters, messageLevel)      );
  eventLoop->addTask( new ParticleAnalyzer("Wide",   acWide, event, eventFilter, nParticleFilters, particleFiltersWide, messageLevel)  );
  eventLoop->run();
}
