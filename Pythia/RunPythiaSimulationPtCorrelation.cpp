//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "Timer.hpp"
#include "Event.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaEventGenerator.hpp"
#include "PythiaEventReader.hpp"
#include "PTCorrelator.hpp"


int main()
{
  Timer * t = new Timer();
  t->start();
  cout << "<INFO> PYTHIA Model Analysis - Single Particle Histograms" << endl;
  EventLoop * eventLoop = new EventLoop("RunPythiaSimulationSingleParticle");
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;
  eventLoop->setNEventRequested(10000000);
  eventLoop->setNEventReported(100);
  eventLoop->setReportLevel(messageLevel);
  eventLoop->setNEventPartialSave(1000000);
  eventLoop->setPartialSave(true);
  eventLoop->setSubsampleAnalysis(false);
  Event * event = Event::getEvent();

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Heavy Ion and Analysis Configuration Parameters
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  TransverseMomentumConfiguration * ac = new TransverseMomentumConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms  = false;
  ac->createHistograms  = true;
  ac->scaleHistograms  = true;
  ac->calculateDerivedHistograms  = true;
  ac->saveHistograms  = true;
  ac->resetHistograms  = false;
  ac->clearHistograms  = false;
  ac->forceHistogramsRewrite  = true;
  ac->inputPath = getenv("WAC_INPUT_PATH");;
  ac->rootInputFileName = "";
  ac->outputPath        = getenv("WAC_OUTPUT_PATH");
  ac->rootOuputFileName =  "/Pythia";
  ac->nBins_mult   = 100;
  ac->min_mult     = 0.0;
  ac->max_mult     = 400.0;
  ac->nBins_cent   = 100;
  ac->min_cent     = 0.0;
  ac->max_cent     = 400.0;
  ac->ptCorrelatorVsCent     = false;
  ac->ptCorrelatorVsMult     = true;
  ac->totEvents = eventLoop->getNEventRequested();
  ac->maxOrder = 4;
  ac->numTypes = 4;

  double minPt  = 0.2;
  double maxPt  = 2.0;
  double minEta = -2.0;
  double maxEta =  2.0;
  double minY   = -2.0;
  double maxY   =  2.0;

  EventFilter     * eventFilter      = new EventFilter(EventFilter::MinBias,0.0,0.0);
  ParticleFilter  * particleFilters[8];
  int nParticleFilters = 0;
  ParticleFilter  * particleFilter     = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged, minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive,minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative,minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive,minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative,minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive,minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative,minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive,minPt, maxPt, minEta, maxEta, minY, maxY );
  particleFilters[nParticleFilters++]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative,minPt, maxPt, minEta, maxEta, minY, maxY );


  int nOptions = 0;
  TString ** pythiaOptions = new TString* [50];
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
  pc->dataOutputUsed = false;
  pc->dataConversionToWac = true;
  
  //eventLoop->addTask( new PythiaEventGenerator("PYTHIA",pc, event,eventFilter,particleFilter, messageLevel) );
  
  pc->dataInputFileName = "Pythia_pp_7000.root";
  pc->dataInputTreeName = "PythiaTree";
  pc->dataInputPath     = getenv("WAC_OUTPUT_DATA_PATH");
  eventLoop->addTask( new PythiaEventReader("PYTHIA",pc, event,eventFilter,particleFilter, messageLevel) );

  eventLoop->addTask( new PTCorrelator("PYTHIA_PTCorrelator_HPHMPiPPiM", ac, event, eventFilter, particleFilters, messageLevel) );
  eventLoop->run();
  t->stop();
  t->print(cout);
}



