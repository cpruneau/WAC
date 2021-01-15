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
#include "PTCorrelator.hpp"
#include "AACollisionGenerator.hpp"
#include "NuDynConfiguration.hpp"


int main()
{
  cout << "<INFO> PYTHIA Model Analysis - Single Particle Histograms" << endl;

  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  EventLoop * eventLoop = new EventLoop("RunPythiaSimulationSingleParticle");
  eventLoop->setNEventRequested(10000000);
  eventLoop->setNEventReported(100000);
  eventLoop->setReportLevel(messageLevel);
  eventLoop->setNEventPartialSave(100000);
  eventLoop->setPartialSave(true);
  eventLoop->setSubsampleAnalysis(true);
  Event * event = Event::getEvent();

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Heavy Ion and Analysis Configuration Parameters
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  NuDynConfiguration * ac = new NuDynConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms  = false;
  ac->createHistograms  = true;
  ac->scaleHistograms  = true;
  ac->calculateDerivedHistograms  = true;
  ac->saveHistograms  = true;
  ac->resetHistograms  = false;
  ac->clearHistograms  = false;
  ac->forceHistogramsRewrite  = true;
  ac->inputPath = "./";
  ac->rootInputFileName = "";
  ac->outputPath        = getenv("WAC_OUTPUT_PATH");
  ac->rootOuputFileName =  "/Pythia6";
  ac->multiplicityType = NuDynConfiguration::TotalMultiplicity;
  ac->nBins_mult   = 100;
  ac->min_mult     = 0.0;
  ac->max_mult     = 400.0;

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

  eventLoop->addTask( new PythiaEventGenerator("PYTHIA",0, event,eventFilter,particleFilter, messageLevel) );
  eventLoop->addTask( new PTCorrelator("PYTHIA_PTCorrelator_HPHMPiPPiM", ac, event, eventFilter, particleFilters, messageLevel) );
  eventLoop->run();
}



