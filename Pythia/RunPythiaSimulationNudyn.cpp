//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "Timer.hpp"
#include "Event.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"
#include "NuDynTask.hpp"

int main()
{
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  EventLoop * eventLoop = new EventLoop("RunPythiaSimulationNuDyn");
  eventLoop->setNEventRequested(1000000);
  eventLoop->setNEventReported(100000);
  eventLoop->setReportLevel(messageLevel);
  eventLoop->setNEventPartialSave(100000);
  eventLoop->setPartialSave(true);
  eventLoop->setSubsampleAnalysis(true);
  Event * event = Event::getEvent();

  // ==========================
  // Generator Section
  // ==========================
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
                                                     14000.0, /* energy in GeV */
                                                     nOptions,
                                                     pythiaOptions);
  EventFilter     * eventFilterGen    = new EventFilter(EventFilter::MinBias,0.0,0.0);
  ParticleFilter  * particleFilterGen = new ParticleFilter(ParticleFilter::Hadron,
                                                           ParticleFilter::Charged,
                                                           0.2,100.0,
                                                           -1.0,1.0,
                                                           -5.0,5.0);
  eventLoop->addTask( new PythiaEventGenerator("PYTHIA",pc, event,eventFilterGen,particleFilterGen, messageLevel) );
  // ==========================
  // Analysis Section
  // ==========================
  NuDynConfiguration * ac = new NuDynConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms         = false;
  ac->createHistograms       = true;
  ac->scaleHistograms        = true;
  ac->calculateDerivedHistograms  = true;
  ac->saveHistograms         = true;
  ac->resetHistograms        = false;
  ac->clearHistograms        = false;
  ac->forceHistogramsRewrite = true;
  ac->inputPath              = getenv("WAC_INPUT_PATH");
  ac->outputPath             = getenv("WAC_OUTPUT_PATH");
  ac->dataSourceName          = "PYTHIA8.0";
  ac->collisionSystemName     = "pp";
  ac->collisionEnergyName     = "14TeV";
  ac->triggerName             = "MB";
  ac->taskTypeName            = "NuDyn";
  ac->eventFilterName         = "MB";
  ac->particleFilterName      = "Varia";
  ac->otherName               = "";
  ac->dataSourceTitle         = "PYTHIA 8.0";
  ac->collisionSystemTitle    = "p-p";
  ac->collisionEnergyTitle    = "14 TeV";
  ac->triggerTitle            = "Min-Bias";
  ac->taskTypeTitle           = "NuDyn";
  ac->eventFilterTitle        = "Min-Bias";
  ac->particleFilterTitle     = "Varia";
  ac->setRootOutputFileName(2);

  ac->multiplicityType = NuDynConfiguration::TotalMultiplicity;
  ac->nBins_mult   = 100;
  ac->min_mult     = 0.0;
  ac->max_mult     = 200.0;

  double minPt  = 0.2;
  double maxPt  = 2.0;
  double minEta = -2.0;
  double maxEta =  2.0;
  double minY   = -2.0;
  double maxY   =  2.0;

  EventFilter     * eventFilter        = new EventFilter(EventFilter::MinBias,0.0,0.0);
  ParticleFilter  * particleFilter_HP  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive, minPt, maxPt, minEta, maxEta, minY, maxY );
  ParticleFilter  * particleFilter_HM  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative, minPt, maxPt, minEta, maxEta, minY, maxY );
  ParticleFilter  * particleFilter_PiP = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive, minPt, maxPt, minEta, maxEta, minY, maxY );
  ParticleFilter  * particleFilter_PiM = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative, minPt, maxPt, minEta, maxEta, minY, maxY );
  ParticleFilter  * particleFilter_KP  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive, minPt, maxPt, minEta, maxEta, minY, maxY );
  ParticleFilter  * particleFilter_KM  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative, minPt, maxPt, minEta, maxEta, minY, maxY );
  ParticleFilter  * particleFilter_PP  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive, minPt, maxPt, minEta, maxEta, minY, maxY );
  ParticleFilter  * particleFilter_PM  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative, minPt, maxPt, minEta, maxEta, minY, maxY );

  eventLoop->addTask(  new NuDynTask("_HPHPHPHP", ac, event, eventFilter,particleFilter_HP,particleFilter_HP,particleFilter_HP,particleFilter_HP, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_HPHPHPHM", ac, event, eventFilter,particleFilter_HP,particleFilter_HP,particleFilter_HP,particleFilter_HM, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_HPHPHMHM", ac, event, eventFilter,particleFilter_HP,particleFilter_HP,particleFilter_HM,particleFilter_HM, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_HPHMHMHM", ac, event, eventFilter,particleFilter_HP,particleFilter_HM,particleFilter_HM,particleFilter_HM, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_HMHMHMHM", ac, event, eventFilter,particleFilter_HM,particleFilter_HM,particleFilter_HM,particleFilter_HM, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_KPKPKPKP", ac, event, eventFilter,particleFilter_KP,particleFilter_KP,particleFilter_KP,particleFilter_KP, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_KPKPKPKM", ac, event, eventFilter,particleFilter_KP,particleFilter_KP,particleFilter_KP,particleFilter_KM, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_KPKPKMKM", ac, event, eventFilter,particleFilter_KP,particleFilter_KP,particleFilter_KM,particleFilter_KM, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_KPKMKMKM", ac, event, eventFilter,particleFilter_KP,particleFilter_KM,particleFilter_KM,particleFilter_KM, messageLevel) );
  eventLoop->addTask(  new NuDynTask("_KMKMKMKM", ac, event, eventFilter,particleFilter_KM,particleFilter_KM,particleFilter_KM,particleFilter_KM, messageLevel) );
  eventLoop->run();
}
