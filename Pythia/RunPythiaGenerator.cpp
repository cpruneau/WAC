//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "Timer.hpp"
#include "Event.hpp"
#include "TwoPartCorrelationAnalyzer.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"

int main()
{
  EventLoop * eventLoop = new EventLoop("RunPythiaGenerator");
  eventLoop->setNEventRequested(1000000);
  eventLoop->setNEventReported(10000);
  //eventLoop->setReportLevel(MessageLogger::Debug);
  Event * event = nullptr; //Event::getEvent();
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
                                                     7000.0, /* energy in GeV */
                                                     nOptions,
                                                     pythiaOptions);
  pc->dataOutputUsed = true;
  pc->dataConversionToWac = false;
  pc->dataOutputFileName = "Pythia_pp_7000_1million.root";
  pc->dataOutputTreeName = "PythiaTree";
  pc->dataOutputPath     = getenv("WAC_OUTPUT_DATA_PATH");
  EventFilter     * eventFilterGen    = nullptr;// new EventFilter(EventFilter::MinBias,0.0,0.0);
  ParticleFilter  * particleFilterGen = nullptr; //new ParticleFilter(ParticleFilter::Hadron,
//                                                           ParticleFilter::Charged,
//                                                           0.2,100.0,
//                                                           -1.0,1.0,
//                                                           -5.0,5.0);

  Task * generator = new PythiaEventGenerator("PYTHIA",pc, event,eventFilterGen,particleFilterGen, MessageLogger::Info);
  //generator->setReportLevel(MessageLogger::Debug);

  eventLoop->addTask(generator);
  eventLoop->run();
}
