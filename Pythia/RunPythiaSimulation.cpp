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
#include "PythiaEventGenerator.hpp"
#include "NuDynTask.hpp"

int main()
{
  cout << "<INFO> PYTHIA Model Analysis - Starting" << endl;

//  long nEventsRequested = 100;
  long nEventsRequested = 50000000;
  int  nEventsReport    = 100000;

  AnalysisConfiguration * ac = new AnalysisConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms  = false;
  ac->createHistograms  = true;
  ac->scaleHistograms  = true;
  ac->calculateDerivedHistograms  = true;
  ac->saveHistograms  = true;
  ac->resetHistograms  = false;
  ac->clearHistograms  = false;
  ac->forceHistogramsRewrite  = true;
  ac->inputPath = "./";
  ac->configurationFileName = "configuration";
  ac->rootInputFileName = "";
  ac->outputPath = getenv("OUTPUT_PATH");
  ac->rootOuputFileName =  "results";
  ac->histoBaseName =  "pythia";

  ac->nBins_pt    = 40;
  ac->min_pt      = 0.2;
  ac->max_pt      = 3.0;
  ac->nBins_eta   = 20;
  ac->min_eta     = -2;
  ac->max_eta     = 2;
  ac->nBins_y     = 20;
  ac->min_y       = -2;
  ac->max_y       = 2;
  ac->nBins_phi   = 36;
  ac->min_phi     = 0.0;
  ac->max_phi     = 2.0*3.1415927;

  ac->nBins_DeltaPlong = 40;
  ac->min_DeltaPlong   = -1.0;
  ac->max_DeltaPlong   =  1.0;
  ac->nBins_DeltaPside = 40;
  ac->min_DeltaPside   = -1.0;
  ac->max_DeltaPside   =  1.0;
  ac->nBins_DeltaPout  = 40;
  ac->min_DeltaPout    = -1.0;
  ac->max_DeltaPout    =  1.0;

  ac->fillPairs        = true;
  ac->fill3D           = false;
  ac->fill6D           = false;
  ac->fillQ3D          = false;
  ac->fillY            = true;

  ac->nuDynVsMult     = true;
  ac->nuDynVsCent     = false;
  ac->nBins_mult   = 100;
  ac->min_mult     = 0.0;
  ac->max_mult     = 400.0;
  ac->nBins_cent   = 20;
  ac->min_cent     = 0.0;
  ac->max_cent     = 100.0;

  EventFilter     * eventFilter      = new EventFilter(EventFilter::MinBias,0.0,0.0);

  ParticleFilter  * particleFilter     = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_HP  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_HM  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_PiP = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_PiM = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_KP  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_KM  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only ParticleFilter  * particleFilter_KP  = new ParticleFilter(ParticleFilter::Pion, ParticleFilter::Positive,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_PP  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter_PM  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative,ac->min_pt+0.001,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only

  int order = 4;// order <= particleFilters1 length
  ParticleFilter * particleFilters1[8];
  particleFilters1[0] = particleFilter_HP;
  particleFilters1[1] = particleFilter_HM;
  particleFilters1[2] = particleFilter_PiP;
  particleFilters1[3] = particleFilter_PiM;
  particleFilters1[4] = particleFilter_KP;
  particleFilters1[5] = particleFilter_KM;
  particleFilters1[6] = particleFilter_PP;
  particleFilters1[7] = particleFilter_PM;

  Event * event = Event::getEvent();
  EventLoop * eventLoop = new EventLoop();
  //eventLoop->addTask( new PythiaEventGenerator("PYTHIA",0, event,eventFilter,particleFilter) );
  eventLoop->addTask( new AACollisionGenerator("PYTHIA_PbPbEventGenerator",ac, event,eventFilter,particleFilter, 208) ); // lead has 208 nucleons

  //t->setReportLevel(MessageLogger::Debug);
  eventLoop->addTask( new PTCorrelator("PYTHIA_PTCorrelator_HPHMPiPPiM", ac, event, eventFilter, particleFilters1, order, nEventsRequested) );
  /*
NuDynTask * t = new NuDynTask("PYTHIA_NuDyn_HPHPHPHP", ac, event, eventFilter,particleFilter_HP,particleFilter_HP,particleFilter_HP,particleFilter_HP);
  //t->setReportLevel(MessageLogger::Debug);
  
  eventLoop->addTask( t );
  eventLoop->addTask( new NuDynTask("PYTHIA_NuDyn_HPHPHPHM", ac, event, eventFilter,particleFilter_HP,particleFilter_HP,particleFilter_HP,particleFilter_HM) );
  eventLoop->addTask( new NuDynTask("PYTHIA_NuDyn_HPHPHMHM", ac, event, eventFilter,particleFilter_HP,particleFilter_HP,particleFilter_HM,particleFilter_HM) );
  eventLoop->addTask( new NuDynTask("PYTHIA_NuDyn_HPHMHMHM", ac, event, eventFilter,particleFilter_HP,particleFilter_HM,particleFilter_HM,particleFilter_HM) );
  eventLoop->addTask( new NuDynTask("PYTHIA_NuDyn_HMHMHMHM", ac, event, eventFilter,particleFilter_HM,particleFilter_HM,particleFilter_HM,particleFilter_HM) );
  */
  //NuDynTask * nudyn_PiPPiPPiPPiP = new NuDynTask("PYTHIA_NuDyn_PiPi",   ac,  event, eventFilter,particleFilter_PiP,particleFilter_PiM);
  eventLoop->run(nEventsRequested,nEventsReport);

  cout << "<INFO> PYTHIA Analysis - Completed" << endl;

}

//  TwoPartCorrelationAnalyzer  * ana1 = new TwoPartCorrelationAnalyzer ("PYTHIA_TPCA_ALL",  ac,  event, eventFilter,particleFilter_HP,particleFilter_HM);
//  TwoPartCorrelationAnalyzer  * ana2 = new TwoPartCorrelationAnalyzer ("PYTHIA_TPCA_PiPi", ac,  event, eventFilter,particleFilter_PiP,particleFilter_PiM);
//  TwoPartCorrelationAnalyzer  * ana3 = new TwoPartCorrelationAnalyzer ("PYTHIA_TPCA_KK",   ac,  event, eventFilter,particleFilter_KP,particleFilter_KM);
//  TwoPartCorrelationAnalyzer  * ana4 = new TwoPartCorrelationAnalyzer ("PYTHIA_TPCA_PP",   ac,  event, eventFilter,particleFilter_PP,particleFilter_PM);


