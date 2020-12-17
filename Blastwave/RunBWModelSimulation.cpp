//
//  cm_figure9_7.c
//  MyMC
//
//  Created by Claude Pruneau on 4/10/17.
//  Copyright © 2017 Claude Pruneau. All rights reserved.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "TRandom.h"

#include "BWModelConfiguration.hpp"
#include "BWModelEventGenerator.hpp"
#include "ParticlePairAnalyzerConfiguration.hpp"
#include "ParticleFilter.hpp"
#include "TwoPartCorrelationAnalyzer.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Run an analysis with Blastwave Model
/////////////////////////////////////////////////////////////////////////////////////////
int main()
{
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  EventLoop * eventLoop = new EventLoop("RunBWModelSimulation");
  eventLoop->setNEventRequested(1000000);
  eventLoop->setNEventReported(100000);
  eventLoop->setReportLevel(messageLevel);
  Event * event = Event::getEvent();

  ParticlePairAnalyzerConfiguration * ac = new ParticlePairAnalyzerConfiguration("TwoPartWithBW","TwoPartWithBW","1.0");
  ac->loadHistograms         = false;
  ac->createHistograms       = true;
  ac->scaleHistograms        = true;
  ac->calculateDerivedHistograms  = false;
  ac->saveHistograms         = true;
  ac->resetHistograms        = false;
  ac->clearHistograms        = false;
  ac->forceHistogramsRewrite = true;
  ac->inputPath              = getenv("WAC_INPUT_PATH");;
  ac->rootInputFileName      = "";
  ac->outputPath             = getenv("WAC_OUTPUT_PATH");
  ac->rootOuputFileName      =  "PYTHIA_softOnHardOff_Pairs_";
  ac->histoBaseName          = "TEST";

  ac->nBins_pt    = 100;
  ac->min_pt      = 0.2;
  ac->max_pt      = 100.0;
  ac->nBins_eta   = 40;
  ac->min_eta     = -2;
  ac->max_eta     = 2;
  ac->nBins_y     = 40;
  ac->min_y       = -2;
  ac->max_y       = 2;
  ac->nBins_phi   = 36;
  ac->validate();  // makes sure the difference boundaries are correct...

  // ==================================================================================
  // Blast Wave Model Generator Parameters
  // ==================================================================================
  BWModelConfiguration *bwc = new BWModelConfiguration("BWModel","BWModel",7);
  bwc->totalMultiplicityMean = 426.0;  /* ALICE 30-40% central events */
  bwc->totalMultiplicitySigma = 15.0;  /* ALICE 30-40% central events */
  bwc->pT_Min = 0.0;
  bwc->pT_Max = 3.0;
  bwc->useAllKinds = true;
  bwc->useBGBlastWave = true;
  bwc->temperatureAllKinds = 0.106;    /* ALICE 30-40% central events */
  bwc->betaAvgBGBW = 0.604;            /* ALICE 30-40% central events */
  bwc->betaProfileBGBW = 0.841;        /* ALICE 30-40% central events */
  bwc->longitudinalFlow = 0.0;         /* a guess */
  /* azimuthal flow */
  bwc->useFlow = true;
  bwc->ellipticFlowAllKinds = 0.09;
  bwc->triangularFlowAllKinds = 0.03;

  EventFilter     * eventFilter      = new EventFilter(EventFilter::MinBias,0.0,0.0);
  ParticleFilter  * particleFilter   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // -ve only
  ParticleFilter  * particleFilter1  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter2  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter3  = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter4  = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter5  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter6  = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only ParticleFilter  * particleFilter5  = new ParticleFilter(ParticleFilter::Pion, ParticleFilter::Positive,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter7  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  ParticleFilter  * particleFilter8  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative,ac->min_pt,ac->max_pt,ac->min_eta,ac->max_eta, ac->min_y,ac->max_y); // +ve only
  eventLoop->addTask( new BWModelEventGenerator("BW",bwc, event,messageLevel) );
  eventLoop->addTask( new TwoPartCorrelationAnalyzer ("BW_TPCA_ALL",  ac,  event, eventFilter,particleFilter1,particleFilter2,messageLevel) );
  eventLoop->addTask( new TwoPartCorrelationAnalyzer ("BW_TPCA_PiPi", ac,  event, eventFilter,particleFilter3,particleFilter4,messageLevel) );
  eventLoop->addTask( new TwoPartCorrelationAnalyzer ("BW_TPCA_KK",   ac,  event, eventFilter,particleFilter5,particleFilter6,messageLevel) );
  eventLoop->addTask( new TwoPartCorrelationAnalyzer ("BW_TPCA_PP",   ac,  event, eventFilter,particleFilter7,particleFilter8,messageLevel) );
  eventLoop->run();
}
