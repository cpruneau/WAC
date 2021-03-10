#include <iostream>
#include <chrono>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>

#include "Event.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "ParticleAnalyzer.hpp"
#include "AACollisionPythiaGenerator.hpp"
#include "AACollisionReader.hpp"
#include "CollisionGeometryGenerator.hpp"
#include "CollisionGeometryAnalyzer.hpp"
#include "CollisionGeometry.hpp"
#include "PTCorrelator.hpp"
#include "RadialBoostTask.hpp"
#include "RadialBoostConfiguration.hpp"
#include "Factory.hpp"
#include "Particle.hpp"
#include "TMath.h"


int main()
{
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  EventLoop * eventLoop = new EventLoop("RunAACollisionPythiaPtSimulation");
  eventLoop->setNEventRequested(1000000);
  eventLoop->setNEventReported(10000);
  eventLoop->setReportLevel(messageLevel);
  eventLoop->setNEventPartialSave(-1);
  eventLoop->setPartialSave(false);
  eventLoop->setSubsampleAnalysis(false);


  ///////////////////////////////////////////////////////////
  // Pythia Generator Configuration
  ////////////////////////////////////////////////////////////
  int nOptions = 0;
  TString ** pythiaOptions  = new TString* [50];
  pythiaOptions[nOptions++] = new TString("Init:showChangedSettings = on");      // list changed settings
  pythiaOptions[nOptions++] = new TString("Init:showChangedParticleData = off"); // list changed particle data
                                                                                 // Pick new random number seed for each run, based on clock.
  pythiaOptions[nOptions++] = new TString("Random:setSeed = on");
  pythiaOptions[nOptions++] = new TString("Random:seed = 0");
  pythiaOptions[nOptions++] = new TString("Next:numberCount = 10000");            // print message every n events
  pythiaOptions[nOptions++] = new TString("Next:numberShowInfo = 1");            // print event information n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowProcess = 0");         // print process record n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowEvent = 0");
  //pythiaOptions[nOptions++] = new TString("SoftQCD:inelastic = on");             //Setting for Minumum bias INEL

  pythiaOptions[nOptions++] = new TString("SoftQCD:all = on");                   // Allow total sigma = elastic/SD/DD/ND
  //pythiaOptions[nOptions++] = new TString("HardQCD:all = on");

  // hard process -- do not turh on with SoftQCD:inelastic = on
  // otherwise, there will be double counting.
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("PhaseSpace:pTHatMin = 60.");
  //pythia.readString("PhaseSpace:pTHatMax = 1000.");

  PythiaConfiguration * pc = new PythiaConfiguration(2212, /* p */
                                                     2212, /* p */
                                                     2760.0, /* energy in GeV */
                                                     nOptions,
                                                     pythiaOptions,
                                                     true,     // only pp
                                                     true,     // remove photons
                                                     10000);
  pc->dataOutputUsed = false;
  pc->dataConversionToWac = true;
  
  //eventLoop->addTask( new PythiaEventGenerator("PYTHIA",pc, event,eventFilter,particleFilter, messageLevel) );
  
  pc->dataInputFileName = "Pythia_pp_2760_10million.root";
  pc->ppdataInputTreeName[0] = "PythiaTree";
  pc->dataInputPath     = getenv("WAC_OUTPUT_DATA_PATH");



  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Transverse Momentum Correlation Configuration Parameters
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
  ac->inputPath = "./";
  ac->rootInputFileName = "";
  ac->outputPath        = getenv("WAC_OUTPUT_PATH");
  ac->rootOuputFileName =  "/Pythia";
  ac->nBins_mult   = 100;
  ac->min_mult     = 0.0;
  ac->max_mult     = 80000.0;
  ac->nBins_cent   = 100;
  ac->min_cent     = 0.0;
  ac->max_cent     = 1.0;
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

  ///////////////////////////////////////////////////////////////////////////////
  // Radial Boost Configuration
  ///////////////////////////////////////////////////////////////////////////////
  RadialBoostConfiguration * rc = new RadialBoostConfiguration("PYTHIARadialBoost", "PYTHIARadialBoost", "1.0");
  rc->param_a = 0.1;//0.1 hardboost, 0.05 soft boost
  rc->param_b = 1.0;
  rc->nBins_phi = 72;
  rc->min_phi = 0.0;
  rc->max_phi = TMath::TwoPi();
  rc->nBins_r = 100;
  rc->min_r = 0.0;
  rc->max_r = 1.0;
  rc->nBins_beta = 100;
  rc->min_beta = 0.0;
  rc->max_beta = 1.0;




  //////////////////////////////////////////////////////////////////////////////
  // Particle and Event Filters
  //////////////////////////////////////////////////////////////////////////////

  EventFilter     * eventFilter        = new EventFilter(EventFilter::MinBias,0.0,0.0);
  int nParticleFilters = 12;
  ParticleFilter  *  particleFilter  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,   minPt, maxPt, minEta, maxEta, minY, maxY);
  ParticleFilter  ** particleFilters = new ParticleFilter*[nParticleFilters];
  particleFilters[8]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,   minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[0]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Positive,  minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[1]   = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Negative,  minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[9]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Charged,   minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[2]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Positive,  minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[3]   = new ParticleFilter(ParticleFilter::Pion,   ParticleFilter::Negative,  minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[10]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Charged,   minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[4]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Positive,  minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[5]   = new ParticleFilter(ParticleFilter::Kaon,   ParticleFilter::Negative,  minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[11]   = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Charged,   minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[6]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Positive,  minPt, maxPt, minEta, maxEta, minY, maxY);
  particleFilters[7]  = new ParticleFilter(ParticleFilter::Proton, ParticleFilter::Negative,  minPt, maxPt, minEta, maxEta, minY, maxY);


 
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Collision Geometry Configuration Parameters
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  CollisionGeometryConfiguration * geometryConfiguration = new CollisionGeometryConfiguration("PYTHIACollisionGeometry","PYTHIACollisionGeometry","1.0");
  geometryConfiguration->nProtonsA  = 82;
  geometryConfiguration->nNeutronsA = 208-82;
  geometryConfiguration->nProtonsB  = 82;
  geometryConfiguration->nNeutronsB = 208-82;
  geometryConfiguration->inputPath  = getenv("WAC_INPUT_PATH");
  geometryConfiguration->outputPath = getenv("WAC_OUTPUT_PATH");
  geometryConfiguration->rootOuputFileName =  "/CollisionGeometry";
  geometryConfiguration->histoBaseName     =  "geom";
  geometryConfiguration->minB = 0.0;
  geometryConfiguration->maxB = 18.0;

  // NucleusGenerator::WoodsSaxon, 6.62 , 0.546, 0.0, 11000,0.0,11.0);
  // NucleusGenerator::WoodsSaxon, 6.62 , 0.546, 0.0, 11000,0.0,11.0);
  geometryConfiguration->aParA = 7.1;
  geometryConfiguration->aParB = 0.535;
  geometryConfiguration->aParC = 0.0;
  geometryConfiguration->aNR   = 10000;
  geometryConfiguration->aMinR = 0.0;
  geometryConfiguration->aMaxR = 8.0;

  geometryConfiguration->bParA = 7.1;
  geometryConfiguration->bParB = 0.535;
  geometryConfiguration->bParC = 0.0;
  geometryConfiguration->bNR   = 10000;
  geometryConfiguration->bMinR = 0.0;
  geometryConfiguration->bMaxR = 8.0;


  geometryConfiguration->nnCrossSection = 4.5;  // in fm^2 -- Config C
  geometryConfiguration->nBins_b = 120;
  geometryConfiguration->min_b   = 0.0;
  geometryConfiguration->max_b   = 18.0;
  geometryConfiguration->nBins_nPart = 100;
  geometryConfiguration->min_nPart   = 0;
  geometryConfiguration->max_nPart   = 500;
  geometryConfiguration->nBins_nBinary = 600;
  geometryConfiguration->min_nBinary   = 0;
  geometryConfiguration->max_nBinary   = 2000;
  geometryConfiguration->calculateDerivedHistograms = true;

  // ========================================================================================================
  Event * event = Event::getEvent();
  CollisionGeometryGenerator * collisionGeometryGenerator = new CollisionGeometryGenerator("PbPbWSGen",      geometryConfiguration, messageLevel);
  CollisionGeometry * collisionGeometry = collisionGeometryGenerator->getCollisionGeometry();
  CollisionGeometryAnalyzer  * collisionGeometryAnalyzer  = new CollisionGeometryAnalyzer ("PbPbWS-ConfigC", geometryConfiguration, collisionGeometry, messageLevel );

  eventLoop->addTask( collisionGeometryGenerator );
  eventLoop->addTask( collisionGeometryAnalyzer );
  //eventLoop->addTask( new AACollisionPythiaGenerator("AAPYTHIA",pc, collisionGeometry, event,eventFilter,particleFilter, messageLevel) );
  eventLoop->addTask( new AACollisionReader("PYTHIA",pc, event,eventFilter,particleFilter, messageLevel, collisionGeometry) );
  //eventLoop->addTask( new RadialBoostTask("PYTHIA_RADIALBOOST",rc, collisionGeometry, event, messageLevel) );
  eventLoop->addTask( new PTCorrelator("PYTHIA_PTCorrelator_HPHMPiPPiM_PbPb2760_NoBoost", ac, event, eventFilter, particleFilters, messageLevel) );
  eventLoop->run();
}



