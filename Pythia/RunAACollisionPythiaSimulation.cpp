//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>

#include "ParticleAnalyzerConfiguration.hpp"
#include "Event.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "ParticleAnalyzer.hpp"
#include "AACollisionPythiaGenerator.hpp"
#include "CollisionGeometryGenerator.hpp"
#include "CollisionGeometryAnalyzer.hpp"
#include "CollisionGeometry.hpp"

int main()
{
  MessageLogger::LogLevel messageLevel = MessageLogger::Info;

  EventLoop * eventLoop = new EventLoop("RunCollisionGeometrySimulation");
  eventLoop->setNEventRequested(1000000);
  eventLoop->setNEventReported(100000);
  eventLoop->setReportLevel(messageLevel);
  eventLoop->setNEventPartialSave(-1);
  eventLoop->setPartialSave(false);
  eventLoop->setSubsampleAnalysis(false);
  // ==========================
  // Generator Section
  // ==========================
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
  pythiaOptions[nOptions++] = new TString("SoftQCD:inelastic = on");             //Setting for Minumum bias INEL

  //pythiaOptions[nOptions++] = new TString("SoftQCD:all = on");                   // Allow total sigma = elastic/SD/DD/ND
                                                                                 //pythiaOptions[nOptions++] = new TString("HardQCD:all = on");

  // hard process -- do not turh on with SoftQCD:inelastic = on
  // otherwise, there will be double counting.
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("PhaseSpace:pTHatMin = 60.");
  //pythia.readString("PhaseSpace:pTHatMax = 1000.");

  PythiaConfiguration * pc = new PythiaConfiguration(2212, /* p */
                                                     2212, /* p */
                                                     14000.0, /* energy in GeV */
                                                     nOptions,
                                                     pythiaOptions,
                                                     true,     // only pp
                                                     true,     // remove photons
                                                     10000);


  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Heavy Ion and Analysis Configuration Parameters
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ParticleAnalyzerConfiguration * ac = new ParticleAnalyzerConfiguration("PYTHIA","PYTHIA","1.0");
  ac->loadHistograms    = false;
  ac->createHistograms  = true;
  ac->scaleHistograms   = true;
  ac->calculateDerivedHistograms  = true;
  ac->saveHistograms    = true;
  ac->inputPath         = getenv("WAC_INPUT_PATH");
  ac->rootInputFileName = "";
  ac->outputPath        = getenv("WAC_OUTPUT_PATH");
  ac->rootOuputFileName =  "/Pythia";
  ac->histoBaseName     =  "pythia";

  ac->nBins_pt    = 40;
  ac->min_pt      = 0.2;
  ac->max_pt      = 2.0;
  ac->nBins_eta   = 20;
  ac->min_eta     = -2;
  ac->max_eta     = 2;
  ac->nBins_y     = 20;
  ac->min_y       = -2;
  ac->max_y       = 2;
  ac->nBins_phi   = 36;
  ac->min_phi     = 0.0;
  ac->max_phi     = 2.0*3.1415927;

  ParticleAnalyzerConfiguration * acWide = new ParticleAnalyzerConfiguration(*ac);
  ac->nBins_eta   = 120;
  ac->min_eta     = -6;
  ac->max_eta     = 6;
  ac->nBins_y     = 120;
  ac->min_y       = -6;
  ac->max_y       = 6;

  EventFilter     * eventFilter        = new EventFilter(EventFilter::MinBias,0.0,0.0);
  int nParticleFilters = 12;
  ParticleFilter  *  particleFilter  = new ParticleFilter(ParticleFilter::Hadron, ParticleFilter::Charged,  acWide->min_pt,acWide->max_pt,acWide->min_eta,acWide->max_eta, acWide->min_y,acWide->max_y);
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
  geometryConfiguration->maxB = 21.0;

  // NucleusGenerator::WoodsSaxon, 6.62 , 0.546, 0.0, 11000,0.0,11.0);
  // NucleusGenerator::WoodsSaxon, 6.62 , 0.546, 0.0, 11000,0.0,11.0);
  geometryConfiguration->aParA = 6.62;
  geometryConfiguration->aParB = 0.546;
  geometryConfiguration->aParC = 0.0;
  geometryConfiguration->aNR   = 11000;
  geometryConfiguration->aMinR = 0.0;
  geometryConfiguration->aMaxR = 11.0;

  geometryConfiguration->bParA = 6.62;
  geometryConfiguration->bParB = 0.546;
  geometryConfiguration->bParC = 0.0;
  geometryConfiguration->bNR   = 11000;
  geometryConfiguration->bMinR = 0.0;
  geometryConfiguration->bMaxR = 11.0;


  geometryConfiguration->nnCrossSection = 6.9;  // in fm^2 -- Config C
  geometryConfiguration->nBins_b = 120;
  geometryConfiguration->min_b   = 0.0;
  geometryConfiguration->max_b   = 24.0;
  geometryConfiguration->nBins_nPart = 100;
  geometryConfiguration->min_nPart   = 0;
  geometryConfiguration->max_nPart   = 500;
  geometryConfiguration->nBins_nBinary = 600;
  geometryConfiguration->min_nBinary   = 0;
  geometryConfiguration->max_nBinary   = 3000;
  geometryConfiguration->calculateDerivedHistograms = true;

  // ========================================================================================================
  Event * event = Event::getEvent();
  CollisionGeometryGenerator * collisionGeometryGenerator = new CollisionGeometryGenerator("PbPbWSGen",      geometryConfiguration, messageLevel);
  CollisionGeometry * collisionGeometry = collisionGeometryGenerator->getCollisionGeometry();
  CollisionGeometryAnalyzer  * collisionGeometryAnalyzer  = new CollisionGeometryAnalyzer ("PbPbWS-ConfigC", geometryConfiguration, collisionGeometry, messageLevel );
  eventLoop->addTask( collisionGeometryGenerator );
  eventLoop->addTask( collisionGeometryAnalyzer );
  eventLoop->addTask( new AACollisionPythiaGenerator("AAPYTHIA",pc, collisionGeometry, event,eventFilter,particleFilter, messageLevel) );
  eventLoop->addTask( new ParticleAnalyzer("Narrow", ac, event, eventFilter, nParticleFilters, particleFilters, messageLevel) );
  eventLoop->addTask( new ParticleAnalyzer("Wide",   acWide, event, eventFilter, nParticleFilters, particleFiltersWide, messageLevel) );
  eventLoop->run();
}



