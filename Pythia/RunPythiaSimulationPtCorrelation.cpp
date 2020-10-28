//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
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
#include "PTCorrelator.hpp"
#include "AACollisionGenerator.hpp"
#include "CollisionGeometryGenerator.hpp"
#include "CollisionGeometry.hpp"
#include "NucleusGenerator.hpp"
#include "HeavyIonConfiguration.hpp"


int main()
{
  auto start = chrono::high_resolution_clock::now(); 

  cout << "<INFO> PYTHIA Model Analysis - Starting" << endl;

//  long nEventsRequested = 100;
  long nEventsRequested = 1000;
  int  nEventsReport    = 100000;

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Heavy Ion and Analysis Configuration Parameters
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  HeavyIonConfiguration * ac = new HeavyIonConfiguration("PYTHIA","PYTHIA","1.0");
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
  ac->outputPath = getenv("OUTPUT_PATH");            // check this 
  ac->rootOuputFileName =  "/Pythia6";                // and this
  ac->histoBaseName =  "pythia";

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

  ac->hardBoost      = false;
  ac->nCollisionsMax = 1;
  ac->param_a        = ac->hardBoost? 0.1 : 0.05;
  ac->param_b        = 1;
  ac->maxOrder       = 4; // order <= particleFilters1 length
  ac->totEvents      = nEventsRequested;




  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Particle and Event Filter Parameters
  ////////////////////////////////////////////////////////////////////////////////////////////////////




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


  ParticleFilter * particleFilters1[8];
  particleFilters1[0] = particleFilter_HP;
  particleFilters1[1] = particleFilter_HM;
  particleFilters1[2] = particleFilter_PiP;
  particleFilters1[3] = particleFilter_PiM;
  particleFilters1[4] = particleFilter_KP;
  particleFilters1[5] = particleFilter_KM;
  particleFilters1[6] = particleFilter_PP;
  particleFilters1[7] = particleFilter_PM;




  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Collision Geometry Configuration Parameters
  ////////////////////////////////////////////////////////////////////////////////////////////////////




  CollisionGeometryConfiguration * ac1 = new CollisionGeometryConfiguration("PYTHIACollisionGeometry","PYTHIACollisionGeometry","1.0");

  ac1->nProtonsA  = 82;
  ac1->nNeutronsA = 208-82;
  ac1->nProtonsB  = 82;
  ac1->nNeutronsB = 208-82;
  ac1->outputPath = getenv("OUTPUT_PATH"); 
  ac1->rootOuputFileName =  "/CollisionGeometry";
  ac1->histoBaseName =  "geom";
  ac1->minB = 0.0;
  ac1->maxB = 18.0;
  ac1->nnCrossSection = 4.5;  // in fm^2

  ac1->nBins_b = 100;
  ac1->min_b   = 0.0;
  ac1->max_b   = 18.0;
  ac1->nBins_nPart = 450;
  ac1->min_nPart   = 0;
  ac1->max_nPart   = 450;
  ac1->nBins_nBinary = 400;
  ac1->min_nBinary   = 0;
  ac1->max_nBinary   = 2000;

  CollisionGeometry * collisionGeometry = new CollisionGeometry(ac1->nProtonsA,ac1->nNeutronsA, ac1->nProtonsB, ac1->nNeutronsB); //Pb(208) has 82 protons and 126 neutrons 
  NucleusGenerator * nucGenA = new NucleusGenerator("PYTHIA_PbPbNucleusGeneratorA", NucleusGenerator::WoodsSaxon, 7.1, 0.535, 0.0, 10000,0.0,8.0);
  NucleusGenerator * nucGenB = new NucleusGenerator("PYTHIA_PbPbNucleusGeneratorB", NucleusGenerator::WoodsSaxon, 7.1, 0.535, 0.0, 10000,0.0,8.0);




  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Run Events
  ////////////////////////////////////////////////////////////////////////////////////////////////////





  Event * event = Event::getEvent();
  EventLoop * eventLoop = new EventLoop();

  //eventLoop->addTask( new PythiaEventGenerator("PYTHIA",0, event,eventFilter,particleFilter) );
  eventLoop->addTask( new CollisionGeometryGenerator("PYTHIA_PbPbCollisionGeometryGenerator",ac1, collisionGeometry,nucGenA,nucGenB) );
  eventLoop->addTask( new AACollisionGenerator("PYTHIA_PbPbEventGenerator",ac, event,eventFilter,particleFilter, collisionGeometry) );
  eventLoop->addTask( new PTCorrelator("PYTHIA_PTCorrelator_HPHMPiPPiM", ac, event, eventFilter, particleFilters1) ); // Note: make sure all filters are distinct

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




  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Timer
  ////////////////////////////////////////////////////////////////////////////////////////////////////




  auto stop = chrono::high_resolution_clock::now(); 
  auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
  int hours = (int)(duration.count()/3600);
  int minutes = (int)((duration.count() - 3600 * hours)/60);
  double seconds = duration.count() - 60 * minutes - 3600 * hours;
  cout << "<INFO> Total Time elapsed "<< (hours) << ":" << (minutes ) << ":" << seconds << endl; 

}



