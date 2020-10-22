//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include "CollisionGeometry.hpp"
#include "CollisionGeometryConfiguration.hpp"
#include "CollisionGeometryGenerator.hpp"
#include "CollisionGeometryAnalyzer.hpp"
#include "EventLoop.hpp"

int main()
{
  cout << "<I> RunCollisionGeometrySimulation() - Starting" << endl;

//  long nEventsRequested = 100;
  long nEventsRequested = 10000000;
  int  nEventsReport    = 10000;

  CollisionGeometryConfiguration * ac = new CollisionGeometryConfiguration("CollisionGeometry","CollisionGeometry","1.0");

  // Au(197,79)+Au(197,79)
  // Type             ParA     ParB   ParC
  // HardSphere       6.5^3    0.0    0.0
  // WoodsSaxon       6.5      0.535  0.0
  // WoodsSaxonHard   6.5      0.020  0.0
  // Exponential      6.5      0.0    0.0
  // Gaussian         6.5^1/2  0.0    0.0
  // DoubleGaussian   -        -      -



// Pb(208,82)+Pb(208,82)
// Type             ParA       ParB   ParC
// HardSphere       7.1^3      0.0    0.0
// WoodsSaxon       7.1        0.535  0.0
// WoodsSaxonHard   7.1        0.020  0.0
// Exponential      7.1        0.0    0.0
// Gaussian         7.1^1/2    0.0    0.0
// DoubleGaussian   -          -      -

  ac->nProtonsA  = 82;
  ac->nNeutronsA = 208-82;
  ac->nProtonsB  = 82;
  ac->nNeutronsB = 208-82;
  ac->outputPath = "/Users/claudeapruneau/Documents/GitHub/run/GeometryStudies/";
  ac->rootOuputFileName =  "CollisionGeometry";
  ac->histoBaseName =  "geom";
  ac->minB = 0.0;
  ac->maxB = 21.0;
  // RHIC Values
//  //ac->nnCrossSection = 4.5;  // in fm^2 -- nominal value
//  //ac->nnCrossSection = 4.2;  // in fm^2 -- Config B
//  ac->nnCrossSection = 4.8;  // in fm^2 -- Config C
  // LHC Values -- ALICE. PHYSICAL REVIEW C 88, 044909 (2013)
  //ac->nnCrossSection = 6.4;  // in fm^2 -- nominal value  gives a totol PbPb xSect of 761.632
  //ac->nnCrossSection = 5.9;  // in fm^2 -- Config B
    ac->nnCrossSection = 6.9;  // in fm^2 -- Config C
  ac->nBins_b = 120;
  ac->min_b   = 0.0;
  ac->max_b   = 24.0;
  ac->nBins_nPart = 100;
  ac->min_nPart   = 0;
  ac->max_nPart   = 500;
  ac->nBins_nBinary = 600;
  ac->min_nBinary   = 0;
  ac->max_nBinary   = 3000;



  ac->calculateDerivedHistograms = true;

  CollisionGeometry * collisionGeometry = new  CollisionGeometry(ac->nProtonsA,
                                                                 ac->nNeutronsA,
                                                                 ac->nProtonsB,
                                                                 ac->nNeutronsB);
  NucleusGenerator * nucleusGeneratorA = new NucleusGenerator("PbWS1",NucleusGenerator::WoodsSaxon, 6.62 , 0.546, 0.0, 11000,0.0,11.0);
  NucleusGenerator * nucleusGeneratorB = new NucleusGenerator("PbWS2",NucleusGenerator::WoodsSaxon, 6.62 , 0.546, 0.0, 11000,0.0,11.0);
  CollisionGeometryGenerator * collisionGeometryGenerator = new CollisionGeometryGenerator("PbPbWSGen",ac, collisionGeometry, nucleusGeneratorA, nucleusGeneratorB);
  CollisionGeometryAnalyzer  * collisionGeometryAnalyzer  = new CollisionGeometryAnalyzer("PbPbWS-ConfigC", ac, collisionGeometry);

  EventLoop * eventLoop = new EventLoop();
  eventLoop->addTask( collisionGeometryGenerator );
  eventLoop->addTask( collisionGeometryAnalyzer );
  eventLoop->run(nEventsRequested,nEventsReport);

  cout << "<I> RunCollisionGeometrySimulation() - Completed" << endl;
}


