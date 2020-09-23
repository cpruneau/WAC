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
  long nEventsRequested = 100000;
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
  ac->maxB = 18.0;
  ac->nnCrossSection = 4.5;  // in fm^2

  ac->nBins_b = 100;
  ac->min_b   = 0.0;
  ac->max_b   = 18.0;
  ac->nBins_nPart = 450;
  ac->min_nPart   = 0;
  ac->max_nPart   = 450;
  ac->nBins_nBinary = 400;
  ac->min_nBinary   = 0;
  ac->max_nBinary   = 2000;
  CollisionGeometry * collisionGeometry = new  CollisionGeometry(ac->nProtonsA,
                                                                 ac->nNeutronsA,
                                                                 ac->nProtonsB,
                                                                 ac->nNeutronsB);
  NucleusGenerator * nucleusGeneratorA = new NucleusGenerator("PbWS1",NucleusGenerator::WoodsSaxon, 7.1, 0.535, 0.0, 10000,0.0,8.0);
  NucleusGenerator * nucleusGeneratorB = new NucleusGenerator("PbWS2",NucleusGenerator::WoodsSaxon, 7.1, 0.535, 0.0, 10000,0.0,8.0);
  CollisionGeometryGenerator * collisionGeometryGenerator = new CollisionGeometryGenerator("PbPbWSGen",ac, collisionGeometry, nucleusGeneratorA, nucleusGeneratorB);
  CollisionGeometryAnalyzer  * collisionGeometryAnalyzer  = new CollisionGeometryAnalyzer("PbPbWS", ac, collisionGeometry);

  EventLoop * eventLoop = new EventLoop();
  eventLoop->addTask( collisionGeometryGenerator );
  eventLoop->addTask( collisionGeometryAnalyzer );
  eventLoop->run(nEventsRequested,nEventsReport);

  cout << "<I> RunCollisionGeometrySimulation() - Completed" << endl;
}


