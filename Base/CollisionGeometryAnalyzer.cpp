// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */
#include "CollisionGeometryAnalyzer.hpp"
ClassImp(CollisionGeometryAnalyzer);

CollisionGeometryAnalyzer::CollisionGeometryAnalyzer(const TString & name,
                                                     CollisionGeometryConfiguration * configuration,
                                                     CollisionGeometry * _collisionGeometry)
:
Task(name, configuration, nullptr),
collisionGeometry(_collisionGeometry)
{
  if (reportDebug()) cout << "CollisionGeometryAnalyzer::CollisionGeometryAnalyzer(...) No ops" << endl;
}

CollisionGeometryAnalyzer::~CollisionGeometryAnalyzer()
{
  if (reportDebug()) cout << "CollisionGeometryAnalyzer::~CollisionGeometryAnalyzer(...) No ops" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize generator
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CollisionGeometryAnalyzer::initialize()
{
  if (reportDebug()) cout << "CollisionGeometryAnalyzer::initialize() Started" << endl;

  CollisionGeometryConfiguration * config = (CollisionGeometryConfiguration *) getTaskConfiguration();
  collisionGeometryHistograms = new CollisionGeometryHistograms(config->histoBaseName,
                                                                config,
                                                                MessageLogger::Info);

  if (reportDebug()) cout << "CollisionGeometryAnalyzer::initialize() Completed" << endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read an ampt event from file
// Copy the event into Event for convenience...
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CollisionGeometryAnalyzer::execute()
{
  if (reportDebug()) cout << "CollisionGeometryAnalyzer::execute() Started" << endl;
  collisionGeometryHistograms->fill(collisionGeometry,1.0);
  if (reportDebug()) cout << "CollisionGeometryAnalyzer::execute() Completed" << endl;
}

void CollisionGeometryAnalyzer::saveHistograms(TFile * outputFile)
{
  if (reportDebug()) cout << "CollisionGeometryAnalyzer::saveHistograms(...) Saving Event histograms to file." << endl;
  if (!outputFile)
    {
    if (reportError()) cout << "CollisionGeometryAnalyzer::saveHistograms(...) outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();
  collisionGeometryHistograms->saveHistograms(outputFile);
  if (reportDebug()) cout << "CollisionGeometryAnalyzer::saveHistograms(...) Completed." << endl;
}


