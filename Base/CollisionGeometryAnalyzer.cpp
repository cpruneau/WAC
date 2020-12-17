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
                                                     CollisionGeometry * _collisionGeometry,
                                                     LogLevel requiredLevel)
:
Task(name, configuration, nullptr,requiredLevel),
collisionGeometry(_collisionGeometry)
{
  // no ops
}

CollisionGeometryAnalyzer::~CollisionGeometryAnalyzer()
{
  // no ops
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize generator
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CollisionGeometryAnalyzer::initialize()
{
  if (reportStart("ParticleAnalyzer",getTaskName(),"initialize()"))
    ;
  CollisionGeometryConfiguration * config = (CollisionGeometryConfiguration *) getTaskConfiguration();

  if (reportInfo("CollisionGeometryAnalyzer",getTaskName(),"initialize()"))
      {
      cout << endl;
      config->printConfiguration(cout);
      }
  collisionGeometryHistograms = new CollisionGeometryHistograms(config->histoBaseName,
                                                                config,
                                                                getReportLevel());
  collisionGeometryHistograms->createHistograms();
  if (reportEnd("ParticleAnalyzer",getTaskName(),"initialize()"))
    ;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read an ampt event from file
// Copy the event into Event for convenience...
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CollisionGeometryAnalyzer::execute()
{
  //if (reportStart("ParticleAnalyzer",getTaskName(),"execute()")) ;
  collisionGeometryHistograms->fill(collisionGeometry,1.0);
  //if (reportEnd("ParticleAnalyzer",getTaskName(),"execute()")) ;
}

void CollisionGeometryAnalyzer::saveHistograms(TFile * outputFile)
{
  if (reportStart("ParticleAnalyzer",getTaskName(),"execute()"))
    ;
  if (!outputFile)
    {
    if (reportError("ParticleAnalyzer",getTaskName(),"execute()")) cout << "outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();
  collisionGeometryHistograms->saveHistograms(outputFile);
  if (reportEnd("ParticleAnalyzer",getTaskName(),"execute()"))
    ;
}


void CollisionGeometryAnalyzer::calculateDerivedHistograms()
{
  collisionGeometryHistograms->calculateDerivedHistograms();
}
