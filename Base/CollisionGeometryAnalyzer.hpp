// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_CollisionGeometryAnalyzer
#define WAC_CollisionGeometryAnalyzer
#include "Task.hpp"
#include "CollisionGeometry.hpp"
#include "CollisionGeometryHistograms.hpp"
#include "CollisionGeometryConfiguration.hpp"

class CollisionGeometryAnalyzer : public Task
{
public:

  CollisionGeometryAnalyzer(const TString & name,
                            CollisionGeometryConfiguration * configuration,
                            CollisionGeometry * _collisionGeometry,
                            LogLevel requiredLevel);
  virtual ~CollisionGeometryAnalyzer();
  virtual void initialize();
  virtual void saveHistograms(TFile * outputFile);
  void execute();
  virtual void calculateDerivedHistograms();
 
  CollisionGeometry * collisionGeometry;
  CollisionGeometryHistograms * collisionGeometryHistograms;

  ClassDef(CollisionGeometryAnalyzer,0)
};

#endif /* WAC_CollisionGeometryAnalyzer */
