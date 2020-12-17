// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_RadialBoostTask
#define WAC_RadialBoostTask
#include <TParameter.h>
#include "TFile.h"
#include "Task.hpp"
#include "Event.hpp"
#include "CollisionGeometry.hpp"
#include "RadialBoostHistos.hpp"

class RadialBoostTask : public Task
{
public:

  //RadialBoostTask();

  RadialBoostTask(const TString &  name,
                  RadialBoostConfiguration * configuration,
                  CollisionGeometry * collisionGeometry,
                  Event * event,
                  LogLevel requiredLevel);
  virtual ~RadialBoostTask();
  virtual void execute();
  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void saveHistograms(TFile * outputFile);
  virtual void scaleHistograms(double factor);
  CollisionGeometry * collisionGeometry;
  RadialBoostHistos * radialBoostHistos;

  double param_b;
  double param_a;
  double max_beta;

  ClassDef(RadialBoostTask,0)
};


#endif /* WAC_RadialBoostTask */
