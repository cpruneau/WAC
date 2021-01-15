// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Two Particle Correlation Analyzer Task
 */

#include "RadialBoostTask.hpp"

ClassImp(RadialBoostTask);

RadialBoostTask::RadialBoostTask(const TString &  _name,
                                 RadialBoostConfiguration * _configuration,
                                 CollisionGeometry * _collisionGeometry,
                                 Event * _event,
                                 LogLevel requiredLevel)
:
Task(_name,_configuration,_event,requiredLevel),
collisionGeometry(_collisionGeometry)
{
  if (reportNoOps("RadialBoostTask",getTaskName(),"CTOR()"))
    ;
}

RadialBoostTask::~RadialBoostTask()
{
  // no ops
}


void RadialBoostTask::createHistograms()
{
  if (reportStart("RadialBoostTask",getTaskName(),"createHistograms()"))
    ;
  RadialBoostConfiguration * ac = (RadialBoostConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  param_b  = ac->param_b; // exponent of order 1
  param_a  = ac->param_a;
  max_beta = ac->max_beta;
  TString prefixName  = getTaskName(); prefixName += "_";
  radialBoostHistos   = new RadialBoostHistos(prefixName,ac,debugLevel);
  radialBoostHistos->createHistograms();
  if (reportEnd("RadialBoostTask",getTaskName(),"createHistograms()"))
    ;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
void RadialBoostTask::loadHistograms(TFile * inputFile)
{
  if (reportStart("RadialBoostTask",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
  RadialBoostConfiguration * ac = (RadialBoostConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName  = getTaskName(); prefixName += "_";
  radialBoostHistos   = new RadialBoostHistos(prefixName,ac,debugLevel);
  radialBoostHistos->loadHistograms(inputFile);
  if (reportEnd("RadialBoostTask",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void RadialBoostTask::saveHistograms(TFile * outputFile)
{
  if (reportStart("RadialBoostTask",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
  if (!outputFile)
    {
    if (reportError()) cout << "RadialBoostTask::saveHistograms(...) outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();
  radialBoostHistos->saveHistograms(outputFile);
  if (reportEnd("RadialBoostTask",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
}

void RadialBoostTask::execute()
{
  //if (reportStart("RadialBoostTask",getTaskName(),"execute()"));
  int nParticles = event->nParticles;
  int oldIndex = -1;
  int currentIndex = -1;
  double beta, betax, betay;
  double rx, ry, phi, r;

  for (int iParticle=0; iParticle<nParticles; iParticle++)
  {
  Particle & particle = * event->getParticleAt(iParticle);
  currentIndex = particle.getSourceIndex();
  if (currentIndex != oldIndex)
    {
    oldIndex = currentIndex;
    if (currentIndex >= collisionGeometry->nBinary)
      {
      if (reportError("RadialBoostTask",getTaskName(),"saveHistograms(TFile * outputFile)"))  cout << "Bad index: " << currentIndex << endl;
      }
    rx  = collisionGeometry->x[currentIndex];
    ry  = collisionGeometry->y[currentIndex];
    r   = sqrt(rx*rx+ry*ry);
    if (r < 1E-6)
      phi = 0.0;
    else
      phi = atan2(ry,rx);
    beta = param_a * TMath::Power(r, param_b);
    if (beta > max_beta) beta = max_beta;
    betax = beta * cos(phi);
    betay = beta * sin(phi);
    radialBoostHistos->fill(r,phi,beta,1.0);
    }
  particle.boost(betax,betay,0.0);
  }
  //if (reportEnd("RadialBoostTask",getTaskName(),"execute()"));
}


void RadialBoostTask::scaleHistograms(double factor)
{
  if (reportStart("RadialBoostTask",getTaskName(),"scaleHistograms(double factor)"))
    ;
  radialBoostHistos->scale(factor);
  if (reportEnd("RadialBoostTask",getTaskName(),"scaleHistograms(double factor)"))
    ;
}
