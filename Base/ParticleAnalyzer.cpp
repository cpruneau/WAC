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

#include "ParticleAnalyzer.hpp"

ClassImp(ParticleAnalyzer);

ParticleAnalyzer::ParticleAnalyzer(const TString &  name,
                                   ParticleAnalyzerConfiguration * configuration,
                                   Event * event,
                                   EventFilter * _eventFilter,
                                   int _nParticleFilters,
                                   ParticleFilter ** _particleFilters,
                                   LogLevel requiredLevel)
:
Task(name,configuration,event,requiredLevel),
nParticleFilters(_nParticleFilters),
eventFilter(_eventFilter),
particleFilters(_particleFilters),
particleHistos(nullptr),
partNames(nullptr)
{
  if (reportStart("ParticleAnalyzer",getTaskName(),"CTOR()"))
    ;
  particleHistos = new ParticleHistos*[nParticleFilters];
  partNames      = new TString*[nParticleFilters];
  nAccepted      = new double[nParticleFilters];
  totalEnergy    = new double[nParticleFilters];
  if (!eventFilter)
    {
    if (reportWarning("ParticleAnalyzer",getTaskName(),"CTOR()")) cout << "eventFilter is null pointer." << endl;
    postTaskWarning();
    return;
    }

  if (nParticleFilters<1)
    {
    if (reportError("ParticleAnalyzer",getTaskName(),"CTOR()")) cout << "nParticleFilters<1." << endl;
    postTaskError();
    return;
    }
  if (!particleFilters)
    {
    if (reportError("ParticleAnalyzer",getTaskName(),"CTOR()")) cout << "particleFilters is null pointer." << endl;
    postTaskError();
    return;
    }
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
 {
  if (!particleFilters[iFilter])
    {
    if (reportError("ParticleAnalyzer",getTaskName(),"CTOR()")) cout << "particleFilter[" << iFilter << "] is a null pointer." << endl;
    postTaskError();
    return;
    }
  else
    {
    partNames[iFilter] = new TString( particleFilters[iFilter]->getName() );
    }
  }
  TString newName = getTaskName();
  newName += "_";
  newName += eventFilter->getName();
  setTaskName(newName);
  if (reportEnd("ParticleAnalyzer",getTaskName(),"CTOR()"))
    ;
}

//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
ParticleAnalyzer::~ParticleAnalyzer()
{
  if (reportStart("ParticleAnalyzer",getTaskName(),"DTOR()"))
    ;
  if (particleHistos != NULL) delete[] particleHistos;
  if (partNames      != NULL) delete[] partNames;
  if (nAccepted      != NULL) delete[] nAccepted;
  if (totalEnergy    != NULL) delete[] totalEnergy;

  if (reportEnd("ParticleAnalyzer",getTaskName(),"DTOR()"))
    ;
}

void ParticleAnalyzer::createHistograms()
{
  if (reportStart("ParticleAnalyzer",getTaskName(),"createHistograms()"))
    ;
  ParticleAnalyzerConfiguration * ac = (ParticleAnalyzerConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName = getTaskName(); prefixName += "_";
  TString histoName;
  if (reportInfo("ParticleAnalyzer",getTaskName(),"createHistograms()"))  cout << "Creating histograms for nParticleFilters:" << nParticleFilters <<  endl;
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  histoName = prefixName + *partNames[iFilter];
  particleHistos[iFilter] = new ParticleHistos(histoName,ac,debugLevel);
  particleHistos[iFilter]->createHistograms();
  }
  if (reportEnd("ParticleAnalyzer",getTaskName(),"createHistograms()"))
    ;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
void ParticleAnalyzer::loadHistograms(TFile * inputFile)
{
  if (reportStart("ParticleAnalyzer",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
  ParticleAnalyzerConfiguration * ac = (ParticleAnalyzerConfiguration *) getTaskConfiguration();
  TString prefixName = getTaskName(); prefixName += "_";
  TString histoName;
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  histoName = prefixName + *partNames[iFilter];
  particleHistos[iFilter] = new ParticleHistos(histoName,ac,getReportLevel());
  particleHistos[iFilter]->loadHistograms(inputFile);
  }
  if (reportEnd("ParticleAnalyzer",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void ParticleAnalyzer::saveHistograms(TFile * outputFile)
{
  if (reportStart("ParticleAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
  if (!outputFile)
    {
    if (reportError("ParticleAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  particleHistos[iFilter]->saveHistograms(outputFile);
  }
  if (reportEnd("ParticleAnalyzer",getTaskName(),"createHistograms()"))
    ;
}

void ParticleAnalyzer::execute()
{
//  if (reportDebug("ParticleAnalyzer",getTaskName(),"execute()")) cout << " nEventProcessed: " << getNEventProcessed() << endl;
//  if (reportDebug("ParticleAnalyzer",getTaskName(),"execute()")) cout << " nEventAccepted: " << getNEventAccepted() << endl;

  incrementEventProcessed();
  if (!eventFilter->accept(*event)) return;

//  if (reportDebug("ParticleAnalyzer",getTaskName(),"execute()")) cout << " ==========  WTF ========= " << endl;

  incrementEventAccepted(); // count events used to fill histograms and for scaling at the end...
  bool accept;
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  nAccepted[iFilter] = 0;
  totalEnergy[iFilter] = 0.0;
  }

  for (int iParticle=0; iParticle<event->nParticles; iParticle++)
  {
  Particle & particle = * event->getParticleAt(iParticle);
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
    {
    accept = particleFilters[iFilter]->accept(particle);
    if (accept)
      {
      nAccepted[iFilter]++;
      totalEnergy[iFilter]  += particle.e;
      particleHistos[iFilter]->fill(particle,1.0);
      }
    }
  }
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  particleHistos[iFilter]->fillMultiplicity(nAccepted[iFilter],totalEnergy[iFilter],1.0);
  }
}


// =========================================================
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
// =========================================================
void ParticleAnalyzer::scaleHistograms(double factor)
{
  if (reportInfo("ParticleAnalyzer",getTaskName(),"scaleHistograms(double factor)"))  cout << "Scale all primary histograms by " << factor << endl;
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  particleHistos[iFilter]->scale(factor);
  }
  if (reportEnd("ParticleAnalyzer",getTaskName(),"scaleHistograms(double factor)"))
    ;
}


// =========================================================
// Reset histograms associated with this task
// Called after partial saves in subsample analyses.
// =========================================================
void ParticleAnalyzer::resetHistograms()
{
  if (reportInfo("ParticleAnalyzer",getTaskName(),"resetHistograms()")) cout << "Will reset histograms of all particle filters" << endl;
  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  particleHistos[iFilter]->reset();
  }
  if (reportEnd("ParticleAnalyzer",getTaskName(),"scaleHistograms(double factor)"))
    ;
}
