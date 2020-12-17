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

#include "GlobalAnalyzer.hpp"

ClassImp(GlobalAnalyzer);

GlobalAnalyzer::GlobalAnalyzer(const TString &  name,
                               GlobalAnalyzerConfiguration * configuration,
                               Event * event,
                               EventFilter * _eventFilter,
                               int _nParticleFilters,
                               ParticleFilter ** _particleFilters,
                               LogLevel requiredLevel)
:
Task(name,configuration,event,requiredLevel),
eventFilter(_eventFilter),
nParticleFilters(_nParticleFilters),
particleFilters(_particleFilters),
globalHistos(nullptr),
filterNames(nullptr),
n(nullptr),
e(nullptr),
q(nullptr),
b(nullptr)
{
  if (reportStart("GlobalAnalyzer",getTaskName(),"CTOR()"))
    ;
  if (!eventFilter)
    {
    if (reportWarning("GlobalAnalyzer",getTaskName(),"CTOR()")) cout << "eventFilter is null pointer." << endl;
    postTaskWarning();
    return;
    }
  if (nParticleFilters<1)
    {
    if (reportError("GlobalAnalyzer",getTaskName(),"CTOR()")) cout << "nParticleFilters<1." << endl;
    postTaskError();
    return;
    }
  if (!particleFilters)
    {
    if (reportError("GlobalAnalyzer",getTaskName(),"CTOR()")) cout << "particleFilters is null pointer." << endl;
    postTaskError();
    return;
    }
  filterNames = new TString*[nParticleFilters];
  n = new double[nParticleFilters];
  e = new double[nParticleFilters];
  q = new double[nParticleFilters];
  b = new double[nParticleFilters];

  if (!event->nFiltered)
    {
    event->nFiltered = new double[nParticleFilters];
    event->eFiltered = new double[nParticleFilters];
    event->qFiltered = new double[nParticleFilters];
    event->bFiltered = new double[nParticleFilters];
    }

  for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  {
  if (!particleFilters[iFilter])
    {
    if (reportError("GlobalAnalyzer",getTaskName(),"CTOR()")) cout << "particleFilter[" << iFilter << "] is a null pointer." << endl;
    postTaskError();
    return;
    }
  else
    {
    filterNames[iFilter] = new TString( particleFilters[iFilter]->getLongName() );
    }
  }
  TString newName = getTaskName();
  newName += "_";
  newName += eventFilter->getName();
  setTaskName(newName);
  if (reportEnd("GlobalAnalyzer",getTaskName(),"CTOR()"))
    ;
}

//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
GlobalAnalyzer::~GlobalAnalyzer()
{
  if (reportStart("GlobalAnalyzer",getTaskName(),"DTOR()"))
    ;
  if (globalHistos!= NULL) delete globalHistos;
  if (filterNames != NULL) delete[] filterNames;
  if (n != NULL) delete[] n;
  if (e != NULL) delete[] e;
  if (q != NULL) delete[] q;
  if (b != NULL) delete[] b;
  if (reportEnd("GlobalAnalyzer",getTaskName(),"DTOR()"))
    ;
}

void GlobalAnalyzer::createHistograms()
{
  if (reportStart("GlobalAnalyzer",getTaskName(),"createHistograms()"))
    ;
  GlobalAnalyzerConfiguration * ac = (GlobalAnalyzerConfiguration *) getTaskConfiguration();
  globalHistos = new GlobalHistos(getTaskName(),ac,nParticleFilters,filterNames,getReportLevel());
  globalHistos->createHistograms();
  if (reportEnd("GlobalAnalyzer",getTaskName(),"createHistograms()"))
    ;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
void GlobalAnalyzer::loadHistograms(TFile * inputFile)
{
  if (reportStart("GlobalAnalyzer",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
  GlobalAnalyzerConfiguration * ac = (GlobalAnalyzerConfiguration *) getTaskConfiguration();
  globalHistos = new GlobalHistos(getTaskName(),ac,nParticleFilters,filterNames,getReportLevel());
  globalHistos->loadHistograms(inputFile);
  if (reportEnd("GlobalAnalyzer",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void GlobalAnalyzer::saveHistograms(TFile * outputFile)
{
  if (reportStart("GlobalAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
  if (!outputFile)
    {
    if (reportError("GlobalAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();
  globalHistos->saveHistograms(outputFile);
  if (reportEnd("GlobalAnalyzer",getTaskName(),"createHistograms()"))
    ;
}

void GlobalAnalyzer::execute()
{
  incrementEventProcessed();
  if (!eventFilter->accept(*event)) return;
  incrementEventAccepted(); // count events used to fill histograms and for scaling at the end...
  GlobalAnalyzerConfiguration * ac = (GlobalAnalyzerConfiguration *) getTaskConfiguration();

  if (ac->countParticles)
    {
    bool accept;
    for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
      {
      n[iFilter] = 0.0;
      e[iFilter] = 0.0;
      q[iFilter] = 0.0;
      b[iFilter] = 0.0;
      }
    for (int iParticle=0; iParticle<event->multiplicity; iParticle++)
      {
      Particle & particle = * event->getParticleAt(iParticle);
      for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
        {
        accept = particleFilters[iFilter]->accept(particle);
        if (accept)
          {
          n[iFilter]++;
          e[iFilter] += particle.e;
          q[iFilter] += particle.charge;
          b[iFilter] += particle.baryon;
          }
        }
      }
    globalHistos->fill(n,e,q,b,1.0);
    if (ac->setEvent)
      {
      for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
        {
        event->nFiltered[iFilter] = n[iFilter];
        event->eFiltered[iFilter] = e[iFilter];
        event->qFiltered[iFilter] = q[iFilter];
        event->bFiltered[iFilter] = b[iFilter];
        }
      }
    }
  else
    {
    globalHistos->fill(event->nFiltered,
                       event->eFiltered,
                       event->qFiltered,
                       event->bFiltered,1.0);
    }
}


// =========================================================
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
// =========================================================
void GlobalAnalyzer::scaleHistograms(double factor)
{
  if (reportInfo("GlobalAnalyzer",getTaskName(),"scaleHistograms(double factor)"))  cout << "Scale histograms by " << factor << endl;
  globalHistos->scale(factor);
  if (reportEnd("GlobalAnalyzer",getTaskName(),"scaleHistograms(double factor)"))
    ;
}


// =========================================================
// Reset histograms associated with this task
// Called after partial saves in subsample analyses.
// =========================================================
void GlobalAnalyzer::resetHistograms()
{
  if (reportInfo("GlobalAnalyzer",getTaskName(),"resetHistograms()")) cout << "Will reset histograms of all particle filters" << endl;
  globalHistos->reset();
  if (reportEnd("GlobalAnalyzer",getTaskName(),"scaleHistograms(double factor)"))
    ;
}
