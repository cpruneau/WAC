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
#include "NuDynTask.hpp"

ClassImp(NuDynTask);

NuDynTask::NuDynTask(const TString &  _name,
                     NuDynConfiguration * _taskConfiguration,
                     Event          * _event,
                     EventFilter    * _eventFilter,
                     ParticleFilter * particleFilter1,
                     ParticleFilter * particleFilter2,
                     ParticleFilter * particleFilter3,
                     ParticleFilter * particleFilter4,
                     LogLevel selectedLevel)
:
Task(_name,_taskConfiguration,_event,selectedLevel),
nuDynHistos       (nullptr),
nuDynDerivedHistos(nullptr),
eventFilter       (_eventFilter),
nParticleFilters  (4),
particleFilters   (new ParticleFilter*[4]),
partNames         (nullptr)
{
  if (reportStart("NuDynTask",getTaskName(),"CTOR"))
    ;
  if (!eventFilter)
    {
    if (reportError("NuDynTask",getTaskName(),"CTOR")) cout << "Given eventFilter is a null pointer." << endl;
    postTaskError();
    return;
    }
  particleFilters[0] = particleFilter1;
  particleFilters[1] = particleFilter2;
  particleFilters[2] = particleFilter3;
  particleFilters[3] = particleFilter4;
  for (int iParticleFilter=0; iParticleFilter<nParticleFilters; iParticleFilter++)
  {
  if (!particleFilters[iParticleFilter])
    {
    if (reportError("NuDynTask",getTaskName(),"CTOR")) cout << "Given particleFilter : " <<  iParticleFilter << " is null pointer." << endl;
    postTaskError();
    return;
    }
  partNames[iParticleFilter] = new TString( particleFilters[iParticleFilter]->getName() );
  }
  createIdentical();
  if (reportEnd("NuDynTask",getTaskName(),"CTOR"))
    ;
}


//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
NuDynTask::~NuDynTask()
{
  if (reportStart("NuDynTask",getTaskName(),"DTOR"))
    ;
  if (nuDynHistos != NULL) delete nuDynHistos;
  if (nuDynDerivedHistos != NULL) delete nuDynDerivedHistos;
  if (partNames != NULL) delete partNames;
  if (reportEnd("NuDynTask",getTaskName(),"DTOR"))
    ;
}


void NuDynTask::createHistograms()
{
  if (reportStart("NuDynTask",getTaskName(),"createHistograms()"))
    ;
  NuDynConfiguration * ac = (NuDynConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName = getTaskName(); prefixName += "_";
  TString histoName = prefixName;
  for (int iParticleFilter=0; iParticleFilter<nParticleFilters; iParticleFilter++)
  {
  histoName +=  * partNames[iParticleFilter];
  }
  nuDynHistos = new NuDynHistos(histoName,identical,ac,debugLevel);
  nuDynHistos->createHistograms();
  if (ac->calculateDerivedHistograms)
    {
    nuDynDerivedHistos = new NuDynDerivedHistos(histoName,ac,debugLevel);
    nuDynDerivedHistos->createHistograms();
    }
  if (reportEnd("NuDynTask",getTaskName(),"createHistograms()"))
    ;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
void NuDynTask::loadHistograms(TFile * inputFile)
{
  if (reportStart("NuDynTask",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
  NuDynConfiguration * ac = (NuDynConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName = getTaskName(); prefixName += "_";
  TString histoName = prefixName;
  for (int iParticleFilter=0; iParticleFilter<nParticleFilters; iParticleFilter++)
  {
  histoName +=  * partNames[iParticleFilter];
  }

  nuDynHistos = new NuDynHistos(histoName,identical,ac,debugLevel);
  nuDynHistos->loadHistograms(inputFile);
  if (ac->calculateDerivedHistograms)
    {
    nuDynDerivedHistos = new NuDynDerivedHistos(histoName,ac,debugLevel);
    nuDynDerivedHistos->loadHistograms(inputFile);
    }
  if (reportEnd("NuDynTask",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void NuDynTask::saveHistograms(TFile * outputFile)
{
  if (reportStart("NuDynTask",getTaskName(),"saveHistograms(TFile * inputFile)"))
    ;
  if (!outputFile)
    {
    if (reportError("NuDynTask",getTaskName(),"saveHistograms(TFile * inputFile)")) cout << "Given inputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();
  nuDynHistos->saveHistograms(outputFile);
  NuDynConfiguration * ac = (NuDynConfiguration *) getTaskConfiguration();
  if (ac->calculateDerivedHistograms)
    {
    nuDynDerivedHistos->saveHistograms(outputFile);
    }
  if (reportEnd("NuDynTask",getTaskName(),"saveHistograms(TFile * inputFile)"))
    ;
}

void NuDynTask::execute()
{
  incrementEventProcessed();
  if (!eventFilter->accept(*event)) return;
  incrementEventAccepted(); // count events used to fill histograms and for scaling at the end...
  NuDynConfiguration & ac = *(NuDynConfiguration *) getTaskConfiguration();
  double nAccepted[4];
  for (int iParticleFilter=0; iParticleFilter<nParticleFilters; iParticleFilter++) nAccepted[iParticleFilter] = 0;

  for (int iParticle=0; iParticle<event->nParticles; iParticle++)
    {
    Particle & particle = * event->getParticleAt(iParticle);
    for (int iParticleFilter=0; iParticleFilter<nParticleFilters; iParticleFilter++)
      {
      if (particleFilters[iParticleFilter]->accept(particle)) nAccepted[iParticleFilter]++;
      }
    }
  switch ( ac.multiplicityType )
    {
      case NuDynConfiguration::Centrality:           nuDynHistos->fill(event->centrality,nAccepted,1.0); break;
      case NuDynConfiguration::TotalMultiplicity:    nuDynHistos->fill(event->multiplicity,nAccepted,1.0); break;
      case NuDynConfiguration::AcceptedMultiplicity: nuDynHistos->fill(event->multiplicity,nAccepted,1.0); break;
    }
}


void NuDynTask::calculateDerivedHistograms()
{
  if (reportStart("NuDynTask",getTaskName(),"calculateDerivedHistograms()"))
    ;
  nuDynDerivedHistos->calculateDerivedHistograms(nuDynHistos);
  if (reportEnd("NuDynTask",getTaskName(),"calculateDerivedHistograms()"))
    ;
}

//////////////////////////////////////////////////////////////
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
//////////////////////////////////////////////////////////////
void NuDynTask::scaleHistograms(double factor)
{
  if (reportInfo("NuDynTask",getTaskName(),"scaleHistograms(double factor)"))   cout << "Scale all primary histograms by " << factor << endl;
  nuDynHistos->scale(factor);
  if (reportEnd("NuDynTask",getTaskName(),"scaleHistograms(double factor)"))
    ;
}

void NuDynTask::resetHistograms()
{
  nuDynHistos->reset();
  nuDynDerivedHistos->reset();
}



void NuDynTask::createIdentical()
{
  identical[0] = 1; //(particleFilters[0]==particleFilters[0]) ? 1 : 0;
  identical[1] = (particleFilters[0]==particleFilters[1]) ? 1 : 0;
  identical[2] = (particleFilters[0]==particleFilters[2]) ? 1 : 0;
  identical[3] = (particleFilters[0]==particleFilters[3]) ? 1 : 0;
  identical[4] = (particleFilters[1]==particleFilters[0]) ? 1 : 0;
  identical[5] = 1; //(particleFilters[1]==particleFilters[1]) ? 1 : 0;
  identical[6] = (particleFilters[1]==particleFilters[2]) ? 1 : 0;
  identical[7] = (particleFilters[1]==particleFilters[3]) ? 1 : 0;
  identical[8] = (particleFilters[2]==particleFilters[0]) ? 1 : 0;
  identical[9] = (particleFilters[2]==particleFilters[1]) ? 1 : 0;
  identical[10] = 1; //(particleFilters[2]==particleFilters[2]) ? 1 : 0;
  identical[11] = (particleFilters[2]==particleFilters[3]) ? 1 : 0;
  identical[12] = (particleFilters[3]==particleFilters[0]) ? 1 : 0;
  identical[13] = (particleFilters[3]==particleFilters[1]) ? 1 : 0;
  identical[14] = (particleFilters[3]==particleFilters[2]) ? 1 : 0;
  identical[15] = 1; //(particleFilters[3]==particleFilters[3]) ? 1 : 0;
}
