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
#include "AnalysisConfiguration.hpp"


ClassImp(NuDynTask);


//////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////
NuDynTask::NuDynTask(const TString &  name,
                     TaskConfiguration * configuration,
                     Event * event,
                     EventFilter * ef,
                     ParticleFilter * pf1,
                     ParticleFilter * pf2,
                     ParticleFilter * pf3,
                     ParticleFilter * pf4)
:
Task(name,configuration,event),
nuDynHistos(NULL),
nuDynDerivedHistos(NULL),
eventFilter(ef),
particleFilter1(pf1),
particleFilter2(pf2),
particleFilter3(pf3),
particleFilter4(pf4),
partName1("U"),
partName2("U"),
partName3("U"),
partName4("U")
{
  if (reportDebug())  cout << "NuDynTask::CTOR(...) Started." << endl;

  if (!eventFilter)
    {
    if (reportError()) cout << "NuDynTask::CTOR(...) eventFilter is null pointer." << endl;
    postTaskError();
    return;
    }

  if (!particleFilter1)
    {
    if (reportError()) cout << "NuDynTask::CTOR(...) particleFilter1 is null pointer." << endl;
    postTaskError();
    return;
    }
  if (!particleFilter2)
    {
    if (reportError()) cout << "NuDynTask::CTOR(...) particleFilter2 is null pointer." << endl;
    postTaskError();
    return;
    }
  if (!particleFilter3)
     {
     if (reportError()) cout << "NuDynTask::CTOR(...) particleFilter3 is null pointer." << endl;
     postTaskError();
     return;
     }
  if (!particleFilter4)
     {
     if (reportError()) cout << "NuDynTask::CTOR(...) particleFilter4 is null pointer." << endl;
     postTaskError();
     return;
     }

  TString newName = getName();
  newName += "_";
  newName += eventFilter->getName();
  setName(newName);
  partName1 = particleFilter1->getName();
  partName2 = particleFilter2->getName();
  partName3 = particleFilter3->getName();
  partName4 = particleFilter4->getName();
}

//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
NuDynTask::~NuDynTask()
{
  if (reportDebug())  cout << "NuDynTask::DTOR(...) Started" << endl;
  if (nuDynHistos != NULL) delete nuDynHistos;
  if (nuDynDerivedHistos != NULL) delete nuDynDerivedHistos;
  if (reportDebug())  cout << "NuDynTask::DTOR(...) Completed" << endl;
}


void NuDynTask::createHistograms()
{
  if (reportDebug())  cout << "NuDynTask::initialize(...) started"<< endl;
  AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString histoName;
  histoName = partName1;
  histoName += partName2;
  histoName += partName3;
  histoName += partName4;
  nuDynHistos = new NuDynHistos(histoName,ac,debugLevel);
  if (ac->calculateDerivedHistograms)
    {
    nuDynDerivedHistos = new NuDynDerivedHistos(histoName,ac,debugLevel);
    }
  if (reportDebug())  cout << "NuDynTask::createHistograms(...) completed"<< endl;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
void NuDynTask::loadHistograms(TFile * inputFile)
{
  if (reportDebug())  cout << "NuDynTask::loadHistograms(...) Starting." << endl;
  /* first load the number of events as from the  cumulated parameter */
  TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get("NoOfEvents");
  eventsProcessed = par->GetVal();
  delete par;
  AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString histoName;
  histoName = partName1;
  histoName += partName2;
  histoName += partName3;
  histoName += partName4;
  nuDynHistos = new NuDynHistos(inputFile,histoName,ac,debugLevel);
  if (ac->calculateDerivedHistograms)
    {
    nuDynDerivedHistos = new NuDynDerivedHistos(inputFile,histoName,ac,debugLevel);
    }
  if (reportDebug())  cout << "NuDynTask::loadHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void NuDynTask::saveHistograms(TFile * outputFile)
{
  if (reportDebug()) cout << "NuDynTask::saveHistograms(...) Saving Event histograms to file." << endl;
  if (!outputFile)
    {
    if (reportError()) cout << "NuDynTask::saveHistograms(...) outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();

  /* first save the number of events as a cumulated parameter */
  TParameter<Long64_t>("NoOfEvents",eventsProcessed,'+').Write();
  AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
  nuDynHistos->saveHistograms(outputFile);
  if (ac->calculateDerivedHistograms)
    {
    nuDynDerivedHistos->saveHistograms(outputFile);
    }
  if (reportDebug()) cout << "NuDynTask::saveHistograms(...) Completed." << endl;
}

void NuDynTask::execute()
{
  //if (reportDebug())  cout << "NuDynTask::execute(...) Starting" << endl;
  if (event != NULL)
    {
    if (reportDebug()) cout << "NuDynTask::execute(...) analyzing " << event->nParticles << " particles" << endl;
    }
  else
    {
    if (reportError()) cout << "NuDynTask::execute(...) event pointer is NULL. Abort." << endl;
    postTaskError();
    return;
    }

  //if (reportDebug()) cout <<"NuDynTask::analyze(...) check if event is acceptable" << endl;
  if (!eventFilter->accept(*event)) return;
  //if (reportDebug()) cout <<"NuDynTask::analyze(...) acceptable event" << endl;

  AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
  if (!ac)
    {
    if (reportError()) cout << "NuDynTask::execute(...) analysisConfiguration null pointer" << endl;
    postTaskError();
    return;
    }

  bool accept1;
  bool accept2;
  bool accept3;
  bool accept4;
  double n[4];
  n[0] = 0;
  n[1] = 0;
  n[2] = 0;
  n[3] = 0;

  for (int iParticle=0; iParticle<event->nParticles; iParticle++)
    {
    Particle & particle = * event->getParticleAt(iParticle);
    //if (reportDebug())  particle.printProperties(cout);
    accept1 = particleFilter1->accept(particle);
    accept2 = particleFilter2->accept(particle);
    accept3 = particleFilter3->accept(particle);
    accept4 = particleFilter4->accept(particle);
    //    if (reportDebug())  cout << "  accept1:" << accept1<< endl;
    //    if (reportDebug())  cout << "  accept2:" << accept2<< endl;
    if (accept1)  n[0]++;
    if (accept2)  n[1]++;
    if (accept3)  n[2]++;
    if (accept4)  n[3]++;
    }
  //cout << " mult:" << event->multiplicity << "   cent:" << event->centrality << " n1:" << n1 << " n2:" << n2 << endl;
  nuDynHistos->fill(event->multiplicity,event->centrality,n,1.0);
  eventsProcessed++;
  if (reportDebug()) cout << "NuDyn::execute() Completed" << endl;
}


void NuDynTask::calculateDerivedHistograms()
{
  if (reportDebug()) cout << "NuDynTask::calculateDerivedHistograms() Starting" << endl;
  nuDynDerivedHistos->calculateDerivedHistograms(nuDynHistos);
  if (reportDebug())  cout << "NuDynTask::calculateDerivedHistograms() Completed" << endl;
}

//////////////////////////////////////////////////////////////
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
//////////////////////////////////////////////////////////////
void NuDynTask::scaleHistograms(double factor)
{
  if (reportDebug())  cout << "NuDynTask::scaleHistograms(..) Scale all primary histograms by " << factor << endl;
  nuDynHistos->scale(factor);
  if (reportDebug())  cout << "NuDynTask::scale(..) Completed"  << endl;
}
