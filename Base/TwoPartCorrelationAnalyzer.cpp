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

#include "TwoPartCorrelationAnalyzer.hpp"
#include "AnalysisConfiguration.hpp"


ClassImp(TwoPartCorrelationAnalyzer);


//////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////
TwoPartCorrelationAnalyzer::TwoPartCorrelationAnalyzer(const TString &  name,
                                                       TaskConfiguration * configuration,
                                                       Event * event,
                                                       EventFilter * ef,
                                                       ParticleFilter * pf1,
                                                       ParticleFilter * pf2)
:
Task(name,configuration,event),
eventFilter(ef),
particleFilter1(pf1),
particleFilter2(pf2),
particle1_Histos(NULL),
particle2_Histos(NULL),
pair11_Histos(NULL),
pair22_Histos(NULL),
pair12_Histos(NULL),
pair11_DerivedHistos(NULL),
pair22_DerivedHistos(NULL),
pair12_DerivedHistos(NULL),
pair12_CIHistos(NULL),
pair12_CDHistos(NULL),
partName1("U"),
partName2("U")
{
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::CTOR(...) Started." << endl;

if (!eventFilter)
  {
   if (reportError()) cout << "TwoPartCorrelationAnalyzer::CTOR(...) eventFilter is null pointer." << endl;
   postTaskError();
   return;
   }

  if (!particleFilter1)
    {
    if (reportError()) cout << "TwoPartCorrelationAnalyzer::CTOR(...) particleFilter1 is null pointer." << endl;
    postTaskError();
    return;
    }
  if (!particleFilter2)
    {
    if (reportError()) cout << "TwoPartCorrelationAnalyzer::CTOR(...) particleFilter2 is null pointer." << endl;
    postTaskError();
    return;
    }

  TString newName = getName();
  newName += "_";
  newName += eventFilter->getName();
  newName += "_";
  setName(newName);
  partName1 = particleFilter1->getName();
  partName2 = particleFilter2->getName();

}

//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
TwoPartCorrelationAnalyzer::~TwoPartCorrelationAnalyzer()
{
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::DTOR(...) Started" << endl;
  if (particle1_Histos != NULL) delete particle1_Histos;
  if (particle2_Histos != NULL) delete particle2_Histos;
  if (pair11_Histos != NULL) delete pair11_Histos;
  if (pair22_Histos != NULL) delete pair22_Histos;
  if (pair12_Histos != NULL) delete pair12_Histos;
  if (pair11_DerivedHistos != NULL) delete pair11_DerivedHistos;
  if (pair22_DerivedHistos != NULL) delete pair22_DerivedHistos;
  if (pair12_DerivedHistos != NULL) delete pair12_DerivedHistos;
  if (pair12_CIHistos != NULL) delete pair12_CIHistos;
  if (pair12_CDHistos != NULL) delete pair12_CDHistos;
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::DTOR(...) Completed" << endl;
}


void TwoPartCorrelationAnalyzer::createHistograms()
{
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::initialize(...) started"<< endl;
  AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();

  TString prefixName = getName();
  particle1_Histos  = new ParticleHistos(prefixName+partName1,ac,debugLevel);
  particle2_Histos  = new ParticleHistos(prefixName+partName2,ac,debugLevel);
  if (ac->fillPairs)
    {
    pair11_Histos = new ParticlePairHistos(prefixName+partName1+partName1,ac,debugLevel);
    pair22_Histos = new ParticlePairHistos(prefixName+partName2+partName2,ac,debugLevel);
    pair12_Histos = new ParticlePairHistos(prefixName+partName1+partName2,ac,debugLevel);
    if (ac->calculateDerivedHistograms)
      {
      pair11_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName1,ac,debugLevel);
      pair22_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName2+partName2,ac,debugLevel);
      pair12_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName2,ac,debugLevel);
      pair12_CIHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CI",ac,debugLevel);
      pair12_CDHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CD",ac,debugLevel);
      }
    }
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::createHistograms(...) completed"<< endl;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::loadHistograms(TFile * inputFile)
{
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::loadHistograms(...) Starting." << endl;
  /* first load the number of events as from the  cumulated parameter */
  TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get("NoOfEvents");
  eventsProcessed = par->GetVal();
  delete par;
  AnalysisConfiguration * analysisConfiguration = (AnalysisConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName = getName();

  particle1_Histos         = new ParticleHistos(inputFile,prefixName+partName1, analysisConfiguration,debugLevel);
  particle2_Histos         = new ParticleHistos(inputFile,prefixName+partName2,analysisConfiguration,debugLevel);
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos        = new ParticlePairHistos(inputFile,prefixName+partName1+partName1,analysisConfiguration,debugLevel);
    pair22_Histos        = new ParticlePairHistos(inputFile,prefixName+partName2+partName2,analysisConfiguration,debugLevel);
    pair12_Histos        = new ParticlePairHistos(inputFile,prefixName+partName1+partName2,analysisConfiguration,debugLevel);
    if (analysisConfiguration->calculateDerivedHistograms)
      {
      pair11_DerivedHistos = new ParticlePairDerivedHistos(inputFile, prefixName+partName1+partName1,analysisConfiguration,debugLevel);
      pair22_DerivedHistos = new ParticlePairDerivedHistos(inputFile, prefixName+partName2+partName2,analysisConfiguration,debugLevel);
      pair12_DerivedHistos = new ParticlePairDerivedHistos(inputFile, prefixName+partName1+partName2,analysisConfiguration,debugLevel);
      pair12_CIHistos      = new ParticlePairCombinedHistos(inputFile,prefixName+partName1+partName2+"CI",analysisConfiguration,debugLevel);
      pair12_CDHistos      = new ParticlePairCombinedHistos(inputFile,prefixName+partName1+partName2+"CD",analysisConfiguration,debugLevel);
      }
    }
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::loadHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// load the base histograms from given file
//////////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::loadBaseHistograms(TFile * inputFile)
{
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::loadHistograms(...) Starting." << endl;

  /* first load the number of events as from the  cumulated parameter */
  TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get("NoOfEvents");
  eventsProcessed = par->GetVal();
  delete par;
  AnalysisConfiguration * analysisConfiguration = (AnalysisConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  
  TString prefixName = getName();
  particle1_Histos         = new ParticleHistos(inputFile,prefixName+partName1, analysisConfiguration,debugLevel);
  particle2_Histos         = new ParticleHistos(inputFile,prefixName+partName2,analysisConfiguration,debugLevel);
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos        = new ParticlePairHistos(inputFile,prefixName+partName1+partName1,analysisConfiguration,debugLevel);
    pair22_Histos        = new ParticlePairHistos(inputFile,prefixName+partName2+partName2,analysisConfiguration,debugLevel);
    pair12_Histos        = new ParticlePairHistos(inputFile,prefixName+partName1+partName2,analysisConfiguration,debugLevel);
    if (analysisConfiguration->calculateDerivedHistograms)
      {
      pair11_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName1,analysisConfiguration,debugLevel);
      pair22_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName2+partName2,analysisConfiguration,debugLevel);
      pair12_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName2,analysisConfiguration,debugLevel);
      pair12_CIHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CI",analysisConfiguration,debugLevel);
      pair12_CDHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CD",analysisConfiguration,debugLevel);
      }
    }

  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::loadHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::saveHistograms(TFile * outputFile)
{
  if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::saveHistograms(...) Saving Event histograms to file." << endl;
  if (!outputFile)
    {
    if (reportError()) cout << "TwoPartCorrelationAnalyzer::saveHistograms(...) outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();

  /* first save the number of events as a cumulated parameter */
  TParameter<Long64_t>("NoOfEvents",eventsProcessed,'+').Write();
  AnalysisConfiguration * analysisConfiguration = (AnalysisConfiguration *) getTaskConfiguration();
  if (!analysisConfiguration)
     {
     if (reportError()) cout << "TwoPartCorrelationAnalyzer::saveHistograms(...) analysisConfiguration is a null  pointer." << endl;
     postTaskError();
     return;
     }

  if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::saveHistograms(...) saving singles." << endl;

  particle1_Histos  ->saveHistograms(outputFile);
  particle2_Histos  ->saveHistograms(outputFile);

  if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::saveHistograms(...) saving pairs -- perhaps." << endl;

  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos ->saveHistograms(outputFile);
    pair22_Histos ->saveHistograms(outputFile);
    pair12_Histos ->saveHistograms(outputFile);
    if (analysisConfiguration->calculateDerivedHistograms)
      {
      if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::saveHistograms(...) saving calculated histograms." << endl;

      pair11_DerivedHistos->saveHistograms(outputFile);
      pair22_DerivedHistos->saveHistograms(outputFile);
      pair12_DerivedHistos->saveHistograms(outputFile);
      pair12_CIHistos     ->saveHistograms(outputFile);
      pair12_CDHistos     ->saveHistograms(outputFile);
      }
    }
if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::saveHistograms(...) Completed." << endl;
}


//////////////////////////////////////////////////////////////
// add histograms to an external list
//////////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::addHistogramsToExtList(TList *list, bool all)
{
  if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::addHistogramsToExtList(...) Saving Event histograms to file." << endl;

  /* first add the number of events as a cumulated parameter */
  list->Add(new TParameter<Long64_t>("NoOfEvents",eventsProcessed,'+'));
  AnalysisConfiguration * analysisConfiguration = (AnalysisConfiguration *) getTaskConfiguration();

  particle1_Histos  ->addHistogramsToExtList(list, all);
  particle2_Histos  ->addHistogramsToExtList(list, all);
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos ->addHistogramsToExtList(list, all);
    pair22_Histos ->addHistogramsToExtList(list, all);
    pair12_Histos ->addHistogramsToExtList(list, all);
    if (analysisConfiguration->calculateDerivedHistograms)
      {
      pair11_DerivedHistos->addHistogramsToExtList(list, all);
      pair22_DerivedHistos->addHistogramsToExtList(list, all);
      pair12_DerivedHistos->addHistogramsToExtList(list, all);
      pair12_CIHistos     ->addHistogramsToExtList(list, all);
      pair12_CDHistos     ->addHistogramsToExtList(list, all);
      }
    }

  if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::addHistogramsToExtList(...) Completed." << endl;
}


void TwoPartCorrelationAnalyzer::execute()
{
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::analyze(...) Starting" << endl;
  if (event != NULL)
    {
    if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::analyze(...) analyzing " << event->nParticles << " particles" << endl;
    }
  else
    {
    if (reportError()) cout << "TwoPartCorrelationAnalyzer::analyze(...) event pointer is NULL. Abort." << endl;
    postTaskError();
    return;
    }

  if (reportDebug()) cout <<"TwoPartCorrelationAnalyzer::analyze(...) check if event is acceptable" << endl;
  if (!eventFilter->accept(*event)) return;
  if (reportDebug()) cout <<"TwoPartCorrelationAnalyzer::analyze(...) acceptable event" << endl;

  AnalysisConfiguration * analysisConfiguration = (AnalysisConfiguration *) getTaskConfiguration();

if (!analysisConfiguration)
  {
  if (reportError()) cout << "TwoPartCorrelationAnalyzer::analyze(...) analysisConfiguration null pointer" << endl;
  postTaskError();
  return;
  }

  bool accept11;
  bool accept21;
  bool accept12;
  bool accept22;

  /* before filtering let's build the particle indexes to hurry up the process */
  for (int iParticle = 0; iParticle < event->nParticles; iParticle++)
    {
    if (reportDebug())cout << "TwoPartCorrelationAnalyzer::analyze(...) particle: " << iParticle << endl;
    Particle *particle = event->getParticleAt(iParticle);
    if (reportDebug()) particle->printProperties(cout);
    particle->ixEtaPhi = analysisConfiguration->getIxEtaPhi(particle->eta,particle->phi);
    particle->ixYPhi = analysisConfiguration->getIxYPhi(particle->y,particle->phi);
    }

  for (int iParticle1=0; iParticle1<event->nParticles; iParticle1++)
    {
    Particle & particle1 = * event->getParticleAt(iParticle1);
    accept11 = particleFilter1->accept(particle1);
    accept21 = particleFilter2->accept(particle1);
    if (reportDebug())  cout << "  accept11:" << accept11<< endl;
    if (reportDebug())  cout << "  accept21:" << accept21<< endl;

    if (accept11)  particle1_Histos->fill(particle1, 1.0);
    if (accept21)  particle2_Histos->fill(particle1, 1.0);
    if (analysisConfiguration->fillPairs)
      {
      for (int iParticle2=0; iParticle2<event->nParticles; iParticle2++)
        {
        if (iParticle1==iParticle2) continue;
        Particle & particle2 = * event->getParticleAt(iParticle2);
        accept12 = particleFilter1->accept(particle2);
        accept22 = particleFilter2->accept(particle2);
        //if (reportDebug())  cout << "  accept12:" << accept12<< endl;
        //if (reportDebug())  cout << "  accept22:" << accept22<< endl;
        if (accept11 && accept12)  pair11_Histos->fill(particle1, particle2, 1.0);
        if (accept21 && accept22)  pair22_Histos->fill(particle1, particle2, 1.0);
        if (accept11 && accept22)  pair12_Histos->fill(particle1, particle2, 1.0);
        }
      }
    }
  eventsProcessed++;
  if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::execute() Completed" << endl;
}


//////////////////////////////////////////////////////////////
// calculate Derived Histograms
//////////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::calculateDerivedHistograms()
{
  if (reportDebug()) cout << "TwoPartCorrelationAnalyzer::calculateDerivedHistograms() Starting" << endl;
  AnalysisConfiguration * analysisConfiguration = (AnalysisConfiguration *) getTaskConfiguration();
  particle1_Histos->completeFill();
  particle2_Histos->completeFill();
  if (analysisConfiguration->fillPairs)
    {
    particle1_Histos->calculateAverages();
    particle2_Histos->calculateAverages();
    pair11_Histos->completeFill();
    pair22_Histos->completeFill();
    pair12_Histos->completeFill();
    pair11_DerivedHistos->calculateDerivedHistograms(particle1_Histos,particle1_Histos,pair11_Histos,analysisConfiguration->binCorrPP);
    pair22_DerivedHistos->calculateDerivedHistograms(particle2_Histos,particle2_Histos,pair22_Histos,analysisConfiguration->binCorrMM);
    pair12_DerivedHistos->calculateDerivedHistograms(particle1_Histos,particle2_Histos,pair12_Histos,analysisConfiguration->binCorrPM);
    pair12_CIHistos->calculate(pair11_DerivedHistos,pair22_DerivedHistos,pair12_DerivedHistos, 0.25, 0.25,0.5);
    pair12_CDHistos->calculate(pair11_DerivedHistos,pair22_DerivedHistos,pair12_DerivedHistos,-0.25,-0.25,0.5);
    }
  else {
    if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::calculateDerivedHistograms() Nothing to do, fill pairs not ordered" << endl;
  }
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::calculateDerivedHistograms() Completed" << endl;
}



//////////////////////////////////////////////////////////////
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
//////////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::scaleHistograms(double factor)
{
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::scaleHistograms(..) Scale all primary histograms by " << factor << endl;
  particle1_Histos->scale(factor);
  particle2_Histos->scale(factor);
  AnalysisConfiguration * analysisConfiguration = (AnalysisConfiguration *) getTaskConfiguration();
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos->scale(factor);
    pair22_Histos->scale(factor);
    pair12_Histos->scale(factor);
    }
  if (reportDebug())  cout << "TwoPartCorrelationAnalyzer::scale(..) Completed"  << endl;
}



