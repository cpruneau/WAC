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

ClassImp(TwoPartCorrelationAnalyzer);

TwoPartCorrelationAnalyzer::TwoPartCorrelationAnalyzer(const TString &  name,
                                                       ParticlePairAnalyzerConfiguration * configuration,
                                                       Event * event,
                                                       EventFilter * ef,
                                                       ParticleFilter * pf1,
                                                       ParticleFilter * pf2,
                                                       LogLevel selectedLevel)
:
Task(name,configuration,event,selectedLevel),
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
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"CTOR"))
    ;
  if (!eventFilter)
  {
   if (reportError("TwoPartCorrelationAnalyzer",getTaskName(),"CTOR")) cout << "Given eventFilter is a null pointer." << endl;
   postTaskError();
   return;
   }

  if (!particleFilter1)
    {
    if (reportError("TwoPartCorrelationAnalyzer",getTaskName(),"CTOR")) cout << "Given particleFilter1 is a null pointer." << endl;
    postTaskError();
    return;
    }
  if (!particleFilter2)
    {
    if (reportError("TwoPartCorrelationAnalyzer",getTaskName(),"CTOR")) cout << "Given particleFilter2 is a null pointer." << endl;
    postTaskError();
    return;
    }
  TString newName = getTaskName();
  newName += "_";
  newName += eventFilter->getName();
  newName += "_";
  setTaskName(newName);
  partName1 = particleFilter1->getName();
  partName2 = particleFilter2->getName();
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"CTOR"))
    ;
}

//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
TwoPartCorrelationAnalyzer::~TwoPartCorrelationAnalyzer()
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"DTOR"))
    ;
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
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"DTOR"))
    ;
}


void TwoPartCorrelationAnalyzer::createHistograms()
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"createHistograms()"))
    ;
  ParticlePairAnalyzerConfiguration * ac = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName = getTaskName();
  particle1_Histos  = new ParticleHistos(prefixName+partName1,ac,debugLevel);
  particle2_Histos  = new ParticleHistos(prefixName+partName2,ac,debugLevel);
  particle1_Histos->createHistograms();
  particle2_Histos->createHistograms();
  if (ac->fillPairs)
    {
    pair11_Histos = new ParticlePairHistos(prefixName+partName1+partName1,ac,debugLevel);
    pair22_Histos = new ParticlePairHistos(prefixName+partName2+partName2,ac,debugLevel);
    pair12_Histos = new ParticlePairHistos(prefixName+partName1+partName2,ac,debugLevel);
    pair11_Histos->createHistograms();
    pair22_Histos->createHistograms();
    pair12_Histos->createHistograms();
    if (ac->calculateDerivedHistograms)
      {
      pair11_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName1,ac,debugLevel);
      pair22_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName2+partName2,ac,debugLevel);
      pair12_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName2,ac,debugLevel);
      pair12_CIHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CI",ac,debugLevel);
      pair12_CDHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CD",ac,debugLevel);
      pair11_DerivedHistos->createHistograms();
      pair22_DerivedHistos->createHistograms();
      pair12_DerivedHistos->createHistograms();
      pair12_CIHistos->createHistograms();
      pair12_CDHistos->createHistograms();
      }
    }
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"createHistograms()"))
    ;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::loadHistograms(TFile * inputFile)
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName  = getTaskName();
  particle1_Histos    = new ParticleHistos(prefixName+partName1,analysisConfiguration,debugLevel);
  particle2_Histos    = new ParticleHistos(prefixName+partName2,analysisConfiguration,debugLevel);
  particle1_Histos->loadHistograms(inputFile);
  particle2_Histos->loadHistograms(inputFile);
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos     = new ParticlePairHistos(prefixName+partName1+partName1,analysisConfiguration,debugLevel);
    pair22_Histos     = new ParticlePairHistos(prefixName+partName2+partName2,analysisConfiguration,debugLevel);
    pair12_Histos     = new ParticlePairHistos(prefixName+partName1+partName2,analysisConfiguration,debugLevel);
    pair11_Histos->loadHistograms(inputFile);
    pair22_Histos->loadHistograms(inputFile);
    pair12_Histos->loadHistograms(inputFile);
    if (analysisConfiguration->calculateDerivedHistograms)
      {
      pair11_DerivedHistos = new ParticlePairDerivedHistos(prefixName+partName1+partName1,analysisConfiguration,debugLevel);
      pair22_DerivedHistos = new ParticlePairDerivedHistos(prefixName+partName2+partName2,analysisConfiguration,debugLevel);
      pair12_DerivedHistos = new ParticlePairDerivedHistos(prefixName+partName1+partName2,analysisConfiguration,debugLevel);
      pair12_CIHistos      = new ParticlePairCombinedHistos(prefixName+partName1+partName2+"CI",analysisConfiguration,debugLevel);
      pair12_CDHistos      = new ParticlePairCombinedHistos(prefixName+partName1+partName2+"CD",analysisConfiguration,debugLevel);
      pair11_DerivedHistos->loadHistograms(inputFile);
      pair22_DerivedHistos->loadHistograms(inputFile);
      pair12_DerivedHistos->loadHistograms(inputFile);
      pair12_CIHistos->loadHistograms(inputFile);
      pair12_CDHistos->loadHistograms(inputFile);
      }
    }
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
}

// ///////////////////////////////////////////////////////////
// load the base histograms from given file
// ///////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::loadBaseHistograms(TFile * inputFile)
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"loadBaseHistograms(TFile * inputFile)"))
    ;
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString prefixName = getTaskName();
  particle1_Histos         = new ParticleHistos(prefixName+partName1, analysisConfiguration,debugLevel);
  particle2_Histos         = new ParticleHistos(prefixName+partName2,analysisConfiguration,debugLevel);
  particle1_Histos->loadHistograms(inputFile);
  particle2_Histos->loadHistograms(inputFile);
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos        = new ParticlePairHistos(prefixName+partName1+partName1,analysisConfiguration,debugLevel);
    pair22_Histos        = new ParticlePairHistos(prefixName+partName2+partName2,analysisConfiguration,debugLevel);
    pair12_Histos        = new ParticlePairHistos(prefixName+partName1+partName2,analysisConfiguration,debugLevel);
    pair11_Histos->loadHistograms(inputFile);
    pair22_Histos->loadHistograms(inputFile);
    pair12_Histos->loadHistograms(inputFile);
    if (analysisConfiguration->calculateDerivedHistograms)
      {
      pair11_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName1,analysisConfiguration,debugLevel);
      pair22_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName2+partName2,analysisConfiguration,debugLevel);
      pair12_DerivedHistos = new ParticlePairDerivedHistos(  prefixName+partName1+partName2,analysisConfiguration,debugLevel);
      pair12_CIHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CI",analysisConfiguration,debugLevel);
      pair12_CDHistos      = new ParticlePairCombinedHistos( prefixName+partName1+partName2+"CD",analysisConfiguration,debugLevel);
      pair11_DerivedHistos->loadHistograms(inputFile);
      pair22_DerivedHistos->loadHistograms(inputFile);
      pair12_DerivedHistos->loadHistograms(inputFile);
      pair12_CIHistos->loadHistograms(inputFile);
      pair12_CDHistos->loadHistograms(inputFile);
      }
    }
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"loadBaseHistograms(TFile * inputFile)"))
    ;
}

// ///////////////////////////////////////////////////////////
// save histograms to given files
// ///////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::saveHistograms(TFile * outputFile)
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
  if (!outputFile)
    {
    if (reportError("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Given outputFile is a null  pointer." << endl;
    postTaskError();
    return;
    }
  outputFile->cd();
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  if (!analysisConfiguration)
     {
     if (reportError("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Task configuration is a null  pointer." << endl;
     postTaskError();
     return;
     }
  if (reportDebug("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Saving singles starting." << endl;
  particle1_Histos  ->saveHistograms(outputFile);
  particle2_Histos  ->saveHistograms(outputFile);
  if (reportDebug("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Saving singles completed." << endl;
  if (analysisConfiguration->fillPairs)
    {
    if (reportDebug("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Saving pairs starting." << endl;
    pair11_Histos ->saveHistograms(outputFile);
    pair22_Histos ->saveHistograms(outputFile);
    pair12_Histos ->saveHistograms(outputFile);
    if (reportDebug("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Saving pairs completed." << endl;
    if (analysisConfiguration->calculateDerivedHistograms)
      {
      if (reportDebug("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Saving derived histograms starting." << endl;
      pair11_DerivedHistos->saveHistograms(outputFile);
      pair22_DerivedHistos->saveHistograms(outputFile);
      pair12_DerivedHistos->saveHistograms(outputFile);
      pair12_CIHistos     ->saveHistograms(outputFile);
      pair12_CDHistos     ->saveHistograms(outputFile);
      if (reportDebug("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "Saving derived histograms completed." << endl;
      }
    }
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
}


// ///////////////////////////////////////////////////////////
// add histograms to an external list
// ///////////////////////////////////////////////////////////
void TwoPartCorrelationAnalyzer::addHistogramsToExtList(TList *list, bool all)
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"addHistogramsToExtList(TList *list, bool all)"))
    ;
  /* first add the number of events as a cumulated parameter */
  list->Add(new TParameter<Long64_t>("nEventProcessed",getNEventProcessed(),'+'));
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();

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
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"addHistogramsToExtList(TList *list, bool all)"))
    ;
}


void TwoPartCorrelationAnalyzer::execute()
{
  incrementEventProcessed();
  if (!eventFilter->accept(*event)) return;
  incrementEventAccepted(); // count events used to fill histograms and for scaling at the end...
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  bool accept11;
  bool accept21;
  bool accept12;
  bool accept22;

  /* before filtering let's build the particle indexes to hurry up the process */
  for (int iParticle = 0; iParticle < event->multiplicity; iParticle++)
  {
  Particle *particle = event->getParticleAt(iParticle);
  particle->ixEtaPhi = analysisConfiguration->getIxEtaPhi(particle->eta,particle->phi);
  particle->ixYPhi = analysisConfiguration->getIxYPhi(particle->y,particle->phi);
  }

  int nAccepted1 = 0;
  int nAccepted2 = 0;
  double totalEnergy1 = 0.0;
  double totalEnergy2 = 0.0;

  for (int iParticle1=0; iParticle1<event->multiplicity; iParticle1++)
  {
  Particle & particle1 = * event->getParticleAt(iParticle1);
  accept11 = particleFilter1->accept(particle1);
  accept21 = particleFilter2->accept(particle1);
  if (accept11)
    {
    nAccepted1++;
    totalEnergy1 += particle1.e;
    particle1_Histos->fill(particle1, 1.0);
    }

  if (accept21)
    {
    nAccepted2++;
    totalEnergy2 += particle1.e;
    particle2_Histos->fill(particle1, 1.0);
    }

  if (analysisConfiguration->fillPairs)
    {
    for (int iParticle2=0; iParticle2<event->nParticles; iParticle2++)
      {
      if (iParticle1==iParticle2) continue;
      Particle & particle2 = * event->getParticleAt(iParticle2);
      accept12 = particleFilter1->accept(particle2);
      accept22 = particleFilter2->accept(particle2);
      if (accept11 && accept12)  pair11_Histos->fill(particle1, particle2, 1.0);
      if (accept21 && accept22)  pair22_Histos->fill(particle1, particle2, 1.0);
      if (accept11 && accept22)  pair12_Histos->fill(particle1, particle2, 1.0);
     }
    }
  }
  particle1_Histos->fillMultiplicity(double(nAccepted1),totalEnergy1,1.0);
  particle2_Histos->fillMultiplicity(double(nAccepted2),totalEnergy2,1.0);

}


// ====================================
// calculate Derived Histograms
// ====================================
void TwoPartCorrelationAnalyzer::calculateDerivedHistograms()
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"calculateDerivedHistograms()"))
    ;
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  particle1_Histos->completeFill();
  particle2_Histos->completeFill();
  if (reportInfo("TwoPartCorrelationAnalyzer",getTaskName(),"calculateDerivedHistograms()"))
    cout << "Fill of singles completed" << endl;
  if (analysisConfiguration->fillPairs)
    {
    if (reportInfo("TwoPartCorrelationAnalyzer",getTaskName(),"calculateDerivedHistograms()"))
      cout << "Pairs -- Calculate averages " << endl;
    particle1_Histos->calculateAverages();
    particle2_Histos->calculateAverages();
    if (reportInfo("TwoPartCorrelationAnalyzer",getTaskName(),"calculateDerivedHistograms()"))
      cout << "Pairs -- complete fill " << endl;
    pair11_Histos->completeFill();
    pair22_Histos->completeFill();
    pair12_Histos->completeFill();
    if (reportInfo("TwoPartCorrelationAnalyzer",getTaskName(),"calculateDerivedHistograms()"))
      cout << "Pairs -- calculate  derived histograms " << endl;
    pair11_DerivedHistos->calculateDerivedHistograms(particle1_Histos,particle1_Histos,pair11_Histos,analysisConfiguration->binCorrPP);
    pair22_DerivedHistos->calculateDerivedHistograms(particle2_Histos,particle2_Histos,pair22_Histos,analysisConfiguration->binCorrMM);
    pair12_DerivedHistos->calculateDerivedHistograms(particle1_Histos,particle2_Histos,pair12_Histos,analysisConfiguration->binCorrPM);
    if (reportInfo("TwoPartCorrelationAnalyzer",getTaskName(),"calculateDerivedHistograms()"))
      cout << "Pairs -- calculate CI and CD combinations " << endl;
    pair12_CIHistos->calculate(pair11_DerivedHistos,pair22_DerivedHistos,pair12_DerivedHistos, 0.25, 0.25,0.5);
    pair12_CDHistos->calculate(pair11_DerivedHistos,pair22_DerivedHistos,pair12_DerivedHistos,-0.25,-0.25,0.5);
    }
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"calculateDerivedHistograms()"))
    ;
}

// ========================================================================
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
// ========================================================================
void TwoPartCorrelationAnalyzer::scaleHistograms(double factor)
{
  if (reportStart("TwoPartCorrelationAnalyzer",getTaskName(),"scaleHistograms(double factor)"))
    ;
  particle1_Histos->scale(factor);
  particle2_Histos->scale(factor);
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos->scale(factor);
    pair22_Histos->scale(factor);
    pair12_Histos->scale(factor);
    }
  if (reportEnd("TwoPartCorrelationAnalyzer",getTaskName(),"scaleHistograms(double factor)"))
    ;
}

// =========================================================
// Reset histograms associated with this task
// Called after partial saves in subsample analyses.
// =========================================================
void TwoPartCorrelationAnalyzer::resetHistograms()
{
  if (reportInfo("TwoPartCorrelationAnalyzer",getTaskName(),"resetHistograms()")) cout << "Will reset histograms of all particle filters" << endl;
  particle1_Histos->reset();
  particle2_Histos->reset();
  ParticlePairAnalyzerConfiguration * analysisConfiguration = (ParticlePairAnalyzerConfiguration *) getTaskConfiguration();
  if (analysisConfiguration->fillPairs)
    {
    pair11_Histos->reset();
    pair22_Histos->reset();
    pair12_Histos->reset();
    }

  if (pair11_DerivedHistos)
    {
    pair11_DerivedHistos->reset();
    pair22_DerivedHistos->reset();
    pair12_DerivedHistos->reset();
    }
  if (pair12_CIHistos)
    {
    pair12_CIHistos->reset();
    pair12_CDHistos->reset();
    }

  if (reportEnd("ParticleAnalyzer",getTaskName(),"scaleHistograms(double factor)"))
    ;
}
