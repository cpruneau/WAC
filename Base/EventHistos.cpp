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

 Class defining Event Histograms
 */

#include "EventHistos.hpp"
ClassImp(EventHistos);

EventHistos::EventHistos(const TString & name,
                         EventAnalyzerConfiguration * configuration,
                         LogLevel  debugLevel)
:
Histograms(name,configuration,5,debugLevel)
{
  initialize();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// CTOR2
////////////////////////////////////////////////////////////////////////////////////////////////////////
EventHistos::EventHistos(TFile * inputFile,
                         const TString & name,
                         EventAnalyzerConfiguration * configuration,
                         LogLevel  debugLevel)
:
Histograms(name,configuration,5,debugLevel)
{
  loadHistograms(inputFile);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
////////////////////////////////////////////////////////////////////////////////////////////////////////
EventHistos::~EventHistos()
{
 /* */
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize the cluster histograms
////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventHistos::initialize()
{
  if (reportDebug()) cout << "EventHistos::initialize() Started." << endl;
  EventAnalyzerConfiguration & ac = *(EventAnalyzerConfiguration*) getConfiguration();
  TString baseName = ac.histoBaseName;
  h_nPartTot = createHistogram(baseName + "EventNpartTot",ac.nBins_nPartTot,ac.min_nPartTot,ac.max_nPartTot,"Event N_{part}","Counts");
  if (reportDebug()) cout << "EventHistos::initialize() Completed." << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// load the cluster histograms from the given file and base name
////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventHistos::loadHistograms(TFile * inputFile)
{
  if (reportDebug()) cout << "EventHistos::loadHistograms() Started." << endl;
  h_nPartTot =  loadH1(inputFile,getConfiguration()->getHistoBaseName() + "EventNpart");
  if (reportDebug()) cout << "EventHistos::loadHistograms() Completed." << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// fill the cluster histograms with the given cluster
////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventHistos::fill(double mult, double weight)
{
 h_nPartTot->Fill(mult, weight);
 }
