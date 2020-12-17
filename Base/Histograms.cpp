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

#include "Histograms.hpp"
ClassImp(Histograms);

Histograms::Histograms(const TString & name,
                       TaskConfiguration * config,
                       int initialCapacity,
                       LogLevel  debugLevel)
:
HistogramCollection(name,initialCapacity,debugLevel),
configuration(config)
{
 /* */
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
////////////////////////////////////////////////////////////////////////////////////////////////////////
Histograms::~Histograms()
{
  /* */
}

void Histograms::reset()
{
  if (reportDebug()) cout << "Histograms::reset()" << endl;
  for (int iHisto=0; iHisto<getNHistograms(); iHisto++)
    {
    getHisto(iHisto)->Reset();
    }
}

// overload this class to create histograms...
void Histograms::createHistograms()
{
  if (reportDebug()) cout << "Histograms::clear()" << endl;
}

// ==============================================================
// load the cluster histograms from the given file and base name
// ==============================================================
void Histograms::loadHistograms(TFile * inputFile)
{
  inputFile = 0; // stop warnings;
  if (reportDebug()) cout << "Histograms::loadHistograms(...) No ops." << endl;
}

TaskConfiguration * Histograms::getConfiguration() const
{
  return configuration;
}

void Histograms::setConfiguration(TaskConfiguration * config)
{
  configuration = config;
}

TString Histograms::getHistoBaseName() const
{
//  TaskConfiguration & ac = *(TaskConfiguration*) getConfiguration();
  TString bn; //ac.histoBaseName;
//  bn += "_";
  bn = getName();
  bn += "_";
  return bn;
}
