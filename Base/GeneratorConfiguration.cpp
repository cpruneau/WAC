// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/
/**
 \class GeneratorConfiguration
 \ingroup WAC

 Event Generator Configuration Class
 */
#include "GeneratorConfiguration.hpp"
ClassImp(GeneratorConfiguration);

GeneratorConfiguration::GeneratorConfiguration(const TString & name, const TString & type, const TString & version)
:
TaskConfiguration(name,type,version),
 dataInputPath(""),
 dataFileMinIndex(0),
 dataFileMaxIndex(1),

impactParameterMinimum( 0.000),
impactParameterMaximum( 1000.0),
ptMinimum    (0.0001),
ptMaximum    (1000.0),
yMinimum     (-4.0),
yMaximum     ( 4.0),
etaMinimum   (-2.0),
etaMaximum   ( 2.0),
chargedOnly  (true),
eventsSource ("")
{
  // no ops
}

GeneratorConfiguration::
GeneratorConfiguration(const GeneratorConfiguration & source)
:
TaskConfiguration (source),
dataInputPath     (  source.dataInputPath ),
dataFileMinIndex  ( source.dataFileMinIndex ),
dataFileMaxIndex  ( source.dataFileMaxIndex ),
impactParameterMinimum (  source.impactParameterMinimum ),
impactParameterMaximum (  source.impactParameterMaximum ),
ptMinimum    ( source.ptMinimum ),
ptMaximum    ( source.ptMaximum ),
yMinimum     ( source.yMinimum ),
yMaximum     ( source.yMaximum ),
etaMinimum   ( source.etaMinimum ),
etaMaximum   ( source.etaMaximum ),
chargedOnly  ( source.chargedOnly ),
eventsSource ( source.eventsSource )
{
  // no ops
}

GeneratorConfiguration & GeneratorConfiguration::operator=(const GeneratorConfiguration & source)
{
  if (this != &source)
    {
    TaskConfiguration::operator=(source);
    dataInputPath  = source.dataInputPath;
    dataFileMinIndex = source.dataFileMinIndex;
    dataFileMaxIndex = source.dataFileMaxIndex;
    impactParameterMinimum = source.impactParameterMinimum;
    impactParameterMaximum = source.impactParameterMaximum;
    ptMinimum = source.ptMinimum;
    ptMaximum = source.ptMaximum;
    yMinimum  = source.yMinimum;
    yMaximum    = source.yMaximum;
    etaMinimum   = source.etaMinimum;
    etaMaximum   = source.etaMaximum;
    chargedOnly  = source.chargedOnly;
    eventsSource = source.eventsSource;
    }
  return *this;
}

void GeneratorConfiguration::printConfiguration(ostream & os)
{
  os
  << " ------------------------------------------------------------------------------------------" << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                           modelName: " << name << endl
  << "                        modelVersion: " << version << endl
  << "              impactParameterMinimum: " << impactParameterMinimum << endl
  << "              impactParameterMaximum: " << impactParameterMaximum << endl
  << "                           ptMinimum: " << ptMinimum << endl
  << "                           ptMaximum: " << ptMaximum << endl
  << "                            yMinimum: " << yMinimum << endl
  << "                            yMaximum: " << yMaximum << endl
  << "                          etaMinimum: " << etaMinimum << endl
  << "                          etaMaximum: " << etaMaximum << endl
  << "                         chargedOnly: " << chargedOnly << endl
  << "                        eventsSource: " << eventsSource << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << " ------------------------------------------------------------------------------------------" << endl;
}
