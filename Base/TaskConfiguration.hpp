// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                 *
 * All rights reserved.                                                *
 * Based on the ROOT package and environment                           *
 *                                                                     *
 * For the licensing terms see LICENSE.                                *
 **********************************************************************/

#ifndef WAC_TaskConfiguration
#define WAC_TaskConfiguration
#include <ostream>
#include "TObject.h"
#include "TString.h"
using std::ostream;
using std::endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////
// TaskConfiguration of a given object or task
////////////////////////////////////////////////////////////////////////////////////////////////////////
class TaskConfiguration : public TObject 
{
public:

  TaskConfiguration();
  TaskConfiguration(const TString & theName, const TString & theType, const TString & theVersion);
  TaskConfiguration(const TaskConfiguration & source);
  virtual ~TaskConfiguration();
  TaskConfiguration & operator=(const TaskConfiguration & source);

  TString getName() const;
  TString getType() const;
  TString getVersion() const;
  virtual void printTaskConfiguration(ostream & os);

  TString getInputRootFileName();
  TString getOutputRootFileName();


  ////////////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////////////
  TString name;
  TString type;
  TString version;
  bool loadHistograms;
  bool createHistograms;
  bool scaleHistograms;
  bool calculateDerivedHistograms;
  bool saveHistograms;
  bool resetHistograms;
  bool clearHistograms;
  bool forceHistogramsRewrite;

  TString inputPath;
//  TString configurationFileName;
  TString outputPath;

  TString dataInputPath;
  TString dataInputFileName;
  TString dataInputTreeName;
  int dataInputFileMinIndex;
  int dataInputFileMaxIndex;



  TString dataSourceName;      // e.g., Pythia, AMPT, ALICE, etc
  TString collisionSystemName; // e.g., pp, PbPb, etc
  TString collisionEnergyName; // e.g., 14TeV
  TString triggerName;         // e.g., MB
  TString taskTypeName;        // e.g., Single, Pair, NuDyn, etc...
  TString eventFilterName;     // e.g., MultGt2Lt40
  TString particleFilterName;  // e.g., HA, HC, HP, PiP, PiM
  TString otherName;           // e.g., Special, WideAccept, etc

  TString rootInputFileName;
  TString rootOuputFileName;

  TString dataSourceTitle;      // e.g., Pythia, AMPT, ALICE, etc
  TString collisionSystemTitle; // e.g., pp, Pb-Pb, etc
  TString collisionEnergyTitle; // e.g., 14 TeV
  TString triggerTitle;         // e.g., Min Bias
  TString taskTypeTitle;        // e.g., Single, Pair, NuDyn, etc...
  TString eventFilterTitle;     // e.g., MultGt2Lt40
  TString particleFilterTitle;  // e.g., HA, HC, HP, PiP, PiM

  TString composeFileName(int option=0);
  void setRootOutputFileName(int option=0);

  ClassDef(TaskConfiguration,0)
  
};

#endif /* TaskConfiguration_hpp */
