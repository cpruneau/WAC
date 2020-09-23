// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_CollisionGeometryConfiguration
#define WAC_CollisionGeometryConfiguration
#include "TaskConfiguration.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class CollisionGeometryConfiguration : public TaskConfiguration
{
public:
  
  CollisionGeometryConfiguration(const TString & name,
                                 const TString & type,
                                 const TString & version);
  virtual ~CollisionGeometryConfiguration(){}
  virtual void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  int nProtonsA, nNeutronsA, nProtonsB, nNeutronsB;

  TString histoBaseName;


  // used for generation
  double minB, maxB;
  double nnCrossSection;

  // used for histograms
  int nBins_b;          double min_b;          double max_b;
  int nBins_nPart;      double min_nPart;      double max_nPart;
  int nBins_nBinary;    double min_nBinary;    double max_nBinary;

ClassDef(CollisionGeometryConfiguration,0)
};

#endif /* WAC_CollisionGeometryConfiguration */
