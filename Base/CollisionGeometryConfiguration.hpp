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
#include "NucleusGenerator.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class CollisionGeometryConfiguration : public TaskConfiguration
{
public:
  
  CollisionGeometryConfiguration(const TString & name,
                                 const TString & type,
                                 const TString & version);
  CollisionGeometryConfiguration(const CollisionGeometryConfiguration & source);
  virtual ~CollisionGeometryConfiguration();
  CollisionGeometryConfiguration & operator=(const CollisionGeometryConfiguration & source);

  virtual void setBXsectPartition(int nBins, double * bValues);
  virtual void setDefaultBXsectPartition();
  virtual void setPbPb2740GeVBXsectPartition();

  virtual void setDefaultConfiguration();
  virtual void setPbPb2740GeVConfiguration();

  virtual void printConfiguration(ostream & os);


  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  TString histoBaseName;

  // Nucleus A Generation
  int nProtonsA, nNeutronsA;
  NucleusGenerator::GeneratorType aGeneratorType;
  double aParA, aParB, aParC;
  int    aNR;
  double aMinR, aMaxR;

  // Nucleus B Generation
  int nProtonsB, nNeutronsB;
  NucleusGenerator::GeneratorType bGeneratorType;
  double bParA, bParB, bParC;
  int    bNR;
  double bMinR, bMaxR;

  // used for collision generation
  double minB, maxB;
  double nnCrossSection;


  // used for histograms
  int nBins_b;          double min_b;          double max_b;
  int nBins_nPart;      double min_nPart;      double max_nPart;
  int nBins_nBinary;    double min_nBinary;    double max_nBinary;

  int nBins_bxSect;
  double * bValues;

ClassDef(CollisionGeometryConfiguration,0)
};

#endif /* WAC_CollisionGeometryConfiguration */
