// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_GlobalHistos
#define WAC_GlobalHistos
#include "Histograms.hpp"
#include "GlobalAnalyzerConfiguration.hpp"

// ==================================================
//
// Histograms of "global" observables.
//
// Observables are computed in the acceptance defined by
// particle filters. The filters can be set to select narrow
// or broad ranges of kinematic variables as well as species
// of particles.
//
// n:  number of selected particles in an  event
// e:  total energy of the selected particles in an event
// q:  total (net) electric charge of the selected particles in an event
// b:  total (net) baryon charge of the selected particles in an event
//
// The fill method of this class must be called at most once per event
// otherwise weird multiple counting will happen...
//
// ==================================================
class GlobalHistos : public Histograms
{
public:

  GlobalHistos(const TString & _name,
               GlobalAnalyzerConfiguration * _configuration,
               int _nFilters,
               TString ** _filterNames,
               LogLevel  _debugLevel);
  virtual ~GlobalHistos();
  void createHistograms();
  void loadHistograms(TFile * inputFile);
  void fill(double *n, double *e, double *q, double *b, double weight);

  TString makeName(const TString & bn,const  TString & filterName1,const  TString & observableName1);
  TString makeName(const TString & bn,const  TString & filterName1,const  TString & observableName1,const  TString & filterName2,const  TString & observableName2);
  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  int nFilters;
  TString ** filterNames;

  TH1 ** h_n;
  TH1 ** h_e;
  TH1 ** h_q;
  TH1 ** h_b;
  TH2 *** h_nVsN;
  TH2 *** h_eVsN;
  TH2 *** h_qVsN;
  TH2 *** h_bVsN;
  TH2 *** h_eVsE;
  TH2 *** h_qVsE;
  TH2 *** h_bVsE;
  TH2 *** h_qVsQ;
  TH2 *** h_bVsQ;
  TH2 *** h_bVsB;

  ClassDef(GlobalHistos,0)

};

#endif /* WAC_GlobalHistos  */



