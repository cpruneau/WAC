//
//  GlobalHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "GlobalHistos.hpp"
ClassImp(GlobalHistos);

GlobalHistos::GlobalHistos(const TString & _name,
                           GlobalAnalyzerConfiguration * _configuration,
                           int _nFilters,
                           TString ** _filterNames,
                           LogLevel  _debugLevel)
:
Histograms(_name,_configuration,100,_debugLevel),
nFilters(_nFilters),
filterNames(_filterNames),
h_n(nullptr),
h_e(nullptr),
h_q(nullptr),
h_b(nullptr),
h_nVsN(nullptr),
h_eVsN(nullptr),
h_qVsN(nullptr),
h_bVsN(nullptr),
h_eVsE(nullptr),
h_qVsE(nullptr),
h_bVsE(nullptr),
h_qVsQ(nullptr),
h_bVsQ(nullptr),
h_bVsB(nullptr)
{
  h_n    = new TH1* [nFilters];
  h_e    = new TH1* [nFilters];
  h_q    = new TH1* [nFilters];
  h_b    = new TH1* [nFilters];
  h_nVsN = new TH2** [nFilters];
  h_eVsN = new TH2** [nFilters];
  h_qVsN = new TH2** [nFilters];
  h_bVsN = new TH2** [nFilters];
  h_eVsE = new TH2** [nFilters];
  h_qVsE = new TH2** [nFilters];
  h_bVsE = new TH2** [nFilters];
  h_qVsQ = new TH2** [nFilters];
  h_bVsQ = new TH2** [nFilters];
  h_bVsB = new TH2** [nFilters];
  if (_configuration->fillCorrelationHistos)
    {
    for (int k1=0; k1<nFilters; k1++)
      {
      int nm1 = nFilters-k1;
      h_nVsN[k1] = new TH2* [nm1];
      h_eVsN[k1] = new TH2* [nm1];
      h_qVsN[k1] = new TH2* [nm1];
      h_bVsN[k1] = new TH2* [nm1];
      h_eVsE[k1] = new TH2* [nm1];
      h_qVsE[k1] = new TH2* [nm1];
      h_bVsE[k1] = new TH2* [nm1];
      h_qVsQ[k1] = new TH2* [nm1];
      h_bVsQ[k1] = new TH2* [nm1];
      h_bVsB[k1] = new TH2* [nm1];
      }
    }

}

GlobalHistos::~GlobalHistos()
{
  //deleteHistograms();
}


TString GlobalHistos::makeName(const TString & bn,const  TString & filterName1,const  TString & observableName1)
{
  TString name = bn;
  name += "_";
  name += filterName1;
  name += "_";
  name += observableName1;
  return name;
}


TString GlobalHistos::makeName(const TString & bn,const  TString & filterName1,const  TString & observableName1,const  TString & filterName2,const  TString & observableName2)
{
  TString name = bn;
  name += "_";
  name += filterName1;
  name += "_";
  name += observableName1;
  name += "_";
  name += filterName2;
  name += "_";
  name += observableName2;
  return name;
}

// for now use the same boundaries for eta and y histogram
void GlobalHistos::createHistograms()
{
  GlobalAnalyzerConfiguration & ac = *(GlobalAnalyzerConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();
  TString name;

  for (int k1=0; k1<nFilters; k1++)
  {
  name = makeName(bn,*filterNames[k1],"n");
  h_n[k1] = createHistogram(name,  ac.nBins_n, ac.min_n,  ac.max_n, "n","N", true,true,true,false);
  name = makeName(bn,*filterNames[k1],"e");
  h_e[k1] = createHistogram(name,  ac.nBins_e, ac.min_e,  ac.max_e, "e","N", true,true,true,false);
  name = makeName(bn,*filterNames[k1],"q");
  h_q[k1] = createHistogram(name,  ac.nBins_q, ac.min_q,  ac.max_q, "q","N", true,true,true,false);
  name = makeName(bn,*filterNames[k1],"b");
  h_b[k1] = createHistogram(name,  ac.nBins_b, ac.min_b,  ac.max_b, "b","N", true,true,true,false);
  int n2 = nFilters-k1;

  if (ac.fillCorrelationHistos)
    {
    for (int k2=k1; k2<n2; k2++)
      {
      if (k1!=k2)
        {
        name = makeName(bn,*filterNames[k1],"n",*filterNames[k2],"n");
        h_nVsN[k1][k2] = createHistogram(name,  ac.nBins_n2, ac.min_n,  ac.max_n, ac.nBins_n2, ac.min_n,  ac.max_n, "n", "n", "N", true,true,true,false);
        name = makeName(bn,*filterNames[k1],"e",*filterNames[k2],"e");
        h_eVsE[k1][k2] = createHistogram(name,  ac.nBins_e2, ac.min_e,  ac.max_e, ac.nBins_e2, ac.min_e,  ac.max_e, "e", "e", "N", true,true,true,false);
        name = makeName(bn,*filterNames[k1],"q",*filterNames[k2],"q");
        h_qVsQ[k1][k2] = createHistogram(name,  ac.nBins_q2, ac.min_q,  ac.max_q, ac.nBins_q2, ac.min_q,  ac.max_q, "q", "q", "N", true,true,true,false);
        name = makeName(bn,*filterNames[k1],"b",*filterNames[k2],"b");
        h_bVsB[k1][k2] = createHistogram(name,  ac.nBins_b2, ac.min_b,  ac.max_b, ac.nBins_b2, ac.min_b,  ac.max_b, "b", "b", "N", true,true,true,false);
        }
      else
        {
        h_nVsN[k1][k2] = nullptr;
        h_eVsE[k1][k2] = nullptr;
        h_qVsQ[k1][k2] = nullptr;
        h_bVsB[k1][k2] = nullptr;
        }
      name = makeName(bn,*filterNames[k1],"n",*filterNames[k2],"e");
      h_eVsN[k1][k2] = createHistogram(name,  ac.nBins_n2, ac.min_n,  ac.max_n, ac.nBins_e2, ac.min_e,  ac.max_e, "n", "e", "N", true,true,true,false);
      name = makeName(bn,*filterNames[k1],"n",*filterNames[k2],"q");
      h_qVsN[k1][k2] = createHistogram(name,  ac.nBins_n2, ac.min_n,  ac.max_n, ac.nBins_q2, ac.min_q,  ac.max_q, "n", "q", "N", true,true,true,false);
      name = makeName(bn,*filterNames[k1],"n",*filterNames[k2],"b");
      h_bVsN[k1][k2] = createHistogram(name,  ac.nBins_n2, ac.min_n,  ac.max_n, ac.nBins_b2, ac.min_b,  ac.max_b, "n", "b", "N", true,true,true,false);
      name = makeName(bn,*filterNames[k1],"e",*filterNames[k2],"q");
      h_qVsE[k1][k2] = createHistogram(name,  ac.nBins_e2, ac.min_e,  ac.max_e, ac.nBins_q2, ac.min_q,  ac.max_q, "e", "q", "N", true,true,true,false);
      name = makeName(bn,*filterNames[k1],"e",*filterNames[k2],"b");
      h_bVsE[k1][k2] = createHistogram(name,  ac.nBins_e2, ac.min_e,  ac.max_e, ac.nBins_b2, ac.min_b,  ac.max_b, "e", "b", "N", true,true,true,false);
      name = makeName(bn,*filterNames[k1],"q",*filterNames[k2],"b");
      h_bVsQ[k1][k2] = createHistogram(name,  ac.nBins_q2, ac.min_q,  ac.max_q, ac.nBins_b2, ac.min_b,  ac.max_b, "q", "b", "N", true,true,true,false);
      }
    }
  }
}


//________________________________________________________________________
void GlobalHistos::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportFatal()) cout << "Attempting to load GlobalHistos from an invalid file pointer" << endl;
    return;
    }
  GlobalAnalyzerConfiguration & ac = *(GlobalAnalyzerConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();

  for (int k1=0; k1<nFilters; k1++)
  {
  h_n[k1] = loadH1(inputFile, makeName(bn,*filterNames[k1],"n"));
  h_e[k1] = loadH1(inputFile, makeName(bn,*filterNames[k1],"e"));
  h_q[k1] = loadH1(inputFile, makeName(bn,*filterNames[k1],"q"));
  h_b[k1] = loadH1(inputFile, makeName(bn,*filterNames[k1],"b"));
  if (ac.fillCorrelationHistos)
    {
    int n2 = nFilters-k1;
    for (int k2=k1; k2<n2; k2++)
      {
      if (k1!=k2)
        {
        h_nVsN[k1][k2] = loadH2(inputFile, makeName(bn,*filterNames[k1],"n",*filterNames[k2],"n"));
        h_eVsE[k1][k2] = loadH2(inputFile, makeName(bn,*filterNames[k1],"e",*filterNames[k2],"e"));
        h_qVsQ[k1][k2] = loadH2(inputFile, makeName(bn,*filterNames[k1],"q",*filterNames[k2],"q"));
        h_bVsB[k1][k2] = loadH2(inputFile, makeName(bn,*filterNames[k1],"b",*filterNames[k2],"b"));
        }
      else
        {
        h_nVsN[k1][k2] = nullptr;
        h_eVsE[k1][k2] = nullptr;
        h_qVsQ[k1][k2] = nullptr;
        h_bVsB[k1][k2] = nullptr;
        }
      h_eVsN[k1][k2] = loadH2(inputFile,  makeName(bn,*filterNames[k1],"n",*filterNames[k2],"e") );
      h_qVsN[k1][k2] = loadH2(inputFile,  makeName(bn,*filterNames[k1],"n",*filterNames[k2],"q") );
      h_bVsN[k1][k2] = loadH2(inputFile,  makeName(bn,*filterNames[k1],"n",*filterNames[k2],"b") );
      h_qVsE[k1][k2] = loadH2(inputFile,  makeName(bn,*filterNames[k1],"e",*filterNames[k2],"q") );
      h_bVsE[k1][k2] = loadH2(inputFile,  makeName(bn,*filterNames[k1],"e",*filterNames[k2],"b") );
      h_bVsQ[k1][k2] = loadH2(inputFile,  makeName(bn,*filterNames[k1],"q",*filterNames[k2],"b") );
      }
    }
  }
}

// =============================================================================
// n : number of particles accepted by filters
// e : total energy accepted by filters
// q : net charge accepted by filters
// b : net baryon number accepted by filters
// =============================================================================
void GlobalHistos::fill(double *n, double *e, double *q, double *b, double weight)
{
  GlobalAnalyzerConfiguration & ac = *(GlobalAnalyzerConfiguration*)getConfiguration();

  for (int k1=0; k1<nFilters; k1++)
  {
  h_n[k1]->Fill(n[k1],weight);
  h_e[k1]->Fill(e[k1],weight);
  h_q[k1]->Fill(q[k1],weight);
  h_b[k1]->Fill(b[k1],weight);
  if (ac.fillCorrelationHistos)
    {
    int n2 = nFilters-k1;
    for (int k2=k1; k2<n2; k2++)
      {
      if (k1!=k2)
        {
        h_nVsN[k1][k2]->Fill(n[k1],n[k2],weight);
        h_eVsE[k1][k2]->Fill(e[k1],e[k2],weight);
        h_qVsQ[k1][k2]->Fill(q[k1],q[k2],weight);
        h_bVsB[k1][k2]->Fill(b[k1],b[k2],weight);
        }
      h_eVsN[k1][k2]->Fill(n[k1],e[k2],weight);
      h_qVsN[k1][k2]->Fill(n[k1],q[k2],weight);
      h_bVsN[k1][k2]->Fill(n[k1],b[k2],weight);
      h_qVsE[k1][k2]->Fill(e[k1],q[k2],weight);
      h_bVsE[k1][k2]->Fill(e[k1],b[k2],weight);
      h_bVsQ[k1][k2]->Fill(q[k1],b[k2],weight);
      }
    }
  }
}
