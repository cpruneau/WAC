//
//  NuDynHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "NuDynHistos.hpp"
ClassImp(NuDynHistos);

NuDynHistos::NuDynHistos(const TString & name,
                         AnalysisConfiguration * configuration,
                         LogLevel  debugLevel)
:
Histograms(name,configuration,400,debugLevel)
{
  initialize();
}

NuDynHistos::NuDynHistos(TFile * inputFile,
                         const TString & name,
                         AnalysisConfiguration * configuration,
                         LogLevel  debugLevel)
:
Histograms(name,configuration,400,debugLevel)
{
  loadHistograms(inputFile);
}

NuDynHistos::~NuDynHistos()
{
  //deleteHistograms();
}

// for now use the same boundaries for eta and y histogram
void NuDynHistos::createHistograms()
{
  AnalysisConfiguration & ac = *getConfiguration();
  TString bn = getHistoBaseName();

  // ================================================================================
  // Naming convention
  // ================================================================================
  // fk_i     : factorial moments of order "k" for particle "i" (i.e., accepted by filter 1)
  // Fk_i..j  : factorial cumulants of order "k" for particles "i, ..., j" (i.e., accepted by filter i..j)
  // rk_i     : ratio of factorial moments of order "k" for particle "i" (i.e., accepted by filter 1)
  //            to product of f1_i
  // 

  // n1 number of particles satisfying filter 1
  // n2 number of particles satisfying filter 2
  // h_f1_1 = <n1>
  // h_f1_2 = <n2>
  // h_f2_11 = <n1(n1-1)>
  // h_f2_22 = <n2(n2-1)>
  // h_f2_12 = <n2(n2-1)>
  // h_r2_11 = <n1(n1-1)>/<n1><n1>
  // h_r2_22 = <n2(n2-1)>/<n2><n2>
  // h_r2_12 = <n1(n2-1)>/<n1><n2>
  // h_nudyn_12 = h_r2_11 + h_r2_22 -2*h_r2_12
  // etc

  int nPart =4;
  h_f1 = new TProfile * [4];
  h_f2 = new TProfile * [10];
  h_f3 = new TProfile * [20];
  h_f4 = new TProfile * [35];

  if (ac.nuDynVsMult)
    {
    h_f1_vsMult = new TProfile * [4];
    h_f2_vsMult = new TProfile * [10];
    h_f3_vsMult = new TProfile * [20];
    h_f4_vsMult = new TProfile * [35];
    }

  if (ac.nuDynVsCent)
    {
    h_f1_vsCent = new TProfile * [4];
    h_f2_vsCent = new TProfile * [10];
    h_f3_vsCent = new TProfile * [20];
    h_f4_vsCent = new TProfile * [35];
    }
  TString histName;
  TString histTitle;

  h_events   = createHistogram(bn+TString("Nevents"),1,ac.min_mult,  ac.max_mult,  "mult","n_{Events}");
  if (ac.nuDynVsMult) h_events_vsMult = createHistogram(bn+TString("Nevents_vsMult"),ac.nBins_mult,ac.min_mult,  ac.max_mult,  "mult","n_{Events}");
  if (ac.nuDynVsCent) h_events_vsCent = createHistogram(bn+TString("Nevents_vsCent"),ac.nBins_cent,ac.min_cent,  ac.max_cent,  "cent","n_{Events}");

//  TString baseName = bn + "f1_";
//  TString baseTitle = "f_{1}";

  for (int i1=0; i1<nPart; i1++)
    {
    histName  = bn + "f1_";
    histName  += i1;
    histTitle = "f_{1}";
    histTitle += "^{";
    histTitle += i1;
    histTitle += "}";
    h_f1[i1]  = createProfile(histName, 1,ac.min_mult,ac.max_mult, "mult", histTitle);
    if (ac.nuDynVsMult)
      {
      histName = bn + "f1_";
      histName += i1;
      histName += "_vsMult";
      h_f1_vsMult[i1]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult",  histTitle);
      }
    if (ac.nuDynVsCent)
      {
      histName = bn + "f1_";
      histName += i1;
      histName += "_vsCent";
      h_f1_vsCent[i1]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent",  histTitle);
      }
    for (int i2=i1; i2<nPart; i2++)
      {
      int index12 = index2(i1,i2);
      histName  = bn + "f2_";
      histName  += i1;
      histName  += i2;
      histTitle = "f_{2}";
      histTitle += "^{";
      histTitle += i1;
      histTitle += i2;
      histTitle += "}";
      h_f2[index12]  = createProfile(histName, 1,ac.min_mult,ac.max_mult,"mult",histTitle);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "f2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_f2_vsMult[index12]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult",histTitle);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "f2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_f2_vsCent[index12]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent",histTitle);
        }
      for (int i3=i2; i3<nPart; i3++)
        {
        int index123 = index3(i1,i2,i3);

        histName  = bn + "f3_";
        histName  += i1;
        histName  += i2;
        histName  += i3;
        histTitle = "f_{3}";
        histTitle += "^{";
        histTitle += i1;
        histTitle += i2;
        histTitle += i3;
        histTitle += "}";
        h_f3[index123]  = createProfile(histName, 1,ac.min_mult,ac.max_mult,"mult",histTitle);
        if (ac.nuDynVsMult)
          {
          histName  = bn + "f3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsMult";
          h_f3_vsMult[index123]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult",histTitle);
          }
        if (ac.nuDynVsCent)
          {
          histName  = bn + "f3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsCent";
          h_f3_vsCent[index123]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent",histTitle);
          }
        for (int i4=i3; i4<nPart; i4++)
          {
          int index1234 = index4(i1,i2,i3,i4);

          histName  = bn + "f4_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName  += i4;
          histTitle = "f_{4}";
          histTitle += "^{";
          histTitle += i1;
          histTitle += i2;
          histTitle += i3;
          histTitle += i4;
          histTitle += "}";
          h_f4[index1234]  = createProfile(histName, 1,ac.min_mult,ac.max_mult,"mult",histTitle);
          if (ac.nuDynVsMult)
            {
            histName  = bn + "f4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsMult";
            h_f4_vsMult[index1234]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult",histTitle);
            }
          if (ac.nuDynVsCent)
            {
            histName  = bn + "f4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsCent";
            h_f4_vsCent[index1234]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent",histTitle);
            }
          }
        }
      }
    }

}

//________________________________________________________________________
void NuDynHistos::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportFatal()) cout << "-Fatal- Attempting to load NuDynHistos from an invalid file pointer" << endl;
    return;
    }
  AnalysisConfiguration & ac = *getConfiguration();
  TString  bn = getHistoBaseName();
  TString histName;
  TString histTitle;

  h_events   = loadH1(inputFile,bn+TString("Nevents"));
  if (ac.nuDynVsMult) h_events_vsMult = loadH1(inputFile,bn+TString("Nevents_vsMult"),ac.nBins_mult,ac.min_mult,  ac.max_mult,  "mult","n_{Events}");
  if (ac.nuDynVsCent) h_events_vsCent = loadH1(inputFile,bn+TString("Nevents_vsCent"),ac.nBins_cent,ac.min_cent,  ac.max_cent,  "cent","n_{Events}");

//  baseName = bn + "f1_";
//  baseTitle = "f_{1}";

  int nPart =4;
  h_f1 = new TProfile * [4];
  h_f2 = new TProfile * [10];
  h_f3 = new TProfile * [20];
  h_f4 = new TProfile * [25];

  if (ac.nuDynVsMult)
     {
     h_f1_vsMult = new TProfile * [4];
     h_f2_vsMult = new TProfile * [10];
     h_f3_vsMult = new TProfile * [20];
     h_f4_vsMult = new TProfile * [35];
     }

   if (ac.nuDynVsCent)
     {
     h_f1_vsCent = new TProfile * [4];
     h_f2_vsCent = new TProfile * [10];
     h_f3_vsCent = new TProfile * [20];
     h_f4_vsCent = new TProfile * [35];
     }


  for (int i1=0; i1<nPart; i1++)
    {
    histName  = bn + "f1_";
    histName  += i1;
    histTitle = "f_{1}";
    histTitle += "^{";
    histTitle += i1;
    histTitle += "}";
    h_f1[i1]  = loadProfile(inputFile,histName);
    if (ac.nuDynVsMult)
      {
      histName = bn + "f1_";
      histName += i1;
      histName += "_vsMult";
      h_f1_vsMult[i1]  = loadProfile(inputFile,histName);
      }
    if (ac.nuDynVsCent)
      {
      histName = bn + "f1_";
      histName += i1;
      histName += "_vsCent";
      h_f1_vsMult[i1]  = loadProfile(inputFile,histName);
      }
    for (int i2=i1; i2<nPart; i2++)
      {
      int index12 = index2(i1,i2);
      histName  = bn + "f2_";
      histName  += i1;
      histName  += i2;
      histTitle = "f_{2}";
      histTitle += "^{";
      histTitle += i1;
      histTitle += i2;
      histTitle += "}";
      h_f2[index12]  = loadProfile(inputFile,histName);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "f2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_f2_vsMult[index12]  = loadProfile(inputFile,histName);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "f2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_f2_vsMult[index12]  = loadProfile(inputFile,histName);
        }
      for (int i3=i2; i3<nPart; i3++)
        {
         int index123 = index3(i1,i2,i3);
        histName  = bn + "f3_";
        histName  += i1;
        histName  += i2;
        histName  += i3;
        histTitle = "f_{3}";
        histTitle += "^{";
        histTitle += i1;
        histTitle += i2;
        histTitle += i3;
        histTitle += "}";
        h_f3[index123]  = loadProfile(inputFile,histName);
        if (ac.nuDynVsMult)
          {
          histName  = bn + "f3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsMult";
          h_f3_vsMult[index123]  = loadProfile(inputFile,histName);
          }
        if (ac.nuDynVsCent)
          {
          histName  = bn + "f3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsCent";
          h_f3_vsMult[index123]  = loadProfile(inputFile,histName);
          }
        for (int i4=i3; i4<nPart; i4++)
          {
          int index1234 = index4(i1,i2,i3,i4);
          histName  = bn + "f4_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName  += i4;
          histTitle = "f_{4}";
          histTitle += "^{";
          histTitle += i1;
          histTitle += i2;
          histTitle += i3;
          histTitle += i4;
          histTitle += "}";
          h_f4[index1234]  = loadProfile(inputFile,histName);
          if (ac.nuDynVsMult)
            {
            histName  = bn + "f3_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsMult";
            h_f4_vsMult[index1234]  = loadProfile(inputFile,histName);
            }
          if (ac.nuDynVsCent)
            {
            histName  = bn + "f3_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsCent";
            h_f4_vsMult[index1234]  = loadProfile(inputFile,histName);
            }
          }
        }
      }
    }
}


void NuDynHistos::fill(double mult, double cent, double * n, double weight)
{
   //if (reportDebug()) cout << "NuDynTask::fill(...) started " << endl;
   AnalysisConfiguration & ac = *getConfiguration();
   h_events->Fill(mult, weight);
  //if (reportDebug()) cout << "NuDynTask::fill(...) 1 " << endl;
  if (ac.nuDynVsMult) h_events_vsMult->Fill(mult, weight);
  if (ac.nuDynVsCent) h_events_vsCent->Fill(cent, weight);
  double fill;
  int nPart = 4;
  //if (reportDebug()) cout << "NuDynTask::fill(...) 2" << endl;
  for (int i1=0; i1<nPart; i1++)
  {

  fill = n[i1];
  h_f1[i1]->Fill(mult,fill,weight);
  if (ac.nuDynVsMult) h_f1_vsMult[i1]->Fill(mult,fill,weight);
  if (ac.nuDynVsCent) h_f1_vsCent[i1]->Fill(cent,fill,weight);
  //if (reportDebug()) cout << "NuDynTask::fill(...) 3" << endl;
  for (int i2=i1; i2<nPart; i2++)
    {
      int index12 = index2(i1,i2);
    //if (reportDebug()) cout << "NuDynHistos::fill(...) i1:" << i1 << "  i2:" << i2 << " index12:" << index12 << endl;
    fill = n[i1]* (n[i2] - (i1==i2?1:0));
     h_f2[index12]->Fill(mult,fill,weight);
     if (ac.nuDynVsMult) h_f2_vsMult[index12]->Fill(mult,fill,weight);
     if (ac.nuDynVsCent) h_f2_vsCent[index12]->Fill(cent,fill,weight);
    for (int i3=i2; i3<nPart; i3++)
      {
      int index123 = index3(i1,i2,i3);
     //if (reportDebug()) cout << "NuDynTask::fill(...) i1:" << i1 << "  i2:" << i2 << "  i3:" << i3 << " index123:" << index123 << endl;
      fill = n[i1]* (n[i2] - (i1==i2?1:0)) * (n[i3] - (i1==i3?1:0) - (i2==i3?1:0));
     h_f3[index123]->Fill(mult,fill,weight);
      if (ac.nuDynVsMult) h_f3_vsMult[index123]->Fill(mult,fill,weight);
      if (ac.nuDynVsCent) h_f3_vsCent[index123]->Fill(cent,fill,weight);
      for (int i4=i3; i4<nPart; i4++)
        {
           int index1234 = index4(i1,i2,i3,i4);
    //if (reportDebug()) cout << "NuDynTask::fill(...) i1:" << i1 << "  i2:" << i2 << "  i3:" << i3 << "  i4:" << i4 << " index1234:" << index1234 << endl;
        fill = n[i1] * (n[i2] - (i1==i2?1:0)) * (n[i3] - (i1==i3?1:0) - (i2==i3?1:0)) * ( n[i4] - (i1==i4?1:0) - (i2==i4?1:0) - (i3==i4?1:0)) ;
        h_f4[index1234]->Fill(mult,fill,weight);
        if (ac.nuDynVsMult) h_f4_vsMult[index1234]->Fill(mult,fill,weight);
        if (ac.nuDynVsCent) h_f4_vsCent[index1234]->Fill(cent,fill,weight);
        }
      }
    }
  }
}

