//
//  NuDynHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//
#include "NuDynHistos.hpp"
ClassImp(NuDynHistos);

NuDynHistos::NuDynHistos(const TString & collectionName,
                         int  * identicalUsed,
                         NuDynConfiguration * analysisConfiguration,
                         LogLevel  debugLevel)
:
Histograms(collectionName,configuration,500,debugLevel),
identical(identicalUsed)
{
  // no ops
}


NuDynHistos::~NuDynHistos()
{
  //deleteHistograms();
}

// for now use the same boundaries for eta and y histogram
void NuDynHistos::createHistograms()
{

  TString bn = getHistoBaseName();
  NuDynConfiguration & ac = *(NuDynConfiguration*) getConfiguration();
  TString suffix = "";
  TString xTitle = "";
  switch ( ac.multiplicityType )
  {
    case NuDynConfiguration::Centrality:           suffix = "_vsCent"; xTitle = "%"; break;
    case NuDynConfiguration::TotalMultiplicity:    suffix = "_vsMult"; xTitle = "mult_{Tot}";  break;
    case NuDynConfiguration::AcceptedMultiplicity: suffix = "_vsMult"; xTitle = "mult_{acc}";  break;
  }


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

  nParticles = 4;
  h_f1 = new TProfile * [4];
  h_f2 = new TProfile * [10];
  h_f3 = new TProfile * [20];
  h_f4 = new TProfile * [35];

  h_f1_vsMult = new TProfile * [4];
  h_f2_vsMult = new TProfile * [10];
  h_f3_vsMult = new TProfile * [20];
  h_f4_vsMult = new TProfile * [35];

  TString histName;
  TString histTitle;

  h_events   = createHistogram(bn+TString("Nevents"),1,ac.min_mult,  ac.max_mult,  xTitle,"n_{Events}");
  h_events_vsMult = createHistogram(bn+TString("Nevents")+suffix,ac.nBins_mult,ac.min_mult,  ac.max_mult,  xTitle,"n_{Events}");

  for (int i1=0; i1<nParticles; i1++)
  {
  histName  = bn + "f1_";
  histName  += i1;
  histTitle = "f_{1}";
  histTitle += "^{";
  histTitle += i1;
  histTitle += "}";
  h_f1[i1]  = createProfile(histName, 1,ac.min_mult,ac.max_mult, xTitle, histTitle);
  histName += suffix;
  h_f1_vsMult[i1]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, xTitle,  histTitle);
  for (int i2=i1; i2<nParticles; i2++)
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
    h_f2[index12]  = createProfile(histName, 1,ac.min_mult,ac.max_mult,xTitle,histTitle);
    histName += suffix;
    h_f2_vsMult[index12]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult,xTitle,histTitle);
    for (int i3=i2; i3<nParticles; i3++)
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
      h_f3[index123]  = createProfile(histName, 1,ac.min_mult,ac.max_mult,xTitle,histTitle);
      histName += suffix;
      h_f3_vsMult[index123]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult,xTitle,histTitle);
      for (int i4=i3; i4<nParticles; i4++)
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
        h_f4[index1234]  = createProfile(histName, 1,ac.min_mult,ac.max_mult,xTitle,histTitle);
        histName += suffix;
        h_f4_vsMult[index1234]  = createProfile(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult,xTitle,histTitle);
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
  TString bn = getHistoBaseName();
  NuDynConfiguration & ac = *(NuDynConfiguration*) getConfiguration();
  TString suffix = "";
  switch ( ac.multiplicityType )
  {
    case NuDynConfiguration::Centrality: suffix = "_vsCent"; break;
    case NuDynConfiguration::TotalMultiplicity: suffix = "_vsMult"; break;
    case NuDynConfiguration::AcceptedMultiplicity: suffix = "_vsMult"; break;
  }
  TString histName;
  TString histTitle;
  h_events        = loadH1(inputFile,bn+TString("Nevents"));
  h_events_vsMult = loadH1(inputFile,bn+TString("Nevents")+suffix);

  int nPart =4;
  h_f1 = new TProfile * [4];
  h_f2 = new TProfile * [10];
  h_f3 = new TProfile * [20];
  h_f4 = new TProfile * [35];

  h_f1_vsMult = new TProfile * [4];
  h_f2_vsMult = new TProfile * [10];
  h_f3_vsMult = new TProfile * [20];
  h_f4_vsMult = new TProfile * [35];

  for (int i1=0; i1<nPart; i1++)
  {
  histName  = bn + "f1_";
  histName  += i1;
  h_f1[i1]  = loadProfile(inputFile,histName);
  histName += suffix;
  h_f1_vsMult[i1]  = loadProfile(inputFile,histName);
  for (int i2=i1; i2<nPart; i2++)
    {
    int index12 = index2(i1,i2);
    histName  = bn + "f2_";
    histName  += i1;
    histName  += i2;
    h_f2[index12]  = loadProfile(inputFile,histName);
    histName += suffix;
    h_f2_vsMult[index12]  = loadProfile(inputFile,histName);
    for (int i3=i2; i3<nPart; i3++)
      {
      int index123 = index3(i1,i2,i3);
      histName  = bn + "f3_";
      histName  += i1;
      histName  += i2;
      histName  += i3;
      h_f3[index123]  = loadProfile(inputFile,histName);
        histName += suffix;
        h_f3_vsMult[index123]  = loadProfile(inputFile,histName);
      for (int i4=i3; i4<nPart; i4++)
        {
        int index1234 = index4(i1,i2,i3,i4);
        histName  = bn + "f4_";
        histName  += i1;
        histName  += i2;
        histName  += i3;
        histName  += i4;
        h_f4[index1234]  = loadProfile(inputFile,histName);
          histName += suffix;
          h_f4_vsMult[index1234]  = loadProfile(inputFile,histName);
        }
      }
    }
  }
}


void NuDynHistos::fill(double mult, double * nAccepted, double weight)
{
  h_events->Fill(mult, weight);
  h_events_vsMult->Fill(mult, weight);
  double fill;
  for (int i1=0; i1<nParticles; i1++)
  {
  fill = nAccepted[i1];
  h_f1[i1]->Fill(mult,fill,weight);
  h_f1_vsMult[i1]->Fill(mult,fill,weight);
  for (int i2=i1; i2<nParticles; i2++)
    {
    int index12 = index2(i1,i2);
    fill = nAccepted[i1]*(nAccepted[i2] - sameFilter(i1,i2));
    h_f2[index12]->Fill(mult,fill,weight);
    h_f2_vsMult[index12]->Fill(mult,fill,weight);
    for (int i3=i2; i3<nParticles; i3++)
      {
      int index123 = index3(i1,i2,i3);
      fill = nAccepted[i1] * (nAccepted[i2] - sameFilter(i1,i2)) * (nAccepted[i3] - sameFilter(i1,i3) - sameFilter(i2,i3));
      h_f3[index123]->Fill(mult,fill,weight);
      h_f3_vsMult[index123]->Fill(mult,fill,weight);
      for (int i4=i3; i4<nParticles; i4++)
        {
        int index1234 = index4(i1,i2,i3,i4);
        fill = nAccepted[i1] * (nAccepted[i2] - sameFilter(i1,i2)) * (nAccepted[i3] - sameFilter(i1,i3) - sameFilter(i2,i3)) * ( nAccepted[i4] - sameFilter(i1,i4) - sameFilter(i2,i4) - sameFilter(i3,i4)) ;
        h_f4[index1234]->Fill(mult,fill,weight);
        h_f4_vsMult[index1234]->Fill(mult,fill,weight);
        }
      }
    }
  }
}

int  NuDynHistos::sameFilter(int i1, int i2)
{
  if (i1>=0 && i1<=3 && i2>=0 && i2<=3)
    return identical[4*i1+i2];
  else
    return 0;
}
