//
//  NuDynDerivedHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include "NuDynDerivedHistos.hpp"
ClassImp(NuDynDerivedHistos);

NuDynDerivedHistos::NuDynDerivedHistos(const TString & name,
                                       AnalysisConfiguration * configuration,
                                       LogLevel  debugLevel)
:
Histograms(name,configuration,600,debugLevel)
{
  h_F2 = new TH1 * [10];
  h_F3 = new TH1 * [20];
  h_F4 = new TH1 * [35];
  h_R2 = new TH1 * [10];
  h_R3 = new TH1 * [20];
  h_R4 = new TH1 * [35];
  h_nudyn = new TH1 * [10];

  if (configuration->nuDynVsMult)
    {
    h_F2_vsMult = new TH1 * [10];
    h_F3_vsMult = new TH1 * [20];
    h_F4_vsMult = new TH1 * [35];
    h_R2_vsMult = new TH1 * [10];
    h_R3_vsMult = new TH1 * [20];
    h_R4_vsMult = new TH1 * [35];
    h_nudyn_vsMult = new TH1 * [10];
    }

  if (configuration->nuDynVsCent)
    {
    h_F2_vsCent = new TH1 * [10];
    h_F3_vsCent = new TH1 * [20];
    h_F4_vsCent = new TH1 * [35];
    h_R2_vsCent = new TH1 * [10];
    h_R3_vsCent = new TH1 * [20];
    h_R4_vsCent = new TH1 * [35];
    h_nudyn_vsCent = new TH1 * [10];
    }

  initialize();
}

NuDynDerivedHistos::NuDynDerivedHistos(TFile * inputFile,
                                       const TString & name,
                                       AnalysisConfiguration * configuration,
                                       LogLevel  debugLevel)
:
Histograms(name,configuration,200,debugLevel)
{
  h_F2 = new TH1 * [10];
  h_F3 = new TH1 * [20];
  h_F4 = new TH1 * [35];
  h_R2 = new TH1 * [10];
  h_R3 = new TH1 * [20];
  h_R4 = new TH1 * [35];
  h_nudyn = new TH1 * [10];

  if (configuration->nuDynVsMult)
    {
    h_F2_vsMult = new TH1 * [10];
    h_F3_vsMult = new TH1 * [20];
    h_F4_vsMult = new TH1 * [35];
    h_R2_vsMult = new TH1 * [10];
    h_R3_vsMult = new TH1 * [20];
    h_R4_vsMult = new TH1 * [35];
    h_nudyn_vsMult = new TH1 * [10];
    }

  if (configuration->nuDynVsCent)
    {
    h_F2_vsCent = new TH1 * [10];
    h_F3_vsCent = new TH1 * [20];
    h_F4_vsCent = new TH1 * [35];
    h_R2_vsCent = new TH1 * [10];
    h_R3_vsCent = new TH1 * [20];
    h_R4_vsCent = new TH1 * [35];
    h_nudyn_vsCent = new TH1 * [10];
    }

  loadHistograms(inputFile);
}

NuDynDerivedHistos::~NuDynDerivedHistos()
{
  cout << "NuDynDerivedHistos::~NuDynDerivedHistos() was called." << endl;

}

// for now use the same boundaries for eta and y histogram
void NuDynDerivedHistos::createHistograms()
{
  AnalysisConfiguration & ac = *getConfiguration();
  TString bn = getHistoBaseName();

  int nPart =4;
  TString histName;
  TString histTitle;
  for (int i1=0; i1<nPart; i1++)
    {
    for (int i2=i1; i2<nPart; i2++)
      {
      int index12 = index2(i1,i2);
      histName  = bn + "F2_";
      histName  += i1;
      histName  += i2;
      histTitle = "F_{2}";
      histTitle += "^{";
      histTitle += i1;
      histTitle += i2;
      histTitle += "}";
      h_F2[index12]  = createHistogram(histName, 1,ac.min_mult,ac.max_mult,"mult", histTitle);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "F2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_F2_vsMult[index12]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult", histTitle);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "F2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_F2_vsCent[index12]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent", histTitle);
        }
      histName  = bn + "R2_";
      histName  += i1;
      histName  += i2;
      histTitle = "R_{2}";
      histTitle += "^{";
      histTitle += i1;
      histTitle += i2;
      histTitle += "}";
      h_R2[index12]  = createHistogram(histName, 1,ac.min_mult,ac.max_mult,"mult", histTitle);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "R2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_R2_vsMult[index12]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult",histTitle);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "R2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_R2_vsCent[index12]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent",histTitle);
        }

      histName  = bn + "nudyn_";
      histName  += i1;
      histName  += i2;
      histTitle = "#nu_{dyn}";
      histTitle += "^{";
      histTitle += i1;
      histTitle += i2;
      histTitle += "}";
      h_nudyn[index12]  = createHistogram(histName, 1,ac.min_mult,ac.max_mult,"mult", histTitle);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "nudyn_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_nudyn_vsMult[index12]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult",histTitle);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "nudyn_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_nudyn_vsCent[index12]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent",histTitle);
        }

      for (int i3=i2; i3<nPart; i3++)
        {
        int index123 = index3(i1,i2,i3);

        histName  = bn + "F3_";
        histName  += i1;
        histName  += i2;
        histName  += i3;
        histTitle = "F_{3}";
        histTitle += "^{";
        histTitle += i1;
        histTitle += i2;
        histTitle += i3;
        histTitle += "}";
        h_F3[index123]  = createHistogram(histName, 1,ac.min_mult,ac.max_mult,"mult", histTitle);
        if (ac.nuDynVsMult)
          {
          histName  = bn + "F3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsMult";
          h_F3_vsMult[index123]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult", histTitle);
          }
        if (ac.nuDynVsCent)
          {
          histName  = bn + "F3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsCent";
          h_F3_vsCent[index123]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent", histTitle);
          }
        histName  = bn + "R3_";
        histName  += i1;
        histName  += i2;
        histName  += i3;
        histTitle = "R_{3}";
        histTitle += "^{";
        histTitle += i1;
        histTitle += i2;
        histTitle += i3;
        histTitle += "}";
        h_R3[index123]  = createHistogram(histName, 1,ac.min_mult,ac.max_mult,"mult", histTitle);
        if (ac.nuDynVsMult)
          {
          histName  = bn + "R3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsMult";
          h_R3_vsMult[index123]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult", histTitle);
          }
        if (ac.nuDynVsCent)
          {
          histName  = bn + "R3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsCent";
          h_R3_vsCent[index123]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent", histTitle);
          }
        for (int i4=i3; i4<nPart; i4++)
          {
          int index1234 = index4(i1,i2,i3,i4);

          histName  = bn + "F4_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName  += i4;
          histTitle = "F_{4}";
          histTitle += "^{";
          histTitle += i1;
          histTitle += i2;
          histTitle += i3;
          histTitle += i4;
          histTitle += "}";
          h_F4[index1234]  = createHistogram(histName, 1,ac.min_mult,ac.max_mult,"mult", histTitle);
          if (ac.nuDynVsMult)
            {
            histName  = bn + "F4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsMult";
            h_F4_vsMult[index1234]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult", histTitle);
            }
          if (ac.nuDynVsCent)
            {
            histName  = bn + "F4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsCent";
            h_F4_vsCent[index1234]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent", histTitle);
            }
          histName  = bn + "R4_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName  += i4;
          histTitle = "R_{4}";
          histTitle += "^{";
          histTitle += i1;
          histTitle += i2;
          histTitle += i3;
          histTitle += i4;
          histTitle += "}";
          h_R4[index1234]  = createHistogram(histName, 1,ac.min_mult,ac.max_mult,"mult", histTitle);
          if (ac.nuDynVsMult)
            {
            histName  = bn + "R4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsMult";
            h_R4_vsMult[index1234]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult,"mult",  histTitle);
            }
          if (ac.nuDynVsCent)
            {
            histName  = bn + "R4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsCent";
            h_R4_vsCent[index1234]  = createHistogram(histName,ac.nBins_mult,ac.min_mult,  ac.max_mult, "cent", histTitle);
            }
          }
        }
      }
    }
}

//________________________________________________________________________
void NuDynDerivedHistos::loadHistograms(TFile * inputFile)
{
  if (!inputFile)
    {
    if (reportFatal()) cout << "-Fatal- Attempting to load NuDynDerivedHistos from an invalid file pointer" << endl;
    return;
    }
  AnalysisConfiguration & ac = *getConfiguration();
  TString  bn = getHistoBaseName();
  int nPart =4;
  TString histName;
  TString histTitle;
  for (int i1=0; i1<nPart; i1++)
    {
    for (int i2=i1; i2<nPart; i2++)
      {
      int index12 = index2(i1,i2);
      histName  = bn + "F2_";
      histName  += i1;
      histName  += i2;
      histTitle = "F_{2}";
      histTitle += "^{";
      histTitle += i1;
      histTitle += i2;
      histTitle += "}";
      h_F2[index12]  = loadH1(inputFile,histName);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "F2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_F2_vsMult[index12]  = loadH1(inputFile,histName);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "F2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_F2_vsMult[index12]  = loadH1(inputFile,histName);
        }
      histName  = bn + "R2_";
      histName  += i1;
      histName  += i2;
      histTitle = "R_{2}";
      histTitle += "^{";
      histTitle += i1;
      histTitle += i2;
      histTitle += "}";
      h_R2[index12]  = loadH1(inputFile,histName);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "R2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_R2_vsMult[index12]  = loadH1(inputFile,histName);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "R2_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_R2_vsCent[index12]  = loadH1(inputFile,histName);
        }

      histName  = bn + "nudyn_";
      histName  += i1;
      histName  += i2;
      histTitle = "#nu_{dyn}";
      h_nudyn[index12]  = loadH1(inputFile,histName);
      if (ac.nuDynVsMult)
        {
        histName  = bn + "nudyn_";
        histName  += i1;
        histName  += i2;
        histName += "_vsMult";
        h_nudyn_vsMult[index12]  = loadH1(inputFile,histName);
        }
      if (ac.nuDynVsCent)
        {
        histName  = bn + "nudyn_";
        histName  += i1;
        histName  += i2;
        histName += "_vsCent";
        h_nudyn_vsCent[index12]  = loadH1(inputFile,histName);
        }

      for (int i3=i2; i3<nPart; i3++)
        {
        int index123 = index3(i1,i2,i3);
        histName  = bn + "F3_";
        histName  += i1;
        histName  += i2;
        histName  += i3;
        histTitle = "F_{3}";
        histTitle += "^{";
        histTitle += i1;
        histTitle += i2;
        histTitle += i3;
        histTitle += "}";
        h_F3[index123]  = loadH1(inputFile,histName);
        if (ac.nuDynVsMult)
          {
          histName  = bn + "F3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsMult";
          h_F3_vsMult[index123]  = loadH1(inputFile,histName);
          }
        if (ac.nuDynVsCent)
          {
          histName  = bn + "F3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsCent";
          h_F3_vsCent[index123]  = loadH1(inputFile,histName);
          }
        histName  = bn + "R3_";
        histName  += i1;
        histName  += i2;
        histName  += i3;
        histTitle = "R_{3}";
        histTitle += "^{";
        histTitle += i1;
        histTitle += i2;
        histTitle += i3;
        histTitle += "}";
        h_R3[index123]  = loadH1(inputFile,histName);
        if (ac.nuDynVsMult)
          {
          histName  = bn + "R3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsMult";
          h_R3_vsMult[index123]  = loadH1(inputFile,histName);
          }
        if (ac.nuDynVsCent)
          {
          histName  = bn + "R3_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName += "_vsCent";
          h_R3_vsCent[index123]  = loadH1(inputFile,histName);
          }
        for (int i4=i3; i4<nPart; i4++)
          {
          int index1234 = index4(i1,i2,i3,i4);
          histTitle = "F_{4}";
          histTitle += "^{";
          histTitle += i1;
          histTitle += i2;
          histTitle += i3;
          histTitle += i4;
          histTitle += "}";

          histName  = bn + "F4_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName  += i4;

          //h_F4[index1234]  = loadH1(inputFile,histName);
          if (ac.nuDynVsMult)
            {
            histName  = bn + "F4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsMult";
            h_F4_vsMult[index1234]  = loadH1(inputFile,histName);
            }
          if (ac.nuDynVsCent)
            {
            histName  = bn + "F4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsCent";
            //h_F4_vsCent[index1234]  = loadH1(inputFile,histName);
            }

          histTitle = "R_{4}";
          histTitle += "^{";
          histTitle += i1;
          histTitle += i2;
          histTitle += i3;
          histTitle += i4;
          histTitle += "}";

          histName  = bn + "R4_";
          histName  += i1;
          histName  += i2;
          histName  += i3;
          histName  += i4;

          //h_R4[index1234]  = loadH1(inputFile,histName);
          if (ac.nuDynVsMult)
            {
            histName  = bn + "R4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsMult";
            h_R4_vsMult[index1234]  = loadH1(inputFile,histName);
            }
          if (ac.nuDynVsCent)
            {
            histName  = bn + "R4_";
            histName  += i1;
            histName  += i2;
            histName  += i3;
            histName  += i4;
            histName += "_vsCent";
            //h_R4_vsCent[index1234]  = loadH1(inputFile,histName);
            }
          }
        }
      }
    }
  if (reportInfo()) cout << "Load completed." << endl;
}

void NuDynDerivedHistos::calculateDerivedHistograms(NuDynHistos * source)
{
  // h_r2_11 = <n1(n1-1)>/<n1><n1>
  // h_r2_22 = <n2(n2-1)>/<n2><n2>
  // h_r2_12 = <n1(n2-1)>/<n1><n2>
  // h_nudyn_12 = h_r2_11 + h_r2_22 -2*h_r2_12

  if (!source)
    {
     if (reportFatal()) cout << "-Fatal- NuDynDerivedHistos::calculateDerivedHistograms(NuDynHistos * source) Source is null." << endl;
     return;
     }

  AnalysisConfiguration & ac = *getConfiguration();
  int nPart = 4;
  for (int i1=0; i1<nPart; i1++)
    {
    for (int i2=i1; i2<nPart; i2++)
      {
      int index12 = index2(i1,i2);
      calculateF2R2(source->h_f1[i1],source->h_f1[i2],source->h_f2[index12], h_F2[index12], h_R2[index12]);
      if (ac.nuDynVsMult) calculateF2R2(source->h_f1_vsMult[i1],source->h_f1_vsMult[i2],source->h_f2_vsMult[index12], h_F2_vsMult[index12], h_R2_vsMult[index12]);
      if (ac.nuDynVsCent) calculateF2R2(source->h_f1_vsCent[i1],source->h_f1_vsCent[i2],source->h_f2_vsCent[index12], h_F2_vsCent[index12], h_R2_vsCent[index12]);
      for (int i3=i2; i3<nPart; i3++)
        {
        int index13 = index2(i1,i3);
        int index23 = index2(i2,i3);
        int index123 = index3(i1,i2,i3);
        calculateF3R3(source->h_f1[i1],source->h_f1[i2],source->h_f1[i3],
                      source->h_f2[index12],source->h_f2[index13],source->h_f2[index23],
                      source->h_f3[index123],
                      h_F3[index123], h_R3[index123]);
        if (ac.nuDynVsMult) calculateF3R3(source->h_f1_vsMult[i1],source->h_f1_vsMult[i2],source->h_f1_vsMult[i3],
                                          source->h_f2_vsMult[index12],source->h_f2_vsMult[index13],source->h_f2_vsMult[index23],
                                          source->h_f3_vsMult[index123],
                                          h_F3_vsMult[index123], h_R3_vsMult[index123]);
        if (ac.nuDynVsCent) calculateF3R3(source->h_f1_vsCent[i1],source->h_f1_vsCent[i2],source->h_f1_vsCent[i3],
                                          source->h_f2_vsCent[index12],source->h_f2_vsCent[index13],source->h_f2_vsCent[index23],
                                          source->h_f3_vsCent[index123],
                                          h_F3_vsCent[index123], h_R3_vsCent[index123]);

        for (int i4=i3; i4<nPart; i4++)
          {
          int index14 = index2(i1,i4);
          int index24 = index2(i2,i4);
          int index34 = index2(i3,i4);
          int index124 = index3(i1,i2,i4);
          int index134 = index3(i1,i3,i4);
          int index234 = index3(i2,i3,i4);
          int index1234 = index4(i1,i2,i3,i4);
          calculateF4R4(source->h_f1[i1],source->h_f1[i2],source->h_f1[i3],source->h_f1[i4],
                        source->h_f2[index12],source->h_f2[index13],source->h_f2[index14],
                        source->h_f2[index23],source->h_f2[index24],source->h_f2[index34],
                        source->h_f3[index123],source->h_f3[index124],source->h_f3[index134],source->h_f3[index234],
                        source->h_f4[index1234],
                        h_F4[index1234], h_R4[index1234]);
          if (ac.nuDynVsMult) calculateF4R4(source->h_f1_vsMult[i1],source->h_f1_vsMult[i2],source->h_f1_vsMult[i3],source->h_f1_vsMult[i4],
                                            source->h_f2_vsMult[index12], source->h_f2_vsMult[index13], source->h_f2_vsMult[index14],
                                            source->h_f2_vsMult[index23], source->h_f2_vsMult[index24], source->h_f2_vsMult[index34],
                                            source->h_f3_vsMult[index123],source->h_f3_vsMult[index124],source->h_f3_vsMult[index134],source->h_f3_vsMult[index234],
                                            source->h_f4_vsMult[index1234],
                                            h_F4_vsMult[index1234], h_R4_vsMult[index1234]);
          if (ac.nuDynVsCent) calculateF4R4(source->h_f1_vsCent[i1],source->h_f1_vsCent[i2],source->h_f1_vsCent[i3],source->h_f1_vsCent[i4],
                                            source->h_f2_vsCent[index12],  source->h_f2_vsCent[index13], source->h_f2_vsCent[index14],
                                            source->h_f2_vsCent[index23],  source->h_f2_vsCent[index24], source->h_f2_vsCent[index34],
                                            source->h_f3_vsCent[index123], source->h_f3_vsCent[index124],source->h_f3_vsCent[index134],source->h_f3_vsCent[index234],
                                            source->h_f4_vsCent[index1234],
                                            h_F4_vsCent[index1234], h_R4_vsCent[index1234]);
          }
        }
      }
    }
  calculateNudyn(h_R2[index2(0,0)],h_R2[index2(0,1)],h_R2[index2(1,1)], h_nudyn[index2(0,1)]);
  calculateNudyn(h_R2[index2(0,0)],h_R2[index2(0,2)],h_R2[index2(2,2)], h_nudyn[index2(0,2)]);
  calculateNudyn(h_R2[index2(0,0)],h_R2[index2(0,3)],h_R2[index2(3,3)], h_nudyn[index2(0,3)]);
  calculateNudyn(h_R2[index2(1,1)],h_R2[index2(1,2)],h_R2[index2(2,2)], h_nudyn[index2(1,2)]);
  calculateNudyn(h_R2[index2(1,1)],h_R2[index2(1,3)],h_R2[index2(3,3)], h_nudyn[index2(1,3)]);
  calculateNudyn(h_R2[index2(2,2)],h_R2[index2(2,3)],h_R2[index2(3,3)], h_nudyn[index2(2,3)]);
}
