
#include "TString.h"
#include "TRandom.h"

int NuDynPlotsHardVsSoft()
{
  TString includesPath = "/Users/claudeapruneau/Documents/GitHub/WAC/Base/";
  gSystem->Load(includesPath+"CanvasCollection.hpp");
  gSystem->Load(includesPath+"CanvasConfiguration.hpp");
  gSystem->Load(includesPath+"TaskConfiguration.hpp");
  gSystem->Load(includesPath+"EnhancedGraph.hpp");
  gSystem->Load(includesPath+"Factory.hpp");
  gSystem->Load(includesPath+"FunctionCollection.hpp");
  gSystem->Load(includesPath+"GraphConfiguration.hpp");
  gSystem->Load(includesPath+"HistogramCollection.hpp");
  gSystem->Load(includesPath+"Histograms.hpp");
  gSystem->Load(includesPath+"RandomGenerators.hpp");
  gSystem->Load(includesPath+"RapidityGenerator.hpp");
  gSystem->Load(includesPath+"Task.hpp");
  gSystem->Load(includesPath+"TaskCollection.hpp");
  gSystem->Load(includesPath+"Property.hpp");
  gSystem->Load(includesPath+"MessageLogger.hpp");
  gSystem->Load(includesPath+"AnalysisConfiguration.hpp");
  gSystem->Load(includesPath+"Event.hpp");
  gSystem->Load(includesPath+"EventFilter.hpp");
  gSystem->Load(includesPath+"EventHistos.hpp");
  gSystem->Load(includesPath+"EventLoop.hpp");
  gSystem->Load(includesPath+"GeneralizedGaussianConfiguration.hpp");
  gSystem->Load(includesPath+"GeneralizedGaussianProfile.hpp");
  gSystem->Load(includesPath+"GeneratorConfiguration.hpp");
  gSystem->Load(includesPath+"TwoPartCorrelationAnalyzer.hpp");
  gSystem->Load(includesPath+"Particle.hpp");
  gSystem->Load(includesPath+"ParticleFilter.hpp");
  gSystem->Load(includesPath+"ParticleHistos.hpp");
  gSystem->Load(includesPath+"ParticlePairCombinedHistos.hpp");
  gSystem->Load(includesPath+"ParticlePairDerivedHistos.hpp");
  gSystem->Load(includesPath+"ParticlePairFilter.hpp");
  gSystem->Load(includesPath+"ParticlePairHistos.hpp");
  gSystem->Load(includesPath+"TrackAndPairConfiguration.hpp");
  gSystem->Load(includesPath+"NuDynTask.hpp");
  gSystem->Load(includesPath+"NuDynHistos.hpp");
  gSystem->Load(includesPath+"NuDynDerivedHistos.hpp");
  gSystem->Load(includesPath+"Plotter.hpp");
  gSystem->Load(includesPath+"NuDynPlotter.hpp");
  gSystem->Load("libBase.dylib");

  HistogramCollection * histogramCollection = new HistogramCollection("Collection",200);
  histogramCollection->setDefaultOptions(1);
  CanvasCollection    * canvasCollection    = new CanvasCollection();
  CanvasConfiguration * canvasConfiguration = new CanvasConfiguration(CanvasConfiguration::Landscape,CanvasConfiguration::Linear);
  GraphConfiguration  ** graphConfigurations = new GraphConfiguration*[40];
  for (int iGraph=0;iGraph<40;iGraph++)
    {
    graphConfigurations[iGraph] = new GraphConfiguration(1,iGraph%10);
    }

  TString inputPath  = "/Users/claudeapruneau/Documents/GitHub/run/NuDynStudies/";
  TString outputPath = "/Users/claudeapruneau/Documents/GitHub/run/NuDynStudies/";


  int nModels            = 2;
  int nSets              = 5;
  int nFiles             = nModels*nSets;
  TString ** fileBaseNames     = new TString*[nSets];
  fileBaseNames[0]  = new TString("resultsPYTHIA_NuDyn_");
  fileBaseNames[1]  = new TString("hardOffPYTHIA_NuDyn_");

  TString ** fileNames         = new TString*[nFiles];
  TFile   ** inputFiles        = new TFile*  [nFiles];
  NuDynHistos ** nuDynHistos   = new NuDynHistos*[nFiles];
  NuDynDerivedHistos ** nuDynDerivedHistos = new NuDynDerivedHistos*[nFiles];

  TString histoName, histoTitle;

  TString ** setNames          = new TString*[nFiles];
  TString ** setTitles         = new TString*[nFiles];
  TString ** setBaseNames    = new TString*[nModels];
  setBaseNames[0] = new TString("HPHPHPHP");
  setBaseNames[1] = new TString("HPHPHPHM");
  setBaseNames[2] = new TString("HPHPHMHM");
  setBaseNames[3] = new TString("HPHMHMHM");
  setBaseNames[4] = new TString("HMHMHMHM");
  TString ** setBaseTitles   = new TString*[nModels];
  setBaseTitles[0] = new TString("h^{+}h^{+}h^{+}h^{+}");
  setBaseTitles[1] = new TString("h^{+}h^{+}h^{+}h^{-}");
  setBaseTitles[2] = new TString("h^{+}h^{+}h^{-}h^{-}");
  setBaseTitles[3] = new TString("h^{+}h^{-}h^{-}h^{-}");
  setBaseTitles[4] = new TString("h^{-}h^{-}h^{-}h^{-}");

  AnalysisConfiguration * ac = new AnalysisConfiguration("","","");
  ac->nuDynVsMult = true;
  ac->nuDynVsCent = false;

  // Load all files of interest
  int iFile;
  for (int iModel=0;iModel<nModels;iModel++)
    {
    for (int iSet=0;iSet<nSets;iSet++)
      {
      iFile = iModel*nSets+iSet;
      TString workString = *fileBaseNames[iModel];
      workString += *setBaseNames[iSet];
      fileNames[iFile] = new TString(workString);
      workString = *fileNames[iFile];
      setNames[iFile]    = new TString(workString);
      setTitles[iFile]   = setBaseTitles[iSet];
      cout << "iModel:" << iModel << " iSet:" << iSet << " Title:" <<  *setTitles[iFile]  << endl;
      workString = inputPath;
      workString += *fileNames[iFile];
      workString += "_MB.root";
      inputFiles[iFile] = new TFile(workString,"OLD");
      if (!inputFiles[iFile])
        {
        cout << "<F> Could not open file: " << workString << endl;
        return 1;
        }
      else
        {
        cout << "<I> Succesfully opened file: " << workString << endl;
        }
      nuDynHistos[iFile]        = new NuDynHistos(inputFiles[iFile], *setBaseNames[iSet], ac, MessageLogger::Debug);
      nuDynDerivedHistos[iFile] = new NuDynDerivedHistos(inputFiles[iFile], *setBaseNames[iSet], ac, MessageLogger::Debug);
      }
    }

//  int index12;
//  index12 = nuDynHistos[0]->index2(0,0);
//  TH1* r2_11 = nuDynDerivedHistos[3]->h_R2_vsMult[index12];
//  if (!r2_11) { cout <<"r2_11 does not exist " << endl; return -1;}
//  index12 = nuDynHistos[0]->index2(0,1);
//  TH1* r2_12 = nuDynDerivedHistos[3]->h_R2_vsMult[index12];
//  if (!r2_12) { cout <<"r2_12 does not exist " << endl; return -1;}
//  index12 = nuDynHistos[0]->index2(1,1);
//  TH1* r2_22 = nuDynDerivedHistos[3]->h_R2_vsMult[index12];
//  if (!r2_22) { cout <<"r2_22 does not exist " << endl; return -1;}
//  index12 = nuDynHistos[0]->index2(0,1);
//  TH1* nuDyn1 = nuDynDerivedHistos[3]->h_nudyn_vsMult[0];
//  if (!nuDyn1) { cout <<"nuDyn does not exist 1" << endl; return -1;}
//  histogramCollection->calculateNudyn(r2_11,r2_12,r2_22,nuDyn1);

  TH1** hh = new TH1*[8];
  Plotter * p;

// f_1
//  hh[0] = nuDynHistos[0]->h_f1_vsMult[0];
//  hh[1] = nuDynHistos[nSets]->h_f1_vsMult[0];
//  if (!hh[0])
//    {
//    cout << "h_f1_vsMult[0] is not loaded" << endl;
//    }
//  if (!hh[1])
//     {
//     cout << "h_f1_vsMult[0] is not loaded" << endl;
//     }
//  TString ** leg_f1   = new TString*[2];
//
//  leg_f1[0] = new TString("f_{1} Hard+Soft");
//  leg_f1[1] = new TString("f_{1} Soft Only");
//  p = new Plotter(false);
//  p->plot(2,"NuDyn_AllvsHardOff_f1",canvasConfiguration, graphConfigurations,
//          "Mult", 0.0, 400.0,
//          "f_{1}", 1.0, -1.0,
//          hh,leg_f1,
//          0.2,0.6, 0.5, 0.85, 0.05);
//  p->printAllCanvas(outputPath);


  // f_2
//  hh[0] = nuDynHistos[0]->h_f2_vsMult[0];
//  hh[1] = nuDynHistos[nSets]->h_f2_vsMult[0];
//  if (!hh[0])
//    {
//    cout << "h_f2_vsMult[0] is not loaded" << endl;
//    }
//  if (!hh[1])
//     {
//     cout << "h_f2_vsMult[0] is not loaded" << endl;
//     }
//  TString ** leg_f2   = new TString*[2];
//
//  leg_f2[0] = new TString("f_{2} Hard+Soft");
//  leg_f2[1] = new TString("f_{2} Soft Only");
//  p = new Plotter(false);
//  p->plot(2,"NuDyn_AllvsHardOff_f2",canvasConfiguration, graphConfigurations,
//          "Mult", 0.0, 400.0,
//          "f_{2}", 1.0, -1.0,
//          hh,leg_f2,
//          0.2,0.6, 0.5, 0.85, 0.05);
//  p->printAllCanvas(outputPath);

  // F_2
//  hh[0] = nuDynDerivedHistos[0]->h_F2_vsMult[0];
//  hh[1] = nuDynDerivedHistos[nSets]->h_F2_vsMult[0];
//  if (!hh[0])
//    {
//    cout << "h_F2_vsMult[0] is not loaded" << endl;
//    }
//  if (!hh[1])
//     {
//     cout << "h_F2_vsMult[0] is not loaded" << endl;
//     }
//  TString ** leg_F2   = new TString*[2];
//
//  leg_F2[0] = new TString("F_{2} Hard+Soft");
//  leg_F2[1] = new TString("F_{2} Soft Only");
//  p = new Plotter(false);
//  p->plot(2,"NuDyn_AllvsHardOff_F2",canvasConfiguration, graphConfigurations,
//          "Mult", 0.0, 400.0,
//          "F_{2}", -6.0, 1.0,
//          hh,leg_F2,
//          0.2,0.2, 0.5, 0.45, 0.05);
//  p->printAllCanvas(outputPath);

  // R_2
   int index12;
   index12 = nuDynHistos[0]->index2(0,1);

   hh[0] = nuDynDerivedHistos[0]->h_R2_vsMult[0];
   hh[1] = nuDynDerivedHistos[3]->h_R2_vsMult[index12];
   hh[2] = nuDynDerivedHistos[nSets]->h_R2_vsMult[0];
   hh[3] = nuDynDerivedHistos[nSets+3]->h_R2_vsMult[index12];
   if (!hh[0])
     {
     cout << "h_R2_vsMult[0] is not loaded" << endl;
     }
   if (!hh[1])
      {
      cout << "h_R2_vsMult[0] is not loaded" << endl;
      }
   TString ** leg_R2   = new TString*[4];
   leg_R2[0] = new TString("R_{2}^{++} Hard+Soft");
   leg_R2[1] = new TString("R_{2}^{+-} Hard+Soft");
   leg_R2[2] = new TString("R_{2}^{++}  Soft Only");
   leg_R2[3] = new TString("R_{2}^{+-}  Soft Only");
   p = new Plotter(false);
   p->plot(4,"NuDyn_AllvsHardOff_R2",canvasConfiguration, graphConfigurations,
           "Mult", 0.0, 400.0,
           "R_{2}", -0.3, 1.5,
           hh,leg_R2,
           0.6,0.5, 0.85, 0.85, 0.05);
   p->printAllCanvas(outputPath);


  // R_4
   int index1234;
   index1234 = nuDynHistos[0]->index4(0,1,2,3);

   hh[0] = nuDynDerivedHistos[0]->h_R4_vsMult[0];
   hh[1] = nuDynDerivedHistos[1]->h_R4_vsMult[index12];
   hh[2] = nuDynDerivedHistos[2]->h_R4_vsMult[index12];
   hh[3] = nuDynDerivedHistos[3]->h_R4_vsMult[index12];
   hh[4] = nuDynDerivedHistos[nSets]->h_R4_vsMult[0];
   hh[5] = nuDynDerivedHistos[nSets+1]->h_R4_vsMult[index12];
   hh[6] = nuDynDerivedHistos[nSets+2]->h_R4_vsMult[index12];
   hh[7] = nuDynDerivedHistos[nSets+3]->h_R4_vsMult[index12];
   if (!hh[0])
     {
     cout << "h_R4_vsMult[0] is not loaded" << endl;
     }
   if (!hh[1])
      {
      cout << "h_R4_vsMult[0] is not loaded" << endl;
      }

//  setBaseNames[0] = new TString("HPHPHPHP");
//  setBaseNames[1] = new TString("HPHPHPHM");
//  setBaseNames[2] = new TString("HPHPHMHM");
//  setBaseNames[3] = new TString("HPHMHMHM");
//  setBaseNames[4] = new TString("HMHMHMHM");

   TString ** leg_R4   = new TString*[8];
   leg_R4[0] = new TString("R_{4}^{++++} Hard+Soft");
   leg_R4[1] = new TString("R_{4}^{+++-} Hard+Soft");
   leg_R4[2] = new TString("R_{4}^{++--} Hard+Soft");
   leg_R4[3] = new TString("R_{4}^{+---} Hard+Soft");
   leg_R4[4] = new TString("R_{4}^{++++}  Soft Only");
   leg_R4[5] = new TString("R_{4}^{+++-}  Soft Only");
   leg_R4[6] = new TString("R_{4}^{++--}  Soft Only");
   leg_R4[7] = new TString("R_{4}^{+---}  Soft Only");
   p = new Plotter(false);
   p->plot(8,"NuDyn_AllvsHardOff_R4",canvasConfiguration, graphConfigurations,
           "Mult", 0.0, 400.0,
           "R_{4}", -1.0, 2.5,
           hh,leg_R4,
           0.6,0.5, 0.85, 0.85, 0.05);
   p->printAllCanvas(outputPath);


return -666;
  //TCanvas * cc = new TCanvas();
  //nuDyn->Draw();
//
//  index12 = nuDynHistos[0]->index2(0,0);
//  r2_11 = nuDynDerivedHistos[4]->h_R2_vsMult[index12];
//  index12 = nuDynHistos[0]->index2(0,1);
//  r2_12 = nuDynDerivedHistos[4]->h_R2_vsMult[index12];
//  index12 = nuDynHistos[0]->index2(1,1);
//  r2_22 = nuDynDerivedHistos[4]->h_R2_vsMult[index12];
//  index12 = nuDynHistos[0]->index2(0,1);
//  TH1* nuDyn2 = nuDynDerivedHistos[4]->h_nudyn_vsMult[index12];
//  if (!nuDyn2) { cout <<"nuDyn does not exist 2" << endl; return -1;}
//  histogramCollection->calculateNudyn(r2_11,r2_12,r2_22,nuDyn2);
//  //nuDyn->Draw("SAME");
//
//  TH1** hh = new TH1*[2];
//  hh[0] = nuDyn1;
//  hh[1] = nuDyn2;
//
//  TString ** leg   = new TString*[2];
//  leg[0] = new TString("PM");
//  leg[1] = new TString("MM");
//  Plotter * p  = new Plotter(false);
//  p->plot(2,"NuDynComp_HardOff",canvasConfiguration, graphConfigurations,
//          "Mult", 0.0, 400.0,
//          "#nu_{dyn}", -1.0, 0.0,
//          hh,leg,
//          0.7,0.6, 0.9, 0.85, 0.05);
//  p->printAllCanvas(outputPath);
//  cout << "<I> NuDynPlots() All done...." << endl;

  return 0;
}

