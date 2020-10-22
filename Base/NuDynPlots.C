
#include "TString.h"
#include "TRandom.h"
//#include "TFile.h"
//#include "AnalysisConfiguration.hpp"
//#include "NuDynHistos.hpp"
//#include "NuDynDerivedHistos.hpp"
//#include "CanvasConfiguration.hpp"
//#include "HistogramCollection.hpp"
//#include "GraphConfiguration.hpp"
//#include "CanvasConfiguration.hpp"
//#include "TRint.h"


//R__LOAD_LIBRARY(/Users/claudeapruneau/opt/WAC/lib/libBase.dylib)

int NuDynPlots()
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
  TString inputPath  = "/Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/";
  TString outputPath = "/Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/";
  int nModels = 5;
  TFile ** inputFiles  = new TFile* [nModels];
  NuDynHistos  ** nuDynHistos = new NuDynHistos*[nModels];
  NuDynDerivedHistos  ** nuDynDerivedHistos = new NuDynDerivedHistos*[nModels];
  NuDynPlotter ** plotter      = new NuDynPlotter*[nModels];
  //TString fileBaseName  = "resultsPYTHIA_NuDyn_";
  TString fileBaseName  = "PYTHIA_softOnHardOff_NuDyn_";
  TString histoName, histoTitle;
  TString fileName;
  TString ** setName  = new TString*[nModels];
  TString ** setTitle = new TString*[nModels];
  TString ** histoBaseNames    = new TString*[nModels];
  TString ** histoBaseTitles   = new TString*[nModels];

  histoBaseNames[0] = new TString("HPHPHPHP");
  histoBaseNames[1] = new TString("HPHPHPHM");
  histoBaseNames[2] = new TString("HPHPHMHM");
  histoBaseNames[3] = new TString("HPHMHMHM");
  histoBaseNames[4] = new TString("HMHMHMHM");

  histoBaseTitles[0] = new TString("h^{+}h^{+}h^{+}h^{+}");
  histoBaseTitles[1] = new TString("h^{+}h^{+}h^{+}h^{-}");
  histoBaseTitles[2] = new TString("h^{+}h^{+}h^{-}h^{-}");
  histoBaseTitles[3] = new TString("h^{+}h^{-}h^{-}h^{-}");
  histoBaseTitles[4] = new TString("h^{-}h^{-}h^{-}h^{-}");

  AnalysisConfiguration * ac = new AnalysisConfiguration("","","");
  ac->nuDynVsMult = true;
  ac->nuDynVsCent = false;

  for (int iModel=0;iModel<nModels;iModel++)
    //for (int iModel=0;iModel<1;iModel++)
    {
    fileName = fileBaseName;
    fileName += *histoBaseNames[iModel];
    setName[iModel]    = new TString(fileName);
    setTitle[iModel]   = histoBaseTitles[iModel];
    cout << "iModel:" << iModel << " Title:" <<  *setTitle[iModel]  << endl;
    inputFiles[iModel] = new TFile(inputPath+fileName+"_MB.root ","OLD");
    if (!inputFiles[iModel])
      {
      cout << "<F> Could not open file: " << inputPath+fileName << endl;
      return 1;
      }
    else
      {
      cout << "<I> Succesfully opened file: " << inputPath+fileName << endl;
      }
    nuDynHistos[iModel]        = new NuDynHistos(inputFiles[iModel], *histoBaseNames[iModel], ac, MessageLogger::Debug);
    nuDynDerivedHistos[iModel] = new NuDynDerivedHistos(inputFiles[iModel], *histoBaseNames[iModel], ac, MessageLogger::Debug);
    plotter[iModel]            = new NuDynPlotter();
    plotter[iModel]->makePlots(*setName[iModel],nuDynHistos[iModel],nuDynDerivedHistos[iModel],
                               canvasConfiguration,graphConfigurations);
    plotter[iModel]->printAllCanvas(outputPath);
    }

//    NuDynPlotter * compPlotter = new NuDynPlotter(false);
//    compPlotter->makeComparisonPlots("ModelComparison_HardOff_", canvasConfiguration,nModels,nuDynHistos,nuDynDerivedHistos,setTitle,graphConfigurations);
//    compPlotter->printAllCanvas(outputPath);
//
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
//
//  //TCanvas * cc = new TCanvas();
//  //nuDyn->Draw();
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

