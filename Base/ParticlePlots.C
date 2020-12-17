
#include "TString.h"
#include "TRandom.h"
#include "ParticlePlotter.hpp"
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


int ParticlePlots()
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
  gSystem->Load(includesPath+"ParticleAnalyzerConfiguration.hpp");
  gSystem->Load(includesPath+"ParticleHistos.hpp");
  gSystem->Load(includesPath+"Plotter.hpp");
  gSystem->Load(includesPath+"ParticlePlotter.hpp");
  gSystem->Load("libBase.dylib");

  HistogramCollection * histogramCollection = new HistogramCollection("Collection",200);
  histogramCollection->setDefaultOptions(1);
  CanvasCollection    * canvasCollection    = new CanvasCollection();
  CanvasConfiguration * canvasConfiguration = new CanvasConfiguration(CanvasConfiguration::Landscape,CanvasConfiguration::Linear);
  GraphConfiguration  ** graphConfigurations = new GraphConfiguration*[40];
  for (int iGraph=0;iGraph<40;iGraph++)
  {
  graphConfigurations[iGraph] = new GraphConfiguration(1,iGraph%10);
  graphConfigurations[iGraph]->markerSize = 0.5;
  graphConfigurations[iGraph]->plotOption = "E0";
 }

  TString inputPath  = "/Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/";
  TString outputPath = "/Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/";
   ///  /Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/PYTHIA_softOnHardOff_Singles_Wide_MB.root

  int nFiles            = 2;
  TString ** fileNames = new TString*[nFiles];
  TString fileName;
  fileNames[1] = new TString("PYTHIA_softOnHardOff_Singles_Narrow_MB.root");
  fileNames[0] = new TString("PYTHIA_softOnHardOff_Singles_Wide_MB.root");
  //fileNames[2] = new TString("hardOffPYTHIA_Singles_MB.root");

//  fileNames[1] = new TString("softOffPYTHIA_Singles_MB.root");
//  fileNames[2] = new TString("softOffPYTHIA_Singles_MB.root");
  TFile ** inputFiles  = new TFile* [nFiles];
  for (int iFile=0;iFile<nFiles;iFile++)
  {
  fileName = inputPath+*fileNames[iFile];
  inputFiles[iFile] = new TFile(fileName,"OLD");
  if (inputFiles[iFile])
    {
    cout << "<I> ParticlePlots() Successfully opened: " << fileName << endl;
    }
  else
    {
    cout << "<E> ParticlePlots() Unable to open: " << fileName << endl;
    return 1;
    }
  }

  int nModels = nFiles;
  TString ** modelNames = new TString*[nModels];
  modelNames[1] = new TString("PYTHIA-MinBiasEta1");
  modelNames[0] = new TString("PYTHIA-MinBiasEta6");
  TString ** modelTitles = new TString*[nModels];
  modelTitles[1] = new TString("PYTHIA MB |#eta|<1");
  modelTitles[0] = new TString("PYTHIA MB |#eta|<6");
  TString ** modelPrefixNames = new TString*[nModels];
  modelPrefixNames[1] = new TString("Narrow_MB_");
  modelPrefixNames[0] = new TString("Wide_MB_");

  int nParticleTypes = 12;
  TString ** particleNames = new TString*[nParticleTypes];
  particleNames[0] = new TString("HC");
  particleNames[1] = new TString("HP");
  particleNames[2] = new TString("HM");
  particleNames[3] = new TString("PiC");
  particleNames[4] = new TString("PiP");
  particleNames[5] = new TString("PiM");
  particleNames[6] = new TString("KC");
  particleNames[7] = new TString("KP");
  particleNames[8] = new TString("KM");
  particleNames[9] = new TString("PC");
  particleNames[10] = new TString("PP");
  particleNames[11] = new TString("PM");
  TString ** particleTitles = new TString*[nParticleTypes];
  particleTitles[0] = new TString("h^{#pm}");
  particleTitles[1] = new TString("h^{+}");
  particleTitles[2] = new TString("h^{-}");
  particleTitles[3] = new TString("#pi^{#pm}");
  particleTitles[4] = new TString("#pi^{+}");
  particleTitles[5] = new TString("#pi^{-}");
  particleTitles[6] = new TString("K^{#pm}");
  particleTitles[7] = new TString("K^{+}");
  particleTitles[8] = new TString("K^{-}");
  particleTitles[9] = new TString("P^{#pm}");
  particleTitles[10] = new TString("P^{+}");
  particleTitles[11] = new TString("P^{-}");

  AnalysisConfiguration * ac = new AnalysisConfiguration("","","");
  ac->fillY  = false;
  ac->fill3D = false;

  ParticlePlotter ** plotters = new ParticlePlotter*[nModels];

  int nSets = nModels*nParticleTypes;
  ParticleHistos  ** particleHistos = new ParticleHistos *[nSets];
  int iSet = 0;
  cout << "<I> ParticlePlots()  Loading histograms" << endl;
  TString prefixName;

  for (int iModel=0;iModel<nModels;iModel++)
    {
    cout << "iModel:" << iModel << " Name:" <<  *modelTitles[iModel]  << endl;
    for (int iParticleType=0;iParticleType<nParticleTypes;iParticleType++)
      {
      cout << "iParticleType:" << iParticleType << " Name:" <<  *particleTitles[iParticleType]  << endl;

      prefixName = *modelPrefixNames[iModel];
      prefixName += *particleNames[iParticleType];
      particleHistos[iSet] = new ParticleHistos(inputFiles[iModel],prefixName,ac,MessageLogger::Debug);
      plotters[iModel] = new ParticlePlotter();

      iSet++;
      }
    TString canvasNameBase = *modelNames[iModel];
    canvasNameBase += "_";
    plotters[iModel]->makePlots(canvasNameBase,12, particleHistos+(iModel*nParticleTypes),particleNames,particleTitles, canvasConfiguration,graphConfigurations);
    plotters[iModel]->printAllCanvas(outputPath);
   }

  Plotter ** compPlotters = new Plotter*[20];
  TString canvasName;
  TH1 ** histograms = new TH1*[20];
  TString ** titles = new TString*[20];


  canvasConfiguration->setAxes(CanvasConfiguration::LogY);
  canvasName = "ModelComp_HC_MultDist";
  histograms[0] = particleHistos[0]->h_n1; titles[0] = modelTitles[0];
  histograms[1] = particleHistos[nParticleTypes]->h_n1; titles[1] = modelTitles[1];
  compPlotters[0] = new Plotter();
  compPlotters[0]->plot(nModels, canvasName,canvasConfiguration,graphConfigurations,
       "N_{h^{#pm}}", 0.0, 200.0,
       "Counts", 1.0E-6, 5.0E-1,
       histograms,titles,0.55, 0.75, 0.8, 0.9, 0.04);
  compPlotters[0]->printAllCanvas(outputPath);


  canvasConfiguration->setAxes(CanvasConfiguration::Linear);
  canvasName = "ModelComp_HC_dndeta";
  histograms[0] = particleHistos[0]->h_n1_eta; titles[0] = modelTitles[0];
  histograms[1] = particleHistos[nParticleTypes]->h_n1_eta; titles[1] = modelTitles[1];
  compPlotters[1] = new Plotter();
  compPlotters[1]->plot(nModels, canvasName,canvasConfiguration,graphConfigurations,
                    "#eta", -6.0, 6.0,
                    "dN/d#eta", 0.0, 2.0,
                    histograms,titles,0.55, 0.17, 0.85, 0.35, 0.04);
  compPlotters[1]->printAllCanvas(outputPath);

  canvasConfiguration->setAxes(CanvasConfiguration::LogY);
  canvasName = "ModelComp_HC_dnpTdpT";
  histograms[0] = particleHistos[0]->h_n1_ptXS; titles[0] = modelTitles[0];
  histograms[1] = particleHistos[nParticleTypes]->h_n1_ptXS; titles[1] = modelTitles[1];
  compPlotters[2] = new Plotter();
  compPlotters[2]->plot(nModels, canvasName,canvasConfiguration,graphConfigurations,
                    "p_{T}", 0.0, 40.0,
                    "dN/p_{T}dp_{T}", 1.0E-9, 3.0E1,
                    histograms,titles,0.55, 0.6, 0.85, 0.9, 0.055);
  compPlotters[2]->printAllCanvas(outputPath);

  return 0;
}

