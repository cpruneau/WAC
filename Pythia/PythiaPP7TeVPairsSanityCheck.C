
#include "TString.h"
#include "TRandom.h"

int PythiaPP7TeVPairsSanityCheck()
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

  cout << "<I> GlobalPlots() Includes loaded." << endl;

  HistogramCollection * histogramCollection = new HistogramCollection("Collection",200);
  histogramCollection->setDefaultOptions(1);
  CanvasCollection    * canvasCollection    = new CanvasCollection();
  CanvasConfiguration * canvasConfigurationLinear = new CanvasConfiguration(CanvasConfiguration::Landscape,CanvasConfiguration::Linear);
  CanvasConfiguration * canvasConfigurationLogY   = new CanvasConfiguration(CanvasConfiguration::Landscape,CanvasConfiguration::LogY);
  CanvasConfiguration * canvasConfigurationLogZ   = new CanvasConfiguration(CanvasConfiguration::LandscapeWide,CanvasConfiguration::LogZ);
  GraphConfiguration  ** graphConfigurations = new GraphConfiguration*[40];
  for (int iGraph=0;iGraph<40;iGraph++)
  {
  graphConfigurations[iGraph] = new GraphConfiguration(2,0);
  //graphConfigurations[iGraph]->markerSize = 0.5;
  graphConfigurations[iGraph]->plotOption = "SURF3";
  }

  cout << "<I> GlobalPlots() Configurations setup." << endl;

  TString inputPath  = "/Users/claudeapruneau/Documents/GitHub/WAC-DATA/InputFiles/PYTHIA_pp_7TeV_softOnHardOff/Pairs/ByCentrality/Partials/";
  TString outputPath = "/Users/claudeapruneau/Documents/GitHub/WAC-Data/InputFiles/PYTHIA_pp_7TeV_softOnHardOff/Pairs/ByCentrality/Partials/";
  ///  /Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/PYTHIA_softOnHardOff_Singles_Wide_MB.root

  int nFiles = 11;
  TFile ** inputFiles  = new TFile*[nFiles];
//  inputFiles[0] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__0.root","OLD");
//  inputFiles[1] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__1.root","OLD");
//  inputFiles[2] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__2.root","OLD");
//  inputFiles[3] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__3.root","OLD");
//  inputFiles[4] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__4.root","OLD");
//  inputFiles[5] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__5.root","OLD");
//  inputFiles[6] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__6.root","OLD");
//  inputFiles[7] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__7.root","OLD");
//  inputFiles[8] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__8.root","OLD");
//  inputFiles[9] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__9.root","OLD");

  inputFiles[0] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__0.root","OLD");
  inputFiles[1] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__1.root","OLD");
  inputFiles[2] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__2.root","OLD");
  inputFiles[3] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__3.root","OLD");
  inputFiles[4] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__4.root","OLD");
  inputFiles[5] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__5.root","OLD");
  inputFiles[6] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__6.root","OLD");
  inputFiles[7] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__7.root","OLD");
  inputFiles[8] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__8.root","OLD");
  inputFiles[9] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__9.root","OLD");

  inputFiles[10] = new TFile("/Users/claudeapruneau/Documents/GitHub/WAC-Data/InputFiles/PYTHIA_pp_7TeV_softOnHardOff/Pairs/ByCentrality/Sums/PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__Sum_0_9.root","OLD");
  if (inputFiles[0])
    {
    cout << "<I> GlobalPlots() Successfully opened" << endl;
    }
  else
    {
    cout << "<E> GlobalPlots() Unable to open first file" << endl;
    return 1;
    }



  TString canvasNameBase = "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft_";

  TString ** histogramNames = new TString*[nFiles];
  histogramNames[0] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[6] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[7] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[8] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[9] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[10] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");

//  histogramNames[0] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[1] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[2] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[3] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[4] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[5] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[6] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[7] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[8] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[9] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  histogramNames[10] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
//  PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10


  TH2 ** histograms = new TH2*[nFiles];
  for (int iFile=0; iFile<nFiles; iFile++)
  {
  histograms[iFile] = (TH2*) inputFiles[iFile]->Get(*histogramNames[iFile]);
  if (histograms[iFile])
    {
    cout << "<I> GlobalPlots() Successfully loaded: " << *histogramNames[iFile] << endl;
    }
  else
    {
    cout << "<E> GlobalPlots() Unable to load: " << *histogramNames[iFile] << endl;
    return 1;
    }
  }



  if (true)
    {
    Plotter * plotter = new Plotter();
    for (int iFile=0; iFile<nFiles; iFile++)
    {
    plotter->plot(canvasNameBase+iFile,canvasConfigurationLinear,graphConfigurations[0],
                  "#Delta#eta", 1.0, -1.0,
                  "#Delta#phi", 1.0, -1.0,
                  "G_{2}", 1.0, -1.0,
                  histograms[iFile]);
    }
    plotter->printAllCanvas(outputPath);

    }


  return 0;
}
