
#include "TString.h"
#include "TRandom.h"

int TwoParticlePlots()
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

  cout << "<I> TwoParticlePlots() Includes loaded." << endl;

  HistogramCollection * histogramCollection = new HistogramCollection("Collection",200);
  histogramCollection->setDefaultOptions(1);
  CanvasCollection    * canvasCollection    = new CanvasCollection();
  CanvasConfiguration * canvasConfigurationLinear = new CanvasConfiguration(CanvasConfiguration::LandscapeWide,CanvasConfiguration::Linear);
  CanvasConfiguration * canvasConfigurationLogY   = new CanvasConfiguration(CanvasConfiguration::LandscapeWide,CanvasConfiguration::LogY);
  CanvasConfiguration * canvasConfigurationLogZ   = new CanvasConfiguration(CanvasConfiguration::LandscapeWide,CanvasConfiguration::LogZ);
  GraphConfiguration  ** graphConfigurations = new GraphConfiguration*[40];
  for (int iGraph=0;iGraph<40;iGraph++)
  {
  graphConfigurations[iGraph] = new GraphConfiguration(1,iGraph%10);
  graphConfigurations[iGraph]->markerSize = 0.5;
  graphConfigurations[iGraph]->plotOption = "E0";
  }

  cout << "<I> TwoParticlePlots() Configurations setup." << endl;

  TString inputPath  = "/Users/claudeapruneau/Documents/GitHub/WAC-DATA/InputFiles/PYTHIA_pp_7TeV_softOnHardOff/Pairs/VsXSect/Sums/";
  TString outputPath = "/Users/claudeapruneau/Documents/GitHub/WAC-DATA/InputFiles/PYTHIA_pp_7TeV_softOnHardOff/Pairs/VsXSect/Sums/";
  ///  /Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/PYTHIA_softOnHardOff_Singles_Wide_MB.root

  TFile ** inputFiles  = new TFile*[12];
  TString ** histogramNames   = new TString*[20];
  TString ** histogramZTitles = new TString*[20];
  TH1 ** histograms = new TH1*[20];
  Plotter * plotter = new Plotter();
  TString canvasNameBase;

  inputFiles[0] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_100_MB__Sum_0_9.root","OLD");
  inputFiles[1] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_5_V0MnGeq66Lt1000__Sum_0_9.root","OLD");
  inputFiles[2] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_5_10_V0MnGeq53Lt66__Sum_0_9.root","OLD");
  inputFiles[3] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_10_20_V0MnGeq38Lt53__Sum_0_9.root","OLD");
  inputFiles[4] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_20_30_V0MnGeq28Lt38__Sum_0_9.root","OLD");
  inputFiles[5] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_30_40_V0MnGeq20Lt28__Sum_0_9.root","OLD");
  inputFiles[6] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_40_50_V0MnGeq16Lt20__Sum_0_9.root","OLD");
  inputFiles[7] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_50_60_V0MnGeq12Lt16__Sum_0_9.root","OLD");
  inputFiles[8] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_60_70_V0MnGeq10Lt12__Sum_0_9.root","OLD");
  inputFiles[9] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_70_80_V0MnGeq7Lt10__Sum_0_9.root","OLD");
  inputFiles[10] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_80_90_V0MnGeq2Lt7__Sum_0_9.root","OLD");
  inputFiles[11] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_90_100_V0MnGeq0Lt2__Sum_0_9.root","OLD");
  if (inputFiles[0])
    {
    cout << "<I> TwoParticlePlots() Successfully opened: PYTHIA_pp_7TeV_softOnHardOff_Pairs_Narrow_HPHM_0_100_MB__Sum_0_9.root" << endl;
    }
  else
    {
    cout << "<E> TwoParticlePlots() Unable to open first file" << endl;
    return 1;
    }

  TString collisionTitle = "pp #sqrt{s_{NN}} = 7 TeV; PYTHIA-8";
  TString multRange;

  histogramZTitles[0] = new TString("R_{2}^{CI}");
  histogramZTitles[1] = new TString("R_{2}^{CD}");
  histogramZTitles[2] = new TString("P_{2}^{CI}");
  histogramZTitles[3] = new TString("P_{2}^{CD}");
  histogramZTitles[4] = new TString("G_{2}^{CI}");
  histogramZTitles[5] = new TString("G_{2}^{CD}");
  canvasNameBase = "PYTHIA_pp_7TeV_softOnHardOff_Pairs_";

  int iFile = 0;
  // ========================================================
  // 0 - 100 %
  // ========================================================
  multRange = "0 - 100%";
   histogramNames[0] = new TString("Narrow_HPHM_0_100_MB_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_0_100_MB_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_0_100_MB_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_0_100_MB_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_0_100_MB_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_0_100_MB_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;
  // ========================================================
  // 0 - 5 %
  // ========================================================
  multRange = "0 - 5%";
  histogramNames[0] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_0_5_V0MnGeq66Lt1000_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 5 - 10 %
  // ========================================================
  multRange = "5 - 10%";
  histogramNames[0] = new TString("Narrow_HPHM_5_10_V0MnGeq53Lt66_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_5_10_V0MnGeq53Lt66_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_5_10_V0MnGeq53Lt66_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_5_10_V0MnGeq53Lt66_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_5_10_V0MnGeq53Lt66_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_5_10_V0MnGeq53Lt66_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 10 - 20 %
  // ========================================================
  multRange = "10 - 20%";
  histogramNames[0] = new TString("Narrow_HPHM_10_20_V0MnGeq38Lt53_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_10_20_V0MnGeq38Lt53_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_10_20_V0MnGeq38Lt53_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_10_20_V0MnGeq38Lt53_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_10_20_V0MnGeq38Lt53_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_10_20_V0MnGeq38Lt53_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 20 - 30 %
  // ========================================================
  multRange = "20 - 30%";
  histogramNames[0] = new TString("Narrow_HPHM_20_30_V0MnGeq28Lt38_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_20_30_V0MnGeq28Lt38_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_20_30_V0MnGeq28Lt38_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_20_30_V0MnGeq28Lt38_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_20_30_V0MnGeq28Lt38_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_20_30_V0MnGeq28Lt38_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 30 - 40 %
  // ========================================================
  multRange = "30 - 40%";
  histogramNames[0] = new TString("Narrow_HPHM_30_40_V0MnGeq20Lt28_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_30_40_V0MnGeq20Lt28_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_30_40_V0MnGeq20Lt28_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_30_40_V0MnGeq20Lt28_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_30_40_V0MnGeq20Lt28_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_30_40_V0MnGeq20Lt28_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 40 - 50 %
  // ========================================================
  multRange = "40 - 50%";
  histogramNames[0] = new TString("Narrow_HPHM_40_50_V0MnGeq16Lt20_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_40_50_V0MnGeq16Lt20_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_40_50_V0MnGeq16Lt20_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_40_50_V0MnGeq16Lt20_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_40_50_V0MnGeq16Lt20_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_40_50_V0MnGeq16Lt20_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 50 - 60 %
  // ========================================================
  multRange = "50 - 60%";
  histogramNames[0] = new TString("Narrow_HPHM_50_60_V0MnGeq12Lt16_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_50_60_V0MnGeq12Lt16_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_50_60_V0MnGeq12Lt16_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_50_60_V0MnGeq12Lt16_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_50_60_V0MnGeq12Lt16_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_50_60_V0MnGeq12Lt16_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 60 - 70 %
  // ========================================================
  multRange = "60 - 70%";
  histogramNames[0] = new TString("Narrow_HPHM_60_70_V0MnGeq10Lt12_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_60_70_V0MnGeq10Lt12_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_60_70_V0MnGeq10Lt12_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_60_70_V0MnGeq10Lt12_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_60_70_V0MnGeq10Lt12_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_60_70_V0MnGeq10Lt12_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 70 - 80 %
  // ========================================================
  multRange = "70 - 80%";
  histogramNames[0] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_70_80_V0MnGeq7Lt10_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 80 - 90 %
  // ========================================================
  multRange = "80 - 90%";
  histogramNames[0] = new TString("Narrow_HPHM_80_90_V0MnGeq2Lt7_HPHMCI_R2_DetaDphi_shft");
  histogramNames[1] = new TString("Narrow_HPHM_80_90_V0MnGeq2Lt7_HPHMCD_R2_DetaDphi_shft");
  histogramNames[2] = new TString("Narrow_HPHM_80_90_V0MnGeq2Lt7_HPHMCI_P2_DetaDphi_shft");
  histogramNames[3] = new TString("Narrow_HPHM_80_90_V0MnGeq2Lt7_HPHMCD_P2_DetaDphi_shft");
  histogramNames[4] = new TString("Narrow_HPHM_80_90_V0MnGeq2Lt7_HPHMCI_G2_DetaDphi_shft");
  histogramNames[5] = new TString("Narrow_HPHM_80_90_V0MnGeq2Lt7_HPHMCD_G2_DetaDphi_shft");
  for (int k=0; k<6; k++)
    {
    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[k])
      {
      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
                                "#Delta#eta", 1.0, -1.0,
                                "#Delta#varphi", 1.0, -1.0,
                                *histogramZTitles[k], 1.0, -1.0,
                                histograms[k],"SURF3");
    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
    }
  iFile++;

  // ========================================================
  // 90 - 100 %
  // ========================================================
//  multRange = "90 - 100%";
//  histogramNames[0] = new TString("Narrow_HPHM_90_100_V0MnGeq0Lt2_HPHMCI_R2_DetaDphi_shft");
//  histogramNames[1] = new TString("Narrow_HPHM_90_100_V0MnGeq0Lt2_HPHMCD_R2_DetaDphi_shft");
//  histogramNames[2] = new TString("Narrow_HPHM_90_100_V0MnGeq0Lt2_HPHMCI_P2_DetaDphi_shft");
//  histogramNames[3] = new TString("Narrow_HPHM_90_100_V0MnGeq0Lt2_HPHMCD_P2_DetaDphi_shft");
//  histogramNames[4] = new TString("Narrow_HPHM_90_100_V0MnGeq0Lt2_HPHMCI_G2_DetaDphi_shft");
//  histogramNames[5] = new TString("Narrow_HPHM_90_100_V0MnGeq0Lt2_HPHMCD_G2_DetaDphi_shft");
//
//  Narrow_HPHM_90_100_V0MnGeq0Lt2
//  iFile = 10;
//  for (int k=0; k<6; k++)
//    {
//    histograms[k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
//    if (histograms[k])
//      {
//      cout << "<I> TwoParticlePlots() Successfully loaded: " << *histogramNames[k] << endl;
//      }
//    else
//      {
//      cout << "<E> TwoParticlePlots() Unable to load: " << *histogramNames[k] << endl;
//      return 1;
//      }
//    plotter->plot(canvasNameBase+*histogramNames[k],canvasConfigurationLinear,graphConfigurations[0],
//                                "#Delta#eta", 1.0, -1.0,
//                                "#Delta#varphi", 1.0, -1.0,
//                                *histogramZTitles[k], 1.0, -1.0,
//                                histograms[k],"SURF3");
//    plotter->createLabel(-0.93,0.95,1,0,0.055,collisionTitle);
//    plotter->createLabel(-0.93,0.75,1,0,0.055,multRange);
//    }
//  iFile++;

  plotter->printAllCanvas(outputPath);

  return 0;
}
