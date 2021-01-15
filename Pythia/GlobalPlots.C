
#include "TString.h"
#include "TRandom.h"

int GlobalPlots()
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
  graphConfigurations[iGraph] = new GraphConfiguration(1,iGraph%10);
  graphConfigurations[iGraph]->markerSize = 0.5;
  graphConfigurations[iGraph]->plotOption = "E0";
  }

  cout << "<I> GlobalPlots() Configurations setup." << endl;

  TString inputPath  = "/Users/claudeapruneau/Documents/GitHub/WAC/InputFiles/";
  TString outputPath = "/Users/claudeapruneau/Documents/GitHub/WAC/OutputFiles/";
  ///  /Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/PYTHIA_softOnHardOff_Singles_Wide_MB.root

  TFile ** inputFiles  = new TFile*[12];
  inputFiles[0] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_MB_0.root","OLD");
  inputFiles[1] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_90_100_V0MnGeq0Lt2_0.root","OLD");
  inputFiles[2] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_80_90_V0MnGeq2Lt7_0.root","OLD");
  inputFiles[3] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_70_80_V0MnGeq7Lt10_0.root","OLD");
  inputFiles[4] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_60_70_V0MnGeq10Lt12_0.root","OLD");
  inputFiles[5] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_50_60_V0MnGeq12Lt16_0.root","OLD");
  inputFiles[6] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_40_50_V0MnGeq16Lt20_0.root","OLD");
  inputFiles[7] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_30_40_V0MnGeq20Lt28_0.root","OLD");
  inputFiles[8] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_20_30_V0MnGeq28Lt38_0.root","OLD");
  inputFiles[9] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_10_20_V0MnGeq38Lt53_0.root","OLD");
  inputFiles[10] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_5_10_V0MnGeq53Lt66_0.root","OLD");
  inputFiles[11] = new TFile(inputPath + "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_0_5_V0MnGeq66Lt1000_0.root","OLD");
  if (inputFiles[0])
    {
    cout << "<I> GlobalPlots() Successfully opened: PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_MB_0.root" << endl;
    }
  else
    {
    cout << "<E> GlobalPlots() Unable to open first file" << endl;
    return 1;
    }



  TString canvasNameBase = "PYTHIA_pp_7TeV_softOnHardOff_GlobalObservables_V0vsTPC_MB_";

  TString ** histogramNames = new TString*[20];
  histogramNames[0] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_n");
  histogramNames[1] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_e");
  histogramNames[2] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_q");
  histogramNames[3] = new TString("V0vsTPC_MB__V0MPtGeq50Lt1000000_n");
  histogramNames[4] = new TString("V0vsTPC_MB__V0MPtGeq50Lt1000000_e");
  histogramNames[5] = new TString("V0vsTPC_MB__V0MPtGeq50Lt1000000_q");

  histogramNames[6] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_n_HCPtGeq200Lt100000EtaGeqM800Lt800_e");
  histogramNames[7] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_n_HCPtGeq200Lt100000EtaGeqM800Lt800_q");
  histogramNames[8] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_n_V0MPtGeq50Lt1000000_n");
  histogramNames[9] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_q_V0MPtGeq50Lt1000000_q");
  histogramNames[10] = new TString("V0vsTPC_MB__HCPtGeq200Lt100000EtaGeqM800Lt800_n_V0MPtGeq50Lt1000000_q");

  TH1 *** histograms = new TH1**[12];
  for (int iFile=0; iFile<1; iFile++)
  {
  histograms[iFile] = new TH1*[11];
  for (int k=0; k<11; k++)
    {
    histograms[iFile][k] = (TH1*) inputFiles[iFile]->Get(*histogramNames[k]);
    if (histograms[iFile][k])
      {
      cout << "<I> GlobalPlots() Successfully loaded: " << *histogramNames[k] << endl;
      }
    else
      {
      cout << "<E> GlobalPlots() Unable to load: " << *histogramNames[k] << endl;
      return 1;
      }
    }
  }



  if (true)
    {


    Plotter * plotter = new Plotter();
    plotter->plot(canvasNameBase+"TPC_n",canvasConfigurationLogY,graphConfigurations[0],
                  "TPC n", 0.0, 100.0,
                  "Yield", 1.0E-8,1.0,
                  histograms[0][0], "TPC Charged Multiplicity",
                  0.5, 0.7, 0.8, 0.9, 0.05);

    plotter->plot(canvasNameBase+"V0M_n",canvasConfigurationLogY,graphConfigurations[0],
                  "V0M n", 0.0, 400.0,
                  "Yield", 1.0E-8,1.0,
                  histograms[0][3], "V0M Charged Multiplicity",
                  0.5, 0.7, 0.8, 0.9, 0.05);

    TH1 ** histos = new TH1*[12];
    TString ** legends = new TString*[12];
    histos[0] = histograms[0][0];
    histos[1] = histograms[0][3];
    legends[0] = new TString("TPC Charged Multiplicity");
    legends[1] = new TString("V0M Charged Multiplicity");
    plotter->plot(2,canvasNameBase+"V0M_n_TPC_n",canvasConfigurationLogY,graphConfigurations,
                  "n", 0.0, 400.0,
                  "Yield", 1.0E-8,1.0,
                  histos,legends,
                  0.5, 0.7, 0.8, 0.9, 0.05);

    histos[0] = histograms[0][1];
    histos[1] = histograms[0][4];
    legends[0] = new TString("TPC Charged Energy");
    legends[1] = new TString("V0M Charged Energy");
    plotter->plot(2,canvasNameBase+"V0M_e_TPC_e",canvasConfigurationLogY,graphConfigurations,
                  "e", 0.0, 1000.0,
                  "Yield", 1.0E-8,1.0,
                  histos,legends,
                  0.5, 0.7, 0.8, 0.9, 0.05);

    histos[0] = histograms[0][2];
    histos[1] = histograms[0][5];
    legends[0] = new TString("TPC Net Charge");
    legends[1] = new TString("V0M Net Charge");
    plotter->plot(2,canvasNameBase+"V0M_q_TPC_q",canvasConfigurationLogY,graphConfigurations,
                  "Q", -25.0, 25.0,
                  "Yield", 1.0E-8,1.0,
                  histos,legends,
                  0.3, 0.2, 0.4, 0.3, 0.05);

    plotter->plot(canvasNameBase+"TPC_n_Vs_V0M_n",canvasConfigurationLogZ,graphConfigurations[0],
                  "TPC n", 0.0, 200.0,
                  "V0M n", 0.0, 200.0,
                  "Yield", 1.0E-7, 2.0,
                  histograms[0][8], "COLZ");

    plotter->plot(canvasNameBase+"TPC_q_Vs_V0M_q",canvasConfigurationLogZ,graphConfigurations[0],
                  "TPC q", -25.0, 25.0,
                  "V0M q",-25.0, 25.0,
                  "Yield", 1.0E-7, 2.0,
                  histograms[0][9], "COLZ");

    histos[0] = histograms[0][3];
    histos[1] = (TH1*) inputFiles[1]->Get("V0vsTPC_90_100_V0MnGeq0Lt2__V0MPtGeq50Lt1000000_n");
    histos[2] = (TH1*) inputFiles[2]->Get("V0vsTPC_80_90_V0MnGeq2Lt7__V0MPtGeq50Lt1000000_n");
    histos[3] = (TH1*) inputFiles[3]->Get("V0vsTPC_70_80_V0MnGeq7Lt10__V0MPtGeq50Lt1000000_n");
    histos[4] = (TH1*) inputFiles[4]->Get("V0vsTPC_60_70_V0MnGeq10Lt12__V0MPtGeq50Lt1000000_n");
    histos[5] = (TH1*) inputFiles[5]->Get("V0vsTPC_50_60_V0MnGeq12Lt16__V0MPtGeq50Lt1000000_n");
    histos[6] = (TH1*) inputFiles[6]->Get("V0vsTPC_40_50_V0MnGeq16Lt20__V0MPtGeq50Lt1000000_n");
    histos[7] = (TH1*) inputFiles[7]->Get("V0vsTPC_30_40_V0MnGeq20Lt28__V0MPtGeq50Lt1000000_n");
    histos[8] = (TH1*) inputFiles[8]->Get("V0vsTPC_20_30_V0MnGeq28Lt38__V0MPtGeq50Lt1000000_n");
    histos[9] = (TH1*) inputFiles[9]->Get("V0vsTPC_10_20_V0MnGeq38Lt53__V0MPtGeq50Lt1000000_n");
    histos[10] = (TH1*) inputFiles[10]->Get("V0vsTPC_5_10_V0MnGeq53Lt66__V0MPtGeq50Lt1000000_n");
    histos[11] = (TH1*) inputFiles[11]->Get("V0vsTPC_0_5_V0MnGeq66Lt1000__V0MPtGeq50Lt1000000_n");

    for (int kk=1;kk<12;kk++)
      {
      histos[kk]->Scale(0.1);
      }
    histos[10]->Scale(0.5);
    histos[11]->Scale(0.5);

    legends[0] = new TString("V0M Net Charge");
    legends[1] = new TString("90 - 100%");
    legends[2] = new TString("80 - 90%");
    legends[3] = new TString("70 - 80%");
    legends[4] = new TString("60 - 70%");
    legends[5] = new TString("50 - 60%");
    legends[6] = new TString("40 - 50%");
    legends[7] = new TString("30 - 40%");
    legends[8] = new TString("20 - 30%");
    legends[9] = new TString("10 - 20%");
    legends[10] = new TString("5 - 10%");
    legends[11] = new TString("0 - 5%");
    plotter->plot(12,canvasNameBase+"V0MvsXset_n",canvasConfigurationLogY,graphConfigurations,
                  "n", 0.0, 200.0,
                  "Yield", 1.0E-8,1.0,
                  histos,legends,
                  0.2, 0.2, 0.4, 0.5, 0.03);
    plotter->printAllCanvas(outputPath);

    }

  if (false)
    {
    // compute x-sect slices...

    double fractions[12];
    fractions[0] = 0.05;
    fractions[1] = 0.1;
    fractions[2] = 0.2;
    fractions[3] = 0.3;
    fractions[4] = 0.4;
    fractions[5] = 0.5;
    fractions[6] = 0.6;
    fractions[7] = 0.7;
    fractions[8] = 0.8;
    fractions[9] = 0.9;
    fractions[10] = 1.00;
    fractions[11] = 2.00;
    double sum = 0.0;
    double count = 0;
    double frac;
    double content;

    TH1* histo = histograms[0][3];

    int n = histo->GetNbinsX();
    for (int i=1; i<=n; i++)
      {
      content = histo->GetBinContent(i);
      sum += content;
      }
    int iFrac = 0;
    for (int i=n-1; i>=0; i--)
      {
      double edge = histo->GetBinLowEdge(i);
      double width = histo->GetBinWidth(i);
      content = histo->GetBinContent(i);
      count += content;
      frac = count/sum;
      //cout << content << "    "  << count << "     " << frac << "    " << edge << endl;
      if (frac>fractions[iFrac])
        {
        double excessFrac = (frac - fractions[iFrac]);
        double adjustedEdge = edge + excessFrac*width;
        cout << frac  << "   " << edge <<  "   " << adjustedEdge << endl;
        iFrac++;
        }
      }
    }
  return 0;
}
