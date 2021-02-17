//==============================================================================
//  calculate the thermodynamic quantities of hadron resonance gas
//
//  Copyright 2016 Chun Shen
//  Email: shen.201@asc.ohio-state.edu
//==============================================================================

#include<sys/time.h>
//#include<gsl/gsl_sf_bessel.h>

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>

#include "Hadron.hpp"
#include "HadronGas.hpp"
#include "HadronGasConfiguration.hpp"
#include "HadronGasHistograms.hpp"
#include "ChemicalPotential.hpp"
#include "TString.h"
#include "TSystem.h"
using namespace std;

//int main(int argc, const char * argv[])


int RunHadronGas()
{
  TString includesPath = "/Users/claudeapruneau/Documents/GitHub/WAC/Base/";
  gSystem->Load(includesPath+"CanvasCollection.hpp");
  gSystem->Load(includesPath+"CanvasConfiguration.hpp");
  gSystem->Load(includesPath+"TaskConfiguration.hpp");
  gSystem->Load(includesPath+"MessageLogger.hpp");
  gSystem->Load(includesPath+"CanvasConfiguration.hpp");
  gSystem->Load(includesPath+"TaskConfiguration.hpp");
  gSystem->Load(includesPath+"Task.hpp");
  gSystem->Load(includesPath+"TaskCollection.hpp");
  gSystem->Load(includesPath+"Plotter.hpp");
  gSystem->Load(includesPath+"GraphConfiguration.hpp");
  gSystem->Load(includesPath+"HistogramCollection.hpp");
  gSystem->Load(includesPath+"Histograms.hpp");
  gSystem->Load("libBase.dylib");

  includesPath = "/Users/claudeapruneau/Documents/GitHub/WAC/HadronGas/";
  gSystem->Load(includesPath+"Hadron.hpp");
  gSystem->Load(includesPath+"HadronGas.hpp");
  gSystem->Load(includesPath+"HadronGasConfiguration.hpp");
  gSystem->Load(includesPath+"HadronGasHistograms.hpp");
  gSystem->Load(includesPath+"HadronGasParticleHistograms.hpp");
  gSystem->Load(includesPath+"ChemicalPotential.hpp");
  gSystem->Load("libHadronGas.dylib");

  std::cout << "==================================================================================" << std::endl;

  TString pdgDataTable = getenv("WAC_SOURCE");

  pdgDataTable += "/EOS/pdg.dat";

  HadronGasConfiguration config("test","test","1");
  config.nT           = 10;
  config.minT         = 0.100;
  config.maxT         = 0.200;
  config.nTypes       = 319;
  config.nStableTypes = 23;
  HadronGas hadronGas(pdgDataTable.Data());
  hadronGas.printParticleBasicProperties(std::cout);

  double temperature;
  double muB = 0.0;
  double muS = 0.0;

  HadronGasHistograms * hgh = new HadronGasHistograms("HG_",&config,MessageLogger::Info);
  hgh->createHistograms();
  HadronGasParticleHistograms ** hgph = new HadronGasParticleHistograms *[config.nT];
  TString ** tempLabels = new TString*[config.nT];
  TString label;

  double dT = (config.maxT-config.minT)/double(config.nT);
  for (int iT=0;iT<config.nT;iT++)
  {
  temperature = config.minT + double(iT)*dT;
  label = "T = ";
  label += int(1000*temperature);
  tempLabels[iT] = new TString(label);

  hadronGas.calculateSystemProperties(temperature,muB,muS);
  hadronGas.printParticleYields(std::cout);
  hgh->fill(hadronGas);
  TString hName = "HGP_";
  hName += iT;
  hgph[iT] = new HadronGasParticleHistograms(hName,&config,MessageLogger::Info);
  hgph[iT] ->createHistograms();
  hgph[iT] ->fill(hadronGas);
  hgph[iT] ->calculateDerivedHistograms();
  }
  CanvasConfiguration * canvasConfiguration = new CanvasConfiguration(CanvasConfiguration::Landscape,CanvasConfiguration::Linear);
  CanvasConfiguration * canvasConfigurationLogY = new CanvasConfiguration(CanvasConfiguration::Landscape,CanvasConfiguration::LogY);
  CanvasConfiguration * canvasConfigurationLogZ = new CanvasConfiguration(CanvasConfiguration::Landscape,CanvasConfiguration::LogZ);
  GraphConfiguration  ** graphConfigurations  = new GraphConfiguration*[40];
  GraphConfiguration  *  graphConfiguration2D = new GraphConfiguration(2,1);
  for (int iGraph=0;iGraph<40;iGraph++)
    {
    graphConfigurations[iGraph] = new GraphConfiguration(1,iGraph%10);
    graphConfigurations[iGraph]->plotOption = "P0.L";
    }
  //TString inputPath  = "/Users/claudeapruneau/Documents/GitHub/run/PythiaStudies/";
  TString outputPath = "/Users/claudeapruneau/Documents/GitHub/WAC-DATA/OutputFiles/HGM/";
  graphConfigurations[0]->plotOption = "P1";
  Plotter * plotter = new Plotter();
  plotter->setDefaultOptions(1);
  plotter->plot("HGM_SysEnergyVsT", canvasConfigurationLogY, graphConfigurations[0],
                "T (GeV)",  config.minT, config.maxT,
                "#varepsilon (GeV/fm^{3})",0.01, 10.0,
                hgh->sysEnergyVsT,
                "Energy Density",0.22, 0.7, 0.4, 0.76, 0.055);

  plotter->plot("HGM_SysEntropyVsT", canvasConfigurationLogY, graphConfigurations[0],
                "T (GeV)",  config.minT, config.maxT,
                "s",0.1, 100.0,
                hgh->sysEntropyVsT,
                "Entropy Density",0.22, 0.7, 0.4, 0.76, 0.055);

  plotter->plot("HGM_PressureVsT", canvasConfigurationLogY, graphConfigurations[0],
                "T (GeV)",  config.minT, config.maxT,
                "P (-)",0.001, 1.0,
                hgh->sysPressureVsT,
                "Pressure",0.22, 0.7, 0.4, 0.76, 0.055);

  double minYield = 1.0E-7;
  double maxYield = 1.0;

  TString ** partLabel = new TString*[23];
  partLabel[0] = new TString("#pi^{0}");
  partLabel[1] = new TString("#pi^{+}");
  partLabel[2] = new TString("#pi^{-}");
  partLabel[3] = new TString("K^{+}");
  partLabel[4] = new TString("K^{-}");
  partLabel[5] = new TString("K^{0}");
  partLabel[6] = new TString("#bar K^{0}");
  partLabel[7] = new TString("p");
  partLabel[8] = new TString("#bar p");
  partLabel[9] = new TString("n");
  partLabel[10] = new TString("#bar n");
  partLabel[11] = new TString("#Lambda");
  partLabel[12] = new TString("#bar#Lambda");
  partLabel[13] = new TString("#Sigma^{+}");
  partLabel[14] = new TString("#bar#Sigma^{+}");
  partLabel[15] = new TString("#Sigma^{-}");
  partLabel[16] = new TString("#bar#Sigma^{-}");
  partLabel[17] = new TString("#Xi^{0}");
  partLabel[18] = new TString("#bar#Xi^{0}");
  partLabel[19] = new TString("#Xi^{-}");
  partLabel[20] = new TString("#bar#Xi^{-}");
  partLabel[21] = new TString("#Omega");
  partLabel[22] = new TString("#bar#Omega");

  int nX = hgph[0]->h_stableThermalYields->GetNbinsX();

  for (int iX=1; iX<=23; iX++)
  {
  cout << " iX:" << iX << "    Label:" << partLabel[iX-1]->Data() << endl;
  hgph[0]->h_stableThermalYields->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[0]->h_stableDecayYields->GetXaxis()  ->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[0]->h_stableDecayToThermalRatio->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[2]->h_stableThermalYields->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[2]->h_stableDecayYields->GetXaxis()  ->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[2]->h_stableDecayToThermalRatio->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[4]->h_stableThermalYields->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[4]->h_stableDecayYields->GetXaxis()  ->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[4]->h_stableDecayToThermalRatio->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());

  hgph[4]->h_stableDecayCorrelatedPairsYields->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[4]->h_stableDecayCorrelatedPairsYields->GetYaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());

  hgph[4]->h_stableDecayCorrelatedPairsNorm->GetXaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());
  hgph[4]->h_stableDecayCorrelatedPairsNorm->GetYaxis()->SetBinLabel(iX,partLabel[iX-1]->Data());

 }

  plotter->plot("HGM_AllYields", canvasConfigurationLogY, graphConfigurations[0],
                "Species",  0.0, 400.0,
                "Yield",minYield,maxYield,
                hgph[2]->h_allYields,
                "Thermal Yields",0.22, 0.7, 0.4, 0.76, 0.055);
  plotter->plot("HGM_StableThermalYields", canvasConfigurationLogY, graphConfigurations[0],
                "Species",  0.0, 23.0,
                "Yield",minYield,maxYield,
                hgph[2]->h_stableThermalYields,
                "Thermal Yields",0.2, 0.20, 0.45, 0.45, 0.055);
  plotter->plot("HGM_StableDecayYields", canvasConfigurationLogY, graphConfigurations[0],
                "Species",  0.0, 23.0,
                "Yield",minYield,maxYield,
                hgph[2]->h_stableDecayYields,
                "Yields w/ Decays",0.2, 0.20, 0.45, 0.45, 0.055);

  int nGraphs = 5;
  TString ** legends = new TString*[nGraphs];
  TH1 ** histograms = new TH1*[nGraphs];

  histograms[0] = hgph[0]->h_stableThermalYields; legends[0] = tempLabels[0];
  histograms[1] = hgph[2]->h_stableThermalYields; legends[1] = tempLabels[2];
  histograms[2] = hgph[4]->h_stableThermalYields; legends[2] = tempLabels[4];
  histograms[3] = hgph[6]->h_stableThermalYields; legends[3] = tempLabels[6];
  histograms[4] = hgph[8]->h_stableThermalYields; legends[4] = tempLabels[8];

  plotter->plot(nGraphs, "HGM_StableThermalYieldsVsT",canvasConfigurationLogY,graphConfigurations,
                        "Species", 0.0, 23.0,
                        "Thermal Yield", minYield,maxYield,
                        histograms,legends,0.2, 0.20, 0.45, 0.45, 0.04);

  histograms[0] = hgph[0]->h_stableDecayYields; legends[0] = tempLabels[0];
  histograms[1] = hgph[2]->h_stableDecayYields; legends[1] = tempLabels[2];
  histograms[2] = hgph[4]->h_stableDecayYields; legends[2] = tempLabels[4];
  histograms[3] = hgph[6]->h_stableDecayYields; legends[3] = tempLabels[6];
  histograms[4] = hgph[8]->h_stableDecayYields; legends[4] = tempLabels[8];

  plotter->plot(nGraphs, "HGM_StableDecayYieldsVsT",canvasConfigurationLogY,graphConfigurations,
                        "Species", 0.0, 23.0,
                        "Thermal + Decay Yield", minYield,maxYield,
                        histograms,legends,0.2, 0.20, 0.45, 0.45, 0.04);


  histograms[0] = hgph[0]->h_stableDecayToThermalRatio; legends[0] = tempLabels[0];
  histograms[1] = hgph[2]->h_stableDecayToThermalRatio; legends[1] = tempLabels[2];
  histograms[2] = hgph[4]->h_stableDecayToThermalRatio; legends[2] = tempLabels[4];
  histograms[3] = hgph[6]->h_stableDecayToThermalRatio; legends[3] = tempLabels[6];
  histograms[4] = hgph[8]->h_stableDecayToThermalRatio; legends[4] = tempLabels[8];
  plotter->plot(nGraphs, "HGM_StableDecayToThermalVsT",canvasConfigurationLogY,graphConfigurations,
                        "Species", 0.0, 23.0,
                        "Decay/Thermal", 0.9,20.0,
                        histograms,legends,0.7, 0.60, 0.85, 0.85, 0.04);



  histograms[0] = (TH1*) hgph[4]->h_stableThermalYields->Clone(); legends[0] = new TString(*tempLabels[4]+" MeV; Thermal Only");
  histograms[1] = (TH1*) hgph[4]->h_stableDecayYields->Clone();   legends[1] = new TString(*tempLabels[4]+" MeV; Thermal + Decays");

  plotter->plot(2, "HGM_StableDecayVsThermalYields",canvasConfigurationLogY,graphConfigurations,
                        "Species", 0.0, 23.0,
                        "Yields", minYield,maxYield,
                        histograms,legends,0.2, 0.20, 0.45, 0.45, 0.04);
  canvasConfigurationLogZ->theta =  35.0;
  canvasConfigurationLogZ->phi   = -125.0;
  graphConfiguration2D->plotOption = "LEGO2";
  hgph[4]->h_stableDecayCorrelatedPairsYields->GetXaxis()->SetNdivisions(23);
  hgph[4]->h_stableDecayCorrelatedPairsYields->GetYaxis()->SetNdivisions(23);
  plotter->plot("HGM_StableDecayCorrelatedPairYields",canvasConfigurationLogZ,graphConfiguration2D,
               "Species 1", 0.0, 23.0,
               "Species 2", 0.0, 23.0,
               "Yield", 1E-7, 1.0,
               hgph[4]->h_stableDecayCorrelatedPairsYields);


  hgph[4]->h_stableDecayPairsYields->GetXaxis()->SetNdivisions(23);
  hgph[4]->h_stableDecayPairsYields->GetYaxis()->SetNdivisions(23);
  plotter->plot("HGM_StableDecayPairYields",canvasConfigurationLogZ,graphConfiguration2D,
               "Species 1", 0.0, 23.0,
               "Species 2", 0.0, 23.0,
               "Yield", 1E-10, 10.0,
               hgph[4]->h_stableDecayPairsYields);




  hgph[4]->h_stableDecayPairsYields->GetXaxis()->SetNdivisions(23);
  hgph[4]->h_stableDecayPairsYields->GetYaxis()->SetNdivisions(23);
  plotter->plot("HGM_StableDecayCorrelatedPairNorm",canvasConfigurationLogZ,graphConfiguration2D,
               "Trigger Species", 0.0, 23.0,
               "Associate Species", 0.0, 23.0,
               "Conditional Yield", 1E-10, 1.0,
               hgph[4]->h_stableDecayCorrelatedPairsNorm);

  plotter->printAllCanvas(outputPath);
 return 0;
}

