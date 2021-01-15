
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

int ParticleSubsampleAnalysis()
{
  cout << "ParticleSubsampleAnalysis() Starting" << endl;
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
  gSystem->Load(includesPath+"SubSampleStatCalculator.hpp");
  gSystem->Load("libBase.dylib");

  MessageLogger::LogLevel  debugLevel = MessageLogger::Info;

  TString inputPathName  = getenv("WAC_INPUT_PATH");
  TString outputPathName = getenv("WAC_OUTPUT_PATH");
  //TString inputFileNameBase = "PYTHIA_softOnHardOff_Singles_LUMP_Wide_MB_";
  //TString outputFileName = "PYTHIA_softOnHardOff_Singles_LUMP_Wide_MB_Sum0_99";
  TString inputFileBaseName = "PYTHIA_softOnHardOff_Singles_LUMP_Narrow_MB_";
  TString outputFileName    = "PYTHIA_softOnHardOff_Singles_LUMP_Narrow_MB_Sum0_99";
  int nFiles = 100;
  TString ** fileNames = new TString* [nFiles];
  cout << "      input Path Name: " << inputPathName << endl;
  cout << " input File Base Name: " << inputFileBaseName << endl;
  cout << "     output Path Name: " << outputPathName << endl;
  cout << "     output File Name: " << outputFileName << endl;
  for (int iFile=0; iFile<nFiles; iFile++)
  {
  TString  fName = inputFileBaseName;
  fName += iFile;
  fileNames[iFile] = new TString(fName);
  cout << "   file Name [" << iFile << "]: " << *fileNames[iFile] << endl;
  }
  SubSampleStatCalculator * calculator = new SubSampleStatCalculator(inputPathName,nFiles,fileNames,outputPathName,outputFileName,debugLevel);
  calculator->execute();

  return 0;

}
