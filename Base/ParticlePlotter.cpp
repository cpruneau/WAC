#include "ParticlePlotter.hpp"
ClassImp(ParticlePlotter);

ParticlePlotter::ParticlePlotter()
:
Plotter()
{  }

ParticlePlotter::~ParticlePlotter()
{
}



void ParticlePlotter::makePlots(const TString & canvasNameBase,
                                int nSets,
                                ParticleHistos ** particleHistos,
                                TString ** particleNames,
                                TString ** particleTitles,
                                CanvasConfiguration * canvasConfiguration,
                                GraphConfiguration  ** graphConfigurations)
{
  TH1 ** histograms = new TH1* [100];
  TString histNames = "";
  for (int iSet=0;iSet<nSets;iSet++)
  {
  histNames += *particleNames[iSet];
  }

  //canvasConfiguration->setAxes(CanvasConfiguration::LogY);
  //canvasConfiguration->setAxes(CanvasConfiguration::Linear);

  for (int iSet=0;iSet<nSets;iSet++)
  {
  histograms[iSet] = particleHistos[iSet]->h_n1;
  }
  canvasConfiguration->setAxes(CanvasConfiguration::LogY);
  plot(nSets, canvasNameBase+histNames+"_MultDist",canvasConfiguration,graphConfigurations,
       "N", 0.0, 200.0,
       "Counts", 1.0E-6, 5.0E-1,
       histograms,particleTitles,0.75, 0.5, 0.85, 0.9, 0.055);

  for (int iSet=0;iSet<nSets;iSet++)
  {
  histograms[iSet] = particleHistos[iSet]->h_n1_eta;
  }
  canvasConfiguration->setAxes(CanvasConfiguration::Linear);
  plot(nSets, canvasNameBase+histNames+"_dndeta",canvasConfiguration,graphConfigurations,
       "#eta", 2.0, -2.0,
       "dN/d#eta", 0.0, 0.8,
       histograms,particleTitles,0.75, 0.2, 0.85, 0.6, 0.055);

  for (int iSet=0;iSet<nSets;iSet++)
  {
  histograms[iSet] = particleHistos[iSet]->h_n1_pt;
  }
  canvasConfiguration->setAxes(CanvasConfiguration::LogY);
  plot(nSets, canvasNameBase+histNames+"_dndpt",canvasConfiguration,graphConfigurations,
       "p_{T}", 0.0, 20.0,
       "dN/dp_{T}", 1.0E-7, 5.0E0,
       histograms,particleTitles,0.7, 0.4, 0.85, 0.9, 0.055);

  for (int iSet=0;iSet<nSets;iSet++)
  {
  histograms[iSet] = particleHistos[iSet]->h_n1_ptXS;
  }
  canvasConfiguration->setAxes(CanvasConfiguration::LogY);
  plot(nSets, canvasNameBase+histNames+"_dnptdpt",canvasConfiguration,graphConfigurations,
       "p_{T}", 0.0, 20.0,
       "dN/p_{T}dp_{T}", 1.0E-9, 3.0E1,
       histograms,particleTitles,0.7, 0.4, 0.85, 0.9, 0.055);

}

//void ParticlePlotter::makeComparisonPlots(const TString & canvasNameBase,
//                                          CanvasConfiguration * canvasConfiguration,
//                                          int nGraphs,
//                                          ParticleHistos ** particleHistos,
//                                          TString ** histLabels,
//                                          GraphConfiguration  ** graphConfigurations)
//{
//  TH1 ** histograms = new TH1* [100];
//  TString ** histTitles = new TString * [100];
//
//  TString histName;
//  TString histTitleBase;
//  TString histTitle;
//
//  canvasConfiguration->setAxes(CanvasConfiguration::Linear);
//  //canvasConfiguration->setAxes(CanvasConfiguration::LogY);
//int index1234 = nuDynHistos[0]->index4(0,1,2,3);
//  histTitleBase = "f_{4}";
//  for (int i1=0; i1<nGraphs; i1++)
//    {
//    histTitle = histTitleBase;
//    histTitle += *histLabels[i1];
//    //cout << "Inserting: " << histTitle  << endl;
//    histograms[i1] = nuDynHistos[i1]->h_f4_vsMult[index1234];
//    if (!histograms[i1])
//        {
//        cout << "<F> Histograms not loaded!!!! ABORT!!!" << endl;
//        return;
//        }
//    histTitles[i1] = new TString(histTitle);
//
//    }
//
//   plot(nGraphs, canvasNameBase+"f4_vsMult",canvasConfiguration,graphConfigurations,
//        "Mult", 20.0, 400.0,
//        "f_{4}", 0.0, -1.2,
//        histograms,histTitles,0.25, 0.6, 0.45, 0.85, 0.055);
//
//  histTitleBase = "F_{4}";
//   for (int i1=0; i1<nGraphs; i1++)
//     {
//     histTitle = histTitleBase;
//     histTitle += *histLabels[i1];
//     cout << "Inserting: " << histTitle  << endl;
//     histograms[i1] = nuDynDerivedHistos[i1]->h_F4_vsMult[index1234]; histTitles[i1] = new TString(histTitle);
//     if (!histograms[i1])
//       {
//       cout << "<F> Histograms not loaded!!!! ABORT!!!" << endl;
//       return;
//       }
//     }
//    plot(nGraphs, canvasNameBase+"F4_vsMult",canvasConfiguration,graphConfigurations,
//         "Mult", 20.0, 400.0,
//         "F_{4}", -50000, 50000,
//         histograms,histTitles,0.25, 0.6, 0.45, 0.85, 0.055);
//
//  histTitleBase = "R_{4}";
//  for (int i1=0; i1<nGraphs; i1++)
//    {
//    histTitle = histTitleBase;
//    histTitle += *histLabels[i1];
//    cout << "Inserting: " << histTitle  << endl;
//    histograms[i1] = nuDynDerivedHistos[i1]->h_R4_vsMult[index1234]; histTitles[i1] = new TString(histTitle);
//    if (!histograms[i1])
//      {
//      cout << "<F> Histograms not loaded!!!! ABORT!!!" << endl;
//      return;
//      }
//    }
//   plot(nGraphs, canvasNameBase+"R4_vsMult",canvasConfiguration,graphConfigurations,
//        "Mult", 20.0, 400.0,
//        "R_{4}", -0.3, 0.3,
//        histograms,histTitles,0.65, 0.6, 0.85, 0.9, 0.055);
//
//
//
//
//}
