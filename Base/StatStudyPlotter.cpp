#include "StatStudyPlotter.hpp"
ClassImp(StatStudyPlotter);

StatStudyPlotter::StatStudyPlotter()
:
Plotter()
{  }

StatStudyPlotter::~StatStudyPlotter()
{
}



void StatStudyPlotter::makePlots(const TString & canvasNameBase,
                                 StatStudyHistograms * statStudyHistograms,
                                 CanvasConfiguration * canvasConfiguration,
                                 GraphConfiguration  ** graphConfigurations)
{
TH1 ** histograms = new TH1* [10];
TString ** hLegends = new TString * [10];

canvasConfiguration->setAxes(CanvasConfiguration::Linear);
//canvasConfiguration->setAxes(CanvasConfiguration::LogY);

histograms[0] = statStudyHistograms->f1_1;     hLegends[0] = new TString("n_{1}^{(1)}");
histograms[1] = statStudyHistograms->f1_2;     hLegends[1] = new TString("n_{1}^{(2)}");
histograms[2] = statStudyHistograms->f1_1Eff;  hLegends[2] = new TString("n_{1}^{(1,Eff)}");
histograms[3] = statStudyHistograms->f1_2Eff;  hLegends[3] = new TString("n_{1}^{(2,Eff)}");
plot(4, canvasNameBase+"f1",canvasConfiguration,graphConfigurations,
     "n_{1}", 0.0, 0.0,
     "P(n_{1})", 0.0, 0.0,
     histograms,hLegends,0.7, 0.4, 0.9, 0.8, 0.055);

histograms[0] = statStudyHistograms->f1_1Ratio;  hLegends[0] = new TString("f_{1}^{(1,Eff)}/f_{1}^{(1)}");
histograms[1] = statStudyHistograms->f1_2Ratio;  hLegends[1] = new TString("f_{1}^{(2,Eff)}/f_{1}^{(2)}");
plot(2, canvasNameBase+"f1Ratio",canvasConfiguration,graphConfigurations,
     "f_{1}^{Eff}/f_{1}", 0.84,0.88,
     "P(f_{1}^{Eff}/f_{1})", 0.0, 0.0,
     histograms,hLegends,0.16, 0.4, 0.32, 0.8, 0.055);

histograms[0] = statStudyHistograms->f2_11;     hLegends[0] = new TString("f_{2}^{(11)}");
histograms[1] = statStudyHistograms->f2_12;     hLegends[1] = new TString("f_{2}^{(12)}");
histograms[2] = statStudyHistograms->f2_22;     hLegends[2] = new TString("f_{2}^{(22)}");
histograms[3] = statStudyHistograms->f2_11Eff;  hLegends[3] = new TString("f_{2}^{(11,Eff)}");
histograms[4] = statStudyHistograms->f2_12Eff;  hLegends[4] = new TString("f_{2}^{(12,Eff)}");
histograms[5] = statStudyHistograms->f2_22Eff;  hLegends[5] = new TString("f_{2}^{(22,Eff)}");
plot(6, canvasNameBase+"f2",canvasConfiguration,graphConfigurations,
     "n_{2}",  0.0, 0.0,
     "P(n_{2})", 0.0, 0.0,
     histograms,hLegends,0.7, 0.4, 0.9, 0.8, 0.055);


histograms[0] = statStudyHistograms->f2_11Ratio;  hLegends[0] = new TString("f_{2}^{(11,Eff)}/f_{2}^{(11)}");
histograms[1] = statStudyHistograms->f2_12Ratio;  hLegends[1] = new TString("f_{2}^{(12,Eff)}/f_{2}^{(12)}");
histograms[2] = statStudyHistograms->f2_22Ratio;  hLegends[2] = new TString("f_{2}^{(22,Eff)}/f_{2}^{(22)}");
plot(3, canvasNameBase+"f2Ratio",canvasConfiguration,graphConfigurations,
     "f_{2}^{Eff}/f_{2}", 0.70,0.80,
     "P(f_{2}^{Eff}/f_{2})", 0.0, 0.0,
     histograms,hLegends,0.70, 0.4, 0.85, 0.9, 0.055);

histograms[0] = statStudyHistograms->r2_11;     hLegends[0] = new TString("R_{2}^{(11)}");
histograms[1] = statStudyHistograms->r2_12;     hLegends[1] = new TString("R_{2}^{(12)}");
histograms[2] = statStudyHistograms->r2_22;     hLegends[2] = new TString("R_{2}^{(22)}");
histograms[3] = statStudyHistograms->r2_11Eff;  hLegends[3] = new TString("R_{2}^{(11,Eff)}");
histograms[4] = statStudyHistograms->r2_12Eff;  hLegends[4] = new TString("R_{2}^{(12,Eff)}");
histograms[5] = statStudyHistograms->r2_22Eff;  hLegends[5] = new TString("R_{2}^{(22,Eff)}");
plot(6, canvasNameBase+"R2",canvasConfiguration,graphConfigurations,
     "R_{2}",  0.0, 0.003,
     "P(R_{2})", 0.0, 0.0,
     histograms,hLegends,0.16, 0.4, 0.32, 0.8, 0.055);

histograms[0] = statStudyHistograms->r2_11Ratio;  hLegends[0] = new TString("R_{2}^{(11,Eff)}/R_{2}^{(11)}");
histograms[1] = statStudyHistograms->r2_12Ratio;  hLegends[1] = new TString("R_{2}^{(12,Eff)}/R_{2}^{(12)}");
histograms[2] = statStudyHistograms->r2_22Ratio;  hLegends[2] = new TString("R_{2}^{(22,Eff)}/R_{2}^{(22)}");
plot(3, canvasNameBase+"R2Ratio",canvasConfiguration,graphConfigurations,
     "R_{2}^{Eff}/R_{2}", 0.8,1.2,
     "P(R_{2}^{Eff}/R_{2})", 0.0, 0.0,
     histograms,hLegends,0.70, 0.4, 0.85, 0.9, 0.055);
//  statStudyHistograms->r2_11Ratio->Fit("gaus");
//  statStudyHistograms->r2_12Ratio->Fit("gaus");
//  statStudyHistograms->r2_22Ratio->Fit("gaus");


//statStudyHistograms->nudyn_12->Fit("gaus");
//statStudyHistograms->nudyn_12Eff->Fit("gaus");
histograms[0] = statStudyHistograms->nudyn_12;     hLegends[0] = new TString("#nu_{dyn}");
histograms[1] = statStudyHistograms->nudyn_12Eff;  hLegends[1] = new TString("#nu_{dyn}^{Eff}");
plot(2, canvasNameBase+"NuDyn",canvasConfiguration,graphConfigurations,
     "#nu_{dyn}", -0.01,0.01,
     "P(nu_{dyn})", 0.0, 0.0,
     histograms,hLegends,0.70, 0.4, 0.85, 0.9, 0.055);

//statStudyHistograms->nudyn_12Ratio->Fit("gaus");
histograms[0] = statStudyHistograms->nudyn_12Ratio;     hLegends[0] = new TString("#nu_{dyn}^{Eff}/#nu_{dyn}");
plot(1, canvasNameBase+"NuDynRatio",canvasConfiguration,graphConfigurations,
     "#nu_{dyn}^{Eff}/#nu_{dyn}", 0.8,1.2,
     "P(#nu_{dyn}^{Eff}/#nu_{dyn})", 0.0, 0.0,
     histograms,hLegends,0.70, 0.4, 0.85, 0.9, 0.055);

histograms[0] = statStudyHistograms->f1_1Ratio_vsSub;     hLegends[0] = new TString("f_{1}^{Eff}/f_{1}");
plot(1, canvasNameBase+"f1_1Ratio_vsSub",canvasConfiguration,graphConfigurations,
     "sub-sample", 0.0, 1005.0,
     "f_{1}^{Eff}/f_{1}", 0.86, 0.88,
     histograms,hLegends,0.75, 0.8, 0.9, 0.9, 0.055);

histograms[0] = statStudyHistograms->r2_11_vsSub;     hLegends[0] = new TString("R_{2}^{11}");
plot(1, canvasNameBase+"r2_vsSub",canvasConfiguration,graphConfigurations,
    "sub-sample", 900.0, 1005.0,
     "R_{2}",  0.0, 0.01,
     histograms,hLegends,0.75, 0.8, 0.9, 0.9, 0.055);

histograms[0] = statStudyHistograms->r2_11Ratio_vsSub;     hLegends[0] = new TString("R_{2}^{11,Eff}/R_{2}^{11}");
plot(1, canvasNameBase+"R2ratio_vsSub",canvasConfiguration,graphConfigurations,
    "sub-sample", 900.0, 1005.0,
     "R_{2}^{Eff}/R_{2}", 0.8, 1.2,
     histograms,hLegends,0.70, 0.8, 0.85, 0.9, 0.055);

histograms[0] = statStudyHistograms->nudyn_12Ratio_vsSub;     hLegends[0] = new TString("#nu_{\rm Dyn}^{12,Eff}/#nu_{\rm Dyn}^{12}");
plot(1, canvasNameBase+"nudyn_12Ratio_vsSub",canvasConfiguration,graphConfigurations,
    "sub-sample", 900.0, 1005.0,
     "#nu_{\rm Dyn}^{12,Eff}/#nu_{\rm Dyn}^{12}", 0.5, 1.5,
     histograms,hLegends,0.70, 0.8, 0.85, 0.9, 0.055);

canvasConfiguration->setAxes(CanvasConfiguration::LogX);

histograms[0] = statStudyHistograms->ef1_1_vsSub;     hLegends[0] = new TString("#delta(f_{1}^{1})/f_{1}^{1}");
histograms[1] = statStudyHistograms->ef1_2_vsSub;     hLegends[1] = new TString("#delta(f_{1}^{2})/f_{1}^{2}");
histograms[2] = statStudyHistograms->ef1_1Eff_vsSub;  hLegends[2] = new TString("#delta(f_{1}^{Eff,1})/f_{1}^{Eff,1}");
histograms[3] = statStudyHistograms->ef1_2Eff_vsSub;  hLegends[3] = new TString("#delta(f_{1}^{Eff,2})/f_{1}^{Eff,2}");
plot(4, canvasNameBase+"ef1_vsSub",canvasConfiguration,graphConfigurations,
     "sub-sample", 0.0, 1005.0,
     "#delta(f_{1})/f_{1}", 0.0, 0.003,
     histograms,hLegends,0.70, 0.65, 0.85, 0.9, 0.05);

histograms[0] = statStudyHistograms->er2_11_vsSub;     hLegends[0] = new TString("#delta(R_{2}^{11})/R_{2}^{11}");
histograms[1] = statStudyHistograms->er2_12_vsSub;     hLegends[1] = new TString("#delta(R_{2}^{12})/R_{2}^{12}");
histograms[2] = statStudyHistograms->er2_22_vsSub;     hLegends[2] = new TString("#delta(R_{2}^{22})/R_{2}^{22}");
histograms[3] = statStudyHistograms->er2_11Eff_vsSub;  hLegends[3] = new TString("#delta(R_{2}^{Eff,11})/R_{2}^{Eff,11}");
histograms[4] = statStudyHistograms->er2_12Eff_vsSub;  hLegends[4] = new TString("#delta(R_{2}^{Eff,12})/R_{2}^{Eff,12}");
histograms[5] = statStudyHistograms->er2_22Eff_vsSub;  hLegends[5] = new TString("#delta(R_{2}^{Eff,22})/R_{2}^{Eff,22})");
plot(6, canvasNameBase+"er2_vsSub",canvasConfiguration,graphConfigurations,
     "sub-sample", 0.0, 1005.0,
     "#delta(R_{2})/R_{2}", 0.0, 0.3,
     histograms,hLegends,0.70, 0.55, 0.85, 0.9, 0.045);

histograms[0] = statStudyHistograms->er2_11Ratio_vsSub;     hLegends[0] = new TString("#delta(R_{2}^{11,Eff}/R_{2}^{11})/R_{2}^{11,Eff}/R_{2}^{11}");
histograms[1] = statStudyHistograms->er2_12Ratio_vsSub;     hLegends[1] = new TString("#delta(R_{2}^{12,Eff}/R_{2}^{12})/R_{2}^{12,Eff}/R_{2}^{12}");
histograms[2] = statStudyHistograms->er2_22Ratio_vsSub;  hLegends[2] = new TString("#delta(R_{2}^{22,Eff}/R_{2}^{22})/R_{2}^{22,Eff}/R_{2}^{22}");
plot(3, canvasNameBase+"er2Ratio_vsSub",canvasConfiguration,graphConfigurations,
     "sub-sample", 0.0, 1005.0,
     "#delta(R_{2}^{Eff}/R_{2})/R_{2}^{Eff}/R_{2}", 0.0, 0.2,
     histograms,hLegends,0.50, 0.65, 0.8, 0.9, 0.05);
}
