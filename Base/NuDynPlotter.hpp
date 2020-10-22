#ifndef WAC_NuDynPlotter
#define WAC_NuDynPlotter
#include "Plotter.hpp"
#include "NuDynHistos.hpp"
#include "NuDynDerivedHistos.hpp"

class NuDynPlotter : public Plotter
{
public:
  
  NuDynPlotter();
  virtual ~NuDynPlotter();
  void makePlots(const TString & canvasNameBase,
                 NuDynHistos * nuDynHistosUsed,
                 NuDynDerivedHistos * nuDynDerivedHistosUsed,
                 CanvasConfiguration * canvasConfigurationUsed,
                 GraphConfiguration  ** graphConfigurationsUsed);

  void makeComparisonPlots(const TString & canvasNameBase,
                           CanvasConfiguration * canvasConfiguration,
                           int nGraphs,
                           NuDynHistos ** nuDynHistos,
                           NuDynDerivedHistos ** nuDynDerivedHistos,
                           TString ** histLabels,
                           GraphConfiguration  ** graphConfigurations);

  void makeNudynComparisonPlots(const TString & canvasNameBase,
                                CanvasConfiguration * canvasConfiguration,
                                int nGraphs,
                                NuDynHistos ** nuDynHistos,
                                NuDynDerivedHistos ** nuDynDerivedHistos,
                                TString ** histLabels,
                                GraphConfiguration  ** graphConfigurations);

  ClassDef(NuDynPlotter,0)
};

#endif /* WAC_NuDynPlotter */
