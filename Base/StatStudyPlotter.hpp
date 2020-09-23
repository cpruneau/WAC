#ifndef WAC_StatStudyPlotter
#define WAC_StatStudyPlotter
#include "Plotter.hpp"
#include "StatStudyHistograms.hpp"

class StatStudyPlotter : public Plotter
{
public:
  
  StatStudyPlotter();
  virtual ~StatStudyPlotter();
  void makePlots(const TString & canvasNameBase,
                 StatStudyHistograms * statStudyHistogramsUsed,
                 CanvasConfiguration * canvasConfigurationUsed,
                 GraphConfiguration  ** graphConfigurationsUsed);

  void makeComparisonPlots(const TString & canvasNameBase,
                           CanvasConfiguration * canvasConfiguration,
                           int nGraphs,
                           StatStudyHistograms ** statStudyHistograms,
                           TString ** hLegends,
                           GraphConfiguration  ** graphConfigurations);

  ClassDef(StatStudyPlotter,0)
};

#endif /* WAC_StatStudyPlotter */
