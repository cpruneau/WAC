#ifndef WAC_ParticlePlotter
#define WAC_ParticlePlotter
#include "Plotter.hpp"
#include "ParticleHistos.hpp"

class ParticlePlotter : public Plotter
{
public:
  
  ParticlePlotter();
  virtual ~ParticlePlotter();
  void makePlots(const TString & canvasNameBase,
                 int nSets,
                 ParticleHistos ** particleHistos,
                 TString ** particleNames,
                 TString ** particleTitles,
                 CanvasConfiguration * canvasConfiguration,
                 GraphConfiguration  ** graphConfigurations);

  ClassDef(ParticlePlotter,0)
};

#endif /* WAC_ParticlePlotter */
