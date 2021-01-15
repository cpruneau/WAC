//
//  BWModelEventGenerator.hpp
//  Created by Victor Gonzalez on 20190206 based in the work of Claude Pruneau and the ALICE AliPhysics SW framework
//  Copyright © 2016-2019 Claude Pruneau. All rights reserved.
//

#ifndef WAC_BWModelEventGenerator
#define WAC_BWModelEventGenerator
#include <iostream>
#include <TF1.h>
#include "Task.hpp"
#include "BWModelConfiguration.hpp"

class BWModelEventGenerator : public Task
{
public:


  BWModelEventGenerator(const TString & name,
                        BWModelConfiguration * configuration,
                        Event * event,
                        LogLevel selectedLevel);
  virtual ~BWModelEventGenerator();
  virtual void initialize();
  virtual void reset();
  void execute();

  static TF1 *GetExponentialSpectra_FN(const char *name, double mass, double temp, double xmin, double xmax);
  static double BGBlastWaveIntegrand(double *x, double *par);
  static Double_t BGBlastWavedNdPt_FN(double *x, double *par);
  static TF1 *GetBGBWSpectra_FN(const char *name, double mass, double temp, double beta_avg,
                                double n, double xmin, double xmax);

  BWModelConfiguration * generatorConfiguration;

  bool   useAllKinds;
  bool   useFlow;

  double pionPercentage;
  double kaonPercentage;
  double protonPercentage;

  TF1   *ptSpectraAllKinds;
  TF1   *azimuthalAngleAllKinds;
  TF1   *ptSpectraPions;
  TF1   *azimuthalAnglePions;
  TF1   *ptSpectraKaons;
  TF1   *azimuthalAngleKaons;
  TF1   *ptSpectraProtons;
  TF1   *azimuthalAngleProtons;

    ClassDef(BWModelEventGenerator,0)
};

#endif /* WAC_BWModelEventGenerator */

