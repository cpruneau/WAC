// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining ParticlePairHistos
 */

#ifndef WAC_ParticlePairFilter
#define WAC_ParticlePairFilter
#include "Particle.hpp"

//////////////////////////////////////////////////////////////////////////////////////////
// Two Particle Filter
//////////////////////////////////////////////////////////////////////////////////////////
class ParticlePairFilter
{
public:
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////////////////////////////////
  ParticlePairFilter(double minDpt,
                     double maxDpt,
                     double minDphi,
                     double maxDphi,
                     double minDeta,
                     double maxDeta);
  virtual ~ParticlePairFilter();
  bool accept(Particle & particle1, Particle &  particle2);

  //////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////////////////////////////////
  double min_dpt;
  double max_dpt;
  double min_dphi;
  double max_dphi;
  double min_deta;
  double max_deta;

  ClassDef(ParticlePairFilter,0)
};

#endif /* WAC_ParticlePairFilter */
