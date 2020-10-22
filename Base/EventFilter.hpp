// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_EventFilter
#define WAC_EventFilter
#include "TString.h"
#include "Event.hpp"

class EventFilter
{
public:

  enum EventSelection   { MinBias, Centrality, Multiplicity, ImpactParameter, Other };


   EventFilter(EventSelection eventSelectionRequested,
               double minimum,
               double maximum);
  virtual ~EventFilter();
  virtual bool accept(Event & event);

  TString getName();
  TString getTitle();
//  TString getLongName();
//  TString getLongTitle();

  EventSelection eventSelectionRequested;
  double minimum;
  double maximum;

  ClassDef(EventFilter,0)
};

#endif /* WAC_EventFilter */
