//
//  EventFilter.cpp
//  MyMC
//
//  Created by Claude Pruneau on 12/12/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include "EventFilter.hpp"

ClassImp(EventFilter);

//////////////////////////////////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////////////////////////////////
EventFilter::EventFilter(EventSelection selection,
                         double min,
                         double max)
:
eventSelectionRequested(selection),
minimum(min),
maximum(max)
{
  // no ops
}

//////////////////////////////////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////////////////////////////////
EventFilter::~EventFilter()
{
  // no ops
}

//////////////////////////////////////////////////////////////////////////////////////////
// accept/reject the given particle based on filter parameter
//////////////////////////////////////////////////////////////////////////////////////////
bool EventFilter::accept(Event & event)
{
  if (eventSelectionRequested==MinBias) return true;
  double condition = 0.0;
  switch (eventSelectionRequested)
    {
      case MinBias: break;
      case Centrality:      condition = event.centrality; break;
      case Multiplicity:    condition = event.multiplicity; break;
      case ImpactParameter: condition = event.impactParameter; break;
      case Other:           condition = event.other; break;
    }

  return (condition >= minimum) && (condition < maximum);
}

//enum EventSelection   { MinBias, Centrality, Multiplicity, ImpactParameter, Other };

TString EventFilter::getName()
{
  TString name;
  switch (eventSelectionRequested)
  {
    case MinBias:         name = "MB"; break;
    case Centrality:      name = Form("CGeq%gLt%g",minimum,maximum); break;
    case Multiplicity:    name = Form("MGeq%gLt%g",minimum,maximum); break;
    case ImpactParameter: name = Form("bGeq%gLt%g",minimum,maximum); break;
    case Other:           name = Form("XGeq%gLt%g",minimum,maximum); break;
  }
  return name;
}



TString EventFilter::getTitle()
{
  TString name;
  switch (eventSelectionRequested)
  {
    case MinBias:         name = "Min-Bias"; break;
    case Centrality:      name = Form("%g < C < %g",minimum,maximum); break;
    case Multiplicity:    name = Form("%g < M < %g",minimum,maximum); break;
    case ImpactParameter: name = Form("%g < b < %g",minimum,maximum); break;
    case Other:           name = Form("%g < X < %g",minimum,maximum); break;
  }
  return name;
}

//TString EventFilter::getLongName();
//TString EventFilter::getLongTitle();
