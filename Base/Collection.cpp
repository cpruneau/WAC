// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class EventLoop
 \ingroup WAC

 Class used to define a generic object Collection
 */
#include "Collection.hpp"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"

ClassImp(Collection<TF1>);
ClassImp(Collection<TH1>);
ClassImp(Collection<TCanvas>);

