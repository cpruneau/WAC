// ////////////////////////////////////////////////////////////////////////////////////////
//  GeneralizedGaussianConfiguration.hpp
//  Created by Victor Gonzalez on 02/10/19 on an original idea from Claude Pruneau.
//  Copyright © 2017-2019 Claude Pruneau. All rights reserved.
// ////////////////////////////////////////////////////////////////////////////////////////
#ifndef WAC_GeneralizedGaussianConfiguration
#define WAC_GeneralizedGaussianConfiguration
#include <ostream>
#include <TMath.h>
#include "TaskConfiguration.hpp"

// Emulate a generalized gaussian shape on the pair correlation

using std::ostream;
using std::endl;

class GeneralizedGaussianConfiguration : public TaskConfiguration
{
public:

  GeneralizedGaussianConfiguration(const TString & name,
                                   const TString & version);
  GeneralizedGaussianConfiguration(const GeneralizedGaussianConfiguration & source);
  GeneralizedGaussianConfiguration & operator=(const GeneralizedGaussianConfiguration & source);

  virtual ~GeneralizedGaussianConfiguration(){}

  virtual void printConfiguration(ostream & os);

  // ////////////////////////////////////////////
  //  Data Members
  // ////////////////////////////////////////////
  double      amplitude;
  double      gammaeta;
  double      gammaphi;
  double      omegaeta;
  double      omegaphi;

   ClassDef(GeneralizedGaussianConfiguration,0)
};


#endif /* WAC_GeneralizedGaussianConfiguration */
