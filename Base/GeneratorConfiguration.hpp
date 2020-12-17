// ////////////////////////////////////////////////////////////////////////////////////////
// GeneratorConfiguration.hpp
//  Created by Claude Pruneau on 4/7/17.
//  Copyright © 2017 Claude Pruneau. All rights reserved.
// ////////////////////////////////////////////////////////////////////////////////////////
#ifndef WAC_GeneratorConfiguration
#define WAC_GeneratorConfiguration
#include "TaskConfiguration.hpp"

class GeneratorConfiguration : public TaskConfiguration
{
public:

  GeneratorConfiguration(const TString & name, const TString & type, const TString & version);
  GeneratorConfiguration(const GeneratorConfiguration & source);
  virtual ~GeneratorConfiguration(){} 
  GeneratorConfiguration & operator=(const GeneratorConfiguration & source);


  virtual void printConfiguration(ostream & os);

  // ////////////////////////////////////////////
  //  Data Members
  // ////////////////////////////////////////////

  TString dataInputPath;
   TString dataInputFileName; // can use templates...
   int dataFileMinIndex;
  int dataFileMaxIndex;


  double   impactParameterMinimum;
  double   impactParameterMaximum;
  double   ptMinimum;
  double   ptMaximum;
  double   yMinimum;
  double   yMaximum;
  double   etaMinimum;
  double   etaMaximum;
  bool     chargedOnly;
  TString  eventsSource;

   ClassDef(GeneratorConfiguration,0)
};


#endif /*WAC_GeneratorConfiguration*/
