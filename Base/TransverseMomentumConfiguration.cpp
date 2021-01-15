/**
 \class Property
 \ingroup WAC

 Utility class used to TransverseMomentumConfiguration
 */

#include "TransverseMomentumConfiguration.hpp"

ClassImp(TransverseMomentumConfiguration);


TransverseMomentumConfiguration::TransverseMomentumConfiguration(const TString & name,
                                       const TString & type,
                                       const TString & version)
:
TaskConfiguration(name,type,version),
nBins_mult(20),
min_mult(0.0),
max_mult(200.0),
nBins_cent(20),
min_cent(0.0),
max_cent(100.0),
ptCorrelatorVsMult(true),
ptCorrelatorVsCent(false),
totEvents(100000000)
{
}

TransverseMomentumConfiguration::TransverseMomentumConfiguration(const TransverseMomentumConfiguration & source)
:
TaskConfiguration( source ),
nBins_mult( source.nBins_mult ),
min_mult(   source.min_mult ),
max_mult(   source.max_mult ),
nBins_cent( source.nBins_cent ),
max_cent(   source.max_cent ),
min_cent(   source.min_cent ),
ptCorrelatorVsMult(source.ptCorrelatorVsMult),
ptCorrelatorVsCent(source.ptCorrelatorVsCent),
totEvents(source.totEvents),
maxOrder(source.maxOrder)
{
}

TransverseMomentumConfiguration::~TransverseMomentumConfiguration()
{

}

TransverseMomentumConfiguration & TransverseMomentumConfiguration::operator=(const TransverseMomentumConfiguration & source)
{
  if (this!=&source)
    {
    TaskConfiguration::operator=( source );
    nBins_mult  =  source.nBins_mult;
    min_mult  =  source.min_mult;
    max_mult  =  source.max_mult;
    nBins_cent  =  source.nBins_cent;
    min_cent  =  source.min_cent;
    max_cent  =  source.max_cent;
    ptCorrelatorVsMult = source.ptCorrelatorVsMult;
    ptCorrelatorVsCent = source.ptCorrelatorVsCent;
    totEvents = source.totEvents;
    maxOrder = source.maxOrder;
    }
  return *this;
}



////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void TransverseMomentumConfiguration::printConfiguration(ostream & os)
{
  printTaskConfiguration(os);
  os
  << "                  TransverseMomentumConfiguration Parameters: " << endl
  << " ------------------------------------------------------------------------------------------" << endl
  << "                        nBins_mult : " << nBins_mult << endl
  << "                          min_mult : " << min_mult << endl
  << "                          max_mult : " << max_mult << endl
  << "                        nBins_cent : " << nBins_cent << endl
  << "                          min_cent : " << min_cent << endl
  << "                          max_cent : " << max_cent << endl
  << "                ptCorrelatorVsMult : " << ptCorrelatorVsMult << endl
  << "                ptCorrelatorVsCent : " << ptCorrelatorVsCent << endl
  << "                         totEvents : " << totEvents << endl
  << "                          maxOrder : " << maxOrder << endl
  << "                          numTypes : " << numTypes << endl;




}
