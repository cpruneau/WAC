#ifndef WAC_HeavyIonConfiguration
#define WAC_HeavyIonConfiguration
#include "NuDynConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class HeavyIonConfiguration : public NuDynConfiguration
{
public:
  
  HeavyIonConfiguration(const TString & name,
                        const TString & type,
                        const TString & version);
  virtual ~HeavyIonConfiguration(){}
  void printConfiguration(ostream & os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  bool hardBoost; //corresponds to a hard boost scenario. If not hard boost the soft boost

  int nCollisionsMax; //holds the max number of binary collisions of an event

  double param_a; double param_b; //parameters in deciding the boost of particles (beta = a * r^b)

  int maxOrder; //order of correlation functions

  int totEvents; //total number of events

  int nThreads; //number of threads for implicit parallel processing

  // event wise parameters



  ///////////////////////////////////////////////////
  // Parameters for Reading events from files
  //Requirements: 
  //One tree per file
  // treess must have branches with the above data
  //All trees have the same name
  //All files have the same name except for a final numeric subscript to distinguish them 
  //No .root extension on the files name
  TString branchName_eventNo;
  TString branchName_mult;
  TString branchName_px;
  TString branchName_py;
  TString branchName_pz;
  TString branchName_ist;
  TString branchName_pdg;
  TString branchName_pE;
  TString treeFile; // file containing all the trees where the pythia events are stored to be read. The full name of the root files WITHOUT the .root extension
  TString eventTreeName; // name of the tree containing the events
  int numFiles;  //total number of root files. Files are assumed to have the same name except for a numeric designation at the end BEFORE the .root extension
  

ClassDef(HeavyIonConfiguration,0)
};

#endif /* WAC_AnalysisConfiguration */
