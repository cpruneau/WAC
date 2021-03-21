#ifndef WAC_PTHistos
#define WAC_PTHistos
#include "Histograms.hpp"
#include "TransverseMomentumConfiguration.hpp"
#include "TMath.h"
#include <cstdlib>
#include <cmath>
#include <chrono>

class PTHistos : public Histograms
{
public:

  PTHistos(const TString & collectionName,
    TaskConfiguration * analysisConfiguration,
    LogLevel  debugLevel,
    int ord);
  PTHistos(TFile * inputFile,
    const TString & collectionName,
    TaskConfiguration * analysisConfiguration,
    LogLevel  debugLevel,
    int ord);
  virtual ~PTHistos();

  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void fillDerivedHistos(double *** transverseMomentumMoments, double ** yields, double * mults, double * cents, double * numParticles);
  virtual void saveHistograms(TFile * outputFile, bool saveAll=false);
  virtual void createHistogramRec(TString * baseName, TString * baseTitle, int depth, int partIndex);
  virtual void loadHistogramRec(TString * baseName, int depth, int partIndex, TFile * inputFile);
  virtual void fillEventHistos(double mult, double cent, double weight);
  virtual void fillTransverseMomentumHistos(double transverseMomentum, int filter, double mult, double cent, double weight);
  virtual void fillNormalizedPTValues(int depth, int partIndex, double product, double relerrpt, int bin, TH1 **OldHistos, TProfile **pTHistos, TH1 **newHistos);
  //virtual void fillNormalizedPTValues( int depth, int partIndex, double product, TH1 *** values, int* reorder, int*  nBin);
  virtual void calculateCumulants(TH1 ** Shistos, TH1 **CHistos, int nBins, double min, double max);
  virtual void calcRecSum(TH1 **CHistos, int *iHisto, int iBin, int *Subset, int len, int *set, int lenSet, double productC, double *used, int &curInd, int productS, double &sum, double &sumError, double relErr, int depth);
  virtual void calculatePTDeviationMoments(double *** transverseMomentumMoments, int bin, int iEvent, int nParticles, TProfile ** pTHisto);

  ////////////////////////////////////////////////////////////////////////////
  //Helper Functions
  ////////////////////////////////////////////////////////////////////////////

  virtual int convert(int * num, int len);
  virtual int* convert(int num, int & len);
  virtual void convertToBinary(int num, char*str, int len );
  virtual int* getSubset(char* subset, int * set, int len, int& lenSub);
  virtual int* getComplementarySubset(char* subset, int * set, int len, int& lenSub);
  virtual int getSubsetNumber(int * subset, int lenSub, int * mainset, int lenSet);



  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  // S is the pT deviation moments
  // s are the normalized moments
  // s* are the moments normalized by average pT's
  // C is the cumulants
  // c are the normalized cumulants
  // c* are the cumulants normalizd by average pT's

  //maximum order of correlation functions (2,3,4 ...) 
  //order is also equal to number of particle filters in PTCorrelator
  int maxOrder;

  //temporary variable used to help the code
  int histoIndex;

  //calculates the total number of combinations of 1, 2, 3, ... maxOrder filters
  //Ex: if there are 3 filters, the combinations are 1, 2, 3, 11, 12, ... 333
  //size of   TProfile ** h_c array
  int size; 

  //keep track of the orders of the correlation functions in the histograms
  //Ex: the cumulant C_(1,2,1) is order 3
  //Useful because the histograms are not in "logical" order, but in "recursive" order
  int * orders;

  //stores the yields of each possible combination of particles
  //i.e. n1, n2, n3, n4, n1*(n1 - 1), n1 * n2, so on NOT in that order 
  double ** yields;

  //store the average transverse momenta of each type of particle
  // first index is filter index(1-4)
  TProfile ** pT;
  TProfile ** pT_vsMult;
  TProfile ** pT_vsCent;

  // Min bias all included NOT IN ORDER (They are in "recursive order")
  // in the order S, s, s*, C, c, c*
  TH1 * h_events;

  TProfile ** S; 
	TProfile ** S_vsMult;
	TProfile ** S_vsCent;

  //first index is function index(S, s, s*), second index is histo index
  TH1 *** hS;
  TH1 *** hC;



  // vs Mult measured in fiducial NOT IN ORDER (They are in "recursive order")
  // in the order S, s, s*, C, c, 
  TH1 * h_events_vsMult;

  //first index is function index(S, s, s*), second index is histo index
  TH1 *** hS_vsMult;
  TH1 *** hC_vsMult;


  // vs Centrality NOT IN ORDER (They are in "recursive order")
  // in the order S, s, s*, C, c, c*
  TH1 * h_events_vsCent;
  //first index is function index(S, s, s*), second index is histo index
  TH1 *** hS_vsCent;
  TH1 *** hC_vsCent;


  //number of pairs, triples ... NOT IN ORDER (They are in "recursive order")
  TProfile ** h_counts;
  TProfile ** h_counts_vsMult;
  TProfile ** h_counts_vsCent;


  //number of functions we are counting 
  // here it is 3 (S, s, s*) (C, c, c*)
  int numFunc;

  //keep track of which index in "recursive order" corresponds to which index in normal order
  int * reorder;

  //opposite of above, plug in the normal order, and spits out the recursive order.
  int * reorder2;

  int totEvents;

  //store the moments of each combination of particle filters (1, 2 ... 4, 11, 12, ... 44 ...) per event
  //first index is event number
  //calculated in normal order, but then changed into recursive order.
  //temporary used variable to help code
  double ** SValues;

  //Used to store the name of the histogram in the createHistogramRec function
  //Later used to name the historgams
  //Not really necessary, but getting rid of them would also require major modifications to the SaveHistograms function 
  //bc the histograms would be created in a different order.
  TString** names;
  TString* *titles;
  TString** names2;
  TString* *titles2;

  ClassDef(PTHistos,0)
};

#endif /* WAC_PTHistos  */
