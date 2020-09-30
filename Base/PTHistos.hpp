#ifndef WAC_PTHistos
#define WAC_PTHistos
#include "Histograms.hpp"
#include "Event.hpp"
#include "TMath.h"
#include <cstdlib>
#include <cmath>
#include <chrono>

class PTHistos : public Histograms
{
public:

  PTHistos(const TString & collectionName,
    AnalysisConfiguration * analysisConfiguration,
    LogLevel  debugLevel,
    int ord);
  PTHistos(TFile * inputFile,
    const TString & collectionName,
    AnalysisConfiguration * analysisConfiguration,
    LogLevel  debugLevel,
    int ord);
  virtual ~PTHistos();

  virtual void createHistograms();
  virtual void loadHistograms(TFile * inputFile);
  virtual void fillDerivedHistos(bool *** acceptances, double * mults, double * cents, double * avgCounts, double * avgpT, double ** SValues, int ** counts, int totEvents);
  virtual void saveHistograms(TFile * outputFile, bool saveAll=false);
  virtual void createHistogramRec(TString * baseName, TString * baseTitle, int depth, int partIndex);
  virtual void loadHistogramRec(TString * baseName, int depth, int partIndex, TFile * inputFile);
  virtual void fillEventHistos(double mult, double cent, double weight);
  virtual void fillNormalizedPTValues( int depth, int partIndex, double product, double * SValues, double  mult, double  cent, double *avgpT);
  virtual void fillNormalizedPTValues( int depth, int partIndex, double product, TH1 *** values, int* reorder, int*  nBin, int  totEvents, double *avgpT);
  virtual void calculateCumulants(TProfile ** Shistos, TH1 **CHistos, int nBins, double min, double max);
  virtual int convert(int * num, int len);
  virtual int* convert(int num, int & len);
  virtual void calcRecSum(TH1 **CHistos, int iBin, double& absESq, double curRelESq, int* iHisto, int* Subset, int len,  int * set, int lenSet, double productC, double* used, int& curInd, int productS, double& sum);
  virtual void convertToBinary(int num, char*str, int len );
  virtual int* getSubset(char* subset, int * set, int len, int& lenSub);
  virtual int* getComplementarySubset(char* subset, int * set, int len, int& lenSub);
  virtual int getSubsetNumber(int * subset, int lenSub, int * mainset, int lenSet);
  virtual void resetHistoRanges(int n);



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

  int histoIndex;

  //size of   TProfile ** h_c array
  int size; 

  //keep track of the orders of the correlation functions in the histograms
  int * orders;

  // Min bias all included NOT IN ORDER (They are in "recursive order")
  // in the order S, s, s*, C, c, 
  TH1 * h_events;
  //first index is function index, second index is histo index
  TProfile *** hS;
  TH1 *** hC;



  // vs Mult measured in fiducial NOT IN ORDER (They are in "recursive order")
  // in the order S, s, s*, C, c, 
  TH1 * h_events_vsMult;
  //first index is function index, second index is histo index
  TProfile *** hS_vsMult;
  TH1 *** hC_vsMult;


  // vs Centrality NOT IN ORDER (They are in "recursive order")
  // in the order S, s, s*, C, c, c*
  TH1 * h_events_vsCent;
  //first index is function index, second index is histo index
  TProfile *** hS_vsCent;
  TH1 *** hC_vsCent;


  //number of pairs, triples ... NOT IN ORDER (They are in "recursive order")
  TProfile ** h_counts;
  TProfile ** h_counts_vsMult;
  TProfile ** h_counts_vsCent;

  //number of functions we are counting 
  // here it is 6 (S, s, s*, C, c, c*)
  int numFunc;

  //keep track of which index in "recursive order" corresponds to which index in normal order
  int * reorder;

  ClassDef(PTHistos,0)
};

#endif /* WAC_PTHistos  */