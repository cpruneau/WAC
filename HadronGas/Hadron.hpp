// Copyright 2016 Chun Shen
#ifndef WAC_Hadron
#define WAC_Hadron
#include <iostream>
#include <iomanip>
#include <string>
#include "TString.h"

using namespace std;

class Hadron
{
private:
  static double hbarC;     // GeV*fm
  static double pionMass;  // Monte-Carlo number according PDG
  static int    trunOrder; // truncated order in the summation

  int pdgCode;        // Monte-Carlo number according PDG
  std::string name;   // Hadron name
  double mass;        // Hadron mass (GeV)
  double width;       // decay width
  int gSpin;          // spin degeneracy
  int baryon;         // baryon number
  int strange;        // strangeness
  int charm;          // charmness
  int bottom;         // bottomness
  int gIsospin;       // isospin degeneracy
  int charge;         // charge
  int nDecays;         // amount of decays listed for this resonance
  int stable;         // defines whether this Hadron is stable
  double mu;          // chemical potential
  int sign;           // Bosons or Fermions
  int channelIdx;

  double yield;        // Hadron thermal yield at given T and mu
  double stableYield;  // Hadron yield after decays

  double dPoverTdmu;
  double deoverTdmu;
  double dndmu;
  double dsdmu;
  double energyDensity;
  double entropyDensity;
  double pressure;     // thermodynamic quantities at given T and mu

  int     nDecayChannels;        // number of decay channels
  double* decaysBranchingRatio;  // branching ratio of each decay channel
                                 // number of daughter Hadrons of each decay channel
  int*  decaysNPart;
  int** decaysHadron;      // identity of daughter Hadrons

  // array to record Hadron decay probability
  // for the Hadron decays into stable ones
  double* decayProbability;
  double** decayProbabilityPairs;
//  double* decayStableProbability;
//  double** decayStableProbabilityPairs;


public:
  Hadron(int monval_in, std::string name_in, double mass_in,
           double width_in, int gSpin_in, int baryon_in, int strange_in,
           int charm_in, int bottom_in, int gIsospin_in, int charge_in,
           int NdecayChannel_in);
  Hadron(const Hadron & source);
  virtual ~Hadron();

  Hadron & operator=(const Hadron & source);

  void instantiateArrays();
  void addResonanceDecays(double branchratio, int Npart, int* decaysHadron);
  void createDecayProbabilityArray(int nStable);
  void calculateChemicalPotential(double muB, double muS);
  void calculateParticleYield(double temperature, double muB, double muS);
  double calculateEnergyDensity(double temperature);
  double calculatePressure(double temperature);
  double calculateEntropyDensity(double temperature);

  int getAntiHadronPdgCode();

  int getPdgCode() {return(pdgCode);}
  std::string getName() {return(name);}
  double getMass() {return(mass);}
  int getBaryon() {return(baryon);}
  int getCharge() {return(charge);}
  int getSpinfactor() {return(gSpin);}
  double getMu() {return(mu);}
  int getSign() {return(sign);}
  double getParticleYield() {return(yield);}
  void setStableParticleYield(double yield_in) {stableYield = yield_in;}
  double getStableParticleYield() {return(stableYield);}


  int getNdecays() {return(nDecays);}
  int getNdecayChannel() {return(nDecayChannels);}
  int getDecaysNpart(int i) {return(decaysNPart[i]);}
  int getDecaysHadron(int i, int j) {return(decaysHadron[i][j]);}
  double getDecaysBranchingRatio(int i) {return(decaysBranchingRatio[i]);}
  int getStable() {return(stable);}
  bool isStable() {return(stable==1);}
  void setStable(int s) {stable = s;}
  void setMu(double chem) {mu = chem;}
  double getDecayProbability(int i) {return(decayProbability[i]);}
  double getDecayProbability(int i,int j) {return(decayProbabilityPairs[i][j]);}

  void setDecayProbability(int i, double val)
  {
  decayProbability[i] = val;
  }
  void setDecayProbability(int i, int j, double val)
  {
  decayProbabilityPairs[i][j] = val;
  }

  double calculate_dsdmu(double temperature);
  double calculate_deoverTdmu(double temperature);
  double calculate_dPoverTdmu(double temperature);
  double calculate_dndmu(double temperature);

  void getBulkViscosityCoefficients(double Tdec, double* bulkvisCoefficients);

  ostream & printBasicProperties(ostream & os);
  ostream & printCalcProperties(ostream & os);

  ClassDef(Hadron,0)
};

#endif  // SRC_Hadron_H_

