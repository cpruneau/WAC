// Adapted to root envinromemt 2021 C. Pruneau
//
// Copyright 2016 Chun Shen
#ifndef WAC_HadronGas
#define WAC_HadronGas
#include <iostream>
#include "TString.h"
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>

#include "Hadron.hpp"
#include "ChemicalPotential.hpp"

using namespace std;

class HadronGas
{
private:
  string particleListFilename;
  vector<Hadron*> particleList;
  vector<Hadron*> stableParticleList;


  double temperature;
  double mu;
  double muB;
  double muS;
  double muQ;
  double systemEnergyDensity;
  double systemEntropyDensity;
  double systemPressure;
  double netBDensity;
  double netSDensity;
  double netQDensity;

public:
  int nPartTypes;
  int nStablePartTypes;

  double ** stableParticlePairYield;
  double ** stableParticlePairCorrelatedYield;

  HadronGas(string HadronTableName);
  virtual ~HadronGas();
  void readParticleListFromFile(string tableName);
  void sortParticleListByParticleMass();
  void createStableParticleList();
  void calculateParticleDecayProbability();
  void calculateParticleDecays();
  void calculateParticleChemicalPotential2(double temperature, ChemicalPotential* muTable);
  void calculateSystemEOS(double muB = 0.0, double muS = 0.0);
  void calculateParticleMu(double muB, double muS);
  void calculateParticleYield(double temperature, double muB = 0.0, double muS = 0.0);
  void calculateSystemEnergyDensity(double temperature);
  void calculateSystemPressure(double temperature);
  void calculateSystemEntropyDensity(double temperature);
  void calculateSystemNetBaryonDensity(double temperature);
  void calculateParticle_dsdmu(double temperature);
  void calculateSystemProperties(double temperature, double muB, double muS);
  void calculateParticleChemicalPotential(double temperature, ChemicalPotential* muTable);
  void calculateSystemEOS_and_output_in_2D();

  int getNHadrons()         { return(particleList.size()); }
  int getNStableHadrons()   { return(stableParticleList.size()); }
  int findParticleIndexFor(int pdgCode);
  int findStableParticleIndexFor(int pdgCode);
  Hadron * getParticle(int index);
  Hadron * getStableParticle(int index);
  
  double getTemperature() const{ return temperature; }
  double getMu() const { return mu; }
  double getMuB() const{ return muB; }
  double getMuS() const{ return muS; }
  double getMuQ() const{ return muQ; }
  double getEnergyDensity()    const{ return systemEnergyDensity; }
  double getEntropyDensity()   const{ return systemEntropyDensity; }
  double getPressure()         const{ return systemPressure; }
  double getNetBaryonDensity() const{ return netBDensity; }
  double getNetStrangenessDensity() const{ return netSDensity; }
  double getNetChargeDensity() const{ return netQDensity; }

  ostream & printParticleBasicProperties(ostream & os);
  ostream & printParticleYields(ostream & os);
  ostream & printParticleCalcProperties(ostream & os);

  void outputParticleChemicalPotentials(ChemicalPotential* muTable);

  ClassDef(HadronGas,0)
};


#endif  // WAC_HadronGas

