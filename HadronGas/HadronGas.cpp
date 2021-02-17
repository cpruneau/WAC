// Copyright 2016 Chun Shen

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "HadronGas.hpp"

ClassImp(HadronGas);

using namespace std;

HadronGas::HadronGas(string filename)
:
particleListFilename(),
particleList(),
stableParticleList(),
nPartTypes(0),
nStablePartTypes(0),
temperature(0),
mu(0),
muB(0),
muS(0),
muQ(0),
systemEnergyDensity(0),
systemEntropyDensity(0),
systemPressure(0),
netBDensity(0),
netSDensity(0),
netQDensity(0)
{
  particleListFilename = filename;  // filename of pdg data file
  readParticleListFromFile(filename);
  sortParticleListByParticleMass();
  createStableParticleList();
  calculateParticleDecayProbability();
}

HadronGas::~HadronGas()
{
  particleList.clear();
  stableParticleList.clear();
}

void HadronGas::readParticleListFromFile(string inputFileName)
{
  // read in Hadron information from pdg data file
  cout << "<I> HadronGas::readParticleListFromFile() Reading particle list from file:" <<  inputFileName << endl;
  ifstream inputFile(inputFileName.c_str());
  int pdgCode;
  string name;
  double mass, width;
  int gSpin, gIsospin;
  int baryon, strange, charm, bottom, charge;
  int nDecays, decayNpart;
  double decayBranchingRatio;
  int decayPart[5] = {0, 0, 0, 0, 0};

  int dummy_int;
  while (1)
    {
    //cout << "Reading in Hadron resonance decay table... - 2 -" << endl;
    inputFile >> pdgCode;
    if (inputFile.eof())
      break;
    inputFile >> name;
    inputFile >> mass;
    inputFile >> width;
    inputFile >> gSpin;
    inputFile >> baryon;
    inputFile >> strange;
    inputFile >> charm;
    inputFile >> bottom;
    inputFile >> gIsospin;
    inputFile >> charge;
    inputFile >> nDecays;
    //cout << "Reading in Hadron resonance decay table... - 3 - " << endl;

    particleList.push_back(new Hadron(pdgCode, name, mass, width, gSpin, baryon, strange,
                                            charm, bottom, gIsospin, charge, nDecays));
    // If baryon, also create/add entry for anti-baryon.
    if (baryon == 1)
      {
      ostringstream antiname;
      antiname << "Anti-" << name;
      particleList.push_back(new Hadron(-pdgCode, antiname.str(), mass, width, gSpin,
                                              -baryon, -strange, -charm, -bottom, gIsospin,
                                              -charge, nDecays));
      }
    //cout << "Reading in Hadron resonance decay table... - 4-" << endl;

    // read decay information
    for (int j = 0; j < nDecays; j++)
      {
      inputFile >> dummy_int;
      inputFile >> decayNpart;
      inputFile >> decayBranchingRatio;
      inputFile >> decayPart[0];
      inputFile >> decayPart[1];
      inputFile >> decayPart[2];
      inputFile >> decayPart[3];
      inputFile >> decayPart[4];
      decayNpart = abs(decayNpart);

      int* tempptr = new int[decayNpart];
      for (int ipart = 0; ipart < decayNpart; ipart++) tempptr[ipart] = decayPart[ipart];

      if (baryon == 0)
        {
        particleList.back()->addResonanceDecays(decayBranchingRatio,decayNpart, tempptr);
        }
      else
        {
        particleList.at(particleList.size()-2)->addResonanceDecays(decayBranchingRatio, decayNpart, tempptr);
        for (int ipart = 0; ipart < decayNpart; ipart++)
          {
          int HadronId = findParticleIndexFor(tempptr[ipart]);
          tempptr[ipart] = particleList[HadronId]->getAntiHadronPdgCode();
          }
        particleList.back()->addResonanceDecays(decayBranchingRatio, decayNpart, tempptr);
        }
      delete [] tempptr;
      }
    }
  //cout << "Reading in Hadron resonance decay table... - 5- " << endl;

  inputFile.close();
  particleList.erase(particleList.begin());  // delete gamma
  cout << "<I> All done! AntiHadrons are added!" << endl;
  nPartTypes = int( particleList.size() );
  cout << "<I> Total number of Hadrons: " <<  nPartTypes << endl;
}

ostream & HadronGas::printParticleBasicProperties(ostream & os)
{
  os << "==================================================================================" << std::endl;
  os << " Basic Particle Properties" << endl;
  os << "==================================================================================" << std::endl;
  for (int k = 1; k < nPartTypes; k++)
  {
  particleList[k]->printBasicProperties(os) << endl;
  }
  return os;
}


ostream &  HadronGas::printParticleYields(ostream & os)
{
  os << "==================================================================================" << std::endl;
  os << " Basic Particle Properties" << endl;
  os << " Temperature: " << temperature  << endl;
  os << "         muB: " << muB << endl;
  os << "==================================================================================" << std::endl;
  os << setw(15) << "Name" << setw(15) <<  "Yield" << setw(15) << "Thermal Yield" << endl;
  os << "==================================================================================" << std::endl;
  for (int i = 0; i < nStablePartTypes; i++)
  {
  os << setw(15) << stableParticleList[i]->getName()
  <<   scientific << setw(15)<< setprecision(4) << stableParticleList[i]->getStableParticleYield()
  <<   scientific << setw(15)<< setprecision(4) << stableParticleList[i]->getParticleYield()
  << endl;
  }
  return os;
}

ostream & HadronGas::printParticleCalcProperties(ostream & os)
{
  for (int k = 1; k < nPartTypes; k++)
  {
  particleList[k]->printCalcProperties(os) << endl;
  }
  return os;
}


void HadronGas::sortParticleListByParticleMass()
{
  cout << "<D> HadronGas::sortHadronListByHadronMass() List size:"
  << particleList.size() << endl;
  double m1, m2;
  for (int i = 1; i < nPartTypes; i++)
  {
  //cout << " i: " << i << endl;
  int k = i;
  int j = i - 1;
  m1 = particleList[k]->getMass();
  m2 = particleList[j]->getMass();
  //cout << "k:" << k << " m1: " << m1 << " j:" << j << " m2:" << m2 << endl;
  while (j >= 0 && (particleList[k]->getMass() < particleList[j]->getMass()) )
    {
    //cout << "k:" << k << " m1: " << m1 << " j:" << j << " m2:" << m2 << endl;
    Hadron* temp = particleList[j];
    particleList[j] = particleList[k];
    particleList[k] = temp;
    k--;
    j--;
    }
  }
  //   for (int i = 0; i < nPartTypes; i++)
  //     {
  //       cout << particleList[i]->getMass() << endl;
  //     }
}

double * createArray1D(int size)
{
  double  * array  = new double [size];
  for (int j = 0; j < size; j++) array[j] = 0.0;
  return array;
}

double ** createArray2D(int size)
{
  double  ** array  = new double* [size];
  for (int j = 0; j < size; j++)
  {
  array[j] = new double[size];
  for (int k = 0; k < size; k++)
    {
    array[j][k] = 0.0;
    }
  }
  return array;
}

void deleteArray2D(double  ** array, int size)
{
  for (int j = 0; j < size; j++) delete[] array[j];
  delete[] array;
}

void HadronGas::createStableParticleList()
{
  for (int iPart = 0; iPart < nPartTypes; iPart++)
  {
  if (particleList[iPart]->isStable())
    {
    stableParticleList.push_back(particleList[iPart]);
    }
  }
  nStablePartTypes = int( stableParticleList.size() );
  cout << "Total number of stable Hadrons: "  << nStablePartTypes << endl;
}

//void HadronGas::calculateParticleDecayProbability -- Original Version()
//{
//  if (nStablePartTypes<1)
//    {
//    cout << "<F> HadronGas::calculateParticleDecayProbability() Must call createStableParticleList() before calling this function" << endl;
//    return;
//    }
//  for (int iPart = 0; iPart < nPartTypes; iPart++)
//  {
//  particleList[iPart]->createDecayProbabilityArray(nStablePartTypes);
//  for (int jStable = 0; jStable < nStablePartTypes; jStable++)
//    {
//    if (particleList[iPart]->getPdgCode() == stableParticleList[jStable]->getPdgCode())
//      {
//      //particleList[iPart]->setStable(1); // not necessary?
//      particleList[iPart]->setDecayProbability(jStable, 1.0);
//      break;
//      }
//    }
//  }
//
//  for (int iPart = 0; iPart < nPartTypes; iPart++)
//  {
//  if (particleList[iPart]->getStable() == 0)
//    {  // unstable resonances
//      double  * temp  = createArray1D(nStablePartTypes);
//      double ** temp2 = createArray2D(nStablePartTypes);
//      for (int jChannel = 0; jChannel < particleList[iPart]->getNdecayChannel(); jChannel++)
//      {
//      for (int kParticle = 0; kParticle < abs(particleList[iPart]->getDecaysNpart(jChannel)); kParticle++)
//        {
//        // disregard photons...
//        if (particleList[iPart]->getDecaysHadron(jChannel, kParticle) == 22)  continue;
//
//        int pdg1 = particleList[iPart]->getDecaysHadron(jChannel, kParticle);
//        int lPart = findParticleIndexFor(pdg1);
//        if (lPart<0)
//          {
//          cout << "<F> Could not find PDG code: " << pdg1 << endl;
//          cout << "<F> Terminate task." << endl;
//          return;
//          }
//        if (lPart > iPart)
//          {
//          cout << "<F> Particle appears to decay into heavier particle:" << particleList[lPart]->getPdgCode() << "   " << particleList[iPart]->getPdgCode()
//          << "   "
//          << particleList[iPart]->getDecaysHadron(jChannel, kParticle)
//          << endl;
//          return;
//          }
//        for (int mStable = 0; mStable < nStablePartTypes; mStable++)
//        temp[mStable] += (particleList[iPart]->getDecaysBranchingRatio(jChannel) * particleList[lPart]->getDecayProbability(mStable));
//        }
//      }
//      // store the results
//      for (int jStable = 0; jStable < nStablePartTypes; jStable++)
//        {
//        particleList[iPart]->setDecayProbability(jStable, temp[jStable]);
//        for (int kStable = 0; kStable < nStablePartTypes; kStable++)
//          {
//          particleList[iPart]->setDecayProbability(jStable,kStable, temp2[jStable][kStable]);
//          }
//        }
//      delete [] temp;
//      deleteArray2D(temp2,nStablePartTypes);
//    }
//  }
//}

void HadronGas::calculateParticleDecayProbability()
{
  if (nStablePartTypes<1)
    {
    cout << "<F> HadronGas::calculateParticleDecayProbability() Must call createStableParticleList() before calling this function" << endl;
    return;
    }

  cout << "----------------------------------------------------------------------------------------" <<  endl;
  cout << "----------------------------------------------------------------------------------------" <<  endl;
  cout << "<I> HadronGas::calculateParticleDecayProbability() nPartTypes:" << nPartTypes << endl;
  cout << "----------------------------------------------------------------------------------------" <<  endl;

  for (int iPart = 0; iPart < nPartTypes; iPart++) particleList[iPart]->createDecayProbabilityArray(nStablePartTypes);
  stableParticlePairCorrelatedYield = createArray2D(nStablePartTypes);

  // paradoxically, we set stable particle to have a decay probability of one.
  for (int iPart = 0; iPart < nStablePartTypes; iPart++)
  {
  stableParticleList[iPart]->setDecayProbability(iPart, 1.0);
  }

  int nDecayParticles;
  int pdgOfDecayParts[20];
  int indexOfDecayParts[20];
  // For each unstable particle
  for (int iPart = 0; iPart < nPartTypes; iPart++)
    {
    if (particleList[iPart]->isStable()) continue;

    //113
    //if (particleList[iPart]->getPdgCode() != 113) continue;

    //int iPdgCode = particleList[iPart]->getPdgCode();
    double  * temp  = createArray1D(nStablePartTypes);
    double ** temp2 = createArray2D(nStablePartTypes);
    for (int jChannel = 0; jChannel < particleList[iPart]->getNdecayChannel(); jChannel++)
      {
      double branchingRatio = particleList[iPart]->getDecaysBranchingRatio(jChannel);
      nDecayParticles = 0;
      for (int kParticle = 0; kParticle < abs(particleList[iPart]->getDecaysNpart(jChannel)); kParticle++)
        {
        int pdgCode = particleList[iPart]->getDecaysHadron(jChannel, kParticle);
        // disregard photons...
        if (pdgCode == 22)  continue;
        pdgOfDecayParts[nDecayParticles] = pdgCode;
        int index = findParticleIndexFor(pdgCode);
        if (index<0)
          {
          cout << "<F> HadronGas::calculateParticleDecayProbability() Logic error; lPart1=" << index << endl;
          cout << "<F> Could not find PDG code: " << pdgCode << endl;
          cout << "<F> Terminate task." << endl;
          return;
          }
        if (index > iPart)
          {
          cout << "<F> HadronGas::calculateParticleDecayProbability() Logic error" << endl;
          cout << "<F> Particle appears to decay into heavier particle: iPart:" << iPart << " into index:" << index << endl;
          cout << "<F> Terminate task." << endl;
          return;
          }
        indexOfDecayParts[nDecayParticles] = index;
        nDecayParticles++;
        }

     // cout << "    nDecayParticles: " << nDecayParticles << endl;

      // now we have a local list of the pdg ids iPart decays into (w/o photons)
      // must look at singles and then at pairs
      // if there are no particles left, just skip this decay channel...
      double prob1, prob2;
      for (int kParticle1 = 0; kParticle1 < nDecayParticles; kParticle1++)
        {
        int lPart1 = indexOfDecayParts[kParticle1];
        //cout << " kParticle1:" << kParticle1 << "  lPart1: " << lPart1 << endl;
        for (int mStable1 = 0; mStable1 < nStablePartTypes; mStable1++)
          {
          prob1 = particleList[lPart1]->getDecayProbability(mStable1);
          temp[mStable1] += branchingRatio * prob1;


          for (int kParticle2 = 0; kParticle2 < nDecayParticles; kParticle2++)
            {
            if (kParticle1 != kParticle2)
              {
              int lPart2 = indexOfDecayParts[kParticle2];
              //cout <<  " kParticle1:" << kParticle1 << "  lPart1: " << lPart1  << " kParticle2:" << kParticle2 << "  lPart2: " << lPart2 << endl;
              for (int mStable2 = 0; mStable2 < nStablePartTypes; mStable2++)
                {
                prob2 = particleList[lPart2]->getDecayProbability(mStable2);
                temp2[mStable1][mStable2] += branchingRatio*prob1*prob2;
                //if (iPart<20 && mStable1<3 && mStable2<5)  cout <<  "--2--- temp2[" << mStable1 << "," << mStable2 << "] ="  << temp2[mStable1][mStable2] << endl;
                }
              }
            }

          for (int mStable3 = 0; mStable3 < nStablePartTypes; mStable3++)
            {
            prob1 = particleList[lPart1]->getDecayProbability(mStable1,mStable3);
            //if (branchingRatio * prob1>0.0) cout << " NON ZERO HAPPENS lPart1: " << lPart1 << " mStable3:" << mStable3<< endl;
            temp2[mStable1][mStable3] += branchingRatio *prob1;
            }



          }
        }
      }
    // store the results
    for (int jStable1 = 0; jStable1 < nStablePartTypes; jStable1++)
      {
      particleList[iPart]->setDecayProbability(jStable1, temp[jStable1]);
      for (int jStable2 = 0; jStable2 < nStablePartTypes; jStable2++)
        {
        //if (iPart<20 && jStable1<3 && jStable2<5) cout <<  "--3-- temp2[" << jStable1 << "," << jStable2 << "] ="  << temp2[jStable1][jStable2] << endl;
        particleList[iPart]->setDecayProbability(jStable1,jStable2, temp2[jStable1][jStable2]);
        }
      }
    delete [] temp;
    deleteArray2D(temp2,nStablePartTypes);
    }
}



int HadronGas::findParticleIndexFor(int pdgCode)
{
  nPartTypes = int ( particleList.size() ); // must keep updating nPartTypes because the list may be changing.
  for (int iPart = 0; iPart < nPartTypes; iPart++)
  {
  if (pdgCode == particleList[iPart]->getPdgCode()) return iPart;
  }
  return -1;
}

int HadronGas::findStableParticleIndexFor(int pdgCode)
{
  nStablePartTypes  = int ( stableParticleList.size() );
  for (int iPart = 0; iPart < nStablePartTypes; iPart++)
  {
  if (pdgCode == stableParticleList[iPart]->getPdgCode()) return iPart;
  }
  return -1;
}



void HadronGas::calculateParticleDecays()
{
  double stableYield;
  for (int jStable = 0; jStable < nStablePartTypes;jStable++)
  {
  stableYield = 0.0;
  for (int kPart = 0; kPart < nPartTypes;kPart++)
    {
    stableYield += (particleList[kPart]->getDecayProbability(jStable)*particleList[kPart]->getParticleYield());
    }
  stableParticleList[jStable]->setStableParticleYield(stableYield);
  }

  for (int jStable1 = 0; jStable1 < nStablePartTypes; jStable1++)
  {
  stableYield = 0.0;
  for (int jStable2 = 0; jStable2 < nStablePartTypes; jStable2++)
    {
    stableYield = 0.0;
    for (int kPart = 0; kPart < nPartTypes;kPart++)
      {
      stableYield += (particleList[kPart]->getDecayProbability(jStable1,jStable2)*particleList[kPart]->getParticleYield());
      }
    stableParticlePairCorrelatedYield[jStable1][jStable2] = stableYield;
    }
  }

}

void HadronGas::outputParticleChemicalPotentials(ChemicalPotential* muTable)
{
  ofstream Hadron_mu_table("chemical_potentials.dat");
  // output names
  Hadron_mu_table << scientific << setw(18) << setprecision(8)
  << 0.0 << "    ";
  for (int i = 0; i < nPartTypes; i++) {
    Hadron_mu_table << scientific << setw(18) << setprecision(8)
    << particleList[i]->getName() << "    ";
  }
  Hadron_mu_table << endl;
  // output mass
  Hadron_mu_table << scientific << setw(18) << setprecision(8)
  << 0.0 << "    ";
  for (int i = 0; i < nPartTypes; i++) {
    Hadron_mu_table << scientific << setw(18) << setprecision(8)
    << particleList[i]->getMass() << "    ";
  }
  Hadron_mu_table << endl;
  // output baryon number
  Hadron_mu_table << scientific << setw(18) << setprecision(8)
  << 0.0 << "    ";
  for (int i = 0; i < nPartTypes; i++) {
    Hadron_mu_table << scientific << setw(18) << setprecision(8)
    << particleList[i]->getBaryon() << "    ";
  }
  Hadron_mu_table << endl;

  double Ti = 0.1;
  double Tf = 0.2;
  double dT = 0.001;
  int nT = (Tf - Ti)/dT + 1;
  for (int i = 0 ; i < nT; i++) {
    double T_local = Ti + i*dT;
    calculateParticleChemicalPotential(T_local, muTable);
    // calculateParticleChemicalPotential2(T_local, muTable);
    Hadron_mu_table << scientific << setw(18) << setprecision(8)
    << T_local << "    ";
    for (int j = 0; j < nPartTypes; j++)
    Hadron_mu_table << scientific << setw(18) << setprecision(8)
    << particleList[j]->getMu() << "    ";
    Hadron_mu_table << endl;
  }
  Hadron_mu_table.close();
}

void HadronGas::calculateParticleChemicalPotential2(double temperature, ChemicalPotential* muTable)
{
  int nStablePartTypes = muTable->getNStableParticles();
  double *muTableStableHadron = new double[nStablePartTypes];
  muTable->output_stable_mu(temperature, muTableStableHadron);
  for (int i = 0; i < nPartTypes ; i++)
  {
  double mu_temp = 0.0;
  for (int j = 0; j < nStablePartTypes; j++) mu_temp += particleList[i]->getDecayProbability(j)*muTableStableHadron[j];
  particleList[i]->setMu(mu_temp);
  }
  delete [] muTableStableHadron;
}

void HadronGas::calculateParticleChemicalPotential(double temperature, ChemicalPotential* muTable)
{
  int N_mu = muTable->getNStableParticles();
  double *muTableStableHadron = new double[N_mu];
  muTable->output_stable_mu(temperature, muTableStableHadron);

  int Idummy;
  char cdummy[256];
  ifstream Hadrontable("EOS/EOS_Hadrontable.dat");
  Hadrontable >> nStablePartTypes;
  if (N_mu != nStablePartTypes)
    {
    cout << "Chemical potential table is not compatible with EOS_Hadrontable.dat" << endl;
    exit(1);
    }
  double *stable_Hadron_monval = new double[nStablePartTypes];
  for (int i = 0; i < nStablePartTypes; i++) {
    Hadrontable >> Idummy >> stable_Hadron_monval[i];
    Hadrontable.getline(cdummy, 256);
  }
  Hadrontable.close();

  for (int i = 0; i < nStablePartTypes; i++)
  for (int j = 0; j < nPartTypes; j++)
  if (particleList[j]->getPdgCode() == stable_Hadron_monval[i])
    {
    particleList[j]->setStable(1);
    particleList[j]->setMu(muTableStableHadron[i]);
    break;
    }
  for (int i = 0; i < nPartTypes ; i++)
  {
  double mu_temp = 0.0;
  if (particleList[i]->getStable() == 0) {
    for (int j=0; j < particleList[i]->getNdecayChannel(); j++) {
      for (int k=0; k < abs(particleList[i]->getDecaysNpart(j)); k++) {
        for (int l=0; l < nPartTypes; l++) {
          if (particleList[i]->getDecaysHadron(j, k)
              == particleList[l]->getPdgCode()) {
            mu_temp += (particleList[i]->getDecaysBranchingRatio(j)
                        *particleList[l]->getMu());
            break;
          }
          if (l == (nPartTypes - 1)
              && particleList[i]->getDecaysHadron(j, k) != 22)
            cout << "warning: can not find Hadron "
            << particleList[i]->getDecaysHadron(j, k) << endl;
        }
      }
    }
    particleList[i]->setMu(mu_temp);
  }
  }

  delete [] stable_Hadron_monval;
  delete [] muTableStableHadron;
}

Hadron * HadronGas::getParticle(int index)
{
  if (index<0) index = 0;
  if (index>=nPartTypes) index = nPartTypes;
  return particleList[index];
}

Hadron * HadronGas::getStableParticle(int index)
{
  if (index<0) index = 0;
  if (index>=nStablePartTypes) index = nStablePartTypes;
  return stableParticleList[index];
}

////////////////////////////////////////////////////////////////////////////////
// calculate Hadron chemical potentials
// need to add support for partial chemical equilibrium
////////////////////////////////////////////////////////////////////////////////
void HadronGas::calculateParticleMu(double muB, double muS)
{
  for (int i = 0; i < nPartTypes; i++)
  particleList[i]->calculateChemicalPotential(muB, muS);
  return;
}

void HadronGas::calculateParticleYield(double temperature,double muB, double muS)
{
  // calculate Hadron yield
  calculateParticleMu(muB, muS);
  for (int i = 0; i < nPartTypes; i++)
  particleList[i]->calculateParticleYield(temperature, muB, muS);
  calculateParticleDecays();
  return;
}

void HadronGas::calculateSystemNetBaryonDensity(double temperature) {
  double result = 0.0;
  for (int i = 0; i < nPartTypes; i++)
  result += particleList[i]->getBaryon()*particleList[i]->getParticleYield();
  netBDensity = result;
  return;
}


void HadronGas::calculateSystemEnergyDensity(double temperature) {
  // calculate the energy density of the system at given T and mu
  double result = 0.0e0;
  for (int i = 0; i < nPartTypes; i++)
  result += particleList[i]->calculateEnergyDensity(temperature);
  systemEnergyDensity = result;
  return;
}

void HadronGas::calculateSystemPressure(double temperature) {
  // calculate the pressure of the system at given T and mu
  double result = 0.0e0;
  for (int i = 0; i < nPartTypes; i++)
  result += particleList[i]->calculatePressure(temperature);
  systemPressure = result;
}

void HadronGas::calculateSystemEntropyDensity(double temperature)
{
  // calculate the entropy density of the system at given T and mu
  double result = 0.0e0;
  for (int i = 0; i < nPartTypes; i++)
  result += particleList[i]->calculateEntropyDensity(temperature);
  systemEntropyDensity = result;
}

void HadronGas::calculateParticle_dsdmu(double temperature)
{
  // calculate dsdmu
  for (int i = 0; i < nPartTypes; i++)
  particleList[i]->calculate_dsdmu(temperature);
}


void HadronGas::calculateSystemProperties(double _temp, double _muB, double _muS)
{
  temperature = _temp;
  muB = _muB;
  muS = _muS;
  muQ = _muS;
  calculateParticleYield(temperature, muB, muS);
  calculateParticle_dsdmu(temperature);
  calculateSystemNetBaryonDensity(temperature);
  calculateSystemEnergyDensity(temperature);
  calculateSystemPressure(temperature);
  calculateSystemEntropyDensity(temperature);
}


void HadronGas::calculateSystemEOS_and_output_in_2D()
{
  // calculate the EOS of given system, e,p,s as functions of T and muB
  int nT = 200;
  double T_i = 0.02;         // unit: (GeV)
  double T_f = 0.2;          // unit: (GeV)
  double dT = (T_f - T_i)/(nT - 1 + 1e-15);
  int n_mu_B = 200;
  double mu_B_i = 0.0;      // unit: (GeV)
  double mu_B_f = 0.8;      // unit: (GeV)
  double dmu_B = (mu_B_f - mu_B_i)/(n_mu_B - 1 + 1e-15);
  double** temp_ptr = new double* [n_mu_B];
  double** muB_ptr = new double* [n_mu_B];
  double** ed_ptr = new double* [n_mu_B];
  double** sd_ptr = new double* [n_mu_B];
  double** pressure_ptr = new double* [n_mu_B];
  double** net_baryon_ptr = new double* [n_mu_B];
  double** cs2_ptr = new double* [n_mu_B];

  for (int i = 0; i < n_mu_B; i++) {
    temp_ptr[i] = new double[nT];
    muB_ptr[i] = new double[nT];
    ed_ptr[i] = new double[nT];
    sd_ptr[i] = new double[nT];
    pressure_ptr[i] = new double[nT];
    net_baryon_ptr[i] = new double[nT];
    cs2_ptr[i] = new double[nT];
    for (int j = 0; j < nT; j++) {
      temp_ptr[i][j] = 0.0;
      muB_ptr[i][j] = 0.0;
      ed_ptr[i][j] = 0.0;
      sd_ptr[i][j] = 0.0;
      pressure_ptr[i][j] = 0.0;
      net_baryon_ptr[i][j] = 0.0;
      cs2_ptr[i][j] = 0.0;
    }
  }

  double muS = 0.0;
  for (int i = 0; i < n_mu_B; i++) {
    double local_mu_B = mu_B_i + i*dmu_B;
    cout << " calculating EOS for muB = " << local_mu_B << " GeV..."
    << endl;
    for (int j = 0; j < nT; j++) {
      muB_ptr[i][j] = local_mu_B;
      temp_ptr[i][j] = T_i + j*dT;

      calculateParticleYield(temp_ptr[i][j], muB_ptr[i][j], muS);
      calculateSystemNetBaryonDensity(temp_ptr[i][j]);
      calculateSystemEnergyDensity(temp_ptr[i][j]);
      calculateSystemPressure(temp_ptr[i][j]);
      calculateSystemEntropyDensity(temp_ptr[i][j]);

      ed_ptr[i][j] = systemEnergyDensity;
      sd_ptr[i][j] = systemEntropyDensity;
      pressure_ptr[i][j] = systemPressure;
      net_baryon_ptr[i][j] = netBDensity;
    }
  }

  // calculate speed of sound cs^2 = dP/de + n/(e+P) dP/drho_B
  for (int i = 0; i < n_mu_B; i++) {
    for (int j = 0; j < nT; j++) {
      double dPde, dPdrho_B;
      if (j == 0)
        dPde = ((pressure_ptr[i][1] - pressure_ptr[i][0])
                /(ed_ptr[i][1] - ed_ptr[i][0] + 1e-15));
      else if (j == nT - 1)
        dPde = ((pressure_ptr[i][nT-1] - pressure_ptr[i][nT-2])
                /(ed_ptr[i][nT-1] - ed_ptr[i][nT-2] + 1e-15));
      else
        dPde = ((pressure_ptr[i][j+1] - pressure_ptr[i][j-1])
                /(ed_ptr[i][j+1] - ed_ptr[i][j-1] + 1e-15));

      if (i == 0)
        dPdrho_B = ((pressure_ptr[1][j] - pressure_ptr[0][j])
                    /(net_baryon_ptr[1][j] - net_baryon_ptr[0][j]
                      + 1e-15));
      else if (i == n_mu_B - 1)
        dPdrho_B = (
                    (pressure_ptr[n_mu_B-1][j] - pressure_ptr[n_mu_B-2][j])
                    /(net_baryon_ptr[n_mu_B-1][j] - net_baryon_ptr[n_mu_B-2][j]
                      + 1e-15));
      else
        dPdrho_B = ((pressure_ptr[i+1][j] - pressure_ptr[i-1][j])
                    /(net_baryon_ptr[i+1][j] - net_baryon_ptr[i-1][j]
                      + 1e-15));

      double v_sound_sq = (dPde + net_baryon_ptr[i][j]/(ed_ptr[i][j]
                                                        + pressure_ptr[i][j])*dPdrho_B);
      cs2_ptr[i][j] = v_sound_sq;
    }
  }

  // output EOS table
  ostringstream EOSfilename;
  EOSfilename << "./EOS_2D_table.dat";
  ofstream output(EOSfilename.str().c_str());
  for (int i = 0; i < n_mu_B; i++) {
    for (int j = 0; j < nT; j++)
    output << scientific << setw(20) << setprecision(8)
    << temp_ptr[i][j] << "   " << muB_ptr[i][j] << "   "
    << ed_ptr[i][j] << "   " << net_baryon_ptr[i][j] << "   "
    << sd_ptr[i][j] << "   " << pressure_ptr[i][j] << "   "
    << cs2_ptr[i][j] << endl;
  }
  output.close();

  // clean up
  for (int i = 0; i < n_mu_B; i++) {
    delete [] temp_ptr[i];
    delete [] muB_ptr[i];
    delete [] ed_ptr[i];
    delete [] sd_ptr[i];
    delete [] pressure_ptr[i];
    delete [] net_baryon_ptr[i];
    delete [] cs2_ptr[i];
  }
  delete [] temp_ptr;
  delete [] muB_ptr;
  delete [] ed_ptr;
  delete [] sd_ptr;
  delete [] pressure_ptr;
  delete [] net_baryon_ptr;
  delete [] cs2_ptr;
}

void HadronGas::calculateSystemEOS(double muB, double muS) {
  // calculate the EOS of given system, e,p,s as functions of T
  // at given muB and muS
  cout << "calculate the EOS of the system with muB = " << muB
  << " GeV and muS = " << muS << " GeV .... " << endl;
  int nT = 200;
  double T_i = 0.01;        // unit: (GeV)
  double T_f = 0.2;          // unit: (GeV)
  double dT = (T_f - T_i)/(nT - 1);
  double* temp_ptr = new double[nT];
  double* ed_ptr = new double[nT];
  double* sd_ptr = new double[nT];
  double* pressure_ptr = new double[nT];
  double* net_baryon_ptr = new double[nT];
  double* cs2_ptr = new double[nT];
  for (int i = 0; i < nT; i++) {
    temp_ptr[i] = T_i + i*dT;
    calculateParticleYield(temp_ptr[i], muB, muS);
    calculateSystemNetBaryonDensity(temp_ptr[i]);
    calculateSystemEnergyDensity(temp_ptr[i]);
    calculateSystemPressure(temp_ptr[i]);
    calculateSystemEntropyDensity(temp_ptr[i]);
    ed_ptr[i] = systemEnergyDensity;
    sd_ptr[i] = systemEntropyDensity;
    pressure_ptr[i] = systemPressure;
    net_baryon_ptr[i] = netBDensity;
  }
  // calculate speed of sound cs^2 = dP/de
  for (int i = 0; i < nT - 1; i++)
  cs2_ptr[i] = ((pressure_ptr[i+1] - pressure_ptr[i])
                /(ed_ptr[i+1] - ed_ptr[i] + 1e-30));
  cs2_ptr[nT-1] = cs2_ptr[nT-2];

  // output EOS table
  ostringstream EOSfilename;
  EOSfilename << "./EOS_muB_" << muB << "_muS_" << muS << ".dat";
  ofstream output(EOSfilename.str().c_str());
  for (int i = 0; i < nT; i++)
  output << scientific << setw(20) << setprecision(8)
  << net_baryon_ptr[i] << "  "
  << temp_ptr[i] << "   "  << ed_ptr[i] << "   "
  << sd_ptr[i] << "   " << pressure_ptr[i] << "   "
  << cs2_ptr[i] << endl;
  output.close();
  delete [] temp_ptr;
  delete [] ed_ptr;
  delete [] sd_ptr;
  delete [] pressure_ptr;
  delete [] net_baryon_ptr;
  delete [] cs2_ptr;
  return;
}

// development

//void HadronGas::calculateParticleTupletProbability()
//{
//  // Build stable particle list.
//  int nPartTypes = nPartTypes;
//  for (int iPart = 0; iPart < nPartTypes; iPart++)
//    {
//    if (particleList[iPart]->isStable())
//      {
//      stableParticleList.push_back(particleList[iPart]);
//      }
//    }
//  int nStablePartTypes = nStablePartTypes;
//  cout << "Total number of stable Hadrons: "  << nStablePartTypes << endl;
//
//  for (int iPart = 0; iPart < nPartTypes; iPart++)
//    {
//    particleList[iPart]->createDecayProbabilityArray(nStablePartTypes);
//    int pdgCode = particleList[iPart]->getPdgCode();
//    for (int j = 0; j < nStablePartTypes; j++)
//      {
//      if (pdgCode == stableParticleList[j]->getPdgCode())
//        {
//        particleList[iPart]->setStable(1); // not necessary?
//        particleList[iPart]->setDecayProbability(j, 1.0);
//        break;
//        }
//      }
//    }
//
//  double * instantiateAndSet1DArray(int nX, double initialValue)
//  {
//  double* temp = new double[nX];
//  for (int iX = 0; iX < nX; iVal++) temp[iX] = initialValue;
//  return temp;
//  }
//
//  for (int i = 0; i < nPartTypes; i++)
//    {
//    if (particleList[i]->getStable() == 0)
//      {  // unstable resonances
//        double* temp = instantiateAndSet1DArray(nStablePartTypes,0.0);
//        for (int j = 0; j < particleList[i]->getNdecayChannel(); j++)
//          {
//          for (int k = 0; k < abs(particleList[i]->getDecaysNpart(j)); k++)
//            {
//            // disregard photons...
//            if (particleList[i]->getDecaysHadron(j, k) == 22)  continue;
//            for (int l = 0; l < nPartTypes; l++)
//              {
//              if (particleList[i]->getDecaysHadron(j, k) == particleList[l]->getPdgCode())
//                {
//                for (int m = 0; m < nStablePartTypes; m++)
//                  temp[m] += (particleList[i]->getDecaysBranchingRatio(j) * particleList[l]->getDecayProbability(m));
//                break;
//                }
//              if (l > i)
//                {
//                cout << particleList[l]->getPdgCode()
//                << "   " << particleList[i]->getPdgCode()
//                << "   "
//                << particleList[i]->getDecaysHadron(j, k)
//                << endl;
//                }
//              }
//            }
//          }
//        for (int j = 0; j < nStablePartTypes; j++) particleList[i]->setDecayProbability(j, temp[j]);
//        delete [] temp;
//      }
//    }
//}
