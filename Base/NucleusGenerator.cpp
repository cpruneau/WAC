// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class NucleusGenerator
 \ingroup WAC

 Class defining NucleusGenerator
 */

#include "NucleusGenerator.hpp"
#include <iostream>
using namespace std;

ClassImp(NucleusGenerator);

NucleusGenerator::NucleusGenerator(const TString & _generatorName,
                                   GeneratorType _generatorType, double _parA, double _parB, double _parC,
                                   int _nR, double _minR, double _maxR)
:
generatorType(_generatorType),
generatorName(_generatorName),
nR(_nR),
minR(_minR),
maxR(_maxR),
parA(_parA),
parB(_parB),
parC(_parC),
rDensity(0),
rProfile(0),
rProfileGen(0),
useRecentering(1),
useNucleonExclusion(1),
exclusionRadius(0.4), // fm
exclusionRadiusSq(0.4*0.4)
{
  initialize();
}

// create the container but do not assign any positions or properties.
void NucleusGenerator::initialize()
{
  TString generatorTypeName;
  switch (generatorType)
       {
         case Uniform:     generatorTypeName = "_Uniform"; break;
         case WoodsSaxon:  generatorTypeName = "_WoodsSaxon"; break;
         case Exponential: generatorTypeName = "_Exponential"; break;
         case Gaussian:    generatorTypeName = "_Gaussian"; break;
         case DoubleGaussian: generatorTypeName = "_DoubleGaussian"; break;
       }

  TString histoName;
  histoName = generatorName + "_rDensity";
  rDensity = new TH1D(histoName,histoName,nR,minR, maxR);
  histoName = generatorName + "_rProfile";
  rProfile = new TH1D(histoName,histoName,nR,minR, maxR);
  histoName = generatorName + "_rProfileGen";
  rProfileGen  = new TH1D(histoName,histoName,nR,minR, maxR);

  double dr = (maxR-minR)/double(nR);
  double r  = minR + dr/2.0;
  double density, profile;
  for (int iR=0; iR<nR; iR++)
    {
    double r2 = r*r;
    switch (generatorType)
      {
        case Uniform: // uniform hard sphere
        density = (r*r*r<parA)? 1.0 : 0.0;
        break;
        case WoodsSaxon: // Woods-Saxon/Fermi
        density = 1.0/(1.0+exp((r-parA)/parB) );
        break;
        case Exponential: // exponential
        density = exp(-r/parA);
        break;
        case Gaussian: // gaussian
        density = exp(-r2/2.0/parA/parA);
        break;
        case DoubleGaussian: //double-gaussian
        density =  (1.0-parC)*exp(-r2/parA/parA)/parA/parA/parA;
        density += parC*exp(-r2/parB/parB)/parB/parB/parB;
        break;
      }
    profile = r2*density;
    rDensity->SetBinContent(iR+1,density);
    rDensity->SetBinError(iR+1,0.0);
    rProfile->SetBinContent(iR+1, profile);
    rProfile->SetBinError(iR+1,0.0);
    r += dr;
    }
}

void NucleusGenerator::generate(Nucleus * nucleus)
{
  double r, cosTheta, phi;
  int nProtons  = nucleus->nProtons;
  int nNucleons = nucleus->nNucleons;
  double xMean=0.0;
  double yMean=0.0;
  double zMean=0.0;

  int iNucleon = 0;
  while (iNucleon < nNucleons)
    {
    //cout << "iNucleon : " << iNucleon << endl;

    Nucleon * nucleon = nucleus->getNucleon(iNucleon);
    generate(r, cosTheta, phi);
    nucleon->setRCosThetaPhi(r, cosTheta, phi);

    // check if this
    int sanityCheck = 0;
    if (useNucleonExclusion)
      {
      bool reject = false;
      for (int jNucleon=0; jNucleon<iNucleon; jNucleon++)
        {
        //cout << "jNucleon : " << jNucleon << endl;

        Nucleon * otherNucleon = nucleus->getNucleon(jNucleon);
        if (nucleon->distanceXYZSq(otherNucleon) < exclusionRadiusSq)
          {
          reject = true;
          break;
          }
        }
      if (reject)
        {
        sanityCheck++;
        if (sanityCheck>200)
          {
          cout << "Nucleon rejected > 200 times" << endl;
          return;
          }
        //cout << "reject  sanityCheck:" << sanityCheck << endl;
        continue;
        }
      }

    xMean += nucleon->x;
    yMean += nucleon->y;
    zMean += nucleon->z;
    if (iNucleon<nProtons)
      nucleon->setNucleonType(Nucleon::Proton);
    else
      nucleon->setNucleonType(Nucleon::Neutron);
    iNucleon++;
    }

  xMean /= double(nNucleons);
  yMean /= double(nNucleons);
  zMean /= double(nNucleons);
  nucleus->x = xMean;
  nucleus->y = yMean;
  nucleus->z = zMean;

  // center of mass, recenter
  if (useRecentering)
    {
    for (int iNucleon=0; iNucleon<nNucleons; iNucleon++)
      {
      Nucleon * nucleon = nucleus->getNucleon(iNucleon);
      nucleon->shift( -xMean, -yMean, -zMean);
      }
    }

  //cout << "Nucleus complted" << endl;
}

void NucleusGenerator::generate(double & r, double & cosTheta, double & phi)
{
  cosTheta = -1 + 2.0*gRandom->Rndm();
  phi      = 2.0*3.1415927*gRandom->Rndm();
  r        = rProfile->GetRandom();
  rProfileGen->Fill(r);
}

void NucleusGenerator::saveHistorgrams()
{
  rDensity->Write();
  rProfile->Write();
  rProfileGen->Write();
}
