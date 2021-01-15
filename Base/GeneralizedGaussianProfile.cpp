// Author: Victor Gonzalez   02/10/19

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau. Victor Gonzalez.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#include "GeneralizedGaussianProfile.hpp"
ClassImp(GeneralizedGaussianProfile);
#include "Event.hpp"

GeneralizedGaussianProfile::GeneralizedGaussianProfile(const TString & name,
                                                       TaskConfiguration * configuration,
                                                       Event * selectedEvent,
                                                       LogLevel selectedLevel)
: Task(name,configuration,selectedEvent,selectedLevel),
profile(nullptr)
{
 // no ops
}

/////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
/////////////////////////////////////////////////////////////////////////////////////////////
GeneralizedGaussianProfile::~GeneralizedGaussianProfile()
{
  if (reportStart("GeneralizedGaussianProfile",getTaskName(),"DTOR"))
    ;
  if (profile != nullptr) delete profile;
  profile = nullptr;
  if (reportEnd("GeneralizedGaussianProfile",getTaskName(),"DTOR"))
    ;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// initialize this analyzer
/////////////////////////////////////////////////////////////////////////////////////////////
void GeneralizedGaussianProfile::initialize()
{
  if (reportStart("GeneralizedGaussianProfile",getTaskName(),"initialize()"))
    ;
  GeneralizedGaussianConfiguration *conf = (GeneralizedGaussianConfiguration *) getTaskConfiguration();

  profile = new TF2("2DGenGauss","[0]*[1]*[2]/4.0/[3]/[4]/TMath::Gamma(1.0/[1])/TMath::Gamma(1.0/[2])*"
                    "TMath::Exp(-1.0*(TMath::Power(TMath::Abs(x/[3]),[1])+TMath::Power(TMath::Abs(y/[4]),[2])))");
  profile->SetParameters(conf->amplitude,conf->gammaeta,conf->gammaphi,conf->omegaeta,conf->omegaphi);
  profile->SetRange(-4.0,-TMath::PiOver2(),4.0,TMath::PiOver2());
  if (reportEnd("GeneralizedGaussianProfile",getTaskName(),"initialize()"))
    ;

}


/////////////////////////////////////////////////////////////////////////////////////////////
// process the given event
/////////////////////////////////////////////////////////////////////////////////////////////
void GeneralizedGaussianProfile::execute()
{
  if (reportStart("GeneralizedGaussianProfile",getTaskName(),"execute()"))
    ;
  if (event != nullptr)
    {
    if (reportDebug("GeneralizedGaussianProfile",getTaskName(),"execute()")) cout << "Processing " << event->nParticles << " particles" << endl;
    for (int itrack = 0; itrack < event->nParticles; itrack++) {
      Particle *particle1 = event->getParticleAt(itrack);
      for (int jtrack = itrack+1; jtrack < event->nParticles; jtrack++) {
        Particle *particle2 = event->getParticleAt(jtrack);
        double deltaeta = particle1->eta - particle2->eta;
        double deltaphi = TVector2::Phi_mpi_pi(particle1->phi - particle2->phi);
        double weight = profile->Eval(deltaeta,deltaphi) + 1.0;
        particle1->pt *= weight;
        particle2->pt *= weight;
      }
    }
    }
  else
    {
    postTaskFatal();
    if (reportFatal("GeneralizedGaussianProfile",getTaskName(),"execute()")) cout << "event pointer is null." << endl;
    }
  if (reportEnd("GeneralizedGaussianProfile",getTaskName(),"execute()"))
    ;
}


