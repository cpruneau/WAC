// Author: Claude Pruneau   12/12/16

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class CanvasConfiguration
 \ingroup WAC

 Utility class used to define the parameters of a root canvas
 */


#include "Event.hpp"

//Factory<Particle> * particleFactory = 0;
//
//extern Factory<Particle> * getParticleFactory(int nParticleMax=5000)
//{
//  if (!particleFactory)
//    {
//    particleFactory = new Factory<Particle>();
//    particleFactory->initialize(nParticleMax);
//    }
//  return particleFactory;
//}


ClassImp(Event);

Event::Event()
 :
index(0),
eventNumber(0),
nProjectile(0),
nTarget(0),
nParticleTotal(0),
nParticles(0),
centrality(0),
multiplicity(0),
impactParameter(0),
other(0),
nFilters(0),
nFiltered(nullptr),
eFiltered(nullptr),
qFiltered(nullptr),
bFiltered(nullptr)
{
}

 //////////////////////////
 // DTOR
 //////////////////////////
 Event::~Event()
 {
 if (nFiltered) delete[] nFiltered;
 if (eFiltered) delete[] eFiltered;
 if (qFiltered) delete[] qFiltered;
 if (bFiltered) delete[] bFiltered;
 }

 ////////////////////////////////////////////////////
 // Call before restarting simulation
 ////////////////////////////////////////////////////
void Event::clear()
 {
 index      = 0;
 nParticles = 0;
 eventNumber = 0;
 nProjectile = 0;
 nTarget = 0;
 nParticleTotal = 0;
 nParticles = 0;
 centrality = 0;
 multiplicity = 0;
 impactParameter = 0;
 other = 0;
 Particle::getFactory()->reset();
 for (int iFilter=0;iFilter<nFilters;iFilter++)
   {
   nFiltered[iFilter] = 0.0;
   eFiltered[iFilter] = 0.0;
   qFiltered[iFilter] = 0.0;
   bFiltered[iFilter] = 0.0;
   }
 }

 ////////////////////////////////////////////////////
 // Call before generating new event
 ////////////////////////////////////////////////////
  void Event::reset()
 {
 index++;
 nParticles = 0;
 eventNumber = 0;
 nProjectile = 0;
 nTarget = 0;
 nParticleTotal = 0;
 nParticles = 0;
 centrality = 0;
 multiplicity = 0;
 impactParameter = 0;
 other = 0;
 Particle::getFactory()->reset();
 for (int iFilter=0;iFilter<nFilters;iFilter++)
   {
   nFiltered[iFilter] = 0.0;
   eFiltered[iFilter] = 0.0;
   qFiltered[iFilter] = 0.0;
   bFiltered[iFilter] = 0.0;
   }
 }


 //
 // Get the particle at the given index
 //
 Particle * Event::getParticleAt(int index)
 {
 if (index>=0 && index<nParticles)
   {
   //return getParticleFactory()->getObjectAt(index);
   return Particle::getFactory()->getObjectAt(index);
   }
 else
   return 0;
 }


 ///////////////////////////////////////////////////////
 // Print properties of this event at the given output
 ///////////////////////////////////////////////////////
void Event::printProperties(ostream & output)
 {
 output << "                 Event# " << index << endl;
 output << " nParticles(Accepted) : " << nParticles << endl;
 output << "          eventNumber : " << eventNumber << endl;
 output << "          nProjectile : " << nProjectile << endl;
 output << "              nTarget : " << nTarget << endl;
 output << "       nParticleTotal : " << nParticleTotal << endl;
 output << "           centrality : " << centrality << endl;
 output << "         multiplicity : " << multiplicity << endl;
 output << "      impactParameter : " << impactParameter << endl;
 output << "                other : " << other << endl;

 for (int iParticle=0; iParticle<nParticles; iParticle++)
   {
     //getParticleFactory()->getObjectAt(iParticle)->printProperties(output);
   Particle::getFactory()->getObjectAt(iParticle)->printProperties(output);
   }
 }

Event * Event::event = 0;
Event * Event::getEvent()
{
  if (!event)
    {
    event = new Event();
    }
  return event;
}

