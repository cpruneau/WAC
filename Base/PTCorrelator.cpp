/**
 \class Task
 \ingroup WAC
 Class defining Transverse Momentum Correlation Analyzer 
 */

#include <chrono>
#include "PTCorrelator.hpp"
#include "HeavyIonConfiguration.hpp"


ClassImp(PTCorrelator);


//////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////
PTCorrelator::PTCorrelator(const TString &  name,
                           NuDynConfiguration * configuration,
                           Event * event,
                           EventFilter * ef,
                           ParticleFilter ** pf,
                           LogLevel requiredLevel)
:
Task(name,configuration,event,requiredLevel),
histos(NULL),
eventFilter(ef),
particleFilters(pf),
maxEvents(0),maxOrder(0),
correlatorIndex(0)
{
  if (reportStart("PTCorrelator",getTaskName(),"CTOR"))
    ;
  //NuDynConfiguration * hc = (NuDynConfiguration *) getTaskConfiguration();

	//maxEvents      = hc->totEvents;
	//maxOrder      = hc->maxOrder;
	partNames     = new TString[maxOrder];
	pT            = new double * [maxEvents];
	acceptances   = new bool **[maxEvents];
	multiplicity  = new double [maxEvents];
	centrality    = new double [maxEvents];
	counts        = new int * [maxEvents];
	numParticles  = new int [maxEvents];

	if (!eventFilter)
	{
		if (reportError("PTCorrelator",getTaskName(),"CTOR")) cout << "eventFilter is a null pointer." << endl;
		postTaskError();
		return;
	}
	for(int i = 0; i < maxOrder; i++)
	{
		if (!particleFilters[i])
		{
			if (reportError("PTCorrelator",getTaskName(),"CTOR")) cout << "PTCorrelator::CTOR(...) particleFilter" << i +1 << " is a null pointer." << endl;
			postTaskError();
			return;
		}
	}

	TString newName = getTaskName();
	newName += "_";
	newName += eventFilter->getName();
	setTaskName(newName);
	for(int i = 0; i < maxOrder; i++)
	{
		partNames[i] = particleFilters[i]->getName();
	}
  if (reportEnd("PTCorrelator",getTaskName(),"CTOR"))
    ;
}


//////////////////////////////////////////////////////////////
// DTOR needs to be implemented
//////////////////////////////////////////////////////////////
PTCorrelator::~PTCorrelator()
{
  if (reportStart("PTCorrelator",getTaskName(),"DTOR"))
    ;
	if (histos != NULL) delete histos;
  if (reportEnd("PTCorrelator",getTaskName(),"DTOR"))
    ;
}


void PTCorrelator::createHistograms()
{
  if (reportStart("PTCorrelator",getTaskName(),"createHistograms()"))
    ;
  HeavyIonConfiguration * ac = (HeavyIonConfiguration *) getTaskConfiguration();
	//LogLevel debugLevel = getReportLevel();
	TString histoName;
	histoName = partNames[0];
	for(int i = 1; i < maxOrder; i++)
	{
		histoName += partNames[i];
	}
	histos = new PTHistos(histoName,ac,MessageLogger::Error, maxOrder);
  histos->createHistograms();
	for(int i = 0; i < maxEvents; i++)
	{
		counts[i] = new int[histos->size];
	}
  if (reportEnd("PTCorrelator",getTaskName(),"createHistograms()"))
    ;
}


//////////////////////////////////////////////////////////////
// load histograms from given files needs to be fixed
//////////////////////////////////////////////////////////////
void PTCorrelator::loadHistograms(TFile * inputFile)
{
  if (reportStart("PTCorrelator",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
  NuDynConfiguration * ac = (NuDynConfiguration *) getTaskConfiguration();

	TString histoName;
	histoName = partNames[0];
	for(int i = 1; i < maxOrder; i++)
	{
		histoName += partNames[i];
	}

	histos = new PTHistos(histoName,ac,MessageLogger::Error, maxOrder);
  histos->loadHistograms(inputFile);
	for(int i = 0; i < maxEvents; i++)
	{
		counts[i] = new int[histos->size];
	}
  if (reportStart("PTCorrelator",getTaskName(),"loadHistograms(TFile * inputFile)"))
    ;
}


//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void PTCorrelator::saveHistograms(TFile * outputFile)
{
  if (reportStart("PTCorrelator",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
  if (!outputFile)
	{
		if (reportError("PTCorrelator",getTaskName(),"saveHistograms(TFile * outputFile)")) cout << "outputFile is a null  pointer." << endl;
		postTaskError();
		return;
	}
	outputFile->cd();
  if (reportEnd("PTCorrelator",getTaskName(),"saveHistograms(TFile * outputFile)"))
    ;
}


void PTCorrelator::execute()
{
	auto start = chrono::high_resolution_clock::now();
  incrementEventProcessed();
  if (!eventFilter->accept(*event)) return;
  incrementEventAccepted(); // count events used to fill histograms and for scaling at the end...

	HeavyIonConfiguration * ac = (HeavyIonConfiguration *) getTaskConfiguration();
	if (!ac)
	{
		if (reportError()) cout << "PTCorrelator::execute(...) analysisConfiguration null pointer" << endl;
		postTaskError();
		return;
	}
	correlatorIndex = 0;
	//int count = 0;
	//int *filters = new int [maxOrder];
	//int *particles = new int [maxOrder];

	histos->fillEventHistos( event->multiplicity, event->centrality, 1.0); 
	storeEventInfo();
	correlatorIndex = 0;
	auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	if (reportDebug()) cout << "Analyzed event "<< getNEventProcessed() << " Time taken by event: " << duration.count() << " microseconds" << " Multiplicity: " << event->multiplicity << endl;
	if(getNEventProcessed() == maxEvents)
	{
		cout << "Max number of collisions: " << ac->nCollisionsMax << endl;
		histos->fillDerivedHistos(acceptances, multiplicity, centrality, numParticles, pT);
	}
}



//////////////////////////////////////////////////////////////
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
//////////////////////////////////////////////////////////////
void PTCorrelator::scaleHistograms(double factor)
{
  if (reportInfo("PTCorrelator",getTaskName(),"scaleHistograms(double factor)"))
    cout << "Scaling histograms by a factor: " << factor << endl;
  histos->scale(factor);
  if (reportEnd("PTCorrelator",getTaskName(),"scaleHistograms(double factor)"))
    ;
}



//////////////////////////////////////////////////////////
//store the transverse momentum of all the particles
//store whether the particles are accepted by each of the filters
//store the multiplicity, centrality, and number of particles  of each event
////////////////////////////////////////////////////////
void PTCorrelator::storeEventInfo()
{
  if (reportStart("PTCorrelator",getTaskName(),"storeEventInfo()"))
    ;	//the transverse momentum
	pT[getNEventProcessed()] = new double [event->nParticles];
	for(int iParticle = 0; iParticle < event->nParticles; iParticle++)
	{
		Particle & particle = * event->getParticleAt(iParticle);
		pT[getNEventProcessed()][iParticle] = particle.pt;
	}

	//The acceptances
	acceptances[getNEventProcessed()] = new bool *[maxOrder];
	for(int i = 0; i < maxOrder; i++)
	{
		acceptances[getNEventProcessed()][i] = new bool [event->nParticles];
		for(int iParticle = 0; iParticle < event->nParticles; iParticle++)
		{
			Particle & particle = * event->getParticleAt(iParticle);
			acceptances[getNEventProcessed()][i][iParticle] = particleFilters[i]->accept(particle);
		}
	}

	numParticles[getNEventProcessed()] = event->nParticles;
	multiplicity[getNEventProcessed()] = event->multiplicity;
	centrality[getNEventProcessed()] = event->centrality;
  if (reportEnd("PTCorrelator",getTaskName(),"storeEventInfo()"))
    ;
}


void PTCorrelator::calculateDerivedHistograms()
{

}


void PTCorrelator::resetHistograms()
{
  histos->reset();
}
