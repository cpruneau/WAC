/**
 \class Task
 \ingroup WAC
 Class defining Transverse Momentum Correlation Analyzer 
 */

#include <chrono>
#include "PTCorrelator.hpp"
#include "AnalysisConfiguration.hpp"
#include "HeavyIonConfiguration.hpp"


ClassImp(PTCorrelator);


//////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////
PTCorrelator::PTCorrelator(const TString &  name,
	TaskConfiguration * configuration,
	Event * event,
	EventFilter * ef,
	ParticleFilter ** pf)
:
Task(name,configuration,event),
histos(NULL),
eventFilter(ef),
particleFilters(pf),
maxOrder(0),
eventAveragept(0),
correlatorIndex(0),
maxEvents(0)
{
	if (reportDebug())  cout << "PTCorrelator::CTOR(...) Started." << endl;
	AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
	HeavyIonConfiguration * hc = (HeavyIonConfiguration *) ac;


	maxEvents      = hc->totEvents;
	maxOrder      = hc->maxOrder;
	partNames     = new TString[maxOrder];
	pT            = new double * [maxEvents];
	acceptances   = new bool **[maxEvents];
	multiplicity  = new double [maxEvents];
	centrality    = new double [maxEvents];
	counts        = new int * [maxEvents];
	numParticles  = new int [maxEvents];

	if (!eventFilter)
	{
		if (reportError()) cout << "PTCorrelator::CTOR(...) eventFilter is null pointer." << endl;
		postTaskError();
		return;
	}
	for(int i = 0; i < maxOrder; i++)
	{
		if (!particleFilters[i])
		{
			if (reportError()) cout << "PTCorrelator::CTOR(...) particleFilter" << i +1 << " is null pointer." << endl;
			postTaskError();
			return;
		}
	}

	TString newName = getName();
	newName += "_";
	newName += eventFilter->getName();
	setName(newName);
	for(int i = 0; i < maxOrder; i++)
	{
		partNames[i] = particleFilters[i]->getName();
	}



}


//////////////////////////////////////////////////////////////
// DTOR needs to be implemented
//////////////////////////////////////////////////////////////
PTCorrelator::~PTCorrelator()
{
	if (reportDebug())  cout << "PTCorrelator::DTOR(...) Started" << endl;
	if (histos != NULL) delete histos;
	if (reportDebug())  cout << "PTCorrelator::DTOR(...) Completed" << endl;
}


void PTCorrelator::createHistograms()
{
	if (reportDebug())  cout << "PTCorrelator::createHistograms(...) started"<< endl;
	HeavyIonConfiguration * ac = (HeavyIonConfiguration *) getTaskConfiguration();
	LogLevel debugLevel = getReportLevel();

	TString histoName;
	histoName = partNames[0];
	for(int i = 1; i < maxOrder; i++)
	{
		histoName += partNames[i];
	}

	histos = new PTHistos(histoName,ac,MessageLogger::Error, maxOrder);
	


	for(int i = 0; i < maxEvents; i++)
	{
		counts[i] = new int[histos->size];
	}

	if (reportDebug())  cout << "PTCorrelator::createHistograms(...) completed"<< endl;


}


//////////////////////////////////////////////////////////////
// load histograms from given files needs to be fixed
//////////////////////////////////////////////////////////////
void PTCorrelator::loadHistograms(TFile * inputFile)
{
	if (reportDebug())  cout << "PTCorrelator::loadHistograms(...) Starting." << endl;
  /* first load the number of events as from the  cumulated parameter */
	TParameter<Long64_t> *par = (TParameter<Long64_t> *) inputFile->Get("NoOfEvents");
	eventsProcessed = par->GetVal();
	delete par;
	AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
	LogLevel debugLevel = getReportLevel();

	TString histoName;
	histoName = partNames[0];
	for(int i = 1; i < maxOrder; i++)
	{
		histoName += partNames[i];
	}

	histos = new PTHistos(inputFile,histoName,ac,MessageLogger::Error, maxOrder);

	for(int i = 0; i < maxEvents; i++)
	{
		counts[i] = new int[histos->size];
	}
	if (reportDebug())  cout << "PTCorrelator::loadHistograms(...) Completed." << endl;
}


//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void PTCorrelator::saveHistograms(TFile * outputFile)
{
	if (reportDebug()) cout << "PTCorrelator::saveHistograms(...) Saving Event histograms to file." << endl;
	if (!outputFile)
	{
		if (reportError()) cout << "PTCorrelator::saveHistograms(...) outputFile is a null  pointer." << endl;
		postTaskError();
		return;
	}
	outputFile->cd();

  /* first save the number of events as a cumulated parameter */
	TParameter<Long64_t>("NoOfEvents",eventsProcessed,'+').Write();
	AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
	histos->saveHistograms(outputFile);
	if (reportDebug()) cout << "PTCorrelator::saveHistograms(...) Completed." << endl;
}


void PTCorrelator::execute()
{

	auto start = chrono::high_resolution_clock::now(); 
	if (event != NULL)
	{
		if (reportDebug()) cout << "PTCorrelator::execute(...) analyzing " << event->nParticles << " particles" << endl;
	}
	else
	{
		if (reportError()) cout << "PTCorrelator::execute(...) event pointer is NULL. Abort." << endl;
		postTaskError();
		return;
	}

  //if (reportDebug()) cout <<"PTCorrelator::analyze(...) check if event is acceptable" << endl;
	if (!eventFilter->accept(*event)) return;
  //if (reportDebug()) cout <<"PTCorrelator::analyze(...) acceptable event" << endl;

	AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
	if (!ac)
	{
		if (reportError()) cout << "PTCorrelator::execute(...) analysisConfiguration null pointer" << endl;
		postTaskError();
		return;
	}

	calculateAverage();

	correlatorIndex = 0;
	int count = 0;
	int *filters = new int [maxOrder];
	int *particles = new int [maxOrder];

	histos->fillEventHistos( event->multiplicity, event->centrality, 1.0); 
	storeEventInfo();

	correlatorIndex = 0;
	eventsProcessed++;


	auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	if (reportDebug()) cout << "Analyzed event "<< eventsProcessed << " Time taken by event: " << duration.count() << " microseconds" << " Multiplicity: " << event->multiplicity << endl; 

	if(eventsProcessed == maxEvents)
	{
		histos->fillDerivedHistos(acceptances, multiplicity, centrality, numParticles, pT);
	}

	if (reportDebug()) cout << "PTCorrelator::execute() Completed" << endl;
}



//////////////////////////////////////////////////////////////
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
//////////////////////////////////////////////////////////////
void PTCorrelator::scaleHistograms(double factor)
{
	if (reportDebug())  cout << "PTCorrelator::scaleHistograms(..) Scale all primary histograms by " << factor << endl;
	histos->scale(factor);
	if (reportDebug())  cout << "PTCorrelator::scale(..) Completed"  << endl;
}



/////////////////////////////////////////////////////////////
// Calculate average pT for each type of particle
////////////////////////////////////////////////////////////
void PTCorrelator::calculateAverage()
{
	if (reportDebug())  cout << "PTCorrelator::calculateAverage(...) Starting." << endl;
	eventAveragept = 0;
	for (int iParticle=0; iParticle<event->nParticles; iParticle++)
	{
		Particle & particle = * event->getParticleAt(iParticle);
		//if (reportDebug())  particle.printProperties(cout);
		eventAveragept += particle.pt;

	}

	eventAveragept /= event->nParticles;
	if (reportDebug())  cout << "PTCorrelator::calculateAverage(...) Completed." << endl;
}



//////////////////////////////////////////////////////////
//store the transverse momentum of all the particles
//store whether the particles are accepted by each of the filters
//store the multiplicity, centrality, and number of particles  of each event
////////////////////////////////////////////////////////
void PTCorrelator::storeEventInfo()
{
	if (reportDebug())  cout << "PTCorrelator::storeEventInfo(...) Starting." << endl;
	//the transverse momentum
	pT[eventsProcessed] = new double [event->nParticles];
	for(int iParticle = 0; iParticle < event->nParticles; iParticle++)
	{
		Particle & particle = * event->getParticleAt(iParticle);
		pT[eventsProcessed][iParticle] = particle.pt;
	}

	//The acceptances
	acceptances[eventsProcessed] = new bool *[maxOrder];
	for(int i = 0; i < maxOrder; i++)
	{
		acceptances[eventsProcessed][i] = new bool [event->nParticles];
		for(int iParticle = 0; iParticle < event->nParticles; iParticle++)
		{
			Particle & particle = * event->getParticleAt(iParticle);
			acceptances[eventsProcessed][i][iParticle] = particleFilters[i]->accept(particle);
		}
	}

	numParticles[eventsProcessed] = event->nParticles;
	multiplicity[eventsProcessed] = event->multiplicity;
	centrality[eventsProcessed] = event->centrality;
	if (reportDebug())  cout << "PTCorrelator::storeEventInfo(...) Completed." << endl;
}