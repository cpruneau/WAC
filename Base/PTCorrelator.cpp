/**
 \class Task
 \ingroup WAC
 Class defining Transverse Momentum Correlation Analyzer 
 */

#include <chrono>
#include "PTCorrelator.hpp"
#include "AnalysisConfiguration.hpp"


ClassImp(PTCorrelator);


//////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////
PTCorrelator::PTCorrelator(const TString &  name,
	TaskConfiguration * configuration,
	Event * event,
	EventFilter * ef,
	ParticleFilter ** pf,
	int ord,
	int events,
	int * maxCollisions)
:
Task(name,configuration,event),
histos(NULL),
eventFilter(ef),
particleFilters(pf),
partNames(new TString[ord]),
maxOrder(ord),
eventAveragept(0),
correlatorIndex(0),
maxEvents(events),
pT(new double * [events]),
acceptances(new bool **[events]),
multiplicity(new double [events]),
centrality(new double [events]),
S(new double * [events]),
avgpT(new double [ord]()),
counts(new int * [events]),
nCollisionsMax(maxCollisions)
{
	setReportLevel(MessageLogger::Debug);
	if (reportDebug())  cout << "PTCorrelator::CTOR(...) Started." << endl;

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
// DTOR
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
	AnalysisConfiguration * ac = (AnalysisConfiguration *) getTaskConfiguration();
	LogLevel debugLevel = getReportLevel();

	TString histoName;
	histoName = partNames[0];
	for(int i = 1; i < maxOrder; i++)
	{
		histoName += partNames[i];
	}

	histos = new PTHistos(histoName,ac,MessageLogger::Error, maxOrder);
	


	avgCounts = new double[histos->size]();

	for(int i = 0; i < maxEvents; i++)
	{
		S[i] = new double[histos->size];
		counts[i] = new int[histos->size];
	}

	if (reportDebug())  cout << "PTCorrelator::createHistograms(...) completed"<< endl;


}


//////////////////////////////////////////////////////////////
// load histograms from given files
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

	avgCounts = new double[histos->size]();

	for(int i = 0; i < maxEvents; i++)
	{
		S[i] = new double[histos->size];
		counts[i] = new int[histos->size];
	}
	if (reportDebug())  cout << "PTCorrelator::loadHistograms(...) Completed." << endl;
	avgCounts = new double[histos->size];
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
	fillSValues(maxOrder - 1, 0, filters, count, particles);
	storeEventInfo();

	correlatorIndex = 0;
	eventsProcessed++;


	auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	if (reportDebug()) cout << "Analyzed event "<< eventsProcessed << " Time taken by event: " << duration.count() << " microseconds" << " Multiplicity: " << event->multiplicity << endl; 

	if(eventsProcessed == maxEvents)
	{
		histos->resetHistoRanges(*nCollisionsMax);
		histos->fillDerivedHistos(acceptances, multiplicity, centrality, avgCounts, avgpT, S, counts, maxEvents);
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
	int * counter = new int[maxOrder]();
	double * pts = new double[maxOrder]();
	eventAveragept = 0;
	for (int iParticle=0; iParticle<event->nParticles; iParticle++)
	{
		Particle & particle = * event->getParticleAt(iParticle);
		if (reportDebug())  particle.printProperties(cout);
		eventAveragept += particle.pt;

		for(int i = 0; i <maxOrder; i++)
		{
			if(particleFilters[i]->accept(particle))
			{
				counter[i]++;
				pts[i] += particle.pt;
			}
		}

	}

	for(int i = 0; i<maxOrder; i++)
	{
		if(counter[i] != 0)	avgpT[i] = (avgpT[i] * eventsProcessed + pts[i] / counter[i])/(eventsProcessed + 1);
		else avgpT[i] = (avgpT[i] * eventsProcessed)/(eventsProcessed + 1);
	}

	eventAveragept /= event->nParticles;
}



//////////////////////////////////////////////////////////
//store the transverse momentum of all the particles
//store whether the particles are accepted by each of the filters
//store the multiplicity of each event
////////////////////////////////////////////////////////
void PTCorrelator::storeEventInfo()
{
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

	//multiplicity
	multiplicity[eventsProcessed] = event->multiplicity;
	centrality[eventsProcessed] = event->centrality;
}









/////////////////////////////////////////////////////////////
// Calculate and fill the pT deviation correlators for orders from 1 to maxOrder for all combinations
// Calculate the average counts of particles, pairs, triples ...
////////////////////////////////////////////////////////////
void PTCorrelator::fillSValues(int depth, int filterIndex, int * filters, int & count, int *particles)
{
	
	double correlators = 0;
	for(int i = filterIndex; i < maxOrder; i++)
	{
		filters[maxOrder - depth- 1] = i;
		count = 0;
		correlators = calculateS(filters, maxOrder - depth, 0, count, particles);

		S[eventsProcessed][correlatorIndex] = correlators;

		counts[eventsProcessed][correlatorIndex] = count;
		avgCounts[correlatorIndex] = (avgCounts[correlatorIndex] * eventsProcessed + count )/(eventsProcessed + 1);

		count = 0;
		correlatorIndex++;
		
		if(depth != 0)
		{
			fillSValues(depth - 1, i , filters, count, particles);
		}
	}
}


/////////////////////////////////////////////////////////////
// Calculate the pT deviation correlators(S) for orders from 1 to maxOrder for all combinations
// checked for correctness
////////////////////////////////////////////////////////////
double PTCorrelator::calculateS(int * filters, int order, int curFilterIndex, int & count, int * particles)
{
	double sum = 0;
	for(int iParticle=0; iParticle<event->nParticles; iParticle++)
	{
		bool accept = true;
		for(int i = 0; i <curFilterIndex; i++ )
		{
			accept = accept && (iParticle != particles[i]);
		}
		if(accept)
		{
			Particle & particle = * event->getParticleAt(iParticle);
			double tempSum = 1;
			if(particleFilters[filters[curFilterIndex]]->accept(particle))
			{
				double deviation = (particle.pt - eventAveragept);
				if(order != 1)
				{
					particles[curFilterIndex] = iParticle;
					tempSum = calculateS(filters, order - 1, curFilterIndex + 1, count, particles);
				}
				else
				{
					count++;
				}
				sum += tempSum * deviation;
			}
		}
	}


	return sum;

}