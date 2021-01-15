/**
 \class Task
 \ingroup WAC
 Class defining Transverse Momentum Correlation Analyzer 
 */

#include <chrono>
#include "PTCorrelator.hpp"

ClassImp(PTCorrelator);

//////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////
PTCorrelator::PTCorrelator(const TString &name,
						   TaskConfiguration *configuration,
						   Event *event,
						   EventFilter *ef,
						   ParticleFilter **pf,
						   LogLevel selectedLevel)
	: Task(name, configuration, event, selectedLevel),
	  histos(NULL),
	  eventFilter(ef),
	  particleFilters(pf),
	  maxOrder(0),
	  correlatorIndex(0),
	  maxEvents(0)
{
	if (reportDebug())
		cout << "PTCorrelator::CTOR(...) Started." << endl;
	TransverseMomentumConfiguration *hc = (TransverseMomentumConfiguration *)getTaskConfiguration();

	maxEvents = hc->totEvents;
	maxOrder = hc->maxOrder;
	partNames = new TString[hc->numTypes];
	multiplicity = new double[maxEvents];
	centrality = new double[maxEvents];
	transverseMomentumMoments = new double **[maxEvents];
	yields = new double *[maxEvents];
	numParticles = new double[maxEvents];

	if (!eventFilter)
	{
		if (reportError())
			cout << "PTCorrelator::CTOR(...) eventFilter is null pointer." << endl;
		postTaskError();
		return;
	}
	for (int i = 0; i < hc->numTypes; i++)
	{
		if (!particleFilters[i])
		{
			if (reportError())
				cout << "PTCorrelator::CTOR(...) particleFilter" << i + 1 << " is null pointer." << endl;
			postTaskError();
			return;
		}
	}

	TString newName = getTaskName();
	newName += "_";
	newName += eventFilter->getName();
	setTaskName(newName);
	for (int i = 0; i < hc->numTypes; i++)
	{
		partNames[i] = particleFilters[i]->getName();
	}
}

//////////////////////////////////////////////////////////////
// DTOR needs to be implemented
//////////////////////////////////////////////////////////////
PTCorrelator::~PTCorrelator()
{
	if (reportDebug())
		cout << "PTCorrelator::DTOR(...) Started" << endl;
	TransverseMomentumConfiguration *hc = (TransverseMomentumConfiguration *)getTaskConfiguration();
	for (int i = 0; i < nEventProcessed; i++)
	{
		for (int j = 0; j < hc->numTypes; j++)
		{
			delete[] transverseMomentumMoments[i][j];
		}
		delete[] transverseMomentumMoments[i];
		delete[] yields[i];
	}
	delete[] transverseMomentumMoments;
	delete[] yields;
	delete[] multiplicity;
	delete[] centrality;
	delete[] numParticles;
	if (histos != NULL)
		delete histos;
	if (reportDebug())
		cout << "PTCorrelator::DTOR(...) Completed" << endl;
}

void PTCorrelator::createHistograms()
{
	if (reportDebug())
		cout << "PTCorrelator::createHistograms(...) started" << endl;
	TransverseMomentumConfiguration *ac = (TransverseMomentumConfiguration *)getTaskConfiguration();
	//LogLevel debugLevel = getReportLevel();

	TString histoName;
	histoName = partNames[0];
	for (int i = 1; i < ac->numTypes; i++)
	{
		histoName += partNames[i];
	}

	histos = new PTHistos(histoName, ac, getReportLevel(), maxOrder);

	if (reportDebug())
		cout << "PTCorrelator::createHistograms(...) completed" << endl;
}

//////////////////////////////////////////////////////////////
// load histograms from given files needs to be fixed
//////////////////////////////////////////////////////////////
void PTCorrelator::loadHistograms(TFile *inputFile)
{
	if (reportDebug())
		cout << "PTCorrelator::loadHistograms(...) Starting." << endl;
	/* first load the number of events as from the  cumulated parameter */
	TParameter<Long64_t> *par = (TParameter<Long64_t> *)inputFile->Get("NoOfEvents");
	nEventProcessed = par->GetVal();
	delete par;
	TransverseMomentumConfiguration *ac = (TransverseMomentumConfiguration *)getTaskConfiguration();
	//LogLevel debugLevel = getReportLevel();

	TString histoName;
	histoName = partNames[0];
	for (int i = 1; i < ac->numTypes; i++)
	{
		histoName += partNames[i];
	}

	histos = new PTHistos(inputFile, histoName, ac, getReportLevel(), maxOrder);

	if (reportDebug())
		cout << "PTCorrelator::loadHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
void PTCorrelator::saveHistograms(TFile *outputFile)
{
	if (reportDebug())
		cout << "PTCorrelator::saveHistograms(...) Saving Event histograms to file." << endl;
	if (!outputFile)
	{
		if (reportError())
			cout << "PTCorrelator::saveHistograms(...) outputFile is a null  pointer." << endl;
		postTaskError();
		return;
	}
	outputFile->cd();
	/* first save the number of events as a cumulated parameter */
	TParameter<Long64_t>("NoOfEvents", nEventProcessed, '+').Write();
	TransverseMomentumConfiguration *ac = (TransverseMomentumConfiguration *)getTaskConfiguration();
	histos->saveHistograms(outputFile);
	if (reportDebug())
		cout << "PTCorrelator::saveHistograms(...) Completed." << endl;
}

void PTCorrelator::execute()
{

	auto start = chrono::high_resolution_clock::now();
	if (event != NULL)
	{
		if (reportDebug())
			cout << "PTCorrelator::execute(...) analyzing " << event->nParticles << " particles" << endl;
	}
	else
	{
		if (reportError())
			cout << "PTCorrelator::execute(...) event pointer is NULL. Abort." << endl;
		postTaskError();
		return;
	}

	//if (reportDebug()) cout <<"PTCorrelator::analyze(...) check if event is acceptable" << endl;
	if (!eventFilter->accept(*event))
		return;
	//if (reportDebug()) cout <<"PTCorrelator::analyze(...) acceptable event" << endl;

	TransverseMomentumConfiguration *ac = (TransverseMomentumConfiguration *)getTaskConfiguration();
	if (!ac)
	{
		if (reportError())
			cout << "PTCorrelator::execute(...) TransverseMomentumConfiguration null pointer" << endl;
		postTaskError();
		return;
	}

	correlatorIndex = 0;
	histos->fillEventHistos(event->multiplicity, event->centrality, 1.0);
	fillTransverseMomentumValues();
	transverseMomentumMoments[nEventProcessed] = new double *[ac->numTypes];
	calculateTransverseMomentumMoments();
	yields[nEventProcessed] = new double[histos->size];
	fillYieldValues();
	storeEventInfo();

	correlatorIndex = 0;
	nEventProcessed++;
	nEventAccepted++;

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	if (reportDebug())
		cout << "Analyzed event " << nEventProcessed << " Time taken by event: " << duration.count() << " microseconds"
			 << " Multiplicity: " << event->multiplicity << endl;

	if (reportDebug())
		cout << "PTCorrelator::execute() Completed" << endl;
}

//////////////////////////////////////////////////////////////
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
//////////////////////////////////////////////////////////////
void PTCorrelator::scaleHistograms(double factor)
{
	if (reportDebug())
		cout << "PTCorrelator::scaleHistograms(..) Scale all primary histograms by " << factor << endl;
	histos->scale(factor);
	if (reportDebug())
		cout << "PTCorrelator::scale(..) Completed" << endl;
}

//////////////////////////////////////////////////////////
//store the transverse momentum of all the particles
//store whether the particles are accepted by each of the filters
//store the multiplicity, centrality, and number of particles  of each event
////////////////////////////////////////////////////////
void PTCorrelator::storeEventInfo()
{
	if (reportDebug())
		cout << "PTCorrelator::storeEventInfo(...) Starting." << endl;
	multiplicity[nEventProcessed] = event->multiplicity;
	centrality[nEventProcessed] = event->centrality;
	numParticles[nEventProcessed] = event->nParticles;
	if (reportDebug())
		cout << "PTCorrelator::storeEventInfo(...) Completed." << endl;
}

void PTCorrelator::fillTransverseMomentumValues()
{
	if (reportDebug())
		cout << "PTCorrelator::fillTransverseMomentumValues(...) Starting." << endl;
	TransverseMomentumConfiguration *hc = (TransverseMomentumConfiguration *)getTaskConfiguration();
	for (int iParticle = 0; iParticle < event->nParticles; iParticle++)
	{
		Particle &particle = *event->getParticleAt(iParticle);
		for (int i = 0; i < hc->numTypes; i++)
		{
			if (particleFilters[i]->accept(particle))
				histos->fillTransverseMomentumHistos(particle.pt, i, event->multiplicity, event->centrality, 1.0);
		}
	}
	if (reportDebug())
		cout << "PTCorrelator::fillTransverseMomentumValues(...) Completed." << endl;
}

void PTCorrelator::fillYieldValues()
{
	if (reportDebug())
		cout << "PTCorrelator::fillYieldValues(...) Starting." << endl;
	TransverseMomentumConfiguration *hc = (TransverseMomentumConfiguration *)getTaskConfiguration();
	double *n = new double[hc->numTypes]();
	int counter = 0;
	double *tempCounts = new double[histos->size]();
	for (int iFilter = 0; iFilter < hc->numTypes; iFilter++)
	{
		for (int iParticle = 0; iParticle < event->nParticles; iParticle++)
		{
			Particle &particle = *event->getParticleAt(iParticle);
			if (particleFilters[iFilter]->accept(particle))
				n[iFilter]++;
		}
	}

	for (int iFilter1 = 0; iFilter1 < hc->numTypes; iFilter1++)
	{
		if (n[iFilter1] > 0)
			tempCounts[counter] = n[iFilter1];
		counter++;
	}

	if (maxOrder > 1)
		for (int iFilter1 = 0; iFilter1 < hc->numTypes; iFilter1++)
		{
			for (int iFilter2 = iFilter1; iFilter2 < hc->numTypes; iFilter2++)
			{
				int same1 = iFilter2 == iFilter1 ? 1 : 0; //subtract out the counting te same particle twice
				if (n[iFilter1] > 0 && (n[iFilter2] - same1) > 0)
					tempCounts[counter] = n[iFilter1] * (n[iFilter2] - same1);
				counter++;
			}
		}

	if (maxOrder > 2)
		for (int iFilter1 = 0; iFilter1 < hc->numTypes; iFilter1++)
		{
			for (int iFilter2 = iFilter1; iFilter2 < hc->numTypes; iFilter2++)
			{
				int same12 = iFilter2 == iFilter1 ? 1 : 0; //subtract out the counting te same particle twice
				for (int iFilter3 = iFilter2; iFilter3 < hc->numTypes; iFilter3++)
				{
					int same13 = iFilter3 == iFilter1 ? 1 : 0; //subtract out the counting te same particle twice
					int same23 = iFilter2 == iFilter3 ? 1 : 0; //subtract out the counting te same particle twice
					if (n[iFilter1] > 0 && (n[iFilter2] - same12) > 0 && (n[iFilter3] - same13 - same23) > 0)
						tempCounts[counter] = n[iFilter1] * (n[iFilter2] - same12) * (n[iFilter3] - same13 - same23);
					counter++;
				}
			}
		}

	if (maxOrder > 3)
		for (int iFilter1 = 0; iFilter1 < hc->numTypes; iFilter1++)
		{
			for (int iFilter2 = iFilter1; iFilter2 < hc->numTypes; iFilter2++)
			{
				for (int iFilter3 = iFilter2; iFilter3 < hc->numTypes; iFilter3++)
				{
					int same12 = iFilter2 == iFilter1 ? 1 : 0;
					int same13 = iFilter3 == iFilter1 ? 1 : 0;
					int same23 = iFilter2 == iFilter3 ? 1 : 0;
					for (int iFilter4 = iFilter3; iFilter4 < hc->numTypes; iFilter4++)
					{
						int same14 = iFilter4 == iFilter1 ? 1 : 0;
						int same24 = iFilter2 == iFilter4 ? 1 : 0;
						int same34 = iFilter3 == iFilter4 ? 1 : 0;
						if (n[iFilter1] > 0 && (n[iFilter2] - same12) > 0 && (n[iFilter3] - same13 - same23) > 0 && (n[iFilter4] - same14 - same24 - same34) > 0)
							tempCounts[counter] = n[iFilter1] * (n[iFilter2] - same12) * (n[iFilter3] - same13 - same23) * (n[iFilter4] - same14 - same24 - same34);
						counter++;
					}
				}
			}
		}
	delete[] n;

	for (int iHisto = 0; iHisto < histos->size; iHisto++)
	{
		yields[nEventProcessed][iHisto] = tempCounts[histos->reorder[iHisto]];
	}
	delete[] tempCounts;
	if (reportDebug())
		cout << "PTCorrelator::fillYieldValues(...) Completed." << endl;
}

void PTCorrelator::calculateTransverseMomentumMoments()
{
	if (reportDebug())
		cout << "PTCorrelator::calculateTransverseMomentumMoments(...) Starting." << endl;
	TransverseMomentumConfiguration *hc = (TransverseMomentumConfiguration *)getTaskConfiguration();
	for (int iFilter = 0; iFilter < hc->numTypes; ++iFilter)
	{
		transverseMomentumMoments[nEventProcessed][iFilter] = new double[maxOrder]();
		for (int iParticle = 0; iParticle < event->nParticles; iParticle++)
		{
			Particle &particle = *event->getParticleAt(iParticle);
			if (particleFilters[iFilter]->accept(particle))
			{
				transverseMomentumMoments[nEventProcessed][iFilter][0] += particle.pt;
				if (maxOrder > 1)
					transverseMomentumMoments[nEventProcessed][iFilter][1] += particle.pt * particle.pt;
				if (maxOrder > 2)
					transverseMomentumMoments[nEventProcessed][iFilter][2] += particle.pt * particle.pt * particle.pt;
				if (maxOrder > 3)
					transverseMomentumMoments[nEventProcessed][iFilter][3] += particle.pt * particle.pt * particle.pt * particle.pt;
			}
		}
	}
	if (reportDebug())
		cout << "PTCorrelator::calculateTransverseMomentumMoments(...) Completed." << endl;
}

void PTCorrelator::reset()
{
	if (reportStart("Task", getTaskName(), "reset()"))
		;
	TransverseMomentumConfiguration *hc = (TransverseMomentumConfiguration *)getTaskConfiguration();
	for (int i = 0; i < nEventProcessed; i++)
	{
		for (int j = 0; j < hc->numTypes; j++)
		{
			delete[] transverseMomentumMoments[i][j];
		}
		delete[] transverseMomentumMoments[i];
		delete[] yields[i];
	}
	delete[] transverseMomentumMoments;
	delete[] yields;
	delete[] multiplicity;
	delete[] centrality;
	delete[] numParticles;

	//This method will be used during subSample analyses only.
	//This means that before resetting, the histos will be saved.
	//In saveHistograms(), hc->totEvents is set to the # of events in the subSample analyses
	//maxEvents, will be the total # of events over all sub-sample analyses. 
	//It is unecesary to use this as the next sub-sample will only have hc->totEvents events.
	multiplicity = new double[hc->totEvents];  
	centrality = new double[hc->totEvents];
	transverseMomentumMoments = new double **[hc->totEvents];
	yields = new double *[hc->totEvents];
	numParticles = new double[hc->totEvents];

	nEventProcessed = 0;
	nEventAccepted = 0;
	if (isTaskOk())
		histos->reset();
	if (reportEnd("Task", getTaskName(), "reset()"))
		;
}

void PTCorrelator::calculateDerivedHistograms()
{
	TransverseMomentumConfiguration *ac = (TransverseMomentumConfiguration *)getTaskConfiguration();
	ac->totEvents = nEventProcessed;
	histos->fillDerivedHistos(transverseMomentumMoments, yields, multiplicity, centrality, numParticles);
}
