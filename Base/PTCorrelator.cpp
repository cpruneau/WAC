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
	setReportLevel(MessageLogger::Debug);
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
	S             = new double * [maxEvents];
	avgpT         = new double [maxOrder]();
	counts        = new int * [maxEvents];

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
	HeavyIonConfiguration * ac = (HeavyIonConfiguration *) getTaskConfiguration();
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
	if (reportDebug())  cout << "PTCorrelator::calculateAverage(...) Starting." << endl;
	int * counter = new int[maxOrder]();
	double * pts = new double[maxOrder]();
	eventAveragept = 0;
	for (int iParticle=0; iParticle<event->nParticles; iParticle++)
	{
		Particle & particle = * event->getParticleAt(iParticle);
		//if (reportDebug())  particle.printProperties(cout);
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
	if (reportDebug())  cout << "PTCorrelator::calculateAverage(...) Completed." << endl;
}



//////////////////////////////////////////////////////////
//store the transverse momentum of all the particles
//store whether the particles are accepted by each of the filters
//store the multiplicity of each event
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

	//multiplicity
	multiplicity[eventsProcessed] = event->multiplicity;
	centrality[eventsProcessed] = event->centrality;
	if (reportDebug())  cout << "PTCorrelator::storeEventInfo(...) Completed." << endl;
}









/////////////////////////////////////////////////////////////
// Calculate and fill the pT deviation correlators for orders from 1 to maxOrder for all combinations
// Calculate the average counts of particles, pairs, triples ...
////////////////////////////////////////////////////////////
void PTCorrelator::fillSValues(int depth, int filterIndex, int * filters, int & count, int *particles)
{
	if (reportDebug())  cout << "PTCorrelator::fillSValues(...) Starting." << endl;
	
	double correlators = 0;
	for(int i = filterIndex; i < maxOrder; i++)
	{
		filters[maxOrder - depth- 1] = i;
		count = 0;
		switch(maxOrder - depth)
		{
			case 1: correlators = calculateS1(filters, count); break;
			case 2: correlators = calculateS2(filters, count); break;
			case 3: correlators = calculateS3(filters, count); break;
			case 4: correlators = calculateS4(filters, count); break;
		}
		//correlators = calculateS(filters, maxOrder - depth, 0, count, particles);

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
	if (reportDebug())  cout << "PTCorrelator::fillSValues(...) Completed." << endl;
}

double PTCorrelator::calculateS1(int * filters, int & count)
{
	if (reportDebug())  cout << "PTCorrelator::calculateS1(...) Starting." << endl;
	double sum = 0;
	for(int iParticle1 = 0; iParticle1<event->nParticles; iParticle1++)
	{
		Particle & particle1 = * event->getParticleAt(iParticle1);
		double deviation1 = (particle1.pt - eventAveragept);
		if(particleFilters[filters[0]]->accept(particle1))
		{
			sum +=  deviation1;
			count++;
		}
	}
	return sum;
	if (reportDebug())  cout << "PTCorrelator::calculateS1(...) Completed." << endl;
}

double PTCorrelator::calculateS2(int * filters, int & count)
{
	if (reportDebug())  cout << "PTCorrelator::calculateS2(...) Starting." << endl;
	double sum = 0;
	for(int iParticle1 = 0; iParticle1<event->nParticles; iParticle1++)
	{
		Particle & particle1 = * event->getParticleAt(iParticle1);
		double deviation1 = (particle1.pt - eventAveragept);
		if(particleFilters[filters[0]]->accept(particle1))
		{
			for(int iParticle2 = 0; iParticle2<event->nParticles; iParticle2++)
			{
				Particle & particle2 = * event->getParticleAt(iParticle2);
				double deviation2 = (particle2.pt - eventAveragept);
				if(particleFilters[filters[1]]->accept(particle2) && iParticle1 != iParticle2)
				{
					sum +=  deviation1 * deviation2;
					count++;
				}
			}
		}
	}
	return sum;
	if (reportDebug())  cout << "PTCorrelator::calculateS2(...) Completed." << endl;
}

double PTCorrelator::calculateS3(int * filters, int & count)
{
	if (reportDebug())  cout << "PTCorrelator::calculateS3(...) Starting." << endl;
	double sum = 0;
	for(int iParticle1 = 0; iParticle1<event->nParticles; iParticle1++)
	{
		Particle & particle1 = * event->getParticleAt(iParticle1);
		double deviation1 = (particle1.pt - eventAveragept);
		if(particleFilters[filters[0]]->accept(particle1))
		{
			for(int iParticle2 = 0; iParticle2<event->nParticles; iParticle2++)
			{
				Particle & particle2 = * event->getParticleAt(iParticle2);
				double deviation2 = (particle2.pt - eventAveragept);
				if(particleFilters[filters[1]]->accept(particle2) && iParticle1 != iParticle2)
				{
					for(int iParticle3 = 0; iParticle3<event->nParticles; iParticle3++)
					{
						Particle & particle3 = * event->getParticleAt(iParticle3);
						double deviation3 = (particle3.pt - eventAveragept);
						if(particleFilters[filters[2]]->accept(particle3) && iParticle1 != iParticle3 && iParticle2!= iParticle3)
						{
							sum +=  deviation1 * deviation2 * deviation3;
							count++;
						}
					}
				}
			}
		}
	}
	if (reportDebug())  cout << "PTCorrelator::calculateS3(...) Completed." << endl;
	return sum;
}

double PTCorrelator::calculateS4(int * filters, int & count)
{
	if (reportDebug())  cout << "PTCorrelator::calculateS4(...) Starting." << endl;
	double sum = 0;
	for(int iParticle1 = 0; iParticle1<event->nParticles; iParticle1++)
	{
		if(particleFilters[filters[0]]->accept(particle1))
		{
			Particle & particle1 = * event->getParticleAt(iParticle1);
			double deviation1 = (particle1.pt - eventAveragept);
			for(int iParticle2 = 0; iParticle2<event->nParticles; iParticle2++)
			{
				if(particleFilters[filters[1]]->accept(particle2) && iParticle1 != iParticle2)
				{
					Particle & particle2 = * event->getParticleAt(iParticle2);
					double deviation2 = (particle2.pt - eventAveragept);
					if(particleFilters[filters[1]]->accept(particle2) && iParticle1 != iParticle2)
					{
						for(int iParticle3 = 0; iParticle3<event->nParticles; iParticle3++)
						{
							if(particleFilters[filters[2]]->accept(particle3) && iParticle1 != iParticle3 && iParticle2!= iParticle3)
							{
								Particle & particle3 = * event->getParticleAt(iParticle3);
								double deviation3 = (particle3.pt - eventAveragept);
								for(int iParticle4 = 0; iParticle4<event->nParticles; iParticle4++)
								{
									if(particleFilters[filters[3]]->accept(particle4) && iParticle1 != iParticle4 && iParticle2!= iParticle4 && iParticle3!= iParticle4)
									{
										Particle & particle4 = * event->getParticleAt(iParticle4);
										double deviation4 = (particle4.pt - eventAveragept);
										sum +=  deviation1 * deviation2 * deviation3 * deviation4;
										count++;
									}
								}
							}
						}
					}
				}
			}
		}
		if (reportDebug())  cout << "PTCorrelator::calculateS4(...) Completed." << endl;
		return sum;
	}


/////////////////////////////////////////////////////////////
// Calculate the pT deviation correlators(S) for orders from 1 to maxOrder for all combinations
// checked for correctness
////////////////////////////////////////////////////////////
	double PTCorrelator::calculateS(int * filters, int order, int curFilterIndex, int & count, int * particles)
	{
		if (reportDebug())  cout << "PTCorrelator::calculateS(...) Starting." << endl;
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

		if (reportDebug())  cout << "PTCorrelator::calculateS(...) Completed." << endl;
		return sum;

	}