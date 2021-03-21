#include "PTHistos.hpp"
#include "TH1.h"
#include "TProfile.h"
#include "TMath.h"

ClassImp(PTHistos);

PTHistos::PTHistos(const TString &name,
				   TaskConfiguration *configuration,
				   LogLevel debugLevel,
				   int ord)
	: Histograms(name, configuration, (1000000), debugLevel), //note the 1000000 number just needs to be larger than the number of histograms
	  maxOrder(ord),
	  histoIndex(0),
	  size(0),
	  numFunc(3)
{
	TransverseMomentumConfiguration *hc = (TransverseMomentumConfiguration *)getConfiguration();
	size = (TMath::Factorial(maxOrder + hc->numTypes)) / (TMath::Factorial(maxOrder) * TMath::Factorial(hc->numTypes)) - 1;
	//initialize(size);
	createHistograms();
}

PTHistos::PTHistos(TFile *inputFile,
				   const TString &name,
				   TaskConfiguration *configuration,
				   LogLevel debugLevel,
				   int ord)
	: Histograms(name, configuration, (TMath::Factorial(2 * ord) / TMath::Factorial(ord)), debugLevel),
	  maxOrder(ord),
	  histoIndex(0),
	  size(0),
	  numFunc(3)
{
	loadHistograms(inputFile);
}

////////////////////////////////////////////////////
//DTOR
//The commented portions are already deleted in the HistogramCollection Dtor
////////////////////////////////////////////////////
PTHistos::~PTHistos()
{
	//deleteHistograms();
	if (reportDebug())
		cout << "PTHistos::DTOR(...) Started" << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	/*for (int i = 0; i < maxOrder; ++i)
	{
		delete pT[i];
		if (ac.ptCorrelatorVsMult) delete pT_vsMult[i];
		if (ac.ptCorrelatorVsCent) delete pT_vsCent[i];
	}*/
	delete[] pT;
	if (ac.ptCorrelatorVsMult)
		delete[] pT_vsMult;
	if (ac.ptCorrelatorVsCent)
		delete[] pT_vsCent;

	/*delete h_events;
	if (ac.ptCorrelatorVsMult) delete h_events_vsMult;
	if (ac.ptCorrelatorVsCent) delete h_events_vsCent;*/
	for (int i = 0; i < numFunc; ++i)
	{
		for (int j = 0; j < size; j++)
		{
			/*delete hS[i][j];
			if (ac.ptCorrelatorVsMult) delete hS_vsMult[i][j];
			if (ac.ptCorrelatorVsCent) delete hS_vsCent[i][j]; 
			delete hC[i][j];
			if (ac.ptCorrelatorVsMult) delete hC_vsMult[i][j];
			if (ac.ptCorrelatorVsCent) delete hC_vsCent[i][j];*/
		}
		delete[] hS[i];
		if (ac.ptCorrelatorVsMult)
			delete[] hS_vsMult[i];
		if (ac.ptCorrelatorVsCent)
			delete[] hS_vsCent[i];
		delete[] hC[i];
		if (ac.ptCorrelatorVsMult)
			delete[] hC_vsMult[i];
		if (ac.ptCorrelatorVsCent)
			delete[] hC_vsCent[i];
		delete[] names[i];
		delete[] titles[i];
		delete[] names2[i];
		delete[] titles2[i];
	}
	delete[] hS;
	if (ac.ptCorrelatorVsMult)
		delete[] hS_vsMult;
	if (ac.ptCorrelatorVsCent)
		delete[] hS_vsCent;
	delete[] hC;
	if (ac.ptCorrelatorVsMult)
		delete[] hC_vsMult;
	if (ac.ptCorrelatorVsCent)
		delete[] hC_vsCent;
	delete[] names;
	delete[] titles;
	delete[] names2;
	delete[] titles2;

	/*for(int j = 0; j < size; j++)
	{
		delete h_counts[j];
		if (ac.ptCorrelatorVsMult) delete h_counts_vsMult[j];
		if (ac.ptCorrelatorVsCent) delete h_counts_vsCent[j]; 
	}*/
	delete[] h_counts;
	if (ac.ptCorrelatorVsMult)
		delete[] h_counts_vsMult;
	if (ac.ptCorrelatorVsCent)
		delete[] h_counts_vsCent;

	delete[] reorder;
	delete[] orders;

	/*for(int i = 0; i < totEvents; i++)
	{
		delete[] SValues[i];// already deleted
	}*/
	delete[] SValues;
	if (reportDebug())
		cout << "PTHistos::DTOR(...) Completed" << endl;
}

// for now use the same boundaries for eta and y histogram
void PTHistos::createHistograms()
{
	if (reportDebug())
		cout << "PTHistos::createHistograms(...) started" << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	TString bn = getHistoBaseName();
	TH1::SetDefaultBufferSize(ac.totEvents);
	totEvents = ac.totEvents; // Note: totEvents will hold the max # of events that will be done. ac.totEvents will initially hold that # then later hold the number of events processed up to that point for each partial save

	// ================================================================================
	// Naming convention
	// ================================================================================
	// S is the pT deviation moments
	// s are the yield-normalized moments
	// s* are the moments normalized by average transverse momenta
	// C is the cumulants
	// c are the yield-normalized cumulants
	// c* are the cumulants normalizd by average transverse momenta

	size = (TMath::Factorial(maxOrder + ac.numTypes)) / (TMath::Factorial(maxOrder) * TMath::Factorial(ac.numTypes)) - 1;

	S = new TProfile *[size];
	S_vsMult = new TProfile *[size];
	S_vsCent = new TProfile *[size];

	hS = new TH1 **[numFunc];
	hS_vsMult = new TH1 **[numFunc];
	hS_vsCent = new TH1 **[numFunc];

	hC = new TH1 **[numFunc];
	hC_vsMult = new TH1 **[numFunc];
	hC_vsCent = new TH1 **[numFunc];

	names = new TString *[numFunc];
	titles = new TString *[numFunc];
	names2 = new TString *[numFunc];
	titles2 = new TString *[numFunc];

	for (int i = 0; i < numFunc; i++)
	{
		hS[i] = new TH1 *[size];
		hS_vsMult[i] = new TH1 *[size];
		hS_vsCent[i] = new TH1 *[size];

		hC[i] = new TH1 *[size];
		hC_vsMult[i] = new TH1 *[size];
		hC_vsCent[i] = new TH1 *[size];

		names[i] = new TString[size];
		titles[i] = new TString[size];
		names2[i] = new TString[size];
		titles2[i] = new TString[size];
	}

	h_counts = new TProfile *[size];
	h_counts_vsMult = new TProfile *[size];
	h_counts_vsCent = new TProfile *[size];

	orders = new int[size];

	h_events = createHistogram(bn + TString("Nevents"), 1, ac.min_mult, ac.max_mult, "mult", "n_{Events}");
	if (ac.ptCorrelatorVsMult)
		h_events_vsMult = createHistogram(bn + TString("Nevents_vsMult"), ac.nBins_mult, ac.min_mult, ac.max_mult, "mult", "n_{Events}");
	if (ac.ptCorrelatorVsCent)
		h_events_vsCent = createHistogram(bn + TString("Nevents_vsCent"), ac.nBins_cent, ac.min_cent, ac.max_cent, "cent", "n_{Events}");

	pT = new TProfile *[ac.numTypes];
	pT_vsMult = new TProfile *[ac.numTypes];
	pT_vsCent = new TProfile *[ac.numTypes];

	for (int i = 0; i < ac.numTypes; i++)
	{
		pT[i] = createProfile(bn + TString("AverageTransverseMomentum_") + (i + 1), 1, ac.min_mult, ac.max_mult, "mult", TString("AverageTransverseMomentum_") + i);
		if (ac.ptCorrelatorVsMult)
			pT_vsMult[i] = createProfile(bn + TString("AverageTransverseMomentum_") + (i + 1) + TString("_vsMult"), ac.nBins_mult, ac.min_mult, ac.max_mult, "mult", TString("AverageTransverseMomentum_") + i);
		if (ac.ptCorrelatorVsCent)
			pT_vsCent[i] = createProfile(bn + TString("AverageTransverseMomentum_") + (i + 1) + TString("_vsCent"), ac.nBins_cent, ac.min_cent, ac.max_cent, "cent", TString("AverageTransverseMomentum_") + i);
	}

	TString *baseName = new TString[2 * numFunc + 1];
	baseName[0] = bn + "S_";
	baseName[1] = bn + "s_";
	baseName[2] = bn + "s*_";
	baseName[3] = bn + "C_";
	baseName[4] = bn + "c_";
	baseName[5] = bn + "c*_";
	baseName[6] = bn + "Yields_";

	TString *baseTitle = new TString[2 * numFunc + 1];
	baseTitle[0] = "S_{";
	baseTitle[1] = "s_{";
	baseTitle[2] = "s*_{";
	baseTitle[3] = "C_{";
	baseTitle[4] = "c_{";
	baseTitle[5] = "c*_{";
	baseTitle[6] = "Yields_{";

	histoIndex = 0;
	createHistogramRec(baseName, baseTitle, maxOrder - 1, 0);

	reorder = new int[size];
	reorder2 = new int[size];
	int counter = 0;
	for (int iOrd = 1; iOrd <= maxOrder; iOrd++)
	{
		for (int iHisto = 0; iHisto < size; iHisto++)
		{
			if (orders[iHisto] == iOrd)
			{
				reorder[iHisto] = counter;
				reorder2[counter] = iHisto;
				counter++;
			}
		}
	}

	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		for (int iFunc = 0; iFunc < numFunc; ++iFunc)
		{
			hS[iFunc][iHisto] = createHistogram(names[iFunc][iHisto], 1, ac.min_mult, ac.max_mult, "mult", titles[iFunc][iHisto] + "}");
			if (ac.ptCorrelatorVsMult)
				hS_vsMult[iFunc][iHisto] = createHistogram(names[iFunc][iHisto] + "_vsMult", ac.nBins_mult, ac.min_mult, ac.max_mult, "mult", titles[iFunc][iHisto] + "}");
			if (ac.ptCorrelatorVsCent)
				hS_vsCent[iFunc][iHisto] = createHistogram(names[iFunc][iHisto] + "_vsCent", ac.nBins_cent, ac.min_cent, ac.max_cent, "cent", titles[iFunc][iHisto] + "}");
		}
	}
	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		for (int i = 0; i < numFunc; ++i)
		{
			hC[i][iHisto] = createHistogram(names2[i][iHisto], 1, ac.min_mult, ac.max_mult, "mult", titles2[i][iHisto] + "}");
			if (ac.ptCorrelatorVsMult)
				hC_vsMult[i][iHisto] = createHistogram(names2[i][iHisto] + "_vsMult", ac.nBins_mult, ac.min_mult, ac.max_mult, "mult", titles2[i][iHisto] + "}");
			if (ac.ptCorrelatorVsCent)
				hC_vsCent[i][iHisto] = createHistogram(names2[i][iHisto] + "_vsCent", ac.nBins_cent, ac.min_cent, ac.max_cent, "cent", titles2[i][iHisto] + "}");
		}
	}
	TString s = "S";
	TString sm = "S_vsMult";
	TString sc = "S_vsCent";
	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		S[iHisto] = createProfile(s + iHisto, 1, ac.min_mult, ac.max_mult, "mult", titles[0][iHisto] + "}");
		if (ac.ptCorrelatorVsMult)
			S_vsMult[iHisto] = createProfile(sm + iHisto, ac.nBins_mult, ac.min_mult, ac.max_mult, "mult", titles[0][iHisto] + "}");
		if (ac.ptCorrelatorVsCent)
			S_vsCent[iHisto] = createProfile(sc + iHisto, ac.nBins_cent, ac.min_cent, ac.max_cent, "cent", titles[0][iHisto] + "}");
	}
	histoIndex = 0;

	if (reportDebug())
		cout << "PTHistos::createHistograms(...) ended" << endl;
	//h_c123vsMultTest = createProfile("c123Test", ac.nBins_mult,ac.min_mult,ac.max_mult,"mult", "c123Test" );
}

/////////////////////////////////////////////////////////
// need to fix so that it mirrors above
////////////////////////////////////////////////////////
void PTHistos::loadHistograms(TFile *inputFile)
{
	if (!inputFile)
	{
		if (reportFatal())
			cout << "-Fatal- Attempting to load NuDynHistos from an invalid file pointer" << endl;
		return;
	}
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());

	TString bn = getHistoBaseName();

	size = (TMath::Factorial(maxOrder + ac.numTypes)) / (TMath::Factorial(maxOrder) * TMath::Factorial(ac.numTypes)) - 1;

	hS = new TH1 **[numFunc];
	hS_vsMult = new TH1 **[numFunc];
	hS_vsCent = new TH1 **[numFunc];

	hC = new TH1 **[numFunc];
	hC_vsMult = new TH1 **[numFunc];
	hC_vsCent = new TH1 **[numFunc];

	for (int i = 0; i < numFunc; i++)
	{
		hS[i] = new TH1 *[size];
		hS_vsMult[i] = new TH1 *[size];
		hS_vsCent[i] = new TH1 *[size];

		hC[i] = new TH1 *[size];
		hC_vsMult[i] = new TH1 *[size];
		hC_vsCent[i] = new TH1 *[size];
	}

	h_counts = new TProfile *[size];
	h_counts_vsMult = new TProfile *[size];
	h_counts_vsCent = new TProfile *[size];

	orders = new int[size];

	h_events = loadH1(inputFile, bn + TString("Nevents"));
	if (ac.ptCorrelatorVsMult)
		h_events_vsMult = loadH1(inputFile, bn + TString("Nevents_vsMult"));
	if (ac.ptCorrelatorVsCent)
		h_events_vsCent = loadH1(inputFile, bn + TString("Nevents_vsCent"));

	TString *baseName = new TString[numFunc + 1];
	baseName[0] = bn + "S_";
	baseName[1] = bn + "s_";
	baseName[2] = bn + "s*_";
	baseName[3] = bn + "C_";
	baseName[4] = bn + "c_";
	baseName[5] = bn + "c*_";
	baseName[6] = bn + "Counts_";

	TString *baseTitle = new TString[numFunc + 1];
	baseTitle[0] = "S_{";
	baseTitle[1] = "s_{";
	baseTitle[2] = "s*_{";
	baseTitle[3] = "C_{";
	baseTitle[4] = "c_{";
	baseTitle[5] = "c*_{";
	baseTitle[6] = "Counts_{";

	loadHistogramRec(baseName, maxOrder - 1, 0, inputFile);

	histoIndex = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
//overloaded HistogramCollection::saveHistograms to save histograms in sequence of lowest order to highest order
//reorders the histograms from the "recursive" order that they are created in
///////////////////////////////////////////////////////////////////////////////////////////
void PTHistos::saveHistograms(TFile *outputFile, bool saveAll)
{

	if (reportDebug())
		cout << "HistogramCollection::saveHistograms(TFile * outputFile) started." << endl;
	outputFile->cd();

	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());

	int numTypes = 1;
	if (ac.ptCorrelatorVsMult)
		numTypes++;
	if (ac.ptCorrelatorVsCent)
		numTypes++;

	int extra = numTypes * (ac.numTypes + 1);

	//saves event and transverse momentum histos
	//there are "extra" of these at the start of the collection
	for (int k = 0; k < extra; k++)
	{
		if (isSaved(options[k]) || saveAll)
			getObjectAt(k)->Write();
	}

	for (int iFunc = 0; iFunc <= 2 * numFunc; iFunc++)
	{
		for (int i = 1; i <= maxOrder; i++)
		{
			for (int k = numTypes; k < getCollectionSize(); k++)
			{
				if (k < size * numTypes + extra && iFunc == 2 * numFunc)
				{
					//These are the "yields" histograms
					//there are size * numTypes of them(plus the "extra" first ones) hence k < size * numTypes + extra
					//Want them to be saved last hence the iFunc == 2 * numFunc
					int k1 = (k - extra); //We want to know the relative order of the histogram in its group.
					int orderIndex = k1 / numTypes;
					//We want the count and count_vsMult to be saved right next to each other.
					//They are created together, so we are only concerned about the clumps of (numtypes) histograms, hence the "/ numTypes"
					if ((isSaved(options[k]) || saveAll) && (orders[orderIndex] == i))
						getObjectAt(k)->Write(); // want to save them in increasing order
				}
				if (k >= size * numTypes + extra && k < size * (numFunc + 1) * numTypes + extra && iFunc < numFunc)
				{
					int k1 = k - size * numTypes - extra;
					int orderIndex = (k1) / (numTypes * (numFunc));
					int funcIndex = (k1) / numTypes - ((numFunc)) * ((k1) / (numTypes * (numFunc))); // We wish to save S, then s, then s*.
					if ((isSaved(options[k]) || saveAll) && (orders[orderIndex] == i) && ((funcIndex) % (numFunc) == (iFunc % (numFunc))))
						getObjectAt(k)->Write();
				}
				if (k >= size * (numFunc + 1) * numTypes + extra && k <= size * (2 * numFunc + 1) * numTypes + extra && iFunc >= numFunc && iFunc < 2 * numFunc)
				{
					int k1 = k - size * (numFunc + 1) * numTypes - extra;
					int orderIndex = (k1) / (numTypes * (numFunc));
					int funcIndex = (k1) / numTypes - ((numFunc)) * ((k1) / (numTypes * (numFunc)));
					if ((isSaved(options[k]) || saveAll) && (orders[orderIndex] == i) && ((funcIndex) % (numFunc) == (iFunc - (numFunc))))
						getObjectAt(k)->Write();
				}
			}
		}
	}

	if (reportDebug())
		cout << "HistogramCollection::saveHistograms(TFile * outputFile) completed." << endl;
}

void PTHistos::fillEventHistos(double mult, double cent, double weight)
{
	if (reportDebug())
		cout << "PTHistos::fillEventHistos(...) started" << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	h_events->Fill(mult, weight);
	if (ac.ptCorrelatorVsMult)
		h_events_vsMult->Fill(mult, weight);
	if (ac.ptCorrelatorVsCent)
		h_events_vsCent->Fill(cent, weight);
	if (reportDebug())
		cout << "HistogramCollection::fillEventHistos(...)completed." << endl;
}

void PTHistos::fillTransverseMomentumHistos(double transverseMomentum, int filter, double mult, double cent, double weight)
{
	if (reportDebug())
		cout << "PTHistos::fillTransverseMomentumHistos(...) started" << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	pT[filter]->Fill(mult, transverseMomentum, weight);
	if (ac.ptCorrelatorVsMult)
		pT_vsMult[filter]->Fill(mult, transverseMomentum, weight);
	if (ac.ptCorrelatorVsCent)
		pT_vsCent[filter]->Fill(cent, transverseMomentum, weight);
	if (reportDebug())
		cout << "HistogramCollection::fillTransverseMomentumHistos(...) completed." << endl;
}

///////////////////////////////////////////////////////////////////////////////////
// recursively create histograms for correlation functions of order 1 - maxOrder
// Note: these histograms are not in sequence of lowest order to highest order. They are in depth first recursive order: 1, 11, 111, 1111, 1112, 1113, 1114, 112, 1121 ...
//The reording occurs when the histograms are saved with the saveHistograms function.
/////////////////////////////////////////////////////////////////////////////////
void PTHistos::createHistogramRec(TString *baseName, TString *baseTitle, int depth, int partIndex)
{
	if (reportDebug())
		cout << "PTHistos::createHistogramsRec(...) started" << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	TString *histoName = new TString[2 * numFunc + 1];
	TString *histoTitle = new TString[2 * numFunc + 1];

	for (int i = partIndex; i < ac.numTypes; i++)
	{

		for (int iFunc = 0; iFunc < numFunc; iFunc++)
		{
			histoName[iFunc] = baseName[iFunc] + (i + 1);
			histoTitle[iFunc] = (depth == maxOrder - 1) ? baseTitle[iFunc] + (i + 1) : baseTitle[iFunc] + ", " + (i + 1);
			histoName[iFunc + numFunc] = baseName[iFunc + numFunc] + (i + 1);
			histoTitle[iFunc + numFunc] = (depth == maxOrder - 1) ? baseTitle[iFunc + numFunc] + (i + 1) : baseTitle[iFunc + numFunc] + ", " + (i + 1);

			names[iFunc][histoIndex] = histoName[iFunc];
			titles[iFunc][histoIndex] = histoTitle[iFunc];

			/*hS[iFunc][histoIndex] = createProfile( histoName[iFunc], 1, ac.max_mult, ac.min_mult,"mult",histoTitle[iFunc] + "}" );
			if (ac.ptCorrelatorVsMult)	hS_vsMult[iFunc][histoIndex] = createProfile(histoName[iFunc] + "_vsMult",ac.nBins_mult, ac.max_mult, ac.min_mult, "mult", histoTitle[iFunc] + "}");
			if (ac.ptCorrelatorVsCent)	hS_vsCent[iFunc][histoIndex] = createProfile(histoName[iFunc] + "_vsCent",ac.nBins_cent,ac.max_mult, ac.min_mult, "cent", histoTitle[iFunc] + "}");*/

			names2[iFunc][histoIndex] = histoName[iFunc + numFunc];
			titles2[iFunc][histoIndex] = histoTitle[iFunc + numFunc];
			/*hC[iFunc][histoIndex] = createHistogram( histoName[iFunc + numFunc], 1, ac.min_mult, ac.max_mult,"mult",histoTitle[iFunc + numFunc] + "}" );
			if (ac.ptCorrelatorVsMult)	hC_vsMult[iFunc][histoIndex] = createHistogram(histoName[iFunc + numFunc] + "_vsMult",ac.nBins_mult,ac.min_mult, ac.max_mult, "mult", histoTitle[iFunc + numFunc] + "}");
			if (ac.ptCorrelatorVsCent)	hC_vsCent[iFunc][histoIndex] = createHistogram(histoName[iFunc + numFunc] + "_vsCent",ac.nBins_cent,ac.min_mult, ac.max_mult, "cent", histoTitle[iFunc + numFunc] + "}");*/
		}

		histoName[2 * numFunc] = baseName[2 * numFunc] + (i + 1);
		histoTitle[2 * numFunc] = (depth == maxOrder - 1) ? baseTitle[2 * numFunc] + (i + 1) : baseTitle[2 * numFunc] + ", " + (i + 1);

		h_counts[histoIndex] = createProfile(histoName[2 * numFunc], 1, ac.min_mult, ac.max_mult, "mult", histoTitle[2 * numFunc] + "}");
		if (ac.ptCorrelatorVsMult)
			h_counts_vsMult[histoIndex] = createProfile(histoName[2 * numFunc] + "_vsMult", ac.nBins_mult, ac.min_mult, ac.max_mult, "mult", histoTitle[2 * numFunc] + "}");
		if (ac.ptCorrelatorVsCent)
			h_counts_vsCent[histoIndex] = createProfile(histoName[2 * numFunc] + "_vsCent", ac.nBins_cent, ac.min_cent, ac.max_cent, "cent", histoTitle[2 * numFunc] + "}");

		orders[histoIndex] = maxOrder - depth;

		histoIndex++;

		if (depth != 0)
			createHistogramRec(histoName, histoTitle, depth - 1, i);
	}

	delete[] histoName;
	delete[] histoTitle;
	if (reportDebug())
		cout << "PTHistos::createHistogramsRec(...) ended" << endl;
	return;
}

//////////////////////////////////////////////
// need to fix so that it mirrors above
//////////////////////////////////////////////
void PTHistos::loadHistogramRec(TString *baseName, int depth, int partIndex, TFile *inputFile)
{

	if (reportDebug())
		cout << "PTHistos::loadHistogramRec(...) Starting." << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	TString *histoName = new TString[2 * numFunc + 1];

	for (int i = partIndex; i < maxOrder; i++)
	{

		for (int iFunc = 0; iFunc < numFunc; iFunc++)
		{
			histoName[iFunc] = baseName[iFunc] + (i + 1);
			histoName[iFunc + numFunc] = baseName[iFunc + numFunc] + (i + 1);

			hS[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc]);
			if (ac.ptCorrelatorVsMult)
				hS_vsMult[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc] + "_vsMult");
			if (ac.ptCorrelatorVsCent)
				hS_vsCent[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc] + "_vsCent");

			//hC[iFunc][histoIndex] = loadHistogram(inputFile, histoName[iFunc + numFunc]);
			//if (ac.ptCorrelatorVsMult)
			//	hC_vsMult[iFunc][histoIndex] = loadHistogram(inputFile, histoName[iFunc + numFunc] + "_vsMult");
			//if (ac.ptCorrelatorVsCent)
			//	hC_vsCent[iFunc][histoIndex] = loadHistogram(inputFile, histoName[iFunc + numFunc] + "_vsCent");
		}

		histoName[2 * numFunc] = baseName[2 * numFunc] + (i + 1);

		h_counts[histoIndex] = loadProfile(inputFile, histoName[2 * numFunc] + "}");
		if (ac.ptCorrelatorVsMult)
			h_counts_vsMult[histoIndex] = loadProfile(inputFile, histoName[2 * numFunc] + "_vsMult");
		if (ac.ptCorrelatorVsCent)
			h_counts_vsCent[histoIndex] = loadProfile(inputFile, histoName[2 * numFunc] + "_vsCent");

		orders[histoIndex] = maxOrder - depth;

		histoIndex++;

		if (depth != 0)
			loadHistogramRec(histoName, depth - 1, i, inputFile);
	}

	delete[] histoName;
	return;
	if (reportDebug())
		cout << "PTHistos::loadHistogramRec(...) Completed." << endl;
}

////////////////////////////////////////////////////////////////
//Fill the Derived Histograms(Moments and Cumulants of tranverse momentum deviations)
//Uses Histogram addition, division, and multiplication.
///////////////////////////////////////////////////////////////
void PTHistos::fillDerivedHistos(double ***transverseMomentumMoments, double **counts, double *mults, double *cents, double *numParticles)
{
	if (reportDebug())
		cout << "PTHistos::fillDerivedHistos(...) Starting." << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)*getConfiguration();
	yields = counts;
	SValues = new double *[ac.totEvents]();
	//fill yields
	for (int iEvent = 0; iEvent < ac.totEvents; iEvent++)
	{
		for (int iHisto = 0; iHisto < size; iHisto++)
		{
			h_counts[iHisto]->Fill(mults[iEvent], yields[iEvent][iHisto], 1.0); //weight is always 1.0
			if (ac.ptCorrelatorVsMult)
				h_counts_vsMult[iHisto]->Fill(mults[iEvent], yields[iEvent][iHisto], 1.0);
			if (ac.ptCorrelatorVsCent)
				h_counts_vsCent[iHisto]->Fill(cents[iEvent], yields[iEvent][iHisto], 1.0);
		}
	}
	//fill SValues
	//fill SValues normalized by counts and average pT's
	for (int iEvent = 0; iEvent < ac.totEvents; iEvent++)
	{
		SValues[iEvent] = new double[size]();
		int bin = pT[0]->FindBin(mults[iEvent]);
		calculatePTDeviationMoments(transverseMomentumMoments, bin, iEvent, numParticles[iEvent], pT); //calculate the deviation using this bin of the average transverse momentum
		for (int iHisto = 0; iHisto < size; iHisto++)
		{
			if (yields[iEvent][iHisto] != 0)
			{
				S[iHisto]->Fill(mults[iEvent], SValues[iEvent][iHisto], 1.0); //weight is always 1.0
			}
		}

		if (ac.ptCorrelatorVsMult)
		{
			bin = pT_vsMult[0]->FindBin(mults[iEvent]);
			calculatePTDeviationMoments(transverseMomentumMoments, bin, iEvent, numParticles[iEvent], pT_vsMult);
			for (int iHisto = 0; iHisto < size; iHisto++)
			{
				if (yields[iEvent][iHisto] != 0)
				{
					S_vsMult[iHisto]->Fill(mults[iEvent], SValues[iEvent][iHisto], 1.0);
				}
			}
		}

		if (ac.ptCorrelatorVsCent)
		{
			bin = pT_vsCent[0]->FindBin(cents[iEvent]);
			calculatePTDeviationMoments(transverseMomentumMoments, bin, iEvent, numParticles[iEvent], pT_vsCent);
			for (int iHisto = 0; iHisto < size; iHisto++)
			{
				if (yields[iEvent][iHisto] != 0)
				{
					S_vsCent[iHisto]->Fill(cents[iEvent], SValues[iEvent][iHisto], 1.0);
				}
			}
		}
		delete SValues[iEvent]; // not needed after this
	}

	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		for (int i = 1; i <= 1; i++)
		{
			hS[0][iHisto]->SetBinContent(i, S[iHisto]->GetBinContent(i));
			hS[0][iHisto]->SetBinError(i, S[iHisto]->GetBinError(i));
		}
		if (ac.ptCorrelatorVsMult)
		{
			for (int i = 1; i <= ac.nBins_mult; i++)
			{
				hS_vsMult[0][iHisto]->SetBinContent(i, S_vsMult[iHisto]->GetBinContent(i));
				hS_vsMult[0][iHisto]->SetBinError(i, S_vsMult[iHisto]->GetBinError(i));
			}
		}
		if (ac.ptCorrelatorVsCent)
		{
			for (int i = 1; i <= ac.nBins_cent; i++)
			{
				hS_vsCent[0][iHisto]->SetBinContent(i, S_vsCent[iHisto]->GetBinContent(i));
				hS_vsCent[0][iHisto]->SetBinError(i, S_vsCent[iHisto]->GetBinError(i));
			}
		}
	}

	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		for (int i = 1; i <= 1; i++)
		{
			if (h_counts[iHisto]->GetBinContent(i) == 0)
			{
				hS[1][iHisto]->SetBinContent(i, 0);
				hS[1][iHisto]->SetBinError(i, 0);
				continue;
			}
			hS[1][iHisto]->SetBinContent(i, hS[0][iHisto]->GetBinContent(i) / h_counts[iHisto]->GetBinContent(i));
			double absError1 = hS[0][iHisto]->GetBinError(i);
			double binContent1 = hS[0][iHisto]->GetBinContent(i);
			double absError2 = h_counts[iHisto]->GetBinError(i);
			double binContent2 = h_counts[iHisto]->GetBinContent(i);
			double abserr = TMath::Sqrt(absError1 * absError1 * binContent2 * binContent2 + absError2 * absError2 * binContent1 * binContent1) / (binContent2 * binContent2);
			hS[1][iHisto]->SetBinError(i, abserr);
		}
		if (ac.ptCorrelatorVsMult)
		{
			for (int i = 1; i <= ac.nBins_mult; i++)
			{
				if (h_counts_vsMult[iHisto]->GetBinContent(i) == 0)
				{
					hS_vsMult[1][iHisto]->SetBinContent(i, 0);
					hS_vsMult[1][iHisto]->SetBinError(i, 0);
					continue;
				}
				hS_vsMult[1][iHisto]->SetBinContent(i, hS_vsMult[0][iHisto]->GetBinContent(i) / h_counts_vsMult[iHisto]->GetBinContent(i));
				double absError1 = hS_vsMult[0][iHisto]->GetBinError(i);
				double binContent1 = hS_vsMult[0][iHisto]->GetBinContent(i);
				double absError2 = h_counts_vsMult[iHisto]->GetBinError(i);
				double binContent2 = h_counts_vsMult[iHisto]->GetBinContent(i);
				double abserr = TMath::Sqrt(absError1 * absError1 * binContent2 * binContent2 + absError2 * absError2 * binContent1 * binContent1) / (binContent2 * binContent2);
				hS_vsMult[1][iHisto]->SetBinError(i, abserr);
			}
		}
		if (ac.ptCorrelatorVsCent)
		{
			for (int i = 1; i <= ac.nBins_cent; i++)
			{
				if (h_counts_vsCent[iHisto]->GetBinContent(i) == 0)
				{
					hS_vsCent[1][iHisto]->SetBinContent(i, 0);
					hS_vsCent[1][iHisto]->SetBinError(i, 0);
					continue;
				}
				hS_vsCent[1][iHisto]->SetBinContent(i, hS_vsCent[0][iHisto]->GetBinContent(i) / h_counts_vsCent[iHisto]->GetBinContent(i));
				double absError1 = hS_vsCent[0][iHisto]->GetBinError(i);
				double binContent1 = hS_vsCent[0][iHisto]->GetBinContent(i);
				double absError2 = h_counts_vsCent[iHisto]->GetBinError(i);
				double binContent2 = h_counts_vsCent[iHisto]->GetBinContent(i);
				double abserr = TMath::Sqrt(absError1 * absError1 * binContent2 * binContent2 + absError2 * absError2 * binContent1 * binContent1) / (binContent2 * binContent2);
				hS_vsCent[1][iHisto]->SetBinError(i, abserr);
			}
		}
	}

	for (int i = 1; i <= 1; i++)
	{
		histoIndex = 0;
		fillNormalizedPTValues(maxOrder - 1, 0, 1, 0, i, hS[0], pT, hS[2]);
	}
	if (ac.ptCorrelatorVsMult)
	{
		for (int i = 1; i <= ac.nBins_mult; i++)
		{
			histoIndex = 0;
			fillNormalizedPTValues(maxOrder - 1, 0, 1, 0, i, hS_vsMult[0], pT_vsMult, hS_vsMult[2]);
		}
	}
	if (ac.ptCorrelatorVsCent)
	{
		for (int i = 1; i <= ac.nBins_cent; i++)
		{
			histoIndex = 0;
			fillNormalizedPTValues(maxOrder - 1, 0, 1, 0, i, hS_vsCent[0], pT_vsCent, hS_vsCent[2]);
		}
	}

	TH1 ***newhCValues = new TH1 **[3];
	for (int iHisto = 0; iHisto < 3; iHisto++)
	{
		newhCValues[iHisto] = new TH1 *[size];
	}
	calculateCumulants(hS[0], newhCValues[0], 1, ac.min_mult, ac.max_mult);
	if (ac.ptCorrelatorVsMult)
		calculateCumulants(hS_vsMult[0], newhCValues[1], ac.nBins_mult, ac.min_mult, ac.max_mult);
	if (ac.ptCorrelatorVsCent)
		calculateCumulants(hS_vsCent[0], newhCValues[2], ac.nBins_cent, ac.min_cent, ac.max_cent);

	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		for (int i = 1; i <= 1; i++)
		{
			hC[0][iHisto]->SetBinContent(i, newhCValues[0][reorder[iHisto]]->GetBinContent(i));
			hC[0][iHisto]->SetBinError(i, newhCValues[0][reorder[iHisto]]->GetBinError(i));
		}
		if (ac.ptCorrelatorVsMult)
		{
			for (int i = 1; i <= ac.nBins_mult; i++)
			{
				hC_vsMult[0][iHisto]->SetBinContent(i, newhCValues[1][reorder[iHisto]]->GetBinContent(i));
				hC_vsMult[0][iHisto]->SetBinError(i, newhCValues[1][reorder[iHisto]]->GetBinError(i));
			}
		}
		if (ac.ptCorrelatorVsCent)
		{
			for (int i = 1; i <= ac.nBins_cent; i++)
			{
				hC_vsCent[0][iHisto]->SetBinContent(i, newhCValues[2][reorder[iHisto]]->GetBinContent(i));
				hC_vsCent[0][iHisto]->SetBinError(i, newhCValues[2][reorder[iHisto]]->GetBinError(i));
			}
		}
	}

	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		for (int i = 1; i <= 1; i++)
		{
			if (h_counts[iHisto]->GetBinContent(i) == 0)
			{
				hC[1][iHisto]->SetBinContent(i, 0);
				hC[1][iHisto]->SetBinError(i, 0);
				continue;
			}
			hC[1][iHisto]->SetBinContent(i, hC[0][iHisto]->GetBinContent(i) / h_counts[iHisto]->GetBinContent(i));
			double absError1 = hC[0][iHisto]->GetBinError(i);
			double binContent1 = hC[0][iHisto]->GetBinContent(i);
			double absError2 = h_counts[iHisto]->GetBinError(i);
			double binContent2 = h_counts[iHisto]->GetBinContent(i);
			double abserr = TMath::Sqrt(absError1 * absError1 * binContent2 * binContent2 + absError2 * absError2 * binContent1 * binContent1) / (binContent2 * binContent2);
			hC[1][iHisto]->SetBinError(i, abserr);
		}
		if (ac.ptCorrelatorVsMult)
		{
			for (int i = 1; i <= ac.nBins_mult; i++)
			{
				if (h_counts_vsMult[iHisto]->GetBinContent(i) == 0)
				{
					hC_vsMult[1][iHisto]->SetBinContent(i, 0);
					hC_vsMult[1][iHisto]->SetBinError(i, 0);
					continue;
				}
				hC_vsMult[1][iHisto]->SetBinContent(i, hC_vsMult[0][iHisto]->GetBinContent(i) / h_counts_vsMult[iHisto]->GetBinContent(i));
				double absError1 = hC_vsMult[0][iHisto]->GetBinError(i);
				double binContent1 = hC_vsMult[0][iHisto]->GetBinContent(i);
				double absError2 = h_counts_vsMult[iHisto]->GetBinError(i);
				double binContent2 = h_counts_vsMult[iHisto]->GetBinContent(i);
				double abserr = TMath::Sqrt(absError1 * absError1 * binContent2 * binContent2 + absError2 * absError2 * binContent1 * binContent1) / (binContent2 * binContent2);
				hC_vsMult[1][iHisto]->SetBinError(i, abserr);
			}
		}
		if (ac.ptCorrelatorVsCent)
		{
			for (int i = 1; i <= ac.nBins_cent; i++)
			{
				if (h_counts_vsCent[iHisto]->GetBinContent(i) == 0)
				{
					hC_vsCent[1][iHisto]->SetBinContent(i, 0);
					hC_vsCent[1][iHisto]->SetBinError(i, 0);
					continue;
				}
				hC_vsCent[1][iHisto]->SetBinContent(i, hC_vsCent[0][iHisto]->GetBinContent(i) / h_counts_vsCent[iHisto]->GetBinContent(i));
				double absError1 = hC_vsCent[0][iHisto]->GetBinError(i);
				double binContent1 = hC_vsCent[0][iHisto]->GetBinContent(i);
				double absError2 = h_counts_vsCent[iHisto]->GetBinError(i);
				double binContent2 = h_counts_vsCent[iHisto]->GetBinContent(i);
				double abserr = TMath::Sqrt(absError1 * absError1 * binContent2 * binContent2 + absError2 * absError2 * binContent1 * binContent1) / (binContent2 * binContent2);
				hC_vsCent[1][iHisto]->SetBinError(i, abserr);
			}
		}
	}

	for (int i = 1; i <= 1; i++)
	{
		histoIndex = 0;
		fillNormalizedPTValues(maxOrder - 1, 0, 1, 0, i, hC[0], pT, hC[2]);
	}
	if (ac.ptCorrelatorVsMult)
	{
		for (int i = 1; i <= ac.nBins_mult; i++)
		{
			histoIndex = 0;
			fillNormalizedPTValues(maxOrder - 1, 0, 1, 0, i, hC_vsMult[0], pT_vsMult, hC_vsMult[2]);
		}
	}
	if (ac.ptCorrelatorVsCent)
	{
		for (int i = 1; i <= ac.nBins_cent; i++)
		{
			histoIndex = 0;
			fillNormalizedPTValues(maxOrder - 1, 0, 1, 0, i, hC_vsCent[0], pT_vsCent, hC_vsCent[2]);
		}
	}

	for (int iHisto = 0; iHisto < 2; iHisto++)
	{
		for (int i = 0; i < size; i++)
		{
			delete newhCValues[iHisto][i];
		}
		delete[] newhCValues[iHisto];
	}
	delete[] newhCValues;

	if (reportDebug())
		cout << "PTHistos::fillDerivedHistos(...) Completed." << endl;
}

////////////////////////////////////////////
//fill  S or C Values normalized by average transverse momenta
/////////////////////////////////////////////
void PTHistos::fillNormalizedPTValues(int depth, int partIndex, double product, double absErrProduct, int bin, TH1 **OldHistos, TProfile **pTHistos, TH1 **newHistos)
{
	if (reportDebug())
		cout << "PTHistos::fillNormalizedPTValues(...) Starting." << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)*getConfiguration();

	for (int i = partIndex; i < ac.numTypes; i++)
	{
		double newProduct;
		newProduct = product * pTHistos[i]->GetBinContent(bin);
		if (newProduct == 0)
			newHistos[histoIndex]->SetBinContent(bin, 0);
		else
			newHistos[histoIndex]->SetBinContent(bin, OldHistos[histoIndex]->GetBinContent(bin) / newProduct); //newProduct content will likely never be , because the pt is more or less always positive

		double absErrPt = pTHistos[i]->GetBinError(bin); //the content will likely never be 0
		double absErrNewProduct = TMath::Sqrt(absErrPt * absErrPt * product * product + absErrProduct * absErrProduct * pTHistos[i]->GetBinContent(bin) * pTHistos[i]->GetBinContent(bin));
		double absErrS = OldHistos[histoIndex]->GetBinError(bin);
		double absErrQuotient;
		if (newProduct == 0)
			absErrQuotient = 0;
		else
			absErrQuotient = TMath::Sqrt(absErrS * absErrS * newProduct * newProduct + absErrNewProduct * absErrNewProduct * OldHistos[histoIndex]->GetBinContent(bin) * OldHistos[histoIndex]->GetBinContent(bin)) / (newProduct * newProduct);
		newHistos[histoIndex]->SetBinError(bin, absErrQuotient);

		if (histoIndex != size - 1)
			histoIndex++;

		if (depth != 0)
			fillNormalizedPTValues(depth - 1, i, newProduct, absErrNewProduct, bin, OldHistos, pTHistos, newHistos);
	}
	if (reportDebug())
		cout << "PTHistos::fillNormalizedPTValues(...) Completed." << endl;
}

///////////////////////////////////////////
//calculate the cumulants of the moments
//////////////////////////////////////////
void PTHistos::calculateCumulants(TH1 **SHistos, TH1 **CHistos, int nBins, double min, double max)
{
	if (reportDebug())
		cout << "PTHistos::calculateCumulants(...) Starting." << endl;

	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());

	TH1 **newSHistos = new TH1 *[size];
	int counter = 0;
	for (int iOrd = 1; iOrd <= maxOrder; iOrd++)
	{
		for (int iHisto = 0; iHisto < size; iHisto++)
		{
			if (orders[iHisto] == iOrd)
			{
				newSHistos[counter] = SHistos[iHisto];
				counter++;
			}
		}
	}

	for (int iBin = 1; iBin <= nBins; iBin++)
	{
		int maxSize = TMath::Factorial(ac.numTypes);
		double *used = new double[maxSize];
		for (int iHisto = 0; iHisto < size; iHisto++)
		{
			int len = 0;
			int *ind = convert(iHisto, len);

			int curInd = 0;
			int *set = new int[len];
			for (int i = 0; i < len; i++)
			{
				set[i] = i + 1;
			}

			TString name = nBins * size + iHisto;
			if (iBin == 1)
				CHistos[iHisto] = new TH1F(name, name, nBins, min, max);

			double sum = 0;
			double sumErr = 0;

			calcRecSum(CHistos, ind, iBin, set, len, set, len, 1, used, curInd, 1, sum, sumErr, 0, 0);
			CHistos[iHisto]->SetBinContent(iBin, newSHistos[iHisto]->GetBinContent(iBin) - sum);
			CHistos[iHisto]->SetBinError(iBin, TMath::Sqrt(newSHistos[iHisto]->GetBinError(iBin) * newSHistos[iHisto]->GetBinError(iBin) + sumErr * sumErr));

			//ind is now uneccesary
			delete[] ind;
		}
		delete[] used;
	}
	return;
	if (reportDebug())
		cout << "PTHistos::calculateCumulants(...) Completed." << endl;
}

///////////////////////////////////////////////////
//used to calculate a recursive partition sum used in calculation of the cumulants
// ex: 	if the highest level cumulant corresponds to C(1,2,3), this function calculates sum = C(1, 2) * C(3) + C(1, 3) * C(2) + C(2, 3) * C(1) + C(1) * C(2) * C(3)
//	   	At the first level, it starts with subset {3}, and gets C(3) * C(1, 2) and adds it to sum and places it on the used array
//	   	It then proceeds to recurse on the proper subsets of {2, 3} and multiply these by C(3), so it gets C(1) * C(2) * C(3) and adds it to sum and places it on the used array
//	   	Going back to the top level, it goes to subset {2} and gets C(2) * C(1, 3) and adds it to sum and places it on the used array
//		It then proceeds to recurse on the proper subsets of {1, 3} and multiply these by C(2), so it gets C(1) * C(2) * C(3) and finds it on the used array so skips
//		Going back to the top level, it goes to subset {2, 3} and gets C(2, 3) * C(1) and adds it to sum and places it on the used array
//		It then proceeds to recurse on the proper subsets of {1} of which there are none
// 		so on...
//     	After this function concludes, C(1, 2, 3) = S(1, 2, 3) - sum
// note that C(1) = S(1) because sum = 0
///////////////////////////////////////////////////
void PTHistos::calcRecSum(TH1 **CHistos, int *iHisto, int iBin, int *Subset, int len, int *set, int lenSet, double productC, double *used, int &curInd, int productS, double &sum, double &sumError, double absErrProductC, int depth)
{
	if (reportDebug())
		cout << "PTHistos::calcRecSum(...) Starting." << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	int lenSub = 0;
	int lenCompSub = 0;
	if (len != 1)
	{
		char *bin = new char[len];
		int *subset;
		int *comp;
		int *compSubset;

		// loop over all possible proper subsets
		for (int iSubset = 1; iSubset < pow(2, len) - 1; iSubset++)
		{
			//get the subset corresponding to iSubset, multiply it by the product(which contains cumulant terms from previous recursive levels)
			convertToBinary(iSubset, bin, len);
			subset = getSubset(bin, Subset, len, lenSub);
			int *s = getSubset(bin, iHisto, len, lenSub);
			int iSub = convert(s, lenSub);

			double tempProduct = productC * CHistos[iSub]->GetBinContent(iBin);
			double absErrTempProduct = TMath::Sqrt(absErrProductC * absErrProductC * CHistos[iSub]->GetBinContent(iBin) * CHistos[iSub]->GetBinContent(iBin) + CHistos[iSub]->GetBinError(iBin) * CHistos[iSub]->GetBinError(iBin) * productC * productC);

			//One of the terms in the cumulant corresponds to the product of the subset cumulant and the complementary subset cumulant. The remaining terms are further partitions which are counted when recursed
			comp = getComplementarySubset(bin, Subset, len, lenCompSub);
			compSubset = getComplementarySubset(bin, iHisto, len, lenCompSub);
			int iComp = convert(compSubset, lenCompSub);

			double tempProduct2 = productC * CHistos[iComp]->GetBinContent(iBin);
			double absErrTempProduct2 = TMath::Sqrt(absErrProductC * absErrProductC * CHistos[iComp]->GetBinContent(iBin) * CHistos[iComp]->GetBinContent(iBin) + CHistos[iComp]->GetBinError(iBin) * CHistos[iComp]->GetBinError(iBin) * productC * productC);

			int tempProduct3 = productS * getSubsetNumber(subset, lenSub, set, lenSet);
			int tempProduct4 = tempProduct3 * getSubsetNumber(comp, lenCompSub, set, lenSet);

			//check to make sure that no terms are counted more than once( ex. C(1) * C(2,3) can be counted for the subset {1} and {2,3} )
			bool accept = true;
			for (int i = 0; i < curInd; i++)
			{
				accept = accept && !(used[i] == tempProduct4);
			}
			if (accept)
			{
				sum += tempProduct2;
				sumError = TMath::Sqrt(sumError * sumError + absErrTempProduct2 * absErrTempProduct2);
				used[curInd] = tempProduct4;
				curInd++;
			}

			delete[] subset;
			delete[] s;
			//partition the complementary subset into smaller subsets and multiply them by tempProduct and add to the sum
			calcRecSum(CHistos, compSubset, iBin, comp, lenCompSub, set, lenSet, tempProduct, used, curInd, tempProduct3, sum, sumError, absErrTempProduct, depth + 1);

			delete[] comp;
			delete[] compSubset;
		}
	}

	if (reportDebug())
		cout << "PTHistos::calcRecSum(...) Completed." << endl;
	return;
}

//////////////////////////////////////////////////////////////////////////////
//Use the pT moments, stored from each event to calculate the pT deviation moments
//Note this is done on a bin by bin and event by event basis
/////////////////////////////////////////////////////////////////////////////
void PTHistos::calculatePTDeviationMoments(double ***transverseMomentumMoments, int bin, int iEvent, int nParticles, TProfile **pTHisto)
{
	if (reportDebug())
		cout << "PTHistos::calculatePTDeviationMoments(...) Starting." << endl;
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());

	double *n1 = new double[ac.numTypes]();
	double *n2 = new double[ac.numTypes]();
	double *n3 = new double[ac.numTypes]();
	double *n4 = new double[ac.numTypes]();
	int counter = 0;
	double *tempSValues = new double[size]();
	for (int iFilter = 0; iFilter < ac.numTypes; iFilter++)
	{
		//Use the pT moments, stored from each event to calculate the pT deviation moments
		//These equations come simply from expanding the deviation using binomial theorem, and then splitting the sum into separate sums
		nParticles = yields[iEvent][reorder2[iFilter]];
		n1[iFilter] = transverseMomentumMoments[iEvent][iFilter][0] - nParticles * pTHisto[iFilter]->GetBinContent(bin);
		if (maxOrder > 1)
			n2[iFilter] = transverseMomentumMoments[iEvent][iFilter][1] - 2 * transverseMomentumMoments[iEvent][iFilter][0] * pTHisto[iFilter]->GetBinContent(bin) + nParticles * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin);
		if (maxOrder > 2)
			n3[iFilter] = transverseMomentumMoments[iEvent][iFilter][2] - 3 * transverseMomentumMoments[iEvent][iFilter][1] * pTHisto[iFilter]->GetBinContent(bin) + 3 * transverseMomentumMoments[iEvent][iFilter][0] * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin) - nParticles * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin);
		if (maxOrder > 3)
			n4[iFilter] = transverseMomentumMoments[iEvent][iFilter][3] - 4 * transverseMomentumMoments[iEvent][iFilter][2] * pTHisto[iFilter]->GetBinContent(bin) + 6 * transverseMomentumMoments[iEvent][iFilter][1] * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin) - 4 * transverseMomentumMoments[iEvent][iFilter][0] * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin) + nParticles * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin) * pTHisto[iFilter]->GetBinContent(bin);
	}

	for (int iFilter1 = 0; iFilter1 < ac.numTypes; iFilter1++)
	{
		tempSValues[counter] = n1[iFilter1];
		counter++;
	}

	if (maxOrder > 1)
		for (int iFilter1 = 0; iFilter1 < ac.numTypes; iFilter1++)
		{
			for (int iFilter2 = iFilter1; iFilter2 < ac.numTypes; iFilter2++)
			{
				double same12 = iFilter2 == iFilter1 ? n2[iFilter1] : 0;
				tempSValues[counter] = n1[iFilter1] * n1[iFilter2] - same12; //subtract out overcounting the same particle Note: These formulas are well known as Newton sums, to convert moments((x+y+z), (x^3 + y^3 + z^3)) to cyclic sums((x + y + z), (xy + xz + yz))
				counter++;
			}
		}

	if (maxOrder > 2)
		for (int iFilter1 = 0; iFilter1 < ac.numTypes; iFilter1++)
		{
			for (int iFilter2 = iFilter1; iFilter2 < ac.numTypes; iFilter2++)
			{
				for (int iFilter3 = iFilter2; iFilter3 < ac.numTypes; iFilter3++)
				{
					if (iFilter1 == iFilter2)
					{
						if (iFilter1 == iFilter3)
							tempSValues[counter] = n1[iFilter1] * n1[iFilter1] * n1[iFilter1] - 3 * n2[iFilter1] * n1[iFilter1] + 2 * n3[iFilter1];
						else
							tempSValues[counter] = (n1[iFilter1] * n1[iFilter1] - n2[iFilter1]) * n1[iFilter3];
					}
					else
					{
						if (iFilter2 == iFilter3)
							tempSValues[counter] = (n1[iFilter2] * n1[iFilter2] - n2[iFilter2]) * n1[iFilter1];
						else
							tempSValues[counter] = n1[iFilter1] * n1[iFilter2] * n1[iFilter3];
					}
					counter++;
				}
			}
		}

	if (maxOrder > 3)
		for (int iFilter1 = 0; iFilter1 < ac.numTypes; iFilter1++)
		{
			for (int iFilter2 = iFilter1; iFilter2 < ac.numTypes; iFilter2++)
			{
				for (int iFilter3 = iFilter2; iFilter3 < ac.numTypes; iFilter3++)
				{
					for (int iFilter4 = iFilter3; iFilter4 < ac.numTypes; iFilter4++)
					{
						if (iFilter1 == iFilter2)
						{
							if (iFilter1 == iFilter3)
							{
								if (iFilter1 == iFilter4)
									tempSValues[counter] = n1[iFilter2] * n1[iFilter1] * n1[iFilter1] * n1[iFilter1] - 6 * n1[iFilter1] * n1[iFilter1] * n2[iFilter1] + 8 * n1[iFilter1] * n3[iFilter1] + 3 * n2[iFilter1] * n2[iFilter1] - 6 * n4[iFilter1];
								else
									tempSValues[counter] = (n1[iFilter1] * n1[iFilter1] * n1[iFilter1] - 3 * n2[iFilter1] * n1[iFilter1] + 2 * n3[iFilter1]) * n1[iFilter4];
							}
							else
							{
								if (iFilter3 == iFilter4)
									tempSValues[counter] = (n1[iFilter1] * n1[iFilter1] - n2[iFilter1]) * (n1[iFilter3] * n1[iFilter3] - n2[iFilter3]);
								else
									tempSValues[counter] = (n1[iFilter1] * n1[iFilter1] - n2[iFilter1]) * n1[iFilter3] * n1[iFilter4];
							}
						}
						else
						{
							if (iFilter2 == iFilter3)
							{
								if (iFilter2 == iFilter4)
									tempSValues[counter] = (n1[iFilter2] * n1[iFilter2] * n1[iFilter2] - 3 * n2[iFilter2] * n1[iFilter2] + 2 * n3[iFilter2]) * n1[iFilter1];
								else
									tempSValues[counter] = (n1[iFilter2] * n1[iFilter2] - n2[iFilter2]) * n1[iFilter1] * n1[iFilter4];
							}
							else
							{
								if (iFilter3 == iFilter4)
									tempSValues[counter] = (n1[iFilter3] * n1[iFilter3] - n2[iFilter3]) * n1[iFilter1] * n1[iFilter2];
								else
									tempSValues[counter] = n1[iFilter1] * n1[iFilter2] * n1[iFilter3] * n1[iFilter4];
							}
						}
						counter++;
					}
				}
			}
		}
	delete[] n1;
	if (maxOrder > 1)
		delete[] n2;
	if (maxOrder > 2)
		delete[] n3;
	if (maxOrder > 3)
		delete[] n4;

	for (int iHisto = 0; iHisto < size; iHisto++)
	{
		SValues[iEvent][iHisto] = tempSValues[reorder[iHisto]];
	}

	delete[] tempSValues;

	if (reportDebug())
		cout << "PTHistos::calculatePTDeviationMoments(...) Completed." << endl;
}
////////////////////////////////////////////////////////////////////////////
//Helper Functions
////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////
// convert a base ac->numTypes integer (represented by num) into the index of the corresponding moment in the array
// ex: for 4th order correlations, the function {1,2,3} is index 19 in the array (1, 2, 3, 4, 11, 12 ... 44, 111, ... 123)
// checked for correctness
///////////////////////////////////////
int PTHistos::convert(int *num, int len)
{
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	int convert = (TMath::Factorial(ac.numTypes + len - 1)) / (TMath::Factorial(ac.numTypes) * TMath::Factorial(len - 1)) - 1;
	for (int i = 0; i < len; i++)
	{
		for (int j = (i == 0) ? 0 : num[i - 1] - 1; j < num[i] - 1; j++)
		{
			convert += (TMath::Factorial(ac.numTypes + len - j - i - 2)) / (TMath::Factorial(ac.numTypes - j - 1) * TMath::Factorial(len - i - 1));
		}
	}
	return convert;
}

//////////////////////////////////////
// convert the index of the moment (represented by num) into a base ac->numTypes integer
// ex: for 4th order correlations, the function {1,2,3} is index 19 in the array (1, 2, 3, 4, 11, 12 ... 44, 111, ... 123)
//////////////////////////////////////
int *PTHistos::convert(int num, int &len)
{
	TransverseMomentumConfiguration &ac = (TransverseMomentumConfiguration &)(*getConfiguration());
	if (reportDebug())
		cout << "PTHistos::convert(...) started." << endl;
	for (len = 2; len <= maxOrder + 1; len++)
	{
		int temp = (TMath::Factorial(ac.numTypes + len - 1)) / (TMath::Factorial(ac.numTypes) * TMath::Factorial(len - 1)) - 1;
		if (num < temp)
		{
			len--;
			num -= (TMath::Factorial(ac.numTypes + len - 1)) / (TMath::Factorial(ac.numTypes) * TMath::Factorial(len - 1)) - 1;
			break;
		}
	}

	int *convert = new int[len];

	for (int i = 0; i < len; i++)
	{
		int temp = 0;
		for (convert[i] = (i == 0) ? 1 : convert[i - 1]; convert[i] <= ac.numTypes; convert[i]++)
		{
			int temp2 = (TMath::Factorial(ac.numTypes + len - (convert[i] - 1) - i - 2)) / (TMath::Factorial(ac.numTypes - (convert[i] - 1) - 1) * TMath::Factorial(len - i - 1));
			temp += temp2;
			if (num < temp)
			{
				num -= (temp - temp2);
				break;
			}
		}
	}
	if (reportDebug())
		cout << "PTHistos::convert(...) Completed." << endl;
	return convert;
}

//////////////////////////////////////////////
//convert num to a binary string of length len (including leading 0's)
// checked for correctness
////////////////////////////////////////////////
void PTHistos::convertToBinary(int num, char *str, int len)
{
	for (; len > 0; len--)
	{
		str[len - 1] = ((num % 2 == 1) ? '1' : '0');
		num /= 2;
	}
}

///////////////////////////////////////
// get the subset corresponding to the binary string bin, from set
// ex: the subset 001101 of [1, 1, 2, 3, 3, 3] is [2, 3, 3]
// note that because set is guaranteed to be ordered, so is subset
//////////////////////////////////////
int *PTHistos::getSubset(char *subset, int *set, int len, int &lenSub)
{
	lenSub = 0;
	for (int i = 0; i < len; i++)
	{
		if (subset[i] == '1')
			lenSub++;
	}

	int *temp = new int[lenSub];

	lenSub = 0;
	for (int i = 0; i < len; i++)
	{
		if (subset[i] == '1')
		{
			temp[lenSub] = set[i];
			lenSub++;
		}
	}
	return temp;
}

///////////////////////////////////////
// get the complementary subset corresponding to the binary string bin, from set
// ex: the complementary subset 001101 of [1, 1, 2, 3, 3, 3] is [1, 1, 3]
// note that because set is guaranteed to be ordered, so is the complementary subset
//////////////////////////////////////
int *PTHistos::getComplementarySubset(char *subset, int *set, int len, int &lenSub)
{
	lenSub = 0;
	for (int i = 0; i < len; i++)
	{
		if (subset[i] == '0')
			lenSub++;
	}

	int *temp = new int[lenSub];

	lenSub = 0;
	for (int i = 0; i < len; i++)
	{
		if (subset[i] == '0')
		{
			temp[lenSub] = set[i];
			lenSub++;
		}
	}
	return temp;
}

////////////////////////////////////////////////////////////
// get the subset number of subset in mainset
// assumes all elements of mainset are distinct
// checked for correctness
////////////////////////////////////////////////////////////
int PTHistos::getSubsetNumber(int *subset, int lenSub, int *mainset, int lenSet)
{
	int num = 0;
	for (; lenSub > 0; lenSub--)
	{
		for (int i = 0; i < lenSet; i++)
		{
			if (subset[lenSub - 1] == mainset[i])
				num += pow(2, lenSet - i - 1);
		}
	}
	return num;
}
