#include "PTHistos.hpp"

ClassImp(PTHistos);

PTHistos::PTHistos(const TString & name,
	AnalysisConfiguration * configuration,
	LogLevel  debugLevel,
	int ord)
:
Histograms(name,configuration,(TMath::Factorial(2 * ord) / TMath::Factorial(ord )),debugLevel),
maxOrder(ord),
histoIndex(0),
size(0),
numFunc(3)
{
	initialize();
}

PTHistos::PTHistos(TFile * inputFile,
	const TString & name,
	AnalysisConfiguration * configuration,
	LogLevel  debugLevel,
	int ord)
:
Histograms(name,configuration,(TMath::Factorial(2 * ord) / TMath::Factorial(maxOrder )),debugLevel),
maxOrder(ord),
histoIndex(0),
size(0),
numFunc(3)
{
	loadHistograms(inputFile);
}

PTHistos::~PTHistos()
{
	//deleteHistograms();
}

// for now use the same boundaries for eta and y histogram
void PTHistos::createHistograms()
{
	if (reportDebug())  cout << "PTHistos::createHistograms(...) started"<< endl;
	AnalysisConfiguration & ac = *getConfiguration();
	TString bn = getHistoBaseName();

	// ================================================================================
	// Naming convention
	// ================================================================================
  	// S is the pT deviation moments
  	// s are the normalized moments
  	// s* are the moments normalized by average pT's
 	// C is the cumulants
	// c are the normalized cumulants
 	// c* are the cumulants normalizd by average pT's

	size = (TMath::Factorial(2 * maxOrder)) / (TMath::Factorial(maxOrder ) * TMath::Factorial(maxOrder )) - 1;

	hS = new TProfile ** [numFunc];
	hS_vsMult = new TProfile ** [numFunc];
	hS_vsCent = new TProfile ** [numFunc];

	hC = new TH1 ** [numFunc];
	hC_vsMult = new TH1 ** [numFunc];
	hC_vsCent = new TH1 ** [numFunc];

	for(int i = 0; i < numFunc; i++)
	{	
		hS[i] = new TProfile * [size];
		hS_vsMult[i] = new TProfile * [size];
		hS_vsCent[i] = new TProfile * [size];

		hC[i] = new TH1 * [size];
		hC_vsMult[i] = new TH1 * [size];
		hC_vsCent[i] = new TH1 * [size];
	}

	h_counts = new TProfile*  [size];
	h_counts_vsMult = new TProfile * [size];
	h_counts_vsCent = new TProfile * [size];


	orders = new int [size];

	h_events   = createHistogram(bn+TString("Nevents"),1,ac.min_mult,  ac.max_mult,  "mult","n_{Events}");
	if (ac.ptCorrelatorVsMult) h_events_vsMult = createHistogram(bn+TString("Nevents_vsMult"),ac.nBins_mult,ac.min_mult,  ac.max_mult,  "mult","n_{Events}");
	if (ac.ptCorrelatorVsCent) h_events_vsCent = createHistogram(bn+TString("Nevents_vsCent"),ac.nBins_cent,ac.min_cent,  ac.max_cent,  "cent","n_{Events}");

	TString * baseName = new TString[2 * numFunc + 1];
	baseName[0] = bn + "S_";
	baseName[1] = bn + "s_";
	baseName[2] = bn + "s*_";
	baseName[3] = bn + "C_";
	baseName[4] = bn + "c_" ;
	baseName[5] = bn + "c*_";
	baseName[6] = bn + "Counts_";

	TString * baseTitle = new TString[2 * numFunc + 1];
	baseTitle[0] = "S_";
	baseTitle[1] = "s_";
	baseTitle[2] = "s*_";
	baseTitle[3] = "C_";
	baseTitle[4] = "c_";
	baseTitle[5] = "c*_";
	baseTitle[6] = "Counts_";

	histoIndex = 0;
	createHistogramRec(baseName, baseTitle, maxOrder - 1, 0);


	histoIndex = 0;
	if (reportDebug())  cout << "PTHistos::createHistograms(...) ended"<< endl;
	//h_c123vsMultTest = createProfile("c123Test", ac.nBins_mult,ac.min_mult,ac.max_mult,"mult", "c123Test" );
}

/////////////////////////////////////////////////////////
// needs to be fixed
////////////////////////////////////////////////////////
void PTHistos::loadHistograms(TFile * inputFile)
{
	if (!inputFile)
	{
		if (reportFatal()) cout << "-Fatal- Attempting to load NuDynHistos from an invalid file pointer" << endl;
		return;
	}
	AnalysisConfiguration & ac = *getConfiguration();

	TString  bn = getHistoBaseName();

	size = (TMath::Factorial(2 * maxOrder)) / (TMath::Factorial(maxOrder ) * TMath::Factorial(maxOrder )) - 1;

	hS = new TProfile ** [numFunc];
	hS_vsMult = new TProfile ** [numFunc];
	hS_vsCent = new TProfile ** [numFunc];

	hC = new TH1 ** [numFunc];
	hC_vsMult = new TH1 ** [numFunc];
	hC_vsCent = new TH1 ** [numFunc];

	for(int i = 0; i < numFunc; i++)
	{	
		hS[i] = new TProfile * [size];
		hS_vsMult[i] = new TProfile * [size];
		hS_vsCent[i] = new TProfile * [size];

		hC[i] = new TH1 * [size];
		hC_vsMult[i] = new TH1 * [size];
		hC_vsCent[i] = new TH1 * [size];
	}

	h_counts = new TProfile*  [size];
	h_counts_vsMult = new TProfile * [size];
	h_counts_vsCent = new TProfile * [size];


	orders = new int [size];

	h_events   = loadH1(inputFile, bn+TString("Nevents"));
	if (ac.ptCorrelatorVsMult) h_events_vsMult = loadH1(inputFile, bn+TString("Nevents_vsMult"));
	if (ac.ptCorrelatorVsCent) h_events_vsCent = loadH1(inputFile, bn+TString("Nevents_vsCent"));

	TString * baseName = new TString[numFunc + 1];
	baseName[0] = bn + "S_";
	baseName[1] = bn + "s_";
	baseName[2] = bn + "s*_";
	baseName[3] = bn + "C_";
	baseName[4] = bn + "c_" ;
	baseName[5] = bn + "c*_";
	baseName[6] = bn + "Counts_";

	TString * baseTitle = new TString[numFunc + 1];
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

//overloaded saveHistograms to save histograms in sequence of lowest order to highest order
void PTHistos::saveHistograms(TFile * outputFile, bool saveAll)
{

	if (reportDebug()) cout << "HistogramCollection::saveHistograms(TFile * outputFile) started."  << endl;
	outputFile->cd();

	AnalysisConfiguration & ac = *getConfiguration();

	int numTypes = 1;
	if (ac.ptCorrelatorVsMult) numTypes++;
	if (ac.ptCorrelatorVsCent) numTypes++;

	for (int k=0; k<3; k++)
	{
		if ((isSaved[k] || saveAll)) histograms[k]->Write();
	}

	for(int iFunc = 0; iFunc <= 2 * numFunc; iFunc++)
	{
		for (int i = 1; i <= maxOrder; i++)
		{
			for (int k=numTypes; k<nHistograms; k++)
			{
				int orderIndex = (k - numTypes)/(numTypes * (2 * numFunc + 1));
				int funcIndex = (k - numTypes)/numTypes - ( (2 * numFunc + 1)) * ((k - numTypes)/(numTypes * (2 * numFunc + 1)));
				if ((isSaved[k] || saveAll) && (orders[orderIndex] == i) && ((funcIndex) % (2 * numFunc + 1) == iFunc)) histograms[k]->Write();
			}
		}
	}

	if (reportDebug()) cout << "HistogramCollection::saveHistograms(TFile * outputFile) completed."  << endl;
}


void PTHistos::fillEventHistos(double mult, double cent, double weight)
{
	AnalysisConfiguration & ac1 = *getConfiguration();
	h_events->Fill(mult,weight);
	if (ac1.ptCorrelatorVsMult) h_events_vsMult->Fill(mult,weight);
	if (ac1.ptCorrelatorVsCent) h_events_vsCent->Fill(cent,weight);
}


// recursively create histograms for correlation functions of order 1 - maxOrder
// Note: these histograms are not in sequence of lowest order to highest order. The reording occurs when the histograms are saved with the saveHistograms function.
void PTHistos::createHistogramRec(TString * baseName, TString * baseTitle, int depth, int partIndex)
{
	if (reportDebug())  cout << "PTHistos::createHistogramsRec(...) started"<< endl;
	AnalysisConfiguration & ac = *getConfiguration();
	TString *histoName = new TString[2 * numFunc +1];
	TString *histoTitle= new TString[2 * numFunc +1];


	for(int i = partIndex; i < maxOrder; i++)
	{
		
		for(int iFunc = 0; iFunc < numFunc ; iFunc++)
		{
			histoName[iFunc]  = baseName[iFunc] + (i + 1) ;
			histoTitle[iFunc] = (depth == maxOrder - 1) ? baseTitle[iFunc] + (i + 1) : baseTitle[iFunc] + ", " + (i + 1);
			histoName[iFunc + numFunc]  = baseName[iFunc + numFunc] + (i + 1);
			histoTitle[iFunc + numFunc] = (depth == maxOrder - 1) ? baseTitle[iFunc + numFunc] + (i + 1) : baseTitle[iFunc + numFunc] + ", " + (i + 1);
			
			
			hS[iFunc][histoIndex] = createProfile( histoName[iFunc], 1,ac.min_mult,ac.max_mult,"mult",histoTitle[iFunc] + "}" );
			if (ac.ptCorrelatorVsMult)	hS_vsMult[iFunc][histoIndex] = createProfile(histoName[iFunc] + "_vsMult",ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult", histoTitle[iFunc] + "}");
			if (ac.ptCorrelatorVsCent)	hS_vsCent[iFunc][histoIndex] = createProfile(histoName[iFunc] + "_vsCent",ac.nBins_cent,ac.min_cent,  ac.max_cent, "cent", histoTitle[iFunc] + "}");
			
			
			hC[iFunc][histoIndex] = createHistogram( histoName[iFunc + numFunc], 1,ac.min_mult,ac.max_mult,"mult",histoTitle[iFunc + numFunc] + "}" );
			if (ac.ptCorrelatorVsMult)	hC_vsMult[iFunc][histoIndex] = createHistogram(histoName[iFunc + numFunc] + "_vsMult",ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult", histoTitle[iFunc + numFunc] + "}");
			if (ac.ptCorrelatorVsCent)	hC_vsCent[iFunc][histoIndex] = createHistogram(histoName[iFunc + numFunc] + "_vsCent",ac.nBins_cent,ac.min_cent,  ac.max_cent, "cent", histoTitle[iFunc + numFunc] + "}");
		}

		
		histoName[2 * numFunc]  = baseName[2 * numFunc] + (i + 1);
		histoTitle[2 * numFunc] = (depth == maxOrder - 1) ? baseTitle[2 * numFunc] + (i + 1) : baseTitle[2 * numFunc] + ", " + (i + 1);

		h_counts[histoIndex]         = createProfile(histoName[2 * numFunc], 1,ac.min_mult,ac.max_mult,"mult", histoTitle[2 * numFunc] + "}" );
		if (ac.ptCorrelatorVsMult)	h_counts_vsMult[histoIndex]  = createProfile(histoName[2 * numFunc] + "_vsMult",ac.nBins_mult,ac.min_mult,  ac.max_mult, "mult", histoTitle[2 * numFunc] + "}");
		if (ac.ptCorrelatorVsCent)	h_counts_vsCent[histoIndex]  = createProfile(histoName[2 * numFunc] + "_vsCent",ac.nBins_cent,ac.min_cent,  ac.max_cent, "cent",histoTitle[2 * numFunc] + "}");

		orders[histoIndex] = maxOrder - depth;

		histoIndex++;

		if(depth != 0)	createHistogramRec(histoName, histoTitle , depth - 1, i);
	}

	delete [] histoName;
	delete [] histoTitle;
	if (reportDebug())  cout << "PTHistos::createHistogramsRec(...) ended"<< endl; 
	return;
}

//////////////////////////////////////////////
// need to fix
//////////////////////////////////////////////
void PTHistos::loadHistogramRec(TString * baseName, int depth, int partIndex, TFile * inputFile)
{

	if (reportDebug())  cout << "PTHistos::loadHistogramRec(...) Starting." << endl;
	AnalysisConfiguration & ac = *getConfiguration();
	TString *histoName = new TString[2 * numFunc +1];


	for(int i = partIndex; i < maxOrder; i++)
	{
		
		for(int iFunc = 0; iFunc < numFunc ; iFunc++)
		{
			histoName[iFunc]  = baseName[iFunc] + (i + 1) ;
			histoName[iFunc + numFunc]  = baseName[iFunc + numFunc] + (i + 1);
			
			
			hS[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc] );
			if (ac.ptCorrelatorVsMult)	hS_vsMult[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc] + "_vsMult");
			if (ac.ptCorrelatorVsCent)	hS_vsCent[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc] + "_vsCent");
			
			
			hC[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc + numFunc] );
			if (ac.ptCorrelatorVsMult)	hC_vsMult[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc + numFunc] + "_vsMult");
			if (ac.ptCorrelatorVsCent)	hC_vsCent[iFunc][histoIndex] = loadProfile(inputFile, histoName[iFunc + numFunc] + "_vsCent");
		}

		
		histoName[2 * numFunc]  = baseName[2 * numFunc] + (i + 1);

		h_counts[histoIndex]         = loadProfile(inputFile, histoName[2 * numFunc] + "}" );
		if (ac.ptCorrelatorVsMult)	h_counts_vsMult[histoIndex]  = loadProfile(inputFile, histoName[2 * numFunc] + "_vsMult");
		if (ac.ptCorrelatorVsCent)	h_counts_vsCent[histoIndex]  = loadProfile(inputFile, histoName[2 * numFunc] + "_vsCent");

		orders[histoIndex] = maxOrder - depth;

		histoIndex++;

		if(depth != 0)	loadHistogramRec(histoName , depth - 1, i, inputFile);
	}

	delete [] histoName;
	return;
	if (reportDebug())  cout << "PTHistos::loadHistogramRec(...) Completed." << endl;
}

void PTHistos::fillDerivedHistos(bool *** acceptances, double * mults, double * cents, double * avgCounts, double * avgpT, double ** SValues, int ** counts, int totEvents)
{
	if (reportDebug())  cout << "PTHistos::fillDerivedHistos(...) Starting." << endl;
	auto start = chrono::high_resolution_clock::now(); 
	AnalysisConfiguration & ac = *getConfiguration();

	//fill SValues
	//fill SValues normalized by counts and average pT's
	for(int iEvent = 0; iEvent < totEvents;iEvent++)
	{
		for(int iHisto = 0; iHisto <size; iHisto++)
		{
			if(counts[iEvent][iHisto] != 0)
			{
				hS[0][iHisto]->Fill(mults[iEvent], SValues[iEvent][iHisto], 1.0);//weight is always 1.0
				if (ac.ptCorrelatorVsMult)	hS_vsMult[0][iHisto]->Fill(mults[iEvent], SValues[iEvent][iHisto] , 1.0);
				if (ac.ptCorrelatorVsCent)	hS_vsCent[0][iHisto]->Fill(cents[iEvent], SValues[iEvent][iHisto] , 1.0);

				hS[1][iHisto]->Fill(mults[iEvent], (SValues[iEvent][iHisto] / avgCounts[iHisto]), 1.0);
				if (ac.ptCorrelatorVsMult)  hS_vsMult[1][iHisto]->Fill(mults[iEvent], (SValues[iEvent][iHisto] / avgCounts[iHisto]), 1.0);
				if (ac.ptCorrelatorVsCent)	hS_vsCent[1][iHisto]->Fill(cents[iEvent], (SValues[iEvent][iHisto] / avgCounts[iHisto]), 1.0);
			}

			h_counts[iHisto]->Fill(mults[iEvent], counts[iEvent][iHisto], 1.0);//weight is always 1.0
			if (ac.ptCorrelatorVsMult)	h_counts_vsMult[iHisto]->Fill(mults[iEvent], counts[iEvent][iHisto] , 1.0);
			if (ac.ptCorrelatorVsCent)	h_counts_vsCent[iHisto]->Fill(cents[iEvent], counts[iEvent][iHisto] , 1.0);
		}
		histoIndex = 0;
		fillNormalizedPTValues(maxOrder - 1, 0, 1,  SValues[iEvent], mults[iEvent], cents[iEvent], avgpT);
		histoIndex = 0;
	}

	TH1 *** newhCValues = new TH1**[3];
	for(int iHisto = 0; iHisto < 3; iHisto++)
	{
		newhCValues[iHisto] = new TH1*[size];
	}

	reorder = new int[size];

	calculateCumulants(hS[0], newhCValues[0], 1 ,ac.min_mult,  ac.max_mult);
	if (ac.ptCorrelatorVsMult)	calculateCumulants(hS_vsMult[0], newhCValues[1], ac.nBins_mult,ac.min_mult,  ac.max_mult);
	if (ac.ptCorrelatorVsCent)	calculateCumulants(hS_vsCent[0], newhCValues[2], ac.nBins_cent,ac.min_cent,  ac.max_cent);

	for(int iHisto = 0; iHisto <size; iHisto++)
	{
		int nBins = 1;
		for(int iBin = 1; iBin <=nBins; iBin++)
		{
			if(newhCValues[0][reorder[iHisto]]->GetBinContent(iBin) != 0)
			{
				hC[0][iHisto]->SetBinContent(iBin, totEvents * newhCValues[0][reorder[iHisto]]->GetBinContent(iBin));
				hC[0][iHisto]->SetBinError(iBin, totEvents * newhCValues[0][reorder[iHisto]]->GetBinError(iBin));

				hC[1][iHisto]->SetBinContent(iBin,  totEvents * newhCValues[0][reorder[iHisto]]->GetBinContent(iBin)/h_counts[iHisto]->GetBinContent(iBin));
				double err = newhCValues[0][reorder[iHisto]]->GetBinContent(iBin)/h_counts[iHisto]->GetBinContent(iBin) * ( newhCValues[0][reorder[iHisto]]->GetBinError(iBin)/newhCValues[0][reorder[iHisto]]->GetBinContent(iBin) + h_counts[iHisto]->GetBinError(iBin)/h_counts[iHisto]->GetBinContent(iBin) );
				hC[1][iHisto]->SetBinError(iBin, totEvents * err);
			}
		}

		if (ac.ptCorrelatorVsMult)	
		{
			nBins =ac.nBins_mult;
			for(int iBin = 1; iBin <=nBins; iBin++)
			{
				if(newhCValues[1][reorder[iHisto]]->GetBinContent(iBin) != 0)
				{
					hC_vsMult[0][iHisto]->SetBinContent(iBin, totEvents * newhCValues[1][reorder[iHisto]]->GetBinContent(iBin));
					hC_vsMult[0][iHisto]->SetBinError(iBin, totEvents * newhCValues[1][reorder[iHisto]]->GetBinError(iBin));

					hC_vsMult[1][iHisto]->SetBinContent(iBin,  totEvents * newhCValues[1][reorder[iHisto]]->GetBinContent(iBin)/h_counts_vsMult[iHisto]->GetBinContent(iBin));
					double err = newhCValues[1][reorder[iHisto]]->GetBinContent(iBin)/h_counts_vsMult[iHisto]->GetBinContent(iBin) * ( newhCValues[1][reorder[iHisto]]->GetBinError(iBin)/newhCValues[1][reorder[iHisto]]->GetBinContent(iBin) + h_counts_vsMult[iHisto]->GetBinError(iBin)/h_counts_vsMult[iHisto]->GetBinContent(iBin) );
					hC_vsMult[1][iHisto]->SetBinError(iBin, totEvents * err);
				}
			}
		}

		if (ac.ptCorrelatorVsCent)
		{
			nBins = ac.nBins_cent;
			for(int iBin = 1; iBin <=nBins; iBin++)
			{
				if(newhCValues[2][reorder[iHisto]]->GetBinContent(iBin) != 0)
				{
					hC_vsCent[0][iHisto]->SetBinContent(iBin, totEvents * newhCValues[2][reorder[iHisto]]->GetBinContent(iBin));
					hC_vsCent[0][iHisto]->SetBinError(iBin, totEvents * newhCValues[2][reorder[iHisto]]->GetBinError(iBin));

					hC_vsCent[1][iHisto]->SetBinContent(iBin, totEvents * newhCValues[2][reorder[iHisto]]->GetBinContent(iBin)/h_counts_vsCent[iHisto]->GetBinContent(iBin));
					double err = newhCValues[2][reorder[iHisto]]->GetBinContent(iBin)/h_counts_vsCent[iHisto]->GetBinContent(iBin) * ( newhCValues[2][reorder[iHisto]]->GetBinError(iBin)/newhCValues[2][reorder[iHisto]]->GetBinContent(iBin) + h_counts_vsCent[iHisto]->GetBinError(iBin)/h_counts_vsCent[iHisto]->GetBinContent(iBin) );
					hC_vsCent[1][iHisto]->SetBinError(iBin, totEvents * err);
				}
			}
		}

	}

	fillNormalizedPTValues(maxOrder - 1, 0, 1,  newhCValues, reorder, new int[3]{1,ac.nBins_mult, ac.nBins_cent} ,totEvents, avgpT);

	auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
	cout << " Time taken to calculate derived histos: " << duration.count() << " microseconds"<< endl;
	if (reportDebug())  cout << "PTHistos::fillDerivedHistos(...) Completed." << endl;




}

////////////////////////////////////////////
//fill  S or C Values normalized by overall event average transverse momenta
/////////////////////////////////////////////
void PTHistos::fillNormalizedPTValues( int depth, int partIndex, double product, double * values, double  mult, double  cent, double *avgpT)
{
	if (reportDebug())  cout << "PTHistos::fillNormalizedPTValues(...) Starting." << endl;
	AnalysisConfiguration & ac = *getConfiguration();

	for(int i = partIndex; i < maxOrder; i++)
	{

		
		double newProduct = product * avgpT[i];

		if(values[histoIndex] != 0)
		{
			hS[2][histoIndex]->Fill(mult, (values[histoIndex] / newProduct), 1.0);
			if (ac.ptCorrelatorVsMult)	hS_vsMult[2][histoIndex]->Fill(mult, (values[histoIndex] / newProduct), 1.0);
			if (ac.ptCorrelatorVsCent)	hS_vsCent[2][histoIndex]->Fill(cent, (values[histoIndex] / newProduct), 1.0);
		}
		if(histoIndex != size - 1)	histoIndex++;

		if(depth != 0)	fillNormalizedPTValues(depth - 1, i, newProduct, values, mult, cent, avgpT);
	}
	if (reportDebug())  cout << "PTHistos::fillNormalizedPTValues(...) Completed." << endl;
}

void PTHistos::fillNormalizedPTValues( int depth, int partIndex, double product, TH1 *** values, int* reorder, int*  nBin, int  totEvents, double *avgpT)
{
	if (reportDebug())  cout << "PTHistos::fillNormalizedPTValues(...) Starting." << endl;
	AnalysisConfiguration & ac = *getConfiguration();

	for(int i = partIndex; i < maxOrder; i++)
	{

		
		double newProduct = product * avgpT[i];

		for(int iBin = 1; iBin <=nBin[0]; iBin++)
		{
			if(values[0][reorder[histoIndex]]->GetBinContent(iBin) != 0)
			{
				hC[2][histoIndex]->SetBinContent(iBin, totEvents * values[0][reorder[histoIndex]]->GetBinContent(iBin) / newProduct);
				double err = values[0][reorder[histoIndex]]->GetBinContent(iBin)/newProduct * ( values[0][reorder[histoIndex]]->GetBinError(iBin)/values[0][reorder[histoIndex]]->GetBinContent(iBin));
				hC[2][histoIndex]->SetBinError(iBin, totEvents * err);
			}
		}
		if (ac.ptCorrelatorVsMult)	
		{
			for(int iBin = 1; iBin <=nBin[1]; iBin++)
			{
				if(values[1][reorder[histoIndex]]->GetBinContent(iBin) != 0)
				{
					hC_vsMult[2][histoIndex]->SetBinContent(iBin, totEvents * values[1][reorder[histoIndex]]->GetBinContent(iBin) / newProduct);
					double err = values[1][reorder[histoIndex]]->GetBinContent(iBin)/newProduct * ( values[1][reorder[histoIndex]]->GetBinError(iBin)/values[1][reorder[histoIndex]]->GetBinContent(iBin));
					hC_vsMult[2][histoIndex]->SetBinError(iBin, totEvents * err);
				}
			}
		}
		if (ac.ptCorrelatorVsCent)	
		{
			for(int iBin = 1; iBin <=nBin[2]; iBin++)
			{
				if(values[2][reorder[histoIndex]]->GetBinContent(iBin) != 0)
				{
					hC_vsCent[2][histoIndex]->SetBinContent(iBin, totEvents * values[2][reorder[histoIndex]]->GetBinContent(iBin) / newProduct);
					double err = values[1][reorder[histoIndex]]->GetBinContent(iBin)/newProduct * ( values[2][reorder[histoIndex]]->GetBinError(iBin)/values[2][reorder[histoIndex]]->GetBinContent(iBin));
					hC_vsCent[2][histoIndex]->SetBinError(iBin, totEvents * err);
				}
			}
		}
		
		if(histoIndex != size - 1)	histoIndex++;

		if(depth != 0)	fillNormalizedPTValues(depth - 1, i, newProduct, values, reorder, nBin, totEvents, avgpT);
	}
	if (reportDebug())  cout << "PTHistos::fillNormalizedPTValues(...) Completed." << endl;
}

///////////////////////////////////////////
//calculate the cumulants of the moments
//////////////////////////////////////////
void PTHistos::calculateCumulants(TProfile ** SHistos, TH1 **CHistos, int nBins, double min, double max)
{	
	AnalysisConfiguration & ac = *getConfiguration();

	TProfile ** newSHistos = new TProfile *[size];
	int counter = 0;
	for(int iOrd = 1; iOrd <= maxOrder; iOrd++)
	{
		for(int iHisto = 0; iHisto < size; iHisto++)
		{
			if(orders[iHisto] == iOrd)
			{
				newSHistos[counter] = SHistos[iHisto];
				reorder[iHisto] = counter;
				counter++;

			}
		}
	}

	int maxSize = TMath::Factorial( maxOrder);
	double * used = new double[maxSize];
	for(int iBin = 1; iBin <=nBins; iBin++)
	{
		for(int iHisto = 0; iHisto < size; iHisto++)
		{
			double sum = 0;
			double absESq = 0;
			
			int len = 0;
			int * ind = convert(iHisto, len);

			int curInd = 0;
			int * set = new int[len];
			for(int i = 0; i < len; i++)
			{
				set[i] = i+1;
			}

			calcRecSum( CHistos, iBin, absESq, 0, ind, set, len, set, len, 1, used, curInd, 1, sum );

			TString name = nBins * size + iHisto;
			
			if(iBin == 1) CHistos[iHisto] = new TH1F(name, name, nBins, min, max);

			double content = newSHistos[iHisto]->GetBinContent(iBin) - sum;
			double error = sqrt(newSHistos[iHisto]->GetBinError(iBin) * newSHistos[iHisto]->GetBinError(iBin) + absESq);

			if(newSHistos[iHisto]->GetBinContent(iBin) != 0)
			{
				CHistos[iHisto]->SetBinContent(iBin, content); 
				CHistos[iHisto]->SetBinError(iBin, error);
			}			
			//if (abs(CHistos[iHisto]->GetBinContent(iBin)) < pow(10, -15)) CHistos[iHisto]->SetBinContent(iBin, 0); // doubles are accurate to 15 decimal places

			//ind is now uneccesary
			delete[] ind;

		}
	}
	
	return;
	if (reportDebug())  cout << "PTHistos::calculateCumulants(...) Completed." << endl;

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
void PTHistos::calcRecSum(TH1 **CHistos, int iBin, double& absESq, double curRelESq, int* iHisto, int* Subset, int len,  int * set, int lenSet, double productC, double* used, int& curInd, int productS, double& sum)
{
	if (reportDebug())  cout << "PTHistos::calcRecSum(...) Starting." << endl;
	int lenSub = 0;
	int lenCompSub = 0;
	if(len != 1)
	{
		char* bin = new char[len];
		int* subset;
		int * comp;
		int * compSubset;

		// loop over all possible proper subsets
		for(int iSubset = 1; iSubset < pow(2, len) - 1; iSubset++)
		{
			//get the subset corresponding to iSubset, multiply it by the product(which contains cumulant terms from previous recursive levels)
			convertToBinary(iSubset, bin, len);
			subset = getSubset(bin, Subset, len, lenSub);
			int * s = getSubset(bin, iHisto, len, lenSub);
			
			int iSub = convert(s, lenSub);
			double tempProduct = productC * CHistos[iSub]->GetBinContent(iBin);
			double tempRelErr = curRelESq + CHistos[iSub]->GetBinError(iBin) * CHistos[iSub]->GetBinError(iBin) / CHistos[iSub]->GetBinContent(iBin) / CHistos[iSub]->GetBinContent(iBin);

			//One of the terms in the cumulant corresponds to the product of the subset cumulant and the complementary subset cumulant. The remaining terms are further partitions which are counted when recursed
			comp = getComplementarySubset(bin, Subset, len, lenCompSub);
			compSubset = getComplementarySubset(bin, iHisto, len, lenCompSub);
			int iComp = convert(compSubset, lenCompSub);

			double tempProduct2 = tempProduct * CHistos[iComp]->GetBinContent(iBin);
			double tempRelErr2  = tempRelErr  + CHistos[iComp]->GetBinError(iBin) * CHistos[iComp]->GetBinError(iBin) / CHistos[iComp]->GetBinContent(iBin) / CHistos[iComp]->GetBinContent(iBin);

			int tempProduct3 = productS * getSubsetNumber(subset, lenSub, set, lenSet);
			int tempProduct4 = tempProduct3 * getSubsetNumber(comp, lenCompSub, set, lenSet);

			//check to make sure that no terms are counted more than once( ex. C(1) * C(2,3) can be counted for the subset {1} and {2,3} )
			bool accept = true;			
			for(int i = 0; i < curInd; i++)
			{
				accept = accept && !(used[i] == tempProduct4);
			}
			if(accept) 
			{
				sum += tempProduct2;
				absESq += tempProduct2 * tempProduct2 * tempRelErr2;
				used[curInd] = tempProduct4;
				curInd++;
			}

			delete [] subset;
			delete[] s;
			//partition the complementary subset into smaller subsets and multiply them by tempProduct and add to the sum
			calcRecSum(CHistos, absESq, tempRelErr, iBin, compSubset, comp, lenCompSub, set, lenSet, tempProduct, used, curInd, tempProduct3, sum);
			
			delete [] comp;
			delete[] compSubset;



		}
	}

	if (reportDebug())  cout << "PTHistos::calcRecSum(...) Completed." << endl;
	return;
	
}


///////////////////////////////////////
// convert a base maxOrder integer (represented by num) into the index of the corresponding moment in the array
// ex: for 4th order correlations, the function {1,2,3} is index 19 in the array (1, 2, 3, 4, 11, 12 ... 44, 111, ... 123)
// checked for correctness
///////////////////////////////////////
int PTHistos::convert(int * num, int len)
{
	int convert = (TMath::Factorial( maxOrder + len - 1)) / (TMath::Factorial(maxOrder ) * TMath::Factorial(len - 1 )) - 1;
	for(int i = 0; i < len; i++)
	{
		for(int j = (i == 0)? 0: num[i - 1] - 1; j < num[i] - 1; j++)
		{
			convert += (TMath::Factorial( maxOrder + len - j - i - 2)) / (TMath::Factorial(maxOrder - j - 1 ) * TMath::Factorial(len - i - 1 ));
		}
	}
	return convert;
}

//////////////////////////////////////
// convert the index of the moment (represented by num) into a base maxOrder integer
// ex: for 4th order correlations, the function {1,2,3} is index 19 in the array (1, 2, 3, 4, 11, 12 ... 44, 111, ... 123)
//////////////////////////////////////
int * PTHistos::convert(int num, int & len)
{
	if (reportDebug())  cout << "PTHistos::convert(...) started." << endl;
	for(len = 2; len <=maxOrder + 1; len++)
	{
		int temp = (TMath::Factorial( maxOrder + len - 1)) / (TMath::Factorial(maxOrder ) * TMath::Factorial(len - 1 )) - 1 ;
		if(num < temp)
		{
			len--;
			num -= (TMath::Factorial( maxOrder + len - 1)) / (TMath::Factorial(maxOrder ) * TMath::Factorial(len - 1 )) - 1 ;
			break;
		}
	}

	int * convert = new int[len];

	for(int i = 0; i <len; i++)
	{
		int temp = 0;
		for(convert[i] = (i == 0)? 1: convert[i - 1] ; convert[i] <= maxOrder; convert[i]++)
		{
			int temp2 = (TMath::Factorial( maxOrder + len - (convert[i] - 1) - i - 2)) / (TMath::Factorial(maxOrder - (convert[i] - 1) - 1 ) * TMath::Factorial(len - i - 1 ));
			temp += temp2;
			if(num < temp)
			{
				num -=  (temp - temp2);
				break;
			}
		}
	}
	if (reportDebug())  cout << "PTHistos::convert(...) Completed." << endl;
	return convert;
}

//////////////////////////////////////////////
//convert num to a binary string of length len (including leading 0's)
// checked for correctness
////////////////////////////////////////////////
void PTHistos::convertToBinary(int num, char*str, int len )
{
	for(; len > 0; len --)
	{
		str[len - 1] = ((num%2 == 1) ? '1':'0');
		num /= 2;
	}
}	


///////////////////////////////////////
// get the subset corresponding to the binary string bin, from set
// ex: the subset 001101 of [1, 1, 2, 3, 3, 3] is [2, 3, 3]
// note that because set is guaranteed to be ordered, so is subset
//////////////////////////////////////
int * PTHistos::getSubset(char* subset, int * set, int len, int& lenSub)
{
	lenSub = 0;
	for(int i = 0; i < len; i ++)
	{
		if(subset[i] == '1') lenSub++;
	}

	int * temp = new int[lenSub];

	lenSub = 0;
	for(int i = 0; i < len; i ++)
	{
		if(subset[i] == '1') 
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
int * PTHistos::getComplementarySubset(char* subset, int * set, int len, int& lenSub)
{
	lenSub = 0;
	for(int i = 0; i < len; i ++)
	{
		if(subset[i] == '0') lenSub++;
	}

	int * temp = new int[lenSub];

	lenSub = 0;
	for(int i = 0; i < len; i ++)
	{
		if(subset[i] == '0') 
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
int PTHistos::getSubsetNumber(int * subset, int lenSub, int * mainset, int lenSet)
{
	int num = 0;
	for(;lenSub > 0; lenSub--)
	{
		for(int i = 0; i<lenSet; i++)
		{
			if(subset[lenSub - 1] == mainset[i]) num += pow(2, lenSet - i - 1);
		}
	}
	return num ;
}

