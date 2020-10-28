#include <iostream>
#include <chrono>
#include <fstream>
#include "omp.h"
#include <TStyle.h>
#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include <TPythia8.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TString.h>
#include "PythiaConfiguration.hpp"
#include "TMath.h"

using namespace std;

int main()
{
  auto start = chrono::high_resolution_clock::now(); 

  cout << "<INFO> PYTHIA Model Analysis - Starting" << endl;
  long nEventsRequested = 10000;
  int  nThreads    = 5;        // make sure this divides the number of events
  int eventsPerThread = nEventsRequested/nThreads;
  TString outputFolder = getenv("OUTPUT_PATH");
  TFile * file = new TFile(outputFolder + "/PythiaEventTree.root", "RECREATE");
  TTree * tree = new TTree("EventTree", "EventTree" );

  int nOptions = 0;
  TString ** pythiaOptions = new TString* [50];
  pythiaOptions[nOptions++] = new TString("Init:showChangedSettings = on");      // list changed settings
  pythiaOptions[nOptions++] = new TString("Init:showChangedParticleData = off"); // list changed particle data
  pythiaOptions[nOptions++] = new TString("Next:numberCount = 10000");            // print message every n events
  pythiaOptions[nOptions++] = new TString("Next:numberShowInfo = 1");            // print event information n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowProcess = 0");         // print process record n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowEvent = 0");
  pythiaOptions[nOptions++] = new TString("SoftQCD:all = on");                   // Allow total sigma = elastic/SD/DD/ND
  //pythiaOptions[nOptions++] = new TString("HardQCD:all = on");
  pythiaOptions[nOptions++] = new TString("Random:setSeed = on");
  pythiaOptions[nOptions++] = new TString("Init:showMultipartonInteractions = off"); // don't list multiparton interaction initialization
  pythiaOptions[nOptions++] = new TString("Init:showProcesses = off");            // don't list initialization info
  PythiaConfiguration * pc = new PythiaConfiguration(2212 /* p */,
                                                     2212 /* p */,
                                                     14000.0, /* energy in GeV */
  nOptions,
  pythiaOptions);

  file->cd();

  ROOT::EnableThreadSafety();
  #pragma omp parallel num_threads(nThreads)
  {
    TString branchname = "EventsFromThread";
    int nMax = 10000;
    TClonesArray * particles = new TClonesArray("TParticle", nMax );
    TBranch * b = tree->Branch(branchname + (omp_get_thread_num()+ 1), particles, 32000 , 0 );

    TPythia8 * pythia8 = new TPythia8(false); // don't print initialization
    for (int iOption=0; iOption<pc->nOptions; iOption++)
    {
      pythia8->ReadString( *pc->options[iOption]);
    }
    auto time = chrono::high_resolution_clock::now();
    unsigned int seed = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
    seed = ((seed  * (omp_get_thread_num()+ 1)) % 90000000) * ((seed  * (omp_get_thread_num()+ 1)) % 90000000) * ((seed  * (omp_get_thread_num()+ 1)) % 90000000) * ((seed  * (omp_get_thread_num()+ 1)) % 90000000) * ((seed  * (omp_get_thread_num()+ 1)) % 90000000);
    seed = seed % 90000000;
    cout << "Seed from thread " << omp_get_thread_num()+ 1 << " is " << seed << endl;
    TString s = "Random:seed = ";
    pythia8->ReadString(s+ seed);
    pythia8->Initialize(pc->beam,pc->target,pc->energy);
    for(int iEvent = 0; iEvent < eventsPerThread; iEvent ++)
    {
      int nparts;
      bool seekingEvent = true;
      while (seekingEvent)
      {
       pythia8->GenerateEvent();
       pythia8->ImportParticles(particles,"Final");
       nparts = particles->GetEntriesFast();
       if (nparts>2) seekingEvent = false;
     }
     b->Fill();
   }
 }


 tree->Write();

 file->Close();


 cout << "<INFO> PYTHIA Analysis - Completed" << endl;




  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Timer
  ////////////////////////////////////////////////////////////////////////////////////////////////////




 auto stop = chrono::high_resolution_clock::now(); 
 auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
 int hours = (int)(duration.count()/3600);
 int minutes = (int)((duration.count() - 3600 * hours)/60);
 double seconds = duration.count() - 60 * minutes - 3600 * hours;
 cout << "<INFO> Total Time elapsed "<< (hours) << ":" << (minutes ) << ":" << seconds << endl; 

}

