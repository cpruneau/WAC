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
#include "AACollisionReader.hpp"
#include "TParticle.h"

using namespace std;

int main()
{
  auto start = chrono::high_resolution_clock::now(); 

  cout << "<INFO> PYTHIA Parallel Event Generation - Starting" << endl;
  long nEventsRequested = 10000000;
  int  nThreads    = 50;        // make sure this divides the number of events per tree
  int eventsPerThread = nEventsRequested/nThreads;
  cout << "<INFO> Events per Thread: " << eventsPerThread << endl;
  TString outputFolder = getenv("OUTPUT_PATH");

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
    TFile** files = new TFile*[nThreads];
    TTree** trees = new TTree*[nThreads];
    TClonesArray** arrays = new TClonesArray*[nThreads];

    ROOT::EnableThreadSafety();
    #pragma omp parallel num_threads(nThreads)
    {

      TFile * file = new TFile(outputFolder + "/PythiaEventTree" + omp_get_thread_num() + ".root", "RECREATE");
      files[omp_get_thread_num()] = file;
      file->cd();

      TTree * tree = new TTree("PythiaEventTree", "PythiaEventTree" );
      trees[omp_get_thread_num()] = tree;
      TString branchName_eventNo = "eventNo";
      TString branchName_mult = "mult";
      TString branchName_px = "px";
      TString branchName_py = "py";
      TString branchName_pz = "pz";
      TString branchName_ist = "ist";
      TString branchName_pdg = "pdg";
      TString branchName_pE = "pE";

      const int nMax = AACollisionReader::arraySize; //30000
      TClonesArray * particles = new TClonesArray("TParticle", nMax );
      arrays[omp_get_thread_num()] = particles;
      double px [nMax] = {};
      double py  [nMax]= {};
      double pz [nMax]= {};
      double pE [nMax]= {};
      int ist [nMax]= {};
      int pdg  [nMax]= {};
      int eventNo, mult;

      tree->Branch(branchName_eventNo, &eventNo);
      tree->Branch(branchName_mult, &mult );
      tree->Branch(branchName_px, px ,"px[30000]/D" );
      tree->Branch(branchName_py, py,"py[30000]/D"  );
      tree->Branch(branchName_pz, pz ,"pz[30000]/D" );
      tree->Branch(branchName_ist, ist,"ist[30000]/I"  );
      tree->Branch(branchName_pdg, pdg,"pdg[30000]/I"  );
      tree->Branch(branchName_pE, pE ,"pE[30000]/D"  );


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

     eventNo = iEvent;
     mult = nparts;
     for(int iParticle = 0; iParticle < nparts; iParticle++)
     {
      TParticle & part = * (TParticle*) particles->At(iParticle);
      ist[iParticle] = part.GetStatusCode();
      pdg[iParticle] = part.GetPdgCode();
      px[iParticle] = part.Px();
      py[iParticle] = part.Py();
      pz[iParticle] = part.Pz();
      pE[iParticle] = part.Energy();
     }

     tree->Fill();
   }

   tree->Write();
   file->Close();
   //delete particles;
   //delete tree;
   //delete file;
 }

/*for(int i = 0; i < nThreads; i++)
{
 delete files[i];
 delete trees[i];
 delete arrays[i];
  }
  delete[] files;
  delete[]trees;
  delete[] arrays;
 delete pythiaOptions;
*/
 cout << "<INFO> Parallel Event Generation - Completed" << endl;
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

