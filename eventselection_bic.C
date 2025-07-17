// eventselection_bic.C
// Filters events by trigger time/number, writes a new ROOT with same tree structure
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
using namespace std;

int eventselection_bic(const char* inFile  = "Run_60184_Waveform.root",
                       const char* outFile = "Run_60184_Selected.root",
                       double timeMin      = 0,
                       double timeMax      = 1e12,
                       int    numMin       = 0,
                       int    numMax       = INT_MAX)
{
  TFile *fIn = TFile::Open(inFile,"READ");
  if (!fIn||fIn->IsZombie()) { cerr<<"Cannot open "<<inFile<<endl; return 1; }
  TTree *tIn = (TTree*)fIn->Get("T");
  if (!tIn) { cerr<<"TTree T not found in "<<inFile<<endl; return 1; }

  // branches
  vector<long long>* trigger_time   = nullptr;
  vector<int>*       trigger_number = nullptr;
  tIn->SetBranchAddress("trigger_time",   &trigger_time);
  tIn->SetBranchAddress("trigger_number", &trigger_number);

  // clone structure but empty
  TFile *fOut = TFile::Open(outFile,"RECREATE");
  TTree *tOut = tIn->CloneTree(0);

  Long64_t n = tIn->GetEntries();
  cout<<"Total events: "<<n<<endl;
  for (Long64_t i=0;i<n;++i) {
    tIn->GetEntry(i);
    long long tt = (trigger_time && !trigger_time->empty()) ? trigger_time->at(0):0;
    int      tn = (trigger_number&& !trigger_number->empty())? trigger_number->at(0):0;
    // debug
    if (i%10000==0) cout<<"Checking entry "<<i<<" time="<<tt<<" num="<<tn<<endl;
    if (tt<timeMin||tt>timeMax)     continue;
    if (tn<numMin||tn>numMax)       continue;
    tOut->Fill();
  }
  fOut->Write();
  fOut->Close();
  fIn->Close();
  cout<<"Wrote selected events to "<<outFile<<endl;
  return 0;
}