// energy_calibration_bic.C
// Applies geom‐based calibration constants to new data, produces calibrated energy tree
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
using namespace std;

// Hash functor for pair<int,int>
struct pair_hash {
  std::size_t operator()(const std::pair<int,int>& p) const noexcept {
    // Combine the two ints
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

int energy_calibration_bic(const char* dataFile = "Waveform_sample.root",
                           const char* calibRoot = "calibration_bic_output.root",
                           const char* dataMapFile = "data_map_sample.txt",
                           const char* outFile   = "energy_calibrated.root")
{
  // Load calibration constants from ROOT file
  unordered_map<int,double> calib;
  TFile *fCalib = TFile::Open(calibRoot, "READ");
  if (!fCalib || fCalib->IsZombie()) {
    cerr<<"Cannot open calibration ROOT "<<calibRoot<<endl;
    return 1;
  }
  TTree *tCal = (TTree*)fCalib->Get("Calibration");
  if (!tCal) {
    cerr<<"Calibration TTree not found in "<<calibRoot<<endl;
    return 1;
  }
  Int_t geomID;
  Double_t C;
  tCal->SetBranchAddress("GeomID",     &geomID);
  tCal->SetBranchAddress("CalibConst", &C);
  Long64_t nCal = tCal->GetEntries();
  for (Long64_t i=0; i<nCal; ++i) {
    tCal->GetEntry(i);
    calib[geomID] = C;
  }
  fCalib->Close();
  cout<<"Loaded "<<calib.size()<<" calibration constants from "<<calibRoot<<endl;

  // Load channel -> geomID mapping
  unordered_map<pair<int,int>, int, pair_hash> channelToGeom;
  ifstream inMap(dataMapFile);
  if (!inMap) { cerr<<"Cannot open data mapping file "<<dataMapFile<<endl; return 1; }
  string line;
  while (getline(inMap,line)) {
    if (line.empty() || line[0]=='#') continue;
    istringstream ss(line);
    int mid0, ch0, module0;
    string geomLabel;
    ss >> mid0 >> ch0 >> module0 >> geomLabel;
    channelToGeom[{mid0,ch0}] = module0;
  }
  inMap.close();
  cout<<"Loaded "<<channelToGeom.size()<<" channel-to-geom mappings"<<endl;

  // Open new data file
  TFile *fData = TFile::Open(dataFile,"READ");
  if (!fData || fData->IsZombie()) { cerr<<"Cannot open data "<<dataFile<<endl; return 1; }
  TTree *tIn = (TTree*)fData->Get("T");
  if (!tIn) { cerr<<"TTree T not found in "<<dataFile<<endl; return 1; }

  // Branches
  vector<int> *MID = nullptr, *ch = nullptr;
  vector<int> *idx = nullptr;
  vector<short> *wave = nullptr;
  tIn->SetBranchAddress("MID", &MID);
  tIn->SetBranchAddress("ch",  &ch);
  tIn->SetBranchAddress("waveform_idx",   &idx);
  tIn->SetBranchAddress("waveform_total", &wave);

  // Output file & tree
  TFile *fOut = TFile::Open(outFile,"RECREATE");
  TTree *tOut = tIn->CloneTree(0);
  double calibratedEnergy = 0;
  tOut->Branch("calibratedEnergy", &calibratedEnergy, "calibratedEnergy/D");

  Long64_t nEntries = tIn->GetEntries();
  cout<<"Entries to calibrate: "<<nEntries<<endl;
  for (Long64_t i=0; i<nEntries; ++i) {
    tIn->GetEntry(i);
    // sum ADC per event
    calibratedEnergy = 0;
    int nCh = MID->size();
    for (int j=0; j<nCh; ++j) {
      int mid0 = MID->at(j), ch0 = ch->at(j);
      auto itMap = channelToGeom.find({mid0,ch0});
      if (itMap == channelToGeom.end()) {
        cerr<<"Warning: no geom mapping for MID "<<mid0<<" CH "<<ch0<<" entry "<<i<<endl;
        continue;
      }
      int geom = itMap->second;
      auto itC = calib.find(geom);
      if (itC == calib.end()) {
        cerr<<"Warning: no calib for GeomID "<<geom<<" at entry "<<i<<endl;
        continue;
      }
      double C = itC->second;
      int start = idx->at(j);
      int end   = (j+1<nCh? idx->at(j+1): wave->size());
      double sumADC=0;
      for(int k=start; k<end; ++k) sumADC += (*wave)[k];
      calibratedEnergy += sumADC * C;
    }
    // debug
    if (i%10000==0) cout<<"Calibrated entry "<<i<<" energy="<<calibratedEnergy<<endl;
    tOut->Fill();
  }
  fOut->Write();
  fOut->Close();
  fData->Close();
  cout<<"Wrote energy-calibrated file "<<outFile<<endl;
  return 0;
}