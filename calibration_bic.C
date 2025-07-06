// calibration_bic.C
// Usage:
//   root -l
//   root[0] .L calibration_bic.C
//   root[1] calibration_bic("Run_60184_Waveform.root","3x5_5GeV_result_new.root", "data_map.txt", "sim_map.txt", true, false);

#include <TFile.h>
#include <TH1D.h>
#include <TH1.h>
#include <TTree.h>
#include "TKey.h"
#include "TIterator.h"
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <unordered_map>
#include <utility>
#include <fstream>
#include <sstream>
#include <TTree.h>
#include <fstream>

struct pair_hash {
  size_t operator()(const std::pair<int,int>& p) const noexcept {
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

using namespace std;

void calibration_bic(const char* dataFile = "Waveform_sample.root",
                     const char* simFile  = "3x8_3GeV_CERN_hist.root",
                     const char* dataMapFile = "data_map_sample.txt", // columns: <DAQ MID>, <CH>, <BIC MID L/R>, + <Geometry mapping MID>
                     const char* simMapFile  = "sim_map_sample.txt", // columns: <Sim results MID>, <Geometry mappring MID>
                     bool useTriggerTime    = true,
                     bool useTriggerNumber  = false)
{
  // --- Data loading ---
  TFile *fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    cerr << "Error: cannot open data file " << dataFile << endl;
    return;
  }
  // Auto-detect TTree in data file
  fData->ls();
  TTree *tData = nullptr;
  {
    TIter nextKeyData(fData->GetListOfKeys());
    TKey *keyData;
    while ((keyData = (TKey*)nextKeyData())) {
      TObject *obj = keyData->ReadObj();
      if (obj->InheritsFrom("TTree")) {
        tData = (TTree*)obj;
        cout << "Using data TTree: " << tData->GetName() << endl;
        break;
      }
    }
    if (!tData) {
      cerr << "Error: no TTree found in " << dataFile << endl;
      return;
    }
  }

  // read flattened waveform from event-builder output
  std::vector<short>* waveform_total = nullptr;
  std::vector<int>* waveform_idx   = nullptr;
  std::vector<int>* MID            = nullptr;
  std::vector<int>* ch             = nullptr;
  tData->SetBranchAddress("waveform_total", &waveform_total);
  tData->SetBranchAddress("waveform_idx",   &waveform_idx);
  tData->SetBranchAddress("MID",            &MID);
  tData->SetBranchAddress("ch",             &ch);
  std::vector<long long>* trigger_time   = nullptr;
  std::vector<int>*       trigger_number = nullptr;
  tData->SetBranchAddress("trigger_time",   &trigger_time);
  tData->SetBranchAddress("trigger_number", &trigger_number);

  // --- Read mapping files ---

  // --- Prepare per-geom histograms for data and simulation ---
  std::map<int, TH1D*> hDataDist;
  std::map<int, TH1D*> hSimDist;
  const int nBins = 200;
  const double maxADC = 1e6; 
  // Allocate histograms after reading mapping files
  unordered_map<pair<int,int>, int, pair_hash> channelToGeom;
  {
    ifstream mf(dataMapFile);
    if (!mf) { cerr<<"Cannot open data mapping file "<<dataMapFile<<endl; return; }
    string line;
    while (getline(mf, line)) {
      if (line.empty() || line[0]=='#') continue;
      istringstream ss(line);
      int mid0, ch0, module0, geom0;
      ss >> mid0 >> ch0 >> module0 >> geom0;
      channelToGeom[{mid0, ch0}] = geom0;
    }
  }
  // Print how many channel mappings loaded
  cout << "Loaded " << channelToGeom.size() << " channel-to-geom mappings" << endl;
  unordered_map<string,int> simToGeom;
  {
    ifstream sm(simMapFile);
    if (!sm) { cerr<<"Cannot open sim mapping file "<<simMapFile<<endl; return; }
    string line;
    while (getline(sm, line)) {
      if (line.empty() || line[0]=='#') continue;
      istringstream ss(line);
      string simName;
      int geom0;
      ss >> simName >> geom0;
      simToGeom[simName] = geom0;
    }
  }
  // Print how many simulation mappings loaded
  cout << "Loaded " << simToGeom.size() << " sim-to-geom mappings" << endl;

  // Allocate per-geom data histograms
  for (auto &kv : channelToGeom) {
    int geom = kv.second;
    if (hDataDist.count(geom)==0) {
      hDataDist[geom] = new TH1D(Form("hData_G%d", geom),
                                 Form("Data INT ADC Geom %d;INT ADC;Events", geom),
                                 nBins, 0, maxADC);
      hDataDist[geom]->SetDirectory(0);   // detach from file
    }
  }
  // Allocate per-geom sim histograms
  for (auto &kv : simToGeom) {
    int geom = kv.second;
    if (hSimDist.count(geom)==0) {
      hSimDist[geom] = new TH1D(Form("hSim_G%d", geom),
                                Form("Sim INT ADC Geom %d;INT ADC;Events", geom),
                                nBins, 0, maxADC);
      hSimDist[geom]->SetDirectory(0);   // detach from file
    }
  }

  // Module Accumulation (sum, count)
  map<int, double> sum_data;
  map<int, long>   count_data;

  Long64_t nD = tData->GetEntries();
  cout << "Data entries: " << nD << endl;
  for (Long64_t i = 0; i < nD; ++i) {
    tData->GetEntry(i);
    long long evtTime = (useTriggerTime && trigger_time && !trigger_time->empty())
                         ? trigger_time->at(0) : 0;
    int      evtNum  = (useTriggerNumber && trigger_number && !trigger_number->empty())
                         ? trigger_number->at(0) : 0;

    // ←── HERE ──► Apply your event filter
    // e.g., if (evtTime < timeMin || evtTime > timeMax) continue;
    //       if (evtNum != desiredTrigger) continue;
    // … rest of loop: accumulate sums, fill histos …

    // Sum per event per geom
    unordered_map<int,double> eventSum;
    int nCh = MID->size();
    for (int j = 0; j < nCh; ++j) {
      int mid_j = MID->at(j);
      int ch_j  = ch->at(j);

      if (mid_j != 41 && mid_j != 42) continue; // only for 41 and 42

      auto it = channelToGeom.find({mid_j, ch_j});
      if (it == channelToGeom.end()) {
        // mapping not found
        cout << "Warning: no geom mapping for MID " << mid_j << " CH " << ch_j << endl;
        continue;
      }
      int geom = it->second;
      int start = waveform_idx->at(j);
      int end   = (j+1 < nCh ? waveform_idx->at(j+1) : waveform_total->size());
      double sum = 0;
      for (int k = start; k < end; ++k) sum += waveform_total->at(k);
      eventSum[geom] += sum;
    }
    for (auto &p : eventSum) {
      int geom = p.first;
      double val = p.second;
      // fill per-event distribution
      if (hDataDist.count(geom)) {
        hDataDist[geom]->Fill(val);
        // hDataDist[geom]->SetDirectory(0);    // <-- detach from fData (now done at creation)
      }
      sum_data[geom]   += val;
      count_data[geom] += 1;
    }
  }
  fData->Close();

  // --- Simulation histograms processing ---
  TFile *fSim = TFile::Open(simFile, "READ");
  if (!fSim || fSim->IsZombie()) {
    cerr << "Error: cannot open sim file " << simFile << endl;
    return;
  }
  cout << "Processing simulation histograms based on sim_map.txt" << endl;
  map<int, double> sum_sim;
  map<int, long>   count_sim;
  for (auto &kv : simToGeom) {
    const string &histName = kv.first;
    int geom = kv.second;
    TH1 *hSim = (TH1*)fSim->Get(histName.c_str());
    if (!hSim) {
      cout<<"Warning: histogram "<<histName<<" not found in "<<simFile<<endl;
      continue;
    }
    // compute integral or mean as representative simIntADC
    double simIntADC = hSim->Integral();
    // fill per-geom distribution and sums
    if (hSimDist.count(geom)) {
      hSimDist[geom]->Fill(simIntADC);
      // hSimDist[geom]->SetDirectory(0);    // <-- detach from file (now done at creation)
    }
    sum_sim[geom]   += simIntADC;
    count_sim[geom] += 1;
  }
  fSim->Close();

  // --- POST-SIM DEBUG INFO ---
  cout << "Post-simulation histogram summary:" << endl;
  for (auto &kv : hSimDist) {
    int geom = kv.first;
    TH1D *h = kv.second;
    cout << Form("  GeomID %2d: hist=%p, entries=%lld, integral=%.0f",
                 geom, (void*)h, h->GetEntries(), h->Integral()) << endl;
  }

  // Prepare output for calibration constants
  TTree *tCalib = new TTree("Calibration", "Per-geom calibration constants");
  int    calibGeom;
  double calibValue;
  tCalib->Branch("GeomID",     &calibGeom, "GeomID/I");
  tCalib->Branch("CalibConst", &calibValue, "CalibConst/D");

  std::ofstream csvOut("calibration_constants.txt");
  csvOut << "#GeomID,CalibConst\n";

  // --- estimate & print calibration constants ---
  cout << "\nGeomID  mean_data  mean_sim  CalibC(sim/data)" << endl;
  for (auto &p : sum_data) {
    int mod = p.first;
    if (count_data[mod]==0 || count_sim[mod]==0) continue;
    double md = sum_data[mod] / count_data[mod];
    double ms = sum_sim[mod]  / count_sim[mod];
    double C  = ms / md;
    printf("  %2d   %10.3f  %10.3f  %8.5f\n", mod, md, ms, C);
    // write to CSV and TTree
    calibGeom  = mod;
    calibValue = C;
    tCalib->Fill();
    csvOut << calibGeom << "," << calibValue << "\n";
  }

  // --- Write out distributions ---
  TFile fout("calibration_bic_output.root", "RECREATE");
  for (auto &kv : hDataDist) {
    kv.second->Write();
  }
  for (auto &kv : hSimDist) {
    kv.second->Write();
  }
  // Write calibration constants
  csvOut.close();
  tCalib->Write();
  fout.Close();
  cout << "Wrote output file calibration_bic_output.root with per-geom distributions." << endl;

}