// calibration_bic.C
// Usage:
//   root -l
//   root[0] .L calibration_bic.C
//   root[1]
//   calibration_bic("Run_60184_Waveform.root","3x5_5GeV_result_new.root",
//   "data_map.txt", "sim_map.txt", true, false);

#include "TIterator.h"
#include "TKey.h"
#include "caloMap.h"
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TTree.h>
#include <unordered_map>
// #include <fstream>
#include <iostream>
// #include <map>
// #include <sstream>
// #include <unordered_map>
// #include <utility>
#include <TCanvas.h>
#include <TLatex.h>
#include <algorithm>
#include <vector>

struct pair_hash {
  size_t operator()(const std::pair<int, int> &p) const noexcept {
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

using namespace std;

void calibration_bic(const char *dataFile = "Data/Waveform_sample.root",
                     const char *simFile = "Sim/3x8_3GeV_CERN_hist.root",
                     const double beamEnergyGeV = 3.0,
                     bool useTriggerTime = true,
                     bool useTriggerNumber = false,
                     int adcThreshold = 0) {
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
    while ((keyData = (TKey *)nextKeyData())) {
      TObject *obj = keyData->ReadObj();
      if (obj->InheritsFrom("TTree")) {
        tData = (TTree *)obj;
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
  std::vector<short> *waveform_total = nullptr;
  std::vector<int> *waveform_idx = nullptr;
  std::vector<int> *MID = nullptr;
  std::vector<int> *ch = nullptr;
  tData->SetBranchAddress("waveform_total", &waveform_total);
  tData->SetBranchAddress("waveform_idx", &waveform_idx);
  tData->SetBranchAddress("MID", &MID);
  tData->SetBranchAddress("ch", &ch);
  std::vector<long long> *trigger_time = nullptr;
  std::vector<int> *trigger_number = nullptr;
  tData->SetBranchAddress("trigger_time", &trigger_time);
  tData->SetBranchAddress("trigger_number", &trigger_number);

  // --- Read mapping files ---

  // --- Prepare per-geom histograms for data and simulation ---
  std::map<std::pair<int,int>, TH1D*> hDataDistLR; // (geomID, lr)
  std::map<int, TH1D*> hSimDist; // sim remains per geom
  const int nBins = 200;
  const double maxADC = 1e6;
  // Allocate histograms after reading mapping files
  // 1) Data mapping from caloMap.h
  auto dataChMap = GetCaloChMap();
  cout << "Loaded " << dataChMap.size()
       << " channel-to-geom entries from caloMap.h" << endl;
  // Build module → GeomID and module labels dynamically
  unordered_map<int, int> moduleToGeom;
  unordered_map<int, string> geomLabel;
  vector<int> uniqueGeoms;
  for (auto &kv : dataChMap) {
    int lr = kv.second[0];            // 0=L, 1=R
    int mod = kv.second[1];           // actual module number
    int col = kv.second[2];           // 0..7
    int layer = kv.second[3];         // 0..7
    int geomID = layer * 8 + col + 1; // 1..24
    moduleToGeom[mod] = geomID;
    // label "M<mod><L/R>"
    geomLabel[geomID] = Form("M%d%c", mod, (lr == 0 ? 'L' : 'R'));
    if (find(uniqueGeoms.begin(), uniqueGeoms.end(), geomID) ==
        uniqueGeoms.end())
      uniqueGeoms.push_back(geomID);
  }
  cout << "Derived " << uniqueGeoms.size() << " GeomIDs" << endl;

  // Allocate one histogram per (geom, L/R)
  for (auto &kv : dataChMap) {
    int lr   = kv.second[0];  // 0=L, 1=R
    int mod  = kv.second[1];
    int geom = moduleToGeom[mod];
    auto key = std::make_pair(geom, lr);
    if (!hDataDistLR.count(key)) {
      std::string name = Form("hData_G%d_%c", geom, lr ? 'R' : 'L');
      std::string title = Form("Data INT ADC Geom %d %c;INT ADC;Events", geom, lr ? 'R' : 'L');
      hDataDistLR[key] = new TH1D(name.c_str(), title.c_str(), nBins, 0, maxADC);
      hDataDistLR[key]->SetDirectory(0);
    }
  }
  // Simulation histos allocated later as before

  // Module Accumulation (sum, count) for data per (geom, lr)
  std::unordered_map<std::pair<int,int>, double, pair_hash> sum_dataLR;
  std::unordered_map<std::pair<int,int>, long, pair_hash> count_dataLR;

  Long64_t nD = tData->GetEntries();
  cout << "Data entries: " << nD << endl;
  for (Long64_t i = 0; i < nD; ++i) {
    tData->GetEntry(i);
    long long evtTime =
        (useTriggerTime && trigger_time && !trigger_time->empty())
            ? trigger_time->at(0)
            : 0;
    int evtNum =
        (useTriggerNumber && trigger_number && !trigger_number->empty())
            ? trigger_number->at(0)
            : 0;

    // ←── HERE ──► Apply your event filter
    // e.g., if (evtTime < timeMin || evtTime > timeMax) continue;
    //       if (evtNum != desiredTrigger) continue;
    // … rest of loop: accumulate sums, fill histos …

    // Sum per event per (geom, lr)
    std::unordered_map<std::pair<int,int>, double, pair_hash> eventSumLR;
    int nCh = MID->size();
    for (int j = 0; j < nCh; ++j) {
      int mid_j = MID->at(j);
      int ch_j = ch->at(j);

      if (mid_j != 41 && mid_j != 42)
        continue; // only for 41 and 42

      // lookup geom; unmapped channels are silently skipped
      auto it = dataChMap.find({mid_j, ch_j});
      if (it == dataChMap.end()) {
        continue;
      }
      int lr = it->second[0];
      int actualMod = it->second[1];
      int geom = moduleToGeom[actualMod];
      auto key = std::make_pair(geom, lr);
      int start = waveform_idx->at(j);
      int end =
          (j + 1 < nCh ? waveform_idx->at(j + 1) : waveform_total->size());
      double sum = 0;
      for (int k = start; k < end; ++k)
        sum += waveform_total->at(k);
      // only integrate if total ADC exceeds threshold
      if (sum < adcThreshold)
        continue;
      eventSumLR[key] += sum;
    }
    // After summing, fill histograms and accumulate sums/counts
    for (auto &p : eventSumLR) {
      auto key = p.first;      // pair(geom,lr)
      double val = p.second;
      if (hDataDistLR.count(key)) {
        hDataDistLR[key]->Fill(val);
      }
      sum_dataLR[key]   += val;
      count_dataLR[key] += 1;
    }
  }
  fData->Close();

  // --- Simulation histograms processing ---
  TFile *fSim = TFile::Open(simFile, "READ");
  if (!fSim || fSim->IsZombie()) {
    cerr << "Error: cannot open sim file " << simFile << endl;
    return;
  }
  map<int, double> sum_sim;
  map<int, long> count_sim;
  // Process simulation by GeomID 1–24 directly
  for (int geom = 1; geom <= 24; ++geom) {
    TH1 *hOrig = (TH1*)fSim->Get(Form("Edep_M%d", geom));
    if (!hOrig) {
      cout << "Warning: sim hist Edep_M" << geom << " missing" << endl;
      continue;
    }
    // Clone entire sim histogram for QA overlay
    TH1D *hCloned = (TH1D*)hOrig->Clone(Form("hSim_G%d", geom));
    hCloned->Scale(0.5); // 0.5 is for unseparated L/R in sim
    hCloned->SetDirectory(0);
    hSimDist[geom] = hCloned;
    // Store mean value for calibration calculation
    double meanVal = hOrig->GetMean() * 0.5; // 0.5 is for unseparated L/R in sim
    sum_sim[geom] = meanVal;
    count_sim[geom] = 1;
  }
  fSim->Close();

  // --- POST-SIM DEBUG INFO ---
  cout << "Post-simulation histogram summary:" << endl;
  for (auto &kv : hSimDist) {
    int geom = kv.first;
    TH1D *h = kv.second;
    cout << Form("  GeomID %2d: hist=%p, entries=%lld, integral=%.0f", geom,
                 (void *)h, h->GetEntries(), h->Integral())
         << endl;
  }

  // Map (GeomID,LR) → actual module number for labeling
  unordered_map<pair<int,int>, int, pair_hash> geomLRToMod;
  for (auto &kv : dataChMap) {
    int lr   = kv.second[0];  // 0=L, 1=R
    int mod  = kv.second[1];  // actual module number
    int geom = moduleToGeom[mod];
    geomLRToMod[{geom, lr}] = mod;
  }

  // Prepare output for calibration constants
  TTree *tCalib = new TTree("Calibration", "Per-geom calibration constants");
  int calibGeom;
  int calibSide;
  double calibValue;
  tCalib->Branch("GeomID", &calibGeom, "GeomID/I");
  tCalib->Branch("Side",   &calibSide, "Side/I");
  tCalib->Branch("CalibConst", &calibValue, "CalibConst/D");

  std::ofstream csvOut("calibration_constant_output/calibration_constants.txt");
  csvOut << "#GeomID,Side,CalibConst\n";

  // --- estimate & print calibration constants ---
  cout << "\nGeomID Module  mean_data  mean_sim(MeV, %)  CalibC(sim/data)" << endl;
  // Write 48 calibration constants
  for (auto &kv : sum_dataLR) {
    auto key = kv.first;
    int geom = key.first;
    int lr   = key.second;
    double md = sum_dataLR[key] / count_dataLR[key];
    double ms = (sum_sim.count(geom) ? sum_sim[geom] : 0.0);
    // percent of half-beam-energy
    double pct = (beamEnergyGeV > 0) ? (ms / (beamEnergyGeV * 1000.0 * 0.5) * 100.0) : 0.0;
    double C  = (md != 0.0 ? (ms / md) : 0.0);
    // Determine actual module label
    int mod = geomLRToMod[{geom, lr}];
    char side = (lr ? 'R' : 'L');
    string label = Form("M%d%c", mod, side);
    printf("  %2d     %-6s  %10.3f  %10.3f (%.1f%%)  %8.5f\n",
           geom, label.c_str(), md, ms, pct, C);
    calibGeom  = geom;
    calibSide  = lr;
    calibValue = C;
    tCalib->Fill();
    csvOut << geom << "," << (lr ? 'R' : 'L') << "," << C << "\n";
  }

  // --- total simulation energy deposit summary ---
  double totalSimE = 0.0;
  for (int geom = 1; geom <= 24; ++geom) {
    auto it = sum_sim.find(geom);
    if (it != sum_sim.end()) {
      totalSimE += it->second;
    }
  }
  double totalPct = 0.0;
  if (beamEnergyGeV > 0) {
    double beamMeV = beamEnergyGeV * 1000.0;
    totalPct = (totalSimE / beamMeV / 2) * 100.0;
  }
  printf("Total sim Edep: %.1f MeV = %.2f%% of beam energy (%.1f GeV)\n",
         totalSimE, totalPct, beamEnergyGeV);

  // --- Write out distributions ---
  // Create output directory if it doesn't exist
  system("mkdir -p calibration_constant_output");
  TFile fout("calibration_constant_output/calibration_bic_output.root", "RECREATE");
  for (auto &kv : hDataDistLR) {
    kv.second->Write();
  }
  for (auto &kv : hSimDist) {
    kv.second->Write();
  }
  // Write calibration constants
  csvOut.close();
  tCalib->Write();
  fout.Close();
  cout << "Wrote output file calibration_bic_output.root with per-geom distributions."
       << endl;

  // --- QA: overlay Data vs Sim vs Calibration constant per module ---
  // Define grid order: rows are geom 17-24, 9-16, 1-8, for all 24 modules
  vector<int> qaOrder = {17, 18, 19, 20, 21, 22, 23, 24, 9, 10, 11, 12,
                         13, 14, 15, 16, 1,  2,  3,  4,  5, 6,  7,  8};
  // For each geom, draw L and R overlays in one pad per module
  TCanvas *cQA = new TCanvas("cQA", "Calibration QA per Module", 2000, 700);
  cQA->Divide(8, 3); // 8 columns, 3 rows for 24 modules
  // Draw per-module summary: Data L, Data R, Sim mean, CalibConst L/R
  for (size_t idx = 0; idx < qaOrder.size(); ++idx) {
    int geom = qaOrder[idx];
    cQA->cd(idx + 1);
    // Draw data L
    auto keyL = make_pair(geom, 0);
    if (hDataDistLR.count(keyL)) {
      hDataDistLR[keyL]->SetLineColor(kBlue);
      // hDataDistLR[keyL]->GetXaxis()->SetRangeUser(0,200);
      hDataDistLR[keyL]->Draw();
    }
    // Draw data R
    auto keyR = make_pair(geom, 1);
    if (hDataDistLR.count(keyR)) {
      hDataDistLR[keyR]->SetLineColor(kGreen+2);
      // hDataDistLR[keyR]->GetXaxis()->SetRangeUser(0,200);
      hDataDistLR[keyR]->Draw("SAME");
    }
    // // Draw sim
    // if (hSimDist.count(geom)) {
    //   hSimDist[geom]->SetLineColor(kRed);
    //   hSimDist[geom]->Draw("SAME");
    // }
    // Compute means and calibration constants
    double meanL = (count_dataLR[keyL]>0 ? sum_dataLR[keyL]/count_dataLR[keyL] : 0);
    double meanR = (count_dataLR[keyR]>0 ? sum_dataLR[keyR]/count_dataLR[keyR] : 0);
    double meanS = (count_sim[geom]>0 ? sum_sim[geom]/count_sim[geom] : 0);
    double CL = (meanL!=0 ? meanS/meanL : 0);
    double CR = (meanR!=0 ? meanS/meanR : 0);
    // Compute standard deviations for L, R, Sim
    double stdL = 0.0, stdR = 0.0, stdS = 0.0;
    if (hDataDistLR.count(keyL)) stdL = hDataDistLR[keyL]->GetMeanError();
    if (hDataDistLR.count(keyR)) stdR = hDataDistLR[keyR]->GetMeanError();
    if (hSimDist.count(geom))   stdS = hSimDist[geom]->GetMeanError() * 0.5;
    // Annotate
    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.06);
    // Show both GeomID and module label
    int mod = -1;
    char sideChar = '?';
    // Prefer L if available, else R
    if (geomLRToMod.count({geom,0})) { mod = geomLRToMod[{geom,0}]; sideChar='L'; }
    else if (geomLRToMod.count({geom,1})) { mod = geomLRToMod[{geom,1}]; sideChar='R'; }
    tex.SetTextColor(kBlack);
    tex.DrawLatex(0.15, 0.85, Form("Geom %d (%s)", geom, Form("M%d", mod)));
    // Data mean ± stddev
    tex.SetTextColor(kBlue);
    tex.DrawLatex(0.15, 0.75, Form("µ_L=%.2g #pm %.2g", meanL, stdL));
    tex.SetTextColor(kGreen+2);
    tex.DrawLatex(0.15, 0.68, Form("µ_R=%.2g #pm %.2g", meanR, stdR));
    // Sim mean ± stddev and percentage of half-beam-energy
    tex.SetTextColor(kRed);
    double pctS = (beamEnergyGeV > 0) ? (meanS / (beamEnergyGeV * 1000.0 * 0.5) * 100.0) : 0.0;
    tex.DrawLatex(0.15, 0.60,
      Form("µ_Sim=%.2g #pm %.2g MeV (%.2g%%)", meanS, stdS, pctS));
    // Calibration constants
    tex.SetTextColor(kBlack);
    tex.DrawLatex(0.15, 0.52, Form("C_L=%.2f, C_R=%.2f", CL, CR));
  }
  cQA->SaveAs("calibration_constant_output/calibration_QA.png");
  // also write canvas into the output root file
  TFile fout2("calibration_constant_output/calibration_bic_output.root", "UPDATE");
  cQA->Write();
  fout2.Close();
}