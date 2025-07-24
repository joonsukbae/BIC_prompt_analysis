// calibration_bic.C
// Macro: calculate calibration constants by comparing data with simulation
// Usage:
//   root -l -q 'calibration_bic.C("Data/Run_60264_Waveform.root", "Sim/3x8_3GeV_CERN_hist.root", 3.0, 1)'

#include "TIterator.h"
#include "TKey.h"
#include "caloMap.h"
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TTree.h>
#include <unordered_map>
#include <iostream>
#include <TCanvas.h>
#include <TLatex.h>
#include <algorithm>
#include <vector>
#include <regex>
#include <string>
using std::string;

#include <cstring>
#include <cstdio>

struct pair_hash {
  size_t operator()(const std::pair<int, int> &p) const noexcept {
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

using namespace std;

// Extract run tag from filename
const char* extractRunTag(const char* dataFile) {
  static char buf[128];
  const char* fname = strrchr(dataFile, '/');
  fname = (fname ? fname + 1 : dataFile);
  
  const char* runptr = strstr(fname, "Run_");
  if (runptr) {
    int runnum = 0;
    if (sscanf(runptr, "Run_%d", &runnum) == 1) {
      snprintf(buf, sizeof(buf), "Run%d", runnum);
      return buf;
    }
  }
  
  const char* wptr = strstr(fname, "_Waveform");
  if (wptr && wptr > fname) {
    size_t len = wptr - fname;
    if (len > sizeof(buf)-1) len = sizeof(buf)-1;
    strncpy(buf, fname, len);
    buf[len] = '\0';
    return buf;
  }
  
  const char* dot = strrchr(fname, '.');
  if (dot && dot > fname) {
    size_t len = dot - fname;
    if (len > sizeof(buf)-1) len = sizeof(buf)-1;
    strncpy(buf, fname, len);
    buf[len] = '\0';
    return buf;
  }
  strncpy(buf, fname, sizeof(buf)-1);
  buf[sizeof(buf)-1] = '\0';
  return buf;
}

void calibration_bic(const char *dataFile = "Data/Waveform_sample.root",
                     const char *simFile = "Sim/3x8_3GeV_CERN_hist.root",
                     const double beamEnergyGeV = 3.0,
                     int targetLayer = 1,
                     int adcThreshold = 0,
                     bool useTriggerTime = true,
                     bool useTriggerNumber = false,
                     double peakThreshold = 0.0,
                     double xMax = 100000.0) {
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
  std::vector<int> *data_length = nullptr;
  tData->SetBranchAddress("waveform_total", &waveform_total);
  tData->SetBranchAddress("waveform_idx", &waveform_idx);
  tData->SetBranchAddress("MID", &MID);
  tData->SetBranchAddress("ch", &ch);
  tData->SetBranchAddress("data_length", &data_length);
  std::vector<long long> *trigger_time = nullptr;
  std::vector<int> *trigger_number = nullptr;
  tData->SetBranchAddress("trigger_time", &trigger_time);
  tData->SetBranchAddress("trigger_number", &trigger_number);

  // --- Prepare per-geom histograms for data and simulation ---
  std::map<std::pair<int,int>, TH1D*> hDataDistLR; // (geomID, lr)
  std::map<int, TH1D*> hSimDist; // sim remains per geom
  
  // Data mapping from caloMap_old.h
  auto dataChMap = GetCaloChMap();
  cout << "Loaded " << dataChMap.size() << " channel-to-geom entries from caloMap.h" << endl;
  
  // Build module → GeomID and module labels dynamically
  unordered_map<int, int> moduleToGeom;
  unordered_map<int, string> geomLabel;
  vector<int> uniqueGeoms;
  for (auto &kv : dataChMap) {
    int lr = kv.second[0];            // 0=L, 1=R
    int mod = kv.second[1];           // actual module number (GeomID for old mapping)
    int col = kv.second[2];           // 0..7
    int layer = kv.second[3];         // 0..7
    int geomID = mod;                 // For old mapping, mod is already the GeomID
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
    int mod  = kv.second[1];  // GeomID for old mapping
    int geom = mod;            // Use mod directly as GeomID
    auto key = std::make_pair(geom, lr);
    if (!hDataDistLR.count(key)) {
      std::string name = Form("hData_G%d_%c", geom, lr ? 'R' : 'L');
      std::string title = Form("Data INT ADC Geom %d %c;INT ADC;Events", geom, lr ? 'R' : 'L');
      hDataDistLR[key] = new TH1D(name.c_str(), title.c_str(), 100, 0, 100000);
      hDataDistLR[key]->SetDirectory(0);
    }
  }

  // Module Accumulation (sum, count) for data per (geom, lr)
  std::unordered_map<std::pair<int,int>, double, pair_hash> sum_dataLR;
  std::unordered_map<std::pair<int,int>, long, pair_hash> count_dataLR;

  Long64_t nD = tData->GetEntries();
  cout << "Data entries: " << nD << endl;
  for (Long64_t i = 0; i < nD; ++i) {
    tData->GetEntry(i);
    
    // data_length.size()=92가 아닌 경우 이벤트 스킵 (92개 채널이 모두 켜진 이벤트만)
    if (data_length->size() != 92) continue;
    
    long long evtTime =
        (useTriggerTime && trigger_time && !trigger_time->empty())
            ? trigger_time->at(0)
            : 0;
    int evtNum =
        (useTriggerNumber && trigger_number && !trigger_number->empty())
            ? trigger_number->at(0)
            : 0;

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
      int geomID = actualMod; // For old mapping, actualMod is already the GeomID
      int col = it->second[2];
      int layer = it->second[3];

      // 새로운 mapping 방식 사용
      geomID = layer * 8 + col + 1;

      // targetLayer에 해당하는 층만 처리
      if (layer == targetLayer) {
        auto key = std::make_pair(geomID, lr);
        if (j >= waveform_idx->size()) {
          cout << "Warning: channel index " << j << " out of range for waveform_idx (size=" << waveform_idx->size() << ")" << endl;
          continue;
        }
        int start = waveform_idx->at(j) + 100;
        int end = waveform_idx->at(j) + 200;
        double sum = 0;
        // ADC/TDC가 번갈아 들어있으므로 짝수 bin만 읽기 (ADC만)
        for (int k = start; k < end; k += 2) {
          if (k >= 0 && k < waveform_total->size()) {
            sum += waveform_total->at(k);
          }
        }
        // only integrate if total ADC exceeds threshold
        if (sum < adcThreshold)
          continue;
        eventSumLR[key] += sum;
      }
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
  
  // Process simulation by GeomID 9-16 (layer 1, col 0-7) for calibration (using Edep)
  // Always use layer 1 (2nd layer, middle layer) for simulation comparison regardless of targetLayer
  int simLayer = 1; // Always use layer 1 (2nd layer, middle layer) for simulation
  cout << "\n=== Loading simulation data for layer " << simLayer << " (GeomID " << (simLayer*8+1) << "-" << (simLayer*8+8) << ") ===" << endl;
  for (int col = 0; col < 8; ++col) {
    int geom = simLayer * 8 + col + 1; // GeomID 9-16 for layer 1 (2nd layer, middle layer)
    TH1 *hOrig = (TH1*)fSim->Get(Form("Edep_M%d", geom));
    if (!hOrig) {
      cout << "Warning: sim hist Edep_M" << geom << " missing" << endl;
      continue;
    }
    // Store mean value for calibration calculation (using Edep)
    double meanVal = hOrig->GetMean(); // Use full module energy, no 0.5 scaling
    sum_sim[geom] = meanVal;
    count_sim[geom] = 1;
    cout << "Loaded sim GeomID " << geom << ": mean = " << meanVal << " MeV" << endl;
  }
  
  // Load simulation Edep histograms for QA plotting (always layer 1)
  map<int, TH1D*> hSimEdep; // For QA plotting with Edep
  for (int col = 0; col < 8; ++col) {
    int geom = simLayer * 8 + col + 1; // GeomID 9-16 for layer 1 (2nd layer, middle layer)
    TH1 *hEdep = (TH1*)fSim->Get(Form("Edep_M%d", geom));
    if (hEdep) {
      TH1D *hCloned = (TH1D*)hEdep->Clone(Form("hSimEdep_G%d", geom));
      hCloned->SetDirectory(0);
      hSimEdep[geom] = hCloned;
      cout << "Loaded Edep histogram for GeomID " << geom << endl;
    } else {
      cout << "Warning: no Edep histogram found for GeomID " << geom << endl;
    }
  }
  fSim->Close();

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

  string runTag = extractRunTag(dataFile);
  string outTxt  = "calibration_constant_output/calibration_constants_" + runTag + "_layer" + to_string(targetLayer) + ".txt";
  string outRoot = "calibration_constant_output/calibration_bic_output_" + runTag + "_layer" + to_string(targetLayer) + ".root";
  string outQA   = "calibration_constant_output/calibration_QA_" + runTag + "_layer" + to_string(targetLayer) + ".png";

  std::ofstream csvOut(outTxt);
  csvOut << "#GeomID,Side,CalibConst\n";

  // --- estimate & print calibration constants ---
  cout << "\nGeomID Module  mean_data(fit)  mean_sim(fit,MeV, %)  CalibC(sim/data)" << endl;
  // Write 48 calibration constants
  for (auto &kv : sum_dataLR) {
    auto key = kv.first;
    int geom = key.first;
    int lr   = key.second;
    
    // Calculate col position for data GeomID
    int dataLayer = (geom - 1) / 8;  // 0, 1, 2, or 3
    int dataCol = (geom - 1) % 8;    // 0, 1, 2, ..., 7
    
    // Find corresponding simulation GeomID (same col, layer 1)
    int simLayer = 1; // Always use layer 1 (2nd layer, middle layer) for simulation
    int simGeom = simLayer * 8 + dataCol + 1; // Same col position
    
    // Simple mean calculation for calibration constants (no fitting needed)
    double md = sum_dataLR[key] / count_dataLR[key];
    double ms_half = (sum_sim.count(simGeom) ? sum_sim[simGeom] : 0.0) * 0.5;
    // percent of half-beam-energy (for individual L/R channels)
    double pct = (beamEnergyGeV > 0) ? (ms_half / (beamEnergyGeV * 1000.0 * 0.5) * 100.0) : 0.0;
    double C = (md > 0) ? (ms_half / md) : 0.0;
    // Determine actual module label
    int mod = geomLRToMod[{geom, lr}];
    char side = (lr ? 'R' : 'L');
    string label = Form("M%d%c", mod, side);
    printf("  %2d     %-6s  %10.3f  %10.3f (%.1f%%)  %8.5f (simGeom=%d)\n",
           geom, label.c_str(), md, ms_half, pct, C, simGeom);
    calibGeom  = geom;
    calibSide  = lr;
    calibValue = C;
    tCalib->Fill();
    csvOut << geom << "," << (lr ? 'R' : 'L') << "," << C << "\n";
  }

  // --- total simulation energy deposit summary ---
  double totalSimE = 0.0;
  for (int geom = 1; geom <= 32; ++geom) {
    auto it = sum_sim.find(geom);
    if (it != sum_sim.end()) {
      totalSimE += it->second;
    }
  }
  double totalPct = 0.0;
  if (beamEnergyGeV > 0) {
    double beamMeV = beamEnergyGeV * 1000.0;
    totalPct = (totalSimE / beamMeV) * 100.0;  // /2 제거
  }
  printf("Total sim Edep: %.1f MeV = %.2f%% of beam energy (%.1f GeV)\n",
         totalSimE, totalPct, beamEnergyGeV);

  // --- Write out distributions ---
  system("mkdir -p calibration_constant_output");
  TFile fout(outRoot.c_str(), "RECREATE");
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
  // Build dynamic grid order based on available GeomIDs
  vector<int> qaOrder;
  for (int layer = 0; layer < 4; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int geom = layer * 8 + col + 1;
      qaOrder.push_back(geom);
    }
  }
  int nCols = 8;
  int nRows = 4;
  TCanvas *cQA = new TCanvas("cQA", "Calibration QA per Module", 2000, 900);
  cQA->Divide(nCols, nRows);
  for (int layer = 0; layer < 4; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int geom = layer * 8 + col + 1;
      int pad = (3 - layer) * 8 + (col + 1);
      cQA->cd(pad);
      
      // targetLayer에 해당하는 층만 데이터 표시
      if (layer == targetLayer) {
        auto keyL = make_pair(geom, 0);
        auto keyR = make_pair(geom, 1);
        try {
          if (hDataDistLR.count(keyL)) {
            hDataDistLR[keyL]->SetLineColor(kBlue);
            hDataDistLR[keyL]->Draw();
          }
          if (hDataDistLR.count(keyR)) {
            hDataDistLR[keyR]->SetLineColor(kGreen+2);
            hDataDistLR[keyR]->Draw("SAME");
          }
          // Don't draw sim Edep histogram (scaling issues), just use for text annotation
          // Compute simple means and calibration constants
          double meanL = (count_dataLR[keyL]>0 ? sum_dataLR[keyL]/count_dataLR[keyL] : 0);
          double meanR = (count_dataLR[keyR]>0 ? sum_dataLR[keyR]/count_dataLR[keyR] : 0);
          
          // Find corresponding simulation GeomID (same col, layer 1)
          int dataLayer = (geom - 1) / 8;  // 0, 1, 2, or 3
          int dataCol = (geom - 1) % 8;    // 0, 1, 2, ..., 7
          int simLayer = 1; // Always use layer 1 (2nd layer, middle layer) for simulation
          int simGeom = simLayer * 8 + dataCol + 1; // Same col position
          
          double meanS_full = (count_sim.count(simGeom) ? sum_sim[simGeom]/count_sim[simGeom] : 0);
          double meanS_half = meanS_full * 0.5; // Half for individual L/R channels
          double CL = (meanL!=0 ? meanS_half/meanL : 0);
          double CR = (meanR!=0 ? meanS_half/meanR : 0);
          // Compute standard deviations for L, R, Sim
          double stdL = 0.0, stdR = 0.0, stdS = 0.0;
          if (hDataDistLR.count(keyL)) stdL = hDataDistLR[keyL]->GetMeanError();
          if (hDataDistLR.count(keyR)) stdR = hDataDistLR[keyR]->GetMeanError();
          if (hSimEdep.count(simGeom))   stdS = hSimEdep[simGeom]->GetMeanError(); // Use Edep for QA
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
          double pctS = (beamEnergyGeV > 0) ? (meanS_half / (beamEnergyGeV * 1000.0 * 0.5) * 100.0) : 0.0; // Half for individual L/R
          tex.DrawLatex(0.15, 0.60,
            Form("µ_Edep=%.2g #pm %.2g MeV (%.2g%%)", meanS_half, stdS, pctS));
          // Calibration constants
          tex.SetTextColor(kBlack);
          tex.DrawLatex(0.15, 0.44, Form("C_L=%.2f, C_R=%.2f", CL, CR));
        } catch (const std::exception& e) {
          cout << "Error plotting GeomID " << geom << ": " << e.what() << endl;
        }
      } else {
        // 다른 층은 빈 pad로 표시
        TLatex tex;
        tex.SetNDC();
        tex.SetTextSize(0.08);
        tex.SetTextColor(kGray);
        tex.DrawLatex(0.5, 0.5, Form("Layer %d", layer));
      }
    }
  }
  cQA->SaveAs(outQA.c_str());
  // also write canvas into the output root file
  TFile fout2(outRoot.c_str(), "UPDATE");
  cQA->Write();
  fout2.Close();
}