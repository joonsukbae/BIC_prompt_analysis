// intADC_QA.C
// QA macro: plot integrated ADC distributions for each channel
// Usage:
//   root -l -q 'intADC_QA.C("Data/Run_60264_Waveform.root", 1)'

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <string>
#include <regex>
#include "caloMap.h"

// Extract run tag from filename
string extractRunTag(const char* filename) {
  string fname(filename);
  regex runPattern(R"(Run_(\d+))");
  smatch match;
  if (regex_search(fname, match, runPattern)) {
    return match[1].str();
  }
  return "unknown";
}

void intADC_QA(const char* dataFile = "Data/Waveform_sample.root", int targetLayer = 1) {
  // 1. Open data file
  TFile* f = TFile::Open(dataFile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open file " << dataFile << std::endl;
    return;
  }
  
  // Extract run tag for output filename
  string runTag = extractRunTag(dataFile);
  
  // 2. Auto-detect TTree
  TTree* t = nullptr;
  {
    TIter nextKey(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextKey())) {
      TObject* obj = key->ReadObj();
      if (obj->InheritsFrom("TTree")) {
        t = (TTree*)obj;
        std::cout << "Using TTree: " << t->GetName() << std::endl;
        break;
      }
    }
    if (!t) {
      std::cerr << "Error: no TTree found in " << dataFile << std::endl;
      return;
    }
  }
  
  // 3. Set branch addresses
  std::vector<short>* waveform_total = nullptr;
  std::vector<int>* waveform_idx = nullptr;
  std::vector<int>* MID = nullptr;
  std::vector<int>* ch = nullptr;
  std::vector<int>* data_length = nullptr;
  t->SetBranchAddress("waveform_total", &waveform_total);
  t->SetBranchAddress("waveform_idx", &waveform_idx);
  t->SetBranchAddress("MID", &MID);
  t->SetBranchAddress("ch", &ch);
  t->SetBranchAddress("data_length", &data_length);

  // 4. Load mapping
  auto chMap = GetCaloChMap();

  // 5. Prepare histograms for (geom, lr) pairs
  std::map<std::pair<int,int>, TH1D*> hIntADC; // (geom, lr) → hist

  // 6. Event loop
  Long64_t nEvt = t->GetEntries();
  std::cout << "Total events: " << nEvt << std::endl;
  int skippedEvents = 0;
  for (Long64_t i = 0; i < nEvt; ++i) {
    t->GetEntry(i);
    
    // Skip events without all 92 channels
    if (data_length->size() != 92) {
      skippedEvents++;
      continue;
    }
    
    int nCh = MID->size();
    for (int j = 0; j < nCh; ++j) {
      int mid = MID->at(j);
      int chn = ch->at(j);
      auto it = chMap.find({mid, chn});
      if (it == chMap.end()) continue; // Skip unmapped channels
      int lr = it->second[0];
      int actualMod = it->second[1];
      int col = it->second[2];
      int layer = it->second[3];
      
      // Use new mapping: geom = layer * 8 + col + 1
      int geom = layer * 8 + col + 1;
      
      // Process only target layer
      if (layer == targetLayer) {
        int start = waveform_idx->at(j) + 100;
        int end = waveform_idx->at(j) + 200;
        double sum = 0;
        // Read only even bins (ADC only, not TDC)
        for (int k = start; k < end; k += 2) {
          if (k >= 0 && k < waveform_total->size())
            sum += waveform_total->at(k);
        }
        auto key = std::make_pair(geom, lr);
        if (!hIntADC.count(key)) {
          TString hname = Form("hIntADC_Geom%d_%c", geom, lr ? 'R' : 'L');
          hIntADC[key] = new TH1D(hname, hname+";intADC;Events", 100, 0, 70000);
          hIntADC[key]->SetDirectory(0);
        }
        hIntADC[key]->Fill(sum);
      }
    }
  }
  f->Close();
  std::cout << "Skipped " << skippedEvents << " events (not all 92 channels)" << std::endl;

  // 7. Create canvas and draw histograms
  TCanvas* c = new TCanvas("cIntADC", "Integrated ADC QA", 1600, 900);
  c->Divide(8, 4); // 4x8 grid
  
  for (int layer = 0; layer < 4; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int geom = layer * 8 + col + 1;
      int pad = (3 - layer) * 8 + (col + 1); // layer 0 at bottom, layer 3 at top
      c->cd(pad);
      
      // Show data only for target layer
      if (layer == targetLayer) {
        auto keyL = std::make_pair(geom, 0);
        auto keyR = std::make_pair(geom, 1);
        
        bool hasData = false;
        if (hIntADC.count(keyL)) {
          hIntADC[keyL]->SetLineColor(kBlue);
          hIntADC[keyL]->Draw();
          hasData = true;
        }
        if (hIntADC.count(keyR)) {
          hIntADC[keyR]->SetLineColor(kRed);
          if (hasData) {
            hIntADC[keyR]->Draw("SAME");
          } else {
            hIntADC[keyR]->Draw();
          }
        }
        
        // Add legend
        if (hIntADC.count(keyL) && hIntADC.count(keyR)) {
          TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
          leg->AddEntry(hIntADC[keyL], "L", "l");
          leg->AddEntry(hIntADC[keyR], "R", "l");
          leg->Draw();
        }
      } else {
        // Empty pad for other layers
        TLatex tex;
        tex.SetNDC();
        tex.SetTextSize(0.08);
        tex.SetTextColor(kGray);
        tex.DrawLatex(0.5, 0.5, Form("Layer %d", layer));
      }
    }
  }
  
  // Create output directory
  system("mkdir -p intADC_QA_output");
  
  // Save canvas with run tag and layer info
  string outFile = "intADC_QA_output/intADC_QA_" + runTag + "_layer" + to_string(targetLayer) + ".png";
  c->SaveAs(outFile.c_str());
  std::cout << "Saved " << outFile << std::endl;
} 