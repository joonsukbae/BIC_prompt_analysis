// energy_calibration_bic.C
// Macro: calibrate energy using calibration constants, draw per-geomID histograms in geometry order.

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "caloMap.h"
#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <unordered_map>

int energy_calibration_bic(
  const char* waveRoot   = "Data/Waveform_sample.root",
  const char* calibRoot  = "calibration_constant_output/calibration_bic_output.root",
  const char* outRoot    = "energy_calibration_output/energy_calibration_QC.root",
  int adcThreshold = 0
) {
  // 1. Load calibration constants from calibRoot - 수정: GeomID와 Side 쌍으로 불러오기
  std::map<std::pair<int,int>, double> geomSideCal; // (GeomID, Side) -> CalibConst
  {
    TFile* fcal = TFile::Open(calibRoot, "READ");
    if (!fcal || fcal->IsZombie()) {
      std::cerr << "Error: cannot open " << calibRoot << "\n";
      return 1;
    }
    TTree* tcal = dynamic_cast<TTree*>(fcal->Get("Calibration"));
    if (!tcal) {
      std::cerr << "Error: TTree \"Calibration\" not found in " << calibRoot << "\n";
      fcal->Close();
      return 1;
    }
    int geomID, side; double cc;
    tcal->SetBranchAddress("GeomID", &geomID);
    tcal->SetBranchAddress("Side", &side);
    tcal->SetBranchAddress("CalibConst", &cc);
    Long64_t n = tcal->GetEntries();
    std::cout << "Found " << n << " calibration constants in tree\n";
    for (Long64_t i=0; i<n; ++i) {
      tcal->GetEntry(i);
      geomSideCal[std::make_pair(geomID, side)] = cc;
      std::cout << "Loaded: GeomID=" << geomID << ", Side=" << side << ", CalibConst=" << cc << std::endl;
    }
    // Print number of loaded calibration constants
    std::cout << "Loaded calibration constants for " << geomSideCal.size() << " (GeomID,Side) pairs\n";
    fcal->Close();
  }

  // 2. Build channelCal: (MID,CH) -> CalibConst using caloMap - 수정: Side 정보 포함
  auto chMap = GetCaloChMap();
  std::map<std::pair<int,int>, double> channelCal;
  std::cout << "\n=== Channel Mapping Debug ===\n";
  for (const auto& kv : chMap) {
    auto key = kv.first; // {MID,CH}
    auto vec = kv.second; // {side,mod,row,layer}
    int side = vec[0];    // 0=L, 1=R
    int layer = vec[3], row = vec[2];
    int geomID = layer*8 + row + 1;
    auto calibKey = std::make_pair(geomID, side);
    auto it = geomSideCal.find(calibKey);
    if (it != geomSideCal.end()) {
      channelCal[key] = it->second;
      std::cout << "MID=" << key.first << ", CH=" << key.second 
                << " -> GeomID=" << geomID << ", Side=" << side 
                << " -> CalibConst=" << it->second << std::endl;
    } else {
      channelCal[key] = 1.0;
      std::cerr << "Warning: no CalibConst for (GeomID=" << geomID << ", Side=" << side << ")\n";
    }
  }
  std::cout << "Mapped " << channelCal.size() << " channels to calibration constants\n";
  
  // 캘리브레이션 상수 통계 출력
  double minCalib = 1e6, maxCalib = -1e6;
  for (const auto& kv : channelCal) {
    minCalib = std::min(minCalib, kv.second);
    maxCalib = std::max(maxCalib, kv.second);
  }
  std::cout << "Calibration constants range: " << minCalib << " to " << maxCalib << std::endl;

  // 3. Open waveform file, setup TTreeReader and channel readers
  TFile* fw = TFile::Open(waveRoot, "READ");
  if (!fw || fw->IsZombie()) {
    std::cerr << "Error: cannot open " << waveRoot << "\n";
    return 1;
  }
  
  // Auto-detect TTree in waveform file
  TTree* tree = nullptr;
  {
    TIter nextKey(fw->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextKey())) {
      TObject* obj = key->ReadObj();
      if (obj->InheritsFrom("TTree")) {
        tree = (TTree*)obj;
        std::cout << "Using TTree: " << tree->GetName() << std::endl;
        break;
      }
    }
    if (!tree) {
      std::cerr << "Error: no TTree found in " << waveRoot << "\n";
      fw->Close();
      return 1;
    }
  }
  
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<short>> vWave(reader, "waveform_total");
  TTreeReaderValue<std::vector<int>> vIdx(reader, "waveform_idx");
  TTreeReaderValue<std::vector<int>> vMID(reader, "MID");
  TTreeReaderValue<std::vector<int>> vCh(reader, "ch");
  int nEventsProcessed = 0;

  // 4. Create 24 histograms (index 1..24, hCal[0] unused)
  std::vector<TH1D*> hCal(25, nullptr);
  for (int g=1; g<=24; ++g) {
    hCal[g] = new TH1D(Form("hCal_G%d",g),
                       Form("Geom %d Calibrated Energy;E_{cal} [MeV];Entries",g),
                       100, 0, 1000);
    hCal[g]->SetDirectory(0);
  }

  // 4b. Total calibrated energy histogram
  TH1D* hTotal = new TH1D("hTotalCal",
                          "Total calibrated energy per event;E_{tot} [MeV];Events",
                          200, 0, 20000);
  hTotal->SetDirectory(0);

  // 4c. L/R 합산 히스토그램 추가
  std::vector<TH1D*> hCalLR(25, nullptr);
  for (int g=1; g<=24; ++g) {
    hCalLR[g] = new TH1D(Form("hCalLR_G%d",g),
                          Form("Geom %d L+R Calibrated Energy;E_{cal} [MeV];Entries",g),
                          100, 0, 1000);
    hCalLR[g]->SetDirectory(0);
  }

  // 5. Loop over events: fill per-geomID histograms
  int debugEvent = 0;
  double minRaw = 1e6, maxRaw = -1e6;
  double minEcal = 1e6, maxEcal = -1e6;
  while (reader.Next()) {
    double sumTotal = 0;
    ++nEventsProcessed;
    
    // 각 GeomID별로 L/R 에너지 누적
    std::map<int, double> geomEnergyL;
    std::map<int, double> geomEnergyR;
    
    for (size_t i = 0; i < (*vMID).size(); ++i) {
      int mid = (*vMID)[i];
      int ch  = (*vCh)[i];
      
      // MID 41, 42만 처리
      if (mid != 41 && mid != 42) continue;
      
      auto key = std::make_pair(mid, ch);
      auto itMap = chMap.find(key);
      if (itMap == chMap.end()) continue;
      
      auto& vec = itMap->second;
      int side = vec[0];    // 0=L, 1=R
      int layer = vec[3], row = vec[2];
      int geomID = layer*8 + row + 1;
      if (geomID < 1 || geomID > 24) continue;
      
      auto itCal = channelCal.find(key);
      double cc = (itCal != channelCal.end() ? itCal->second : 1.0);
      
      int start = (*vIdx)[i];
      int end   = (i + 1 < (*vMID).size() ? (*vIdx)[i+1] : (*vWave).size());
      double sumRaw = 0;
      for (int idx = start; idx < end; ++idx) {
        sumRaw += (*vWave)[idx];
      }
      // only integrate if total ADC exceeds threshold (same as calibration_bic.C)
      if (sumRaw < adcThreshold)
        continue;
      double ecal = sumRaw * cc;
      
      // Raw ADC와 Ecal 범위 추적
      minRaw = std::min(minRaw, sumRaw);
      maxRaw = std::max(maxRaw, sumRaw);
      minEcal = std::min(minEcal, ecal);
      maxEcal = std::max(maxEcal, ecal);
      
      // 첫 번째 이벤트에서 디버깅 출력
      if (debugEvent < 3) {
        std::cout << "Event " << debugEvent << ": MID=" << mid << ", CH=" << ch 
                  << " -> GeomID=" << geomID << ", Side=" << side 
                  << " -> Raw=" << sumRaw << ", CalibConst=" << cc 
                  << " -> Ecal=" << ecal << " MeV\n";
      }
      
      sumTotal += ecal;
      
      // L/R별로 에너지 누적
      if (side == 0) { // L
        geomEnergyL[geomID] += ecal;
      } else { // R
        geomEnergyR[geomID] += ecal;
      }
    }
    
    if (debugEvent < 3) {
      std::cout << "Event " << debugEvent << " total energy: " << sumTotal << " MeV\n";
      debugEvent++;
    }
    
    // 이벤트별로 GeomID의 L+R 합산 에너지를 히스토그램에 채우기
    for (int g = 1; g <= 24; ++g) {
      double sumLR = geomEnergyL[g] + geomEnergyR[g];
      if (sumLR > 0) {
        hCal[g]->Fill(sumLR);      // 개별 채널이 아닌 이벤트별 L+R 합산
        hCalLR[g]->Fill(sumLR);
      }
    }
    
    hTotal->Fill(sumTotal);
  }
  fw->Close();
  std::cout << "Processed " << nEventsProcessed << " events\n";
  std::cout << "Raw ADC range: " << minRaw << " to " << maxRaw << std::endl;
  std::cout << "Calibrated energy range: " << minEcal << " to " << maxEcal << " MeV\n";

  // 6. Draw histograms in geometry order on canvas (3 rows x 8 cols)
  TCanvas* c = new TCanvas("cCalQC","Energy Calibration QC",1600,800);
  c->Divide(8,3);
  // Pads 1-8: Geom17-24; 9-16: Geom9-16; 17-24: Geom1-8
  for (int pad=1; pad<=24; ++pad) {
    c->cd(pad);
    int geomID = 0;
    if (pad >= 1 && pad <=8) geomID = 16+pad;
    else if (pad >=9 && pad<=16) geomID = pad;
    else if (pad >=17 && pad<=24) geomID = pad-16;
    if (geomID>=1 && geomID<=24 && hCal[geomID]) {
      hCal[geomID]->Draw("hist");
    }
  }
  c->SaveAs("energy_calibration_output/energy_calibration_QC.png");

  // 7. Write all histograms to output root file
  // Create output directory if it doesn't exist
  system("mkdir -p energy_calibration_output");
  TFile* fo = TFile::Open(outRoot, "RECREATE");
  for (int g=1; g<=24; ++g) {
    if (hCal[g]) hCal[g]->Write();
    if (hCalLR[g]) hCalLR[g]->Write();
  }
  // write total calibrated energy histogram
  if (hTotal) hTotal->Write();
  fo->Close();
  long totalCalEntries = 0;
  for (int g=1; g<=24; ++g) {
    if (hCal[g]) totalCalEntries += hCal[g]->GetEntries();
  }
  std::cout << "Total calibrated hits across all modules: " << totalCalEntries << "\n";

  // 8. Print mean, sigma, resolution for each geomID
  std::cout << "\n=== Individual Channel Statistics ===\n";
  for (int g=1; g<=24; ++g) {
    if (!hCal[g]) continue;
    double mean = hCal[g]->GetMean();
    double sigma = hCal[g]->GetStdDev();
    double res = (mean>0) ? (sigma/mean)*100. : 0.;
    std::cout << "GeomID " << g << ": mean = " << mean
              << " MeV, sigma = " << sigma
              << " MeV, resolution = " << res << " %\n";
  }
  
  // 8b. Print L/R 합산 통계
  std::cout << "\n=== L+R Combined Statistics ===\n";
  for (int g=1; g<=24; ++g) {
    if (!hCalLR[g]) continue;
    double mean = hCalLR[g]->GetMean();
    double sigma = hCalLR[g]->GetStdDev();
    double res = (mean>0) ? (sigma/mean)*100. : 0.;
    std::cout << "GeomID " << g << " (L+R): mean = " << mean
              << " MeV, sigma = " << sigma
              << " MeV, resolution = " << res << " %\n";
  }
  
  // Print total per-event energy distribution statistics
  std::cout << "\nTotal energy per event: mean = "
            << hTotal->GetMean() << " MeV, sigma = "
            << hTotal->GetStdDev() << " MeV\n";
  return 0;
}