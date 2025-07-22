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
#include <regex>
#include <string>
#include <cstring>
#include <cstdio>

// input 파일명에서 run number 또는 고유 문자열 추출 (C 스타일, char* 반환)
const char* extractRunTag(const char* dataFile) {
  static char buf[128];
  // 파일명만 추출
  const char* fname = strrchr(dataFile, '/');
  fname = (fname ? fname + 1 : dataFile);
  // Run_60184_Waveform.root → Run60184
  const char* runptr = strstr(fname, "Run_");
  if (runptr) {
    int runnum = 0;
    if (sscanf(runptr, "Run_%d", &runnum) == 1) {
      snprintf(buf, sizeof(buf), "Run%d", runnum);
      return buf;
    }
  }
  // Waveform_sample.root → Waveform_sample
  const char* wptr = strstr(fname, "_Waveform");
  if (wptr && wptr > fname) {
    size_t len = wptr - fname;
    if (len > sizeof(buf)-1) len = sizeof(buf)-1;
    strncpy(buf, fname, len);
    buf[len] = '\0';
    return buf;
  }
  // 확장자 제거
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
  std::vector<TH1D*> hRawADC(25, nullptr); // Raw ADC histograms for fitting
  for (int g=1; g<=24; ++g) {
    hCal[g] = new TH1D(Form("hCal_G%d",g),
                       Form("Geom %d Calibrated Energy;E_{cal} [MeV];Entries",g),
                       100, 0, 1000);
    hCal[g]->SetDirectory(0);
    
    hRawADC[g] = new TH1D(Form("hRawADC_G%d",g),
                           Form("Geom %d Raw ADC;ADC;Entries",g),
                           200, 0, 100000);
    hRawADC[g]->SetDirectory(0);
  }

  // 4b. Total calibrated energy histogram
  TH1D* hTotal = new TH1D("hTotalCal",
                          "Total calibrated energy per event;E_{tot} [MeV];Events",
                          200, 0, 20000);
  hTotal->SetDirectory(0);
  
  // 4d. Total raw ADC histogram for fitting
  TH1D* hTotalRawADC = new TH1D("hTotalRawADC",
                                 "Total raw ADC per event;ADC;Events",
                                 200, 0, 1000000);
  hTotalRawADC->SetDirectory(0);

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
      
      // 첫 번째 이벤트에서만 간단한 디버깅 출력
      if (debugEvent < 1) {
        std::cout << "Event " << debugEvent << ": MID=" << mid << ", CH=" << ch 
                  << " -> GeomID=" << geomID << ", Side=" << side 
                  << " -> Ecal=" << ecal << " MeV\n";
      }
      
      sumTotal += ecal;
      
      // Fill raw ADC histogram for fitting
      hRawADC[geomID]->Fill(sumRaw);
      
      // L/R별로 에너지 누적
      if (side == 0) { // L
        geomEnergyL[geomID] += ecal;
      } else { // R
        geomEnergyR[geomID] += ecal;
      }
    }
    
    // Fill total raw ADC histogram for fitting
    double totalRawADC = 0.0;
    for (size_t i = 0; i < (*vMID).size(); ++i) {
      int mid = (*vMID)[i];
      int ch  = (*vCh)[i];
      
      if (mid != 41 && mid != 42) continue;
      
      auto key = std::make_pair(mid, ch);
      auto itMap = chMap.find(key);
      if (itMap == chMap.end()) continue;
      
      int start = (*vIdx)[i];
      int end   = (i + 1 < (*vMID).size() ? (*vIdx)[i+1] : (*vWave).size());
      double sumRaw = 0;
      for (int idx = start; idx < end; ++idx) {
        sumRaw += (*vWave)[idx];
      }
      if (sumRaw >= adcThreshold) {
        totalRawADC += sumRaw;
      }
    }
    hTotalRawADC->Fill(totalRawADC);
    
    if (debugEvent < 1) {
      std::cout << "Event " << debugEvent << " total energy: " << sumTotal << " MeV\n";
      std::cout << "Event " << debugEvent << " total raw ADC: " << totalRawADC << "\n";
      debugEvent++;
    }
    
    // 이벤트별로 GeomID의 L+R 합산 에너지를 히스토그램에 채우기
    // Note: L/R 각각에 대해 캘리브레이션 상수가 적용되었으므로 합산이 올바름
    for (int g = 1; g <= 24; ++g) {
      double sumLR = geomEnergyL[g] + geomEnergyR[g];
      if (sumLR > 0) {
        hCal[g]->Fill(sumLR);      // L+R 합산 (각각 캘리브레이션 적용됨)
        hCalLR[g]->Fill(sumLR);
      }
    }
    
    hTotal->Fill(sumTotal);
  }
  fw->Close();
  std::cout << "Processed " << nEventsProcessed << " events\n";
  std::cout << "Raw ADC range: " << minRaw << " to " << maxRaw << std::endl;
  std::cout << "Calibrated energy range: " << minEcal << " to " << maxEcal << " MeV\n";

  // output 파일명에 runTag 적용
  const char* runTag = extractRunTag(waveRoot);
  char outRootFile[256], outPngFile[256];
  snprintf(outRootFile, sizeof(outRootFile), "energy_calibration_output/energy_calibration_QC_%s.root", runTag);
  snprintf(outPngFile, sizeof(outPngFile), "energy_calibration_output/energy_calibration_QC_%s.png", runTag);

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
  c->SaveAs(outPngFile);

  // 7. Write all histograms to output root file
  // Create output directory if it doesn't exist
  system("mkdir -p energy_calibration_output");
  TFile* fo = TFile::Open(outRootFile, "RECREATE");
  for (int g=1; g<=24; ++g) {
    if (hCal[g]) hCal[g]->Write();
    if (hCalLR[g]) hCalLR[g]->Write();
    if (hRawADC[g]) hRawADC[g]->Write(); // Save raw ADC histograms for analysis
  }
  // write total calibrated energy histogram
  if (hTotal) hTotal->Write();
  if (hTotalRawADC) hTotalRawADC->Write(); // Save total raw ADC histogram for fitting
  fo->Close();
  long totalCalEntries = 0;
  for (int g=1; g<=24; ++g) {
    if (hCal[g]) totalCalEntries += hCal[g]->GetEntries();
  }
  std::cout << "Total calibrated hits across all modules: " << totalCalEntries << "\n";

  // 8. Print statistics for each GeomID (L+R only, since we're using combined histograms)
  std::cout << "\n=== Per-GeomID Statistics (L+R Combined) ===\n";
  for (int g=1; g<=24; ++g) {
    if (!hCal[g]) continue;
    
    // Get fitted mean from raw ADC histogram
    double meanRawADC_fit = 0.0;
    if (hRawADC[g] && hRawADC[g]->GetEntries() > 10) {
      TF1* fRaw = new TF1("fRaw", "gaus", hRawADC[g]->GetXaxis()->GetXmin(), hRawADC[g]->GetXaxis()->GetXmax());
      fRaw->SetParameters(hRawADC[g]->GetMaximum(), hRawADC[g]->GetMean(), hRawADC[g]->GetRMS());
      hRawADC[g]->Fit(fRaw, "Q");
      meanRawADC_fit = fRaw->GetParameter(1);
      delete fRaw;
    } else {
      meanRawADC_fit = hRawADC[g] ? hRawADC[g]->GetMean() : 0.0;
    }
    
    // Get fitted mean from calibrated energy histogram
    double meanLR_fit = 0.0;
    if (hCal[g]->GetEntries() > 10) {
      TF1* fCal = new TF1("fCal", "gaus", hCal[g]->GetXaxis()->GetXmin(), hCal[g]->GetXaxis()->GetXmax());
      fCal->SetParameters(hCal[g]->GetMaximum(), hCal[g]->GetMean(), hCal[g]->GetRMS());
      hCal[g]->Fit(fCal, "Q");
      meanLR_fit = fCal->GetParameter(1);
      delete fCal;
    } else {
      meanLR_fit = hCal[g]->GetMean();
    }
    
    double sigmaLR = hCal[g]->GetStdDev();
    double resLR = (meanLR_fit>0) ? (sigmaLR/meanLR_fit)*100. : 0.;
    
    std::cout << "GeomID " << std::setw(2) << g << ": "
              << "E(fit)=" << std::setw(6) << std::fixed << std::setprecision(1) << meanLR_fit 
              << " MeV (σ=" << std::setw(5) << std::fixed << std::setprecision(1) << sigmaLR 
              << ", " << std::setw(5) << std::fixed << std::setprecision(1) << resLR << "%)\n";
  }
  
  // Print total per-event energy distribution statistics
  std::cout << "\nTotal energy per event: mean = "
            << hTotal->GetMean() << " MeV, sigma = "
            << hTotal->GetStdDev() << " MeV\n";
  
  // Calculate total energy sum across all events (for comparison with simulation)
  double totalEnergySum = 0.0;
  for (int g = 1; g <= 24; ++g) {
    if (hCal[g]) {
      totalEnergySum += hCal[g]->GetMean() * hCal[g]->GetEntries();
    }
  }
  std::cout << "Total energy sum across all events: " << totalEnergySum << " MeV\n";
  std::cout << "Average energy per hit event: " << (totalEnergySum / nEventsProcessed) << " MeV\n";
  
  // Calculate total energy per event (sum of all GeomIDs)
  // Note: This is the sum of average energies per GeomID, not per-event total
  double totalEnergyPerEvent = 0.0;
  for (int g = 1; g <= 24; ++g) {
    if (hCal[g]) {
      totalEnergyPerEvent += hCal[g]->GetMean();
    }
  }
  std::cout << "Sum of average energies per GeomID: " << totalEnergyPerEvent << " MeV\n";
  std::cout << "Expected simulation energy per event: 2777 MeV\n";
  std::cout << "Ratio (calibrated/simulation): " << (totalEnergyPerEvent / 2777.0) * 100.0 << "%\n";
  
  // Calculate actual per-event total energy from histogram
  double actualTotalEnergyPerEvent = hTotal->GetMean();
  std::cout << "Actual per-event total energy (from histogram): " << actualTotalEnergyPerEvent << " MeV\n";
  std::cout << "Ratio (actual/simulation): " << (actualTotalEnergyPerEvent / 2777.0) * 100.0 << "%\n";
  
  // Fit total energy distribution for energy resolution
  std::cout << "\n=== Total Energy Distribution Fitting ===" << std::endl;
  if (hTotal->GetEntries() > 10) {
    TF1* fTotal = new TF1("fTotal", "gaus", hTotal->GetXaxis()->GetXmin(), hTotal->GetXaxis()->GetXmax());
    fTotal->SetParameters(hTotal->GetMaximum(), hTotal->GetMean(), hTotal->GetStdDev());
    hTotal->Fit(fTotal, "Q");
    double meanTotal_fit = fTotal->GetParameter(1);
    double sigmaTotal_fit = fTotal->GetParameter(2);
    double resolutionTotal = (meanTotal_fit > 0) ? (sigmaTotal_fit / meanTotal_fit) * 100.0 : 0.0;
    std::cout << "Fitted total energy: mean = " << meanTotal_fit << " MeV, sigma = " << sigmaTotal_fit << " MeV\n";
    std::cout << "Energy resolution: " << resolutionTotal << "%\n";
    delete fTotal;
  }
  
  // Fit total raw ADC distribution
  std::cout << "\n=== Total Raw ADC Distribution Fitting ===" << std::endl;
  if (hTotalRawADC->GetEntries() > 10) {
    TF1* fRaw = new TF1("fRaw", "gaus", hTotalRawADC->GetXaxis()->GetXmin(), hTotalRawADC->GetXaxis()->GetXmax());
    fRaw->SetParameters(hTotalRawADC->GetMaximum(), hTotalRawADC->GetMean(), hTotalRawADC->GetStdDev());
    hTotalRawADC->Fit(fRaw, "Q");
    double meanRaw_fit = fRaw->GetParameter(1);
    double sigmaRaw_fit = fRaw->GetParameter(2);
    double resolutionRaw = (meanRaw_fit > 0) ? (sigmaRaw_fit / meanRaw_fit) * 100.0 : 0.0;
    std::cout << "Fitted total raw ADC: mean = " << meanRaw_fit << ", sigma = " << sigmaRaw_fit << "\n";
    std::cout << "ADC resolution: " << resolutionRaw << "%\n";
    delete fRaw;
  }
  
  // Debug: Check a few specific GeomIDs
  std::cout << "\n=== Debug: Specific GeomID Values ===\n";
  for (int g = 10; g <= 12; ++g) {
    if (hCal[g]) {
      std::cout << "GeomID " << g << ": mean = " << hCal[g]->GetMean() 
                << " MeV, entries = " << hCal[g]->GetEntries() << "\n";
    }
  }
  return 0;
}