// energy_calibration_all_layers.C
// Macro: calibrate energy using calibration constants from all layers, draw per-geomID histograms in geometry order.
// This macro processes all layers (4x8 or 3x8) using separate calibration files for each layer.

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "caloMap_old.h"
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

int energy_calibration_all_layers(
  const char* waveRoot   = "Data/Waveform_sample.root",
  const char* calibRoot0 = "calibration_constant_output/calibration_bic_output_layer0.root",
  const char* calibRoot1 = "calibration_constant_output/calibration_bic_output_layer1.root",
  const char* calibRoot2 = "calibration_constant_output/calibration_bic_output_layer2.root",
  const char* calibRoot3 = "calibration_constant_output/calibration_bic_output_layer3.root",
  const char* outRoot    = "energy_calibration_output/energy_calibration_all_layers_QC.root",
  int adcThreshold = 0,
  bool isNewType = false  // false: 4x8 (구형), true: 3x8 (신형)
) {
  // 1. Load calibration constants from all layers
  std::map<std::pair<int,int>, double> geomSideCal; // (GeomID, Side) -> CalibConst
  
  // Load calibration constants for each layer
  const char* calibRoots[4] = {calibRoot0, calibRoot1, calibRoot2, calibRoot3};
  int maxLayers = isNewType ? 3 : 4; // 3x8 or 4x8
  
  for (int layer = 0; layer < maxLayers; ++layer) {
    TFile* fcal = TFile::Open(calibRoots[layer], "READ");
    if (!fcal || fcal->IsZombie()) {
      std::cerr << "Error: cannot open " << calibRoots[layer] << "\n";
      continue;
    }
    TTree* tcal = dynamic_cast<TTree*>(fcal->Get("Calibration"));
    if (!tcal) {
      std::cerr << "Error: TTree \"Calibration\" not found in " << calibRoots[layer] << "\n";
      fcal->Close();
      continue;
    }
    int geomID, side; double cc;
    tcal->SetBranchAddress("GeomID", &geomID);
    tcal->SetBranchAddress("Side", &side);
    tcal->SetBranchAddress("CalibConst", &cc);
    Long64_t n = tcal->GetEntries();
    std::cout << "Found " << n << " calibration constants in " << calibRoots[layer] << "\n";
    for (Long64_t i=0; i<n; ++i) {
      tcal->GetEntry(i);
      geomSideCal[std::make_pair(geomID, side)] = cc;
      std::cout << "Loaded: GeomID=" << geomID << ", Side=" << side << ", CalibConst=" << cc << std::endl;
    }
    fcal->Close();
  }
  std::cout << "Loaded calibration constants for " << geomSideCal.size() << " (GeomID,Side) pairs\n";

  // 2. Build channelCal: (MID,CH) -> CalibConst using caloMap_old.h - 수정: Side 정보 포함
  auto chMap = GetCaloChMapOld();
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
  TTreeReaderValue<std::vector<int>> vDataLength(reader, "data_length");
  int nEventsProcessed = 0;

  // 4. Create histograms for all modules (4x8 or 3x8)
  std::vector<TH1D*> hCal(33, nullptr);
  std::vector<TH1D*> hRawADC(33, nullptr); // Raw ADC histograms for fitting
  for (int g=1; g<=32; ++g) {
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
  std::vector<TH1D*> hCalLR(33, nullptr);
  for (int g=1; g<=32; ++g) {
    hCalLR[g] = new TH1D(Form("hCalLR_G%d",g),
                          Form("Geom %d L+R Calibrated Energy;E_{cal} [MeV];Entries",g),
                          100, 0, 1000);
    hCalLR[g]->SetDirectory(0);
  }
  
  // 4e. Simulation Edep histograms for comparison (all layers)
  std::vector<TH1D*> hSimEdep(33, nullptr);
  TFile* fSim = TFile::Open("Sim/4x8_5GeV_3rd_result_new.root", "READ");
  if (fSim && !fSim->IsZombie()) {
    for (int layer = 0; layer < maxLayers; ++layer) {
      for (int col=0; col<8; ++col) {
        int g = layer * 8 + col + 1; // GeomID for each layer
        TH1* hOrig = (TH1*)fSim->Get(Form("Edep_M%d", g));
        if (hOrig) {
          hSimEdep[g] = (TH1D*)hOrig->Clone(Form("hSimEdep_G%d", g));
          hSimEdep[g]->SetDirectory(0);
          cout << "Loaded simulation Edep for GeomID " << g << endl;
        }
      }
    }
    fSim->Close();
  }

  // 5. Loop over events: fill per-geomID histograms
  int debugEvent = 0;
  double minRaw = 1e6, maxRaw = -1e6;
  double minEcal = 1e6, maxEcal = -1e6;
  while (reader.Next()) {
    // data_length.size()=92가 아닌 경우 이벤트 스킵 (92개 채널이 모두 켜진 이벤트만)
    if ((*vDataLength).size() != 92) continue;
    
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
      int layer = vec[3], col = vec[2];
      // 새로운 mapping 방식 사용
      int geomID = layer * 8 + col + 1;
      
      // 모든 층 처리 (maxLayers까지)
      if (layer < maxLayers && geomID >= 1 && geomID <= 32) {
        auto itCal = channelCal.find(key);
        double cc = (itCal != channelCal.end() ? itCal->second : 1.0);
        
        int start = (*vIdx)[i];
        int end   = (i + 1 < (*vMID).size() ? (*vIdx)[i+1] : (*vWave).size());
        double sumRaw = 0;
        // ADC/TDC가 번갈아 들어있으므로 짝수 bin만 읽기 (ADC만)
        // waveform overlay plot과 동일한 구간 사용
        for (int idx = start; idx < end; idx += 2) {
          if (idx >= 0 && idx < (*vWave).size()) {
            sumRaw += (*vWave)[idx];
          }
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
    }
    
    // Fill total raw ADC histogram for fitting (all layers)
    double totalRawADC = 0.0;
    for (size_t i = 0; i < (*vMID).size(); ++i) {
      int mid = (*vMID)[i];
      int ch  = (*vCh)[i];
      
      if (mid != 41 && mid != 42) continue;
      
      auto key = std::make_pair(mid, ch);
      auto itMap = chMap.find(key);
      if (itMap == chMap.end()) continue;
      
      auto& vec = itMap->second;
      int layer = vec[3];
      
      // 모든 층 처리 (maxLayers까지)
      if (layer < maxLayers) {
        int start = (*vIdx)[i];
        int end   = (i + 1 < (*vMID).size() ? (*vIdx)[i+1] : (*vWave).size());
        double sumRaw = 0;
        // ADC/TDC가 번갈아 들어있으므로 짝수 bin만 읽기 (ADC만)
        for (int idx = start; idx < end; idx += 2) {
          if (idx >= 0 && idx < (*vWave).size()) {
            sumRaw += (*vWave)[idx];
          }
        }
        if (sumRaw >= adcThreshold) {
          totalRawADC += sumRaw;
        }
      }
    }
    hTotalRawADC->Fill(totalRawADC);
    
    if (debugEvent < 1) {
      std::cout << "Event " << debugEvent << " total energy: " << sumTotal << " MeV\n";
      std::cout << "Event " << debugEvent << " total raw ADC: " << totalRawADC << "\n";
      debugEvent++;
    }
    
    // 이벤트별로 GeomID의 L+R 합산 에너지를 히스토그램에 채우기 (all layers)
    // Note: L/R 각각에 대해 캘리브레이션 상수가 적용되었으므로 합산이 올바름
    for (int layer = 0; layer < maxLayers; ++layer) {
      for (int col = 0; col < 8; ++col) {
        int g = layer * 8 + col + 1; // GeomID for each layer
        double sumLR = geomEnergyL[g] + geomEnergyR[g];
        if (sumLR > 0) {
          hCal[g]->Fill(sumLR);      // L+R 합산 (각각 캘리브레이션 적용됨)
          hCalLR[g]->Fill(sumLR);
        }
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
  snprintf(outRootFile, sizeof(outRootFile), "energy_calibration_output/energy_calibration_all_layers_QC_%s.root", runTag);
  snprintf(outPngFile, sizeof(outPngFile), "energy_calibration_output/energy_calibration_all_layers_QC_%s.png", runTag);

  // 6. Draw histograms in geometry order on canvas (all layers)
  TCanvas* c = new TCanvas("cCalQC","Energy Calibration QC (All Layers)",1600,900);
  c->Divide(8,4); // 4x8 grid for all layers
  for (int layer = 0; layer < 4; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int geomID = layer * 8 + col + 1;
      int pad = (3 - layer) * 8 + (col + 1); // layer 0 at bottom, layer 3 at top
      c->cd(pad);
      
      // 해당 층이 maxLayers 범위 내에 있으면 데이터 표시
      if (layer < maxLayers && hCal[geomID]) {
        // Peak normalization - maximum value to 1
        double dataMax = hCal[geomID]->GetMaximum();
        if (dataMax > 0) hCal[geomID]->Scale(1.0/dataMax);
        hCal[geomID]->GetYaxis()->SetRangeUser(0, 1.1); // Fixed Y-axis range
        hCal[geomID]->Draw("hist");
        if (hSimEdep[geomID]) {
          hSimEdep[geomID]->SetLineColor(kRed);
          // Peak normalization for simulation too
          double simMax = hSimEdep[geomID]->GetMaximum();
          if (simMax > 0) hSimEdep[geomID]->Scale(1.0/simMax);
          hSimEdep[geomID]->GetYaxis()->SetRangeUser(0, 1.1); // Fixed Y-axis range
          hSimEdep[geomID]->Draw("SAME");
        }
      } else {
        // 해당 층이 없으면 빈 pad로 표시
        TLatex tex;
        tex.SetNDC();
        tex.SetTextSize(0.08);
        tex.SetTextColor(kGray);
        tex.DrawLatex(0.5, 0.5, Form("Layer %d", layer));
      }
    }
  }
  c->SaveAs(outPngFile);

  // 6b. Create total energy comparison plot (Data vs Simulation)
  TCanvas* cTotal = new TCanvas("cTotal","Total Energy Comparison",1200,800);
  cTotal->Divide(2,1);
  
  // Calculate total energy per event for all modules
  TH1D* hTotalData = new TH1D("hTotalData", "Total Calibrated Energy (Data);E_{tot} [MeV];Events", 200, 0, 5000);
  TH1D* hTotalSim = new TH1D("hTotalSim", "Total Energy Deposit (Simulation);E_{tot} [MeV];Events", 200, 0, 5000);
  
  // Sum up all module energies for each event
  for (int layer = 0; layer < maxLayers; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int g = layer * 8 + col + 1;
      if (hCal[g] && hSimEdep[g]) {
        // Add data histogram
        hTotalData->Add(hCal[g]);
        // Add simulation histogram
        hTotalSim->Add(hSimEdep[g]);
      }
    }
  }
  
  // Plot total energy comparison
  cTotal->cd(1);
  hTotalData->SetLineColor(kBlue);
  hTotalData->SetLineWidth(2);
  hTotalData->Draw("hist");
  hTotalSim->SetLineColor(kRed);
  hTotalSim->SetLineWidth(2);
  hTotalSim->Draw("SAME");
  
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->AddEntry(hTotalData, "Data (Calibrated)", "l");
  leg->AddEntry(hTotalSim, "Simulation", "l");
  leg->Draw();
  
  // Calculate and display resolution
  cTotal->cd(2);
  double dataMean = hTotalData->GetMean();
  double dataSigma = hTotalData->GetStdDev();
  double simMean = hTotalSim->GetMean();
  double simSigma = hTotalSim->GetStdDev();
  
  double dataResolution = (dataMean > 0) ? (dataSigma / dataMean) * 100.0 : 0.0;
  double simResolution = (simMean > 0) ? (simSigma / simMean) * 100.0 : 0.0;
  
  // Create text display
  TLatex tex;
  tex.SetNDC();
  tex.SetTextSize(0.04);
  tex.DrawLatex(0.1, 0.9, Form("Data: Mean = %.1f MeV, #sigma = %.1f MeV", dataMean, dataSigma));
  tex.DrawLatex(0.1, 0.8, Form("Data Resolution: %.1f%%", dataResolution));
  tex.DrawLatex(0.1, 0.7, Form("Simulation: Mean = %.1f MeV, #sigma = %.1f MeV", simMean, simSigma));
  tex.DrawLatex(0.1, 0.6, Form("Sim Resolution: %.1f%%", simResolution));
  tex.DrawLatex(0.1, 0.5, Form("Ratio (Data/Sim): %.2f", dataMean/simMean));
  
  cTotal->SaveAs(Form("energy_calibration_output/total_energy_comparison_all_layers_%s.png", runTag));

  // 7. Write all histograms to output root file
  // Create output directory if it doesn't exist
  system("mkdir -p energy_calibration_output");
  TFile* fo = TFile::Open(outRootFile, "RECREATE");
  for (int layer = 0; layer < maxLayers; ++layer) {
    for (int col=0; col<8; ++col) {
      int g = layer * 8 + col + 1; // GeomID for each layer
      if (hCal[g]) hCal[g]->Write();
      if (hCalLR[g]) hCalLR[g]->Write();
      if (hRawADC[g]) hRawADC[g]->Write(); // Save raw ADC histograms for analysis
      if (hSimEdep[g]) hSimEdep[g]->Write(); // Save simulation Edep histograms
    }
  }
  // write total calibrated energy histogram
  if (hTotal) hTotal->Write();
  if (hTotalRawADC) hTotalRawADC->Write(); // Save total raw ADC histogram for fitting
  if (hTotalData) hTotalData->Write(); // Save total data energy histogram
  if (hTotalSim) hTotalSim->Write(); // Save total simulation energy histogram
  fo->Close();
  long totalCalEntries = 0;
  for (int layer = 0; layer < maxLayers; ++layer) {
    for (int col=0; col<8; ++col) {
      int g = layer * 8 + col + 1; // GeomID for each layer
      if (hCal[g]) totalCalEntries += hCal[g]->GetEntries();
    }
  }
  std::cout << "Total calibrated hits across all modules: " << totalCalEntries << "\n";

  // 8. Print statistics for each GeomID (L+R only, since we're using combined histograms)
  std::cout << "\n=== Per-GeomID Statistics (L+R Combined) ===\n";
  for (int layer = 0; layer < maxLayers; ++layer) {
    for (int col=0; col<8; ++col) {
      int g = layer * 8 + col + 1; // GeomID for each layer
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
      
      // Get simulation Edep mean for comparison
      double meanSimEdep = 0.0;
      if (hSimEdep[g]) {
        meanSimEdep = hSimEdep[g]->GetMean();
      }
      
      std::cout << "GeomID " << std::setw(2) << g << ": "
                << "E(fit)=" << std::setw(6) << std::fixed << std::setprecision(1) << meanLR_fit 
                << " MeV (σ=" << std::setw(5) << std::fixed << std::setprecision(1) << sigmaLR 
                << ", " << std::setw(5) << std::fixed << std::setprecision(1) << resLR << "%)";
      if (meanSimEdep > 0) {
        std::cout << " | Sim Edep=" << std::setw(6) << std::fixed << std::setprecision(1) << meanSimEdep 
                  << " MeV | Ratio=" << std::setw(5) << std::fixed << std::setprecision(1) << (meanLR_fit/meanSimEdep*100.0) << "%";
      }
      std::cout << "\n";
    }
  }
  
  // Print total per-event energy distribution statistics
  std::cout << "\nTotal energy per event: mean = "
            << hTotal->GetMean() << " MeV, sigma = "
            << hTotal->GetStdDev() << " MeV\n";
  
  // Calculate total energy sum across all events (for comparison with simulation)
  double totalEnergySum = 0.0;
  for (int layer = 0; layer < maxLayers; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int g = layer * 8 + col + 1; // GeomID for each layer
      if (hCal[g]) {
        totalEnergySum += hCal[g]->GetMean() * hCal[g]->GetEntries();
      }
    }
  }
  std::cout << "Total energy sum across all events: " << totalEnergySum << " MeV\n";
  std::cout << "Average energy per hit event: " << (totalEnergySum / nEventsProcessed) << " MeV\n";
  
  // Calculate total energy per event (sum of all GeomIDs)
  // Note: This is the sum of average energies per GeomID, not per-event total
  double totalEnergyPerEvent = 0.0;
  for (int layer = 0; layer < maxLayers; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int g = layer * 8 + col + 1; // GeomID for each layer
      if (hCal[g]) {
        totalEnergyPerEvent += hCal[g]->GetMean();
      }
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