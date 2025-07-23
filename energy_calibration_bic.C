// energy_calibration_bic.C
// Macro: calibrate energy using calibration constants, draw per-geomID histograms in geometry order.
// Usage:
//   root -l -q 'energy_calibration_bic.C("Data/Run_60264_Waveform.root", "calibration_constant_output/calibration_bic_output_Run60264_layer2.root", "energy_calibration_output/energy_calibration_QC_Run60264_layer2.root", 2)'

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

int energy_calibration_bic(
  const char* waveRoot   = "Data/Waveform_sample.root",
  const char* calibRoot  = "calibration_constant_output/calibration_bic_output_layer2.root",
  const char* outRoot    = "energy_calibration_output/energy_calibration_QC.root",
  int targetLayer = 2,
  int adcThreshold = 0,
  bool isNewType = false  // false: 4x8 (구형), true: 3x8 (신형)
) {
  // 1. Load calibration constants from single layer
  std::map<std::pair<int,int>, double> geomSideCal; // (GeomID, Side) -> CalibConst
  
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
  std::cout << "Found " << n << " calibration constants in " << calibRoot << "\n";
  for (Long64_t i=0; i<n; ++i) {
    tcal->GetEntry(i);
    geomSideCal[std::make_pair(geomID, side)] = cc;
  }
  fcal->Close();

  // 2. Build channelCal: (MID,CH) -> CalibConst using caloMap_old.h
  auto chMap = GetCaloChMapOld();
  std::map<std::pair<int,int>, double> channelCal;
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
    } else {
      channelCal[key] = 1.0;
    }
  }
  std::cout << "Mapped " << channelCal.size() << " channels to calibration constants\n";

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

  // 4. Create histograms for target layer only
  std::vector<TH1D*> hCal(33, nullptr);
  std::vector<TH1D*> hRawADC(33, nullptr);
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

  // Total calibrated energy histogram
  TH1D* hTotal = new TH1D("hTotalCal",
                          "Total calibrated energy per event;E_{tot} [MeV];Events",
                          200, 0, 20000);
  hTotal->SetDirectory(0);
  
  // Total raw ADC histogram for fitting
  TH1D* hTotalRawADC = new TH1D("hTotalRawADC",
                                 "Total raw ADC per event;ADC;Events",
                                 200, 0, 1000000);
  hTotalRawADC->SetDirectory(0);

  // L/R 합산 히스토그램
  std::vector<TH1D*> hCalLR(33, nullptr);
  for (int g=1; g<=32; ++g) {
    hCalLR[g] = new TH1D(Form("hCalLR_G%d",g),
                          Form("Geom %d L+R Calibrated Energy;E_{cal} [MeV];Entries",g),
                          100, 0, 1000);
    hCalLR[g]->SetDirectory(0);
  }
  
  // Simulation Edep histograms for comparison (target layer only)
  std::vector<TH1D*> hSimEdep(33, nullptr);
  TFile* fSim = TFile::Open("Sim/4x8_5GeV_3rd_result_new.root", "READ");
  if (fSim && !fSim->IsZombie()) {
    // 시뮬레이션은 항상 고정된 층(3층, 코드상 2층) 사용
    int simLayer = 2; // Always use layer 2 (3rd layer) for simulation
    for (int col=0; col<8; ++col) {
      int g = simLayer * 8 + col + 1; // GeomID for fixed simulation layer
      TH1* hOrig = (TH1*)fSim->Get(Form("Edep_M%d", g));
      if (hOrig) {
        hSimEdep[g] = (TH1D*)hOrig->Clone(Form("hSimEdep_G%d", g));
        hSimEdep[g]->SetDirectory(0);
      }
    }
    fSim->Close();
  }

  // 5. Loop over events: fill per-geomID histograms
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
      
      // targetLayer에 해당하는 층만 처리
      if (layer == targetLayer && geomID >= 1 && geomID <= 32) {
        auto itCal = channelCal.find(key);
        double cc = (itCal != channelCal.end() ? itCal->second : 1.0);
      
        int start = (*vIdx)[i];
        int end   = (i + 1 < (*vMID).size() ? (*vIdx)[i+1] : (*vWave).size());
        double sumRaw = 0;
        // ADC/TDC가 번갈아 들어있으므로 짝수 bin만 읽기 (ADC만)
        for (int idx = start; idx < end; idx += 2) {
          if (idx >= 0 && idx < (*vWave).size()) {
            sumRaw += (*vWave)[idx];
          }
        }
        // only integrate if total ADC exceeds threshold
        if (sumRaw < adcThreshold)
          continue;
        double ecal = sumRaw * cc;
        
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
    
    // Fill total raw ADC histogram for fitting (target layer only)
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
      
      // targetLayer에 해당하는 층만 처리
      if (layer == targetLayer) {
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
    
    // 이벤트별로 GeomID의 L+R 합산 에너지를 히스토그램에 채우기 (target layer only)
    for (int col = 0; col < 8; ++col) {
      int g = targetLayer * 8 + col + 1; // GeomID for target layer only
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

  // output 파일명에 runTag 적용
  const char* runTag = extractRunTag(waveRoot);
  char outRootFile[256], outPngFile[256];
  snprintf(outRootFile, sizeof(outRootFile), "energy_calibration_output/energy_calibration_QC_%s.root", runTag);
  snprintf(outPngFile, sizeof(outPngFile), "energy_calibration_output/energy_calibration_QC_%s.png", runTag);

  // 6. Draw histograms in geometry order on canvas (target layer only)
  TCanvas* c = new TCanvas("cCalQC","Energy Calibration QC (Target Layer)",1600,900);
  c->Divide(8,4); // 4x8 grid for all layers
  for (int layer = 0; layer < 4; ++layer) {
    for (int col = 0; col < 8; ++col) {
      int geomID = layer * 8 + col + 1;
      int pad = (3 - layer) * 8 + (col + 1); // layer 0 at bottom, layer 3 at top
      c->cd(pad);
      
      // targetLayer에 해당하는 층만 데이터 표시
      if (layer == targetLayer && hCal[geomID]) {
        // Peak normalization - maximum value to 1
        double dataMax = hCal[geomID]->GetMaximum();
        if (dataMax > 0) hCal[geomID]->Scale(1.0/dataMax);
        hCal[geomID]->GetYaxis()->SetRangeUser(0, 1.1); // Fixed Y-axis range
        hCal[geomID]->Draw("hist");
        
        // 같은 col 위치의 시뮬레이션 히스토그램 찾기
        int simLayer = 2; // Always use layer 2 (3rd layer) for simulation
        int simGeomID = simLayer * 8 + col + 1; // Same col position
        if (hSimEdep[simGeomID]) {
          hSimEdep[simGeomID]->SetLineColor(kRed);
          // Peak normalization for simulation too
          double simMax = hSimEdep[simGeomID]->GetMaximum();
          if (simMax > 0) hSimEdep[simGeomID]->Scale(1.0/simMax);
          hSimEdep[simGeomID]->GetYaxis()->SetRangeUser(0, 1.1); // Fixed Y-axis range
          hSimEdep[simGeomID]->Draw("SAME");
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
  c->SaveAs(outPngFile);

  // 6b. Create total energy comparison plot (Data vs Simulation) - target layer only
  TCanvas* cTotal = new TCanvas("cTotal","Total Energy Comparison",1200,800);
  cTotal->Divide(2,1);
  
  // Use already calculated total energy histogram for data
  TH1D* hTotalData = hTotal; // Use the already calculated total energy histogram
  TH1D* hTotalSim = new TH1D("hTotalSim", "Total Energy Deposit (Simulation);E_{tot} [MeV];Events", 200, 0, 5000);
  
  // Calculate total simulation energy by summing individual module energies
  double totalSimEnergy = 0.0;
  int validModules = 0;
  
  // First calculate the mean total energy
  for (int col = 0; col < 8; ++col) {
    int simGeomID = 2 * 8 + col + 1; // Simulation GeomID for fixed layer 2
    if (hSimEdep[simGeomID]) {
      totalSimEnergy += hSimEdep[simGeomID]->GetMean();
      validModules++;
    }
  }
  
  // Fill simulation histogram with multiple events to show distribution
  if (validModules > 0) {
    // Generate multiple events based on simulation statistics
    int nSimEvents = 10000; // Generate 10000 simulation events
    for (int i = 0; i < nSimEvents; ++i) {
      double eventTotalEnergy = 0.0;
      for (int col = 0; col < 8; ++col) {
        int simGeomID = 2 * 8 + col + 1; // Simulation GeomID for fixed layer 2
        if (hSimEdep[simGeomID]) {
          // Sample from each module's energy distribution
          double moduleEnergy = hSimEdep[simGeomID]->GetRandom();
          eventTotalEnergy += moduleEnergy;
        }
      }
      hTotalSim->Fill(eventTotalEnergy);
    }
  }
  
  // Plot total energy comparison
  cTotal->cd(1);
  
  // Normalize histograms to have maximum value of 1
  double dataMax = hTotalData->GetMaximum();
  double simMax = hTotalSim->GetMaximum();
  
  if (dataMax > 0) hTotalData->Scale(1.0 / dataMax);
  if (simMax > 0) hTotalSim->Scale(1.0 / simMax);
  
  hTotalData->SetLineColor(kBlue);
  hTotalData->SetLineWidth(2);
  hTotalData->GetYaxis()->SetRangeUser(0, 1.1); // Fixed Y-axis range
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
  
  cTotal->SaveAs(Form("energy_calibration_output/total_energy_comparison_%s.png", runTag));

  // 7. Write all histograms to output root file
  system("mkdir -p energy_calibration_output");
  TFile* fo = TFile::Open(outRootFile, "RECREATE");
  for (int col=0; col<8; ++col) {
    int dataGeomID = targetLayer * 8 + col + 1; // Data GeomID for target layer
    int simGeomID = 2 * 8 + col + 1; // Simulation GeomID for fixed layer 2
    if (hCal[dataGeomID]) hCal[dataGeomID]->Write();
    if (hCalLR[dataGeomID]) hCalLR[dataGeomID]->Write();
    if (hRawADC[dataGeomID]) hRawADC[dataGeomID]->Write();
    if (hSimEdep[simGeomID]) hSimEdep[simGeomID]->Write();
  }
  if (hTotal) hTotal->Write();
  if (hTotalRawADC) hTotalRawADC->Write();
  if (hTotalData) hTotalData->Write();
  if (hTotalSim) hTotalSim->Write();
  fo->Close();
  
  long totalCalEntries = 0;
  for (int col=0; col<8; ++col) {
    int g = targetLayer * 8 + col + 1;
    if (hCal[g]) totalCalEntries += hCal[g]->GetEntries();
  }
  std::cout << "Total calibrated hits across target layer modules: " << totalCalEntries << "\n";

  // 8. Print statistics for each GeomID
  std::cout << "\n=== Per-GeomID Statistics (L+R Combined) ===\n";
  for (int col=0; col<8; ++col) {
    int g = targetLayer * 8 + col + 1;
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
    int simGeomID = 2 * 8 + col + 1; // Simulation GeomID for fixed layer 2
    if (hSimEdep[simGeomID]) {
      meanSimEdep = hSimEdep[simGeomID]->GetMean();
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
  
  // Print total per-event energy distribution statistics
  std::cout << "\nTotal energy per event: mean = "
            << hTotal->GetMean() << " MeV, sigma = "
            << hTotal->GetStdDev() << " MeV\n";
  
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
  
  return 0;
}