# PS Calibration Macros

## Usage

### calibration_bic.C
```bash
root -l -q 'calibration_bic.C("Data/<Data_Waveform.root>", "Sim/3x8_3GeV_CERN_hist.root", <beamEnergyGeV energy>, <targetLayer>)'
```
Output: `calibration_constant_output/calibration_bic_output_Run60264_layer1.root`

### energy_calibration_bic.C
```bash
root -l -q 'energy_calibration_bic.C("Data/Run_60264_Waveform.root", "calibration_constant_output/calibration_bic_output_Run60264_layer1.root", "Sim/3x8_3GeV_CERN_hist.root", "energy_calibration_output/energy_calibration_QC_Run60264_layer1.root", <beamEnergyGeV energy>, <targetLayer>)'
```
Output: `energy_calibration_output/energy_calibration_QC_Run60264_layer1.root`

## Parameters

- `targetLayer`: detector layer to analyze (0=bottom, 1=middle, 2=top)
- `beamEnergyGeV`: beam energy in GeV (default: 3.0)
- `dataFile`: input data file path
- `simFile`: simulation file path (for calibration)
- `calibFile`: calibration constants file (for energy calibration)
- `outFile`: output file path

## Input Files

- Data files: `Data/Run_XXXXX_Waveform.root`
- Simulation files: `Sim/3x8_3GeV_CERN_hist.root`
- Calibration files: `calibration_constant_output/calibration_bic_output_RunXXXXX_layerX.root`

## Files

- `calibration_bic.C`: calibration constants
- `caloMap.h`: channel mapping
- `energy_calibration_bic.C`: energy calibration
- `README.md`: this file