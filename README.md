# Preshower (PS) Calibration and QA Framework

이 프로젝트는 Preshower 검출기의 캘리브레이션 및 QA를 위한 ROOT 매크로들의 일관된 분석 프레임워크입니다.

## 주요 특징

- **일관된 GeomID 매핑**: `layer * 8 + col + 1` 방식으로 통일 (new mapping)
- **이벤트 필터링**: 92개 채널이 모두 활성화된 이벤트만 처리
- **ADC 적분**: 짝수 bin만 사용하여 ADC만 처리
- **물리적 검출기 구조**: 3층×8열 그리드로 플롯 배치 (층 0이 아래, 층 2가 위)
- **고정 축 범위**: 시각적 비교를 위한 히스토그램 축 범위 고정
- **시뮬레이션 고정층**: 시뮬레이션은 항상 2층(코드상 1층, 가운데 층) 사용

## 매크로 설명

### 1. intADC_QA.C
**기능**: 각 채널별 적분 ADC 분포를 QA 플롯으로 표시
**사용법**:
```bash
root -l -q 'intADC_QA.C("Data/Run_60264_Waveform.root", 1)'
```
**출력**: `intADC_QA_layer1.png`

### 2. calibration_bic.C
**기능**: 데이터와 시뮬레이션을 비교하여 캘리브레이션 상수 계산
**사용법**:
```bash
root -l -q 'calibration_bic.C("Data/Run_60264_Waveform.root", "Sim/3x8_3GeV_CERN_hist.root", 3.0, 1)'
```
**출력**: 
- `calibration_constant_output/calibration_constants_Run60264_layer1.txt`
- `calibration_constant_output/calibration_bic_output_Run60264_layer1.root`
- `calibration_constant_output/calibration_QA_Run60264_layer1.png`

### 3. energy_calibration_bic.C
**기능**: 캘리브레이션 상수를 사용하여 에너지 캘리브레이션 수행
**사용법**:
```bash
root -l -q 'energy_calibration_bic.C("Data/Run_60264_Waveform.root", "calibration_constant_output/calibration_bic_output_Run60264_layer1.root", "Sim/3x8_3GeV_CERN_hist.root", "energy_calibration_output/energy_calibration_QC_Run60264_layer1.root", 3.0, 1)'
```
**출력**:
- `energy_calibration_output/energy_calibration_QC_Run60264_layer1.root`
- `energy_calibration_output/energy_calibration_QC_Run60264_layer1.png`
- `energy_calibration_output/total_energy_comparison_Run60264_layer1.png`

## 매개변수 설명

- `targetLayer`: 분석할 층 (0-2, 0이 아래층, 1이 가운데 층)
- `adcThreshold`: ADC 임계값 (기본값: 0)
- `beamEnergyGeV`: 빔 에너지 (GeV, 캘리브레이션 계산용, 기본값: 3.0 GeV)

## 파일 구조

```
202507_PS_prompt_analysis/
├── Data/                          # 데이터 파일
│   └── Run_60264_Waveform.root
├── Sim/                           # 시뮬레이션 파일
│   └── 3x8_3GeV_CERN_hist.root
├── calibration_constant_output/    # 캘리브레이션 상수 출력
├── energy_calibration_output/     # 에너지 캘리브레이션 출력
├── intADC_QA.C                   # ADC QA 매크로
├── calibration_bic.C             # 캘리브레이션 상수 계산 매크로
├── energy_calibration_bic.C      # 에너지 캘리브레이션 매크로
├── caloMap.h                     # 채널 매핑 헤더 (new mapping)
└── README.md                     # 이 파일
```

## 분석 워크플로우

1. **QA 단계**: `intADC_QA.C`로 데이터 품질 확인
2. **캘리브레이션 단계**: `calibration_bic.C`로 캘리브레이션 상수 계산
3. **에너지 캘리브레이션 단계**: `energy_calibration_bic.C`로 최종 에너지 캘리브레이션

## 주요 수정사항

- **New Mapping**: `caloMap.h` 사용으로 업데이트된 채널 매핑
- **3x8 구조**: 3층×8열 검출기 구조로 변경
- **Target Layer**: 기본값을 1(가운데 층)로 설정
- **Beam Energy**: 3 GeV로 변경 (3x8 시뮬레이션에 맞춤)
- **시뮬레이션 파일**: `Sim/3x8_3GeV_CERN_hist.root` 사용
- **시뮬레이션 층**: 2층(코드상 1층, 가운데 층) 고정 사용

## 주의사항

- 모든 매크로는 92개 채널이 모두 활성화된 이벤트만 처리
- 시뮬레이션 파일은 `Sim/3x8_3GeV_CERN_hist.root` 형식 필요
- 데이터 파일은 `Data/Run_XXXXX_Waveform.root` 형식 필요
- New mapping을 사용하므로 `caloMap.h` 헤더 파일 필요