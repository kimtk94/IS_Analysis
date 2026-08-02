# 새 Git 저장소 구축 및 데이터 이전 가이드

이 문서는 현재 `IS_Analysis`에서 정리한 UKB-PPP EUR/EAS pQTL workflow를
새 Git 저장소에서 다시 시작할 때 필요한 코드, Google Drive 데이터 위치,
복사 범위 및 실행 순서를 정리한다.

## 1. 코드와 데이터 루트 분리

새 저장소에서는 Colab의 일시적인 코드 clone과 Drive의 영구 데이터를 분리한다.

```bash
export CODE_ROOT="/content/IS_Analysis_V3"
export WORK_ROOT="/content/drive/MyDrive/IS_Analysis_V3"
export SOURCE_WORK_ROOT="/content/drive/MyDrive/IS_Analysis_V2"
```

| 루트 | 역할 | 영구 보존 |
|---|---|---|
| `CODE_ROOT` | 새 Git clone, Python/R/shell 코드, committed fixture | 아니요 |
| `WORK_ROOT` | manifest, reference, test data, QC, 분석 결과 | 예 |
| `SOURCE_WORK_ROOT` | V2에서 재사용할 작은 metadata/reference의 읽기 전용 원본 | 예 |

`CODE_ROOT` 아래에 실제 pQTL/GWAS 원자료나 실행 결과를 두지 않는다. Colab
세션이 종료되면 이 경로가 사라질 수 있다.

## 2. 새 Git 저장소에 가져갈 코드

최소한 다음 tracked 파일과 디렉터리를 새 저장소로 옮긴다.

```text
AGENTS.md
README.md
requirements.txt
scripts/
workflow/
config/
docs/
tests/fixtures/
```

다음 항목은 Git에 commit하지 않는다.

```text
data/rawdata/
results/
*.tar
실제 UKB-PPP/GWAS 전체 summary statistics
Google Drive 인증 정보와 토큰
```

새 저장소를 만든 뒤에는 먼저 현재 저장소의 tracked 파일만 복사하거나 Git
history를 이전하고, `.gitignore`가 대용량 raw/output을 제외하는지 확인한다.

## 3. 현재 workflow 요약

### Stage A — manifest 및 reference setup

1. EUR/EAS source archive metadata로 `ukb_ppp_download_manifest.tsv`를 만든다.
2. manifest의 gene symbol로 GRCh38 gene-coordinate table을 만든다.
3. ancestry, gene, source archive identity와 checksum/size metadata를 보존한다.

필수 입력:

```text
${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv
${WORK_ROOT}/data/reference/gene_coordinates_hg38.tsv
```

### Stage B — 15-gene paired batch 구성

`scripts/ukb_ppp_batch_manifest_runner_fast.py`는 EUR과 EAS가 모두 있는 gene을
기본 15개씩 묶는다. 각 batch에서 필요한 source archive만 staging하고 검증한
뒤 ancestry별 R preparation을 실행한다.

주요 checkpoint:

```text
${WORK_ROOT}/results/qc/batch_pipeline/batch_manifest.tsv
${WORK_ROOT}/results/qc/batch_pipeline/execution_plan.tsv
${WORK_ROOT}/results/qc/batch_pipeline/batch_progress.tsv
${WORK_ROOT}/results/qc/batch_pipeline/downloads/
${WORK_ROOT}/results/qc/batch_pipeline/processing_logs/
${WORK_ROOT}/results/qc/batch_pipeline/raw_cleanup/
```

### Stage C — canonical standardization과 instrument selection

`scripts/01_prepare_exposure_fast.R`은 각 tar의 summary member를 스트리밍으로
읽고 두 종류의 산출물을 만든다.

1. canonical standardized summary
2. GRCh38 cis window, p-value 및 F-statistic 조건을 통과한 instrument candidate

경로:

```text
${WORK_ROOT}/results/standardized/pqtl/{EUR,EAS}/batch_###/
${WORK_ROOT}/results/instrument_candidates/{EUR,EAS}/
${WORK_ROOT}/results/exposure_batches/{EUR,EAS}/
```

`results/exposure_batches`는 legacy 호환용 filtered candidate이며 full summary가
아니다.

### Stage D — downstream analysis

현재 repository의 downstream workflow는 다음을 분리된 단계로 실행하도록
구성되어 있다.

- ancestry-matched LD clumping 및 fine-mapping:
  `workflow/ancestry_ld_finemap.py`
- exposure/outcome harmonization, MR 및 sensitivity checkpoint:
  `workflow/causal_checkpoint_analysis.py`
- brain eQTL colocalization:
  `workflow/annotation/brain_eqtl_colocalization.py`
- GIGASTROKE outcome standardization:
  `workflow/gigastroke_outcome_adapter.py`

MR screening에는 filtered cis instrument가 필요하다. Fine-mapping과 coloc에는
선정 locus의 비유의 variant를 포함한 regional full summary가 필요하다. 저장
공간을 줄이려면 모든 protein의 genome-wide summary를 영구 저장하는 대신 MR로
선정한 locus의 regional summary만 다시 materialize할 수 있다.

## 4. V2에서 V3로 복사할 데이터

### 필수로 복사

다음 두 파일은 작고 새 workflow를 시작하는 데 필수다.

```text
${SOURCE_WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv
${SOURCE_WORK_ROOT}/data/reference/gene_coordinates_hg38.tsv
```

복사 명령:

```bash
mkdir -p \
  "${WORK_ROOT}/data/metadata" \
  "${WORK_ROOT}/data/reference"

cp -n \
  "${SOURCE_WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv" \
  "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv"

cp -n \
  "${SOURCE_WORK_ROOT}/data/reference/gene_coordinates_hg38.tsv" \
  "${WORK_ROOT}/data/reference/gene_coordinates_hg38.tsv"
```

`cp -n`을 사용하여 이미 검토한 V3 파일을 덮어쓰지 않는다. 복사 후에는
coordinate의 `genome_build` 값과 manifest의 EUR/EAS source 수를 다시 확인한다.

### 선택적으로 복사

과거 진행상태를 이어서 사용할 때만 다음 QC 문서를 복사한다.

```text
${SOURCE_WORK_ROOT}/results/qc/batch_pipeline/batch_manifest.tsv
${SOURCE_WORK_ROOT}/results/qc/batch_pipeline/downloads/
${SOURCE_WORK_ROOT}/results/qc/batch_pipeline/raw_cleanup/
${SOURCE_WORK_ROOT}/results/qc/batch_pipeline/variant_extraction_audit/
${SOURCE_WORK_ROOT}/results/qc/batch_pipeline/retained_variant_audit/
${SOURCE_WORK_ROOT}/results/qc/batch_pipeline/multi_gene_filter_audit/
```

새로운 15-gene batch plan은 과거 10-gene batch ID와 의미가 다를 수 있으므로,
기존 `batch_manifest.tsv`를 그대로 재시작 상태로 사용하지 않는다. provenance
참고용 별도 디렉터리에 복사한다.

```bash
mkdir -p "${WORK_ROOT}/results/qc/legacy_v2_reference"
cp -a \
  "${SOURCE_WORK_ROOT}/results/qc/batch_pipeline/." \
  "${WORK_ROOT}/results/qc/legacy_v2_reference/"
```

### 기본적으로 복사하지 않음

다음 대용량 또는 stale 가능성이 있는 항목은 새 workflow에 자동 복사하지 않는다.

```text
${SOURCE_WORK_ROOT}/data/rawdata/pqtl/selected_targets/
${SOURCE_WORK_ROOT}/results/exposure_batches/
${SOURCE_WORK_ROOT}/results/standardized/pqtl/
${SOURCE_WORK_ROOT}/results/instrument_candidates/
```

필요한 source archive는 download manifest를 기준으로 새 batch가 선택한 것만
다운로드한다. 기존 filtered exposure는 audit/reference로 사용할 수 있지만 새
15-gene batch의 완료 증거로 승계하지 않는다.

## 5. Drive V3 test workspace setup

Google Drive를 mount한 뒤 새 Git clone에서 다음을 실행한다.

```bash
export CODE_ROOT="/content/IS_Analysis_V3"
export WORK_ROOT="/content/drive/MyDrive/IS_Analysis_V3"
export SOURCE_WORK_ROOT="/content/drive/MyDrive/IS_Analysis_V2"

cd "${CODE_ROOT}"
bash scripts/setup_ido1_test_drive.sh
source "${WORK_ROOT}/ido1_test.env"
```

setup helper는 V3 디렉터리가 없으면 생성하고, V3에 없는 manifest와 coordinate만
V2에서 복사한다. raw 및 derived output은 복사하지 않는다.

생성되는 주요 test 경로:

```text
${WORK_ROOT}/data/rawdata/pqtl/selected_targets/{EUR,EAS}/
${WORK_ROOT}/results/qc/ido1_test_pipeline/
${WORK_ROOT}/results/test/ido1/exposure_batches/
${WORK_ROOT}/results/test/ido1/standardized/pqtl/
${WORK_ROOT}/results/test/ido1/instrument_candidates/
```

## 6. IDO1 focused test 실행

다음 실행은 IDO1이 포함된 15-gene batch 하나만 선택한다. IDO1은 압축 해제된
summary stream 기준 20 MB, 나머지 gene은 header 포함 1,000줄로 제한한다.

```bash
python3 -u "${CODE_ROOT}/scripts/ukb_ppp_batch_manifest_runner_fast.py" \
  --base "${WORK_ROOT}/data/rawdata/pqtl/selected_targets" \
  --qc-dir "${IDO1_TEST_QC_DIR}" \
  --outdir "${IDO1_TEST_OUTDIR}" \
  --standardized-dir "${IDO1_TEST_STANDARDIZED_DIR}" \
  --instrument-dir "${IDO1_TEST_INSTRUMENT_DIR}" \
  --download-manifest "${WORK_ROOT}/data/metadata/ukb_ppp_download_manifest.tsv" \
  --gene-coordinate-file "${WORK_ROOT}/data/reference/gene_coordinates_hg38.tsv" \
  --batch-size 15 \
  --focus-gene IDO1 \
  --focus-max-bytes 20000000 \
  --other-max-file-lines 1000 \
  --run \
  --stop-on-error
```

tar 전송과 archive integrity 검증에는 전체 tar가 필요하다. 20 MB 제한은 tar
내부의 압축 해제된 IDO1 summary stream에 적용된다. test 검증이 끝나기 전에는
raw 삭제 옵션을 추가하지 않는다.

## 7. 새 저장소에서 확인할 완료 조건

새 Git/Drive 조합이 준비된 뒤 다음을 확인한다.

1. `WORK_ROOT`의 모든 runtime 경로가 Drive 아래인지 확인한다.
2. 15-gene execution plan에서 IDO1 batch가 정확히 하나 선택되는지 확인한다.
3. EUR/EAS 각각의 selected source가 download audit에 기록되는지 확인한다.
4. 비-IDO1 canonical test 파일이 header 포함 1,000줄 이하인지 확인한다.
5. IDO1 canonical 입력 stream limit가 20,000,000 bytes로 로그에 기록되는지
   확인한다.
6. source별 standardization 및 instrument-selection status를 확인한다.
7. test 결과가 `CODE_ROOT/results`가 아니라 `WORK_ROOT/results/test/ido1`에
   생성되는지 확인한다.
8. committed fixture 검증은 실제 raw가 아니라 `tests/fixtures`만 사용하여
   `bash scripts/codex_smoke_test.sh`로 실행한다.

## 8. 재현성을 위해 기록할 값

각 실행에서 다음을 QC manifest 또는 실행 노트에 남긴다.

```text
새 Git remote URL
Git commit SHA
CODE_ROOT
WORK_ROOT
download manifest checksum
gene-coordinate checksum 및 genome build
batch size
focus gene
focus byte/line limits
p-value, F-statistic, cis-window 설정
실행 시작/종료 시간
raw cleanup 여부
```

이 정보를 남기면 새 Git 저장소의 코드와 Drive V3의 산출물을 명확히 연결할 수
있다.
