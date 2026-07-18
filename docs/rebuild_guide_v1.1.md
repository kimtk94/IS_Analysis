# IS_Analysis 재구축 가이드

> 상태: 재구축 설계 기준 v1.1  
> 작성일: 2026-07-18  
> 적용 대상: Google Drive `IS_Analysis`의 코드, 설정, 결과 및 `Backup`  
> 기본 원칙: 기존 폴더는 근거 보존용으로 동결하고, 새 구현은 `IS_Analysis_v2`와 GitHub에서 시작한다.

## 1. 목적

이 문서는 Codex가 `IS_Analysis`를 다시 설계하고 구현할 때 따라야 할 단일 기준이다. 목표는 단순히 현재 스크립트를 고치는 것이 아니라 다음을 만족하는 재현 가능한 cis-pQTL–ischemic stroke 분석 파이프라인을 만드는 것이다.

- 원시 자료, 중간 산출물, 최종 결과의 계보를 추적한다.
- 파일 존재 여부가 아니라 내용과 생물정보학적 조건으로 성공 여부를 판단한다.
- 중단 후 안전하게 재시작할 수 있다.
- 10개 유전자의 EUR/EAS 제한행 smoke test를 통과한 뒤 범위를 확장한다.
- genome build, variant identity, allele direction을 분석 전에 확정한다.
- cis instrument selection, LD clumping, MR, sensitivity analysis, colocalization, replication을 서로 분리된 단계로 실행한다.
- 코드·테스트·설정은 GitHub를, 대용량 데이터·실행 결과는 Google Drive를 기준 저장소로 사용한다.

## 2. 이번 검토 범위

검토에 포함한 자료:

- 최상위 `Backup`, `R`, `scripts`, `config`, `requirements`, `results` 구조
- Python/R 코드와 설정 파일 및 시점별 backup
- `audit_report_20260718_123115.xlsx`
- `qc_integrity_report.xlsx`
- original/recovery/recovery_exact batch 구조와 대표 실행 로그
- 주요 타깃 유전자 상태와 batch별 row count

내용 검토에서 제외한 자료:

- 대용량 `.tar`, `.gz` 등 raw data 본문
- 동일 패턴의 반복 로그 전체 본문
- 수천 개 batch/per-gene 폴더의 동일 구조 파일

대용량 파일은 파일명, 경로, 크기, 상태 메타데이터만 사용했다. 반복 결과는 전체 개수와 상태를 집계하고 대표 성공·실패·header-only 사례를 확인했다.

## 3. 현재 상태 요약

### 3.1 확인된 구조

- original exposure batch: 651개
- recovery exact batch 폴더: 641개
- original batch log: 약 652개
- recovery exact log: 약 641개
- 최신 audit 대상 source: 2,617개
- 최신 audit에서 raw가 실제로 존재한다고 기록된 source: 1,845개
- 처리 후 raw가 삭제되었다고 기록된 source: 772개

### 3.2 최신 audit의 표면적 결과

| 지표 | 값 |
|---|---:|
| expected sources | 2,617 |
| raw download PASS | 0 |
| raw download REVIEW | 2,617 |
| prepare complete | 2,617 |
| final ready | 0 |
| final review | 2,617 |

이 결과는 그대로 신뢰하면 안 된다. `raw download PASS=0`과 `prepare complete=2,617`이 동시에 나온 원인은 아래 결함 때문이다.

### 3.3 batch 실제 row 상태

| 구분 | batch 수 |
|---|---:|
| original output가 1행 이상 | 480 |
| original output가 0행(header-only) | 171 |
| recovery output가 1행 이상 | 0 |
| recovery output가 0행(header-only) | 640 |
| original과 recovery 모두 0행 | 171 |

즉, recovery 실행 로그의 `Completed batch`는 유효한 variant가 생성되었다는 의미가 아니다. 다수 recovery 파일은 24개 컬럼의 header만 가진 258-byte 파일이다.

## 4. 중대한 결함

### P0-1. raw 파일 탐색 경로가 실제 구조와 맞지 않는다

`audit_config.json`의 `raw_dirs`는 `data/rawdata/pqtl/selected_targets`이고, `find_raw_file()`은 여기에 파일명만 붙인다. 실제 파일은 `selected_targets/EUR/<source_file>` 또는 ancestry 하위 폴더에 있다.

재구축 기준:

- raw key는 `(dataset_id, ancestry, source_file)`로 정의한다.
- raw 경로는 config의 ancestry-aware template으로 생성한다.
- `RAW_PRESENT`, `RAW_DELETED_PROCESSED_VERIFIED`, `RAW_MISSING_UNRECOVERED`를 구분한다.

### P0-2. header-only 파일이 성공으로 판정된다

현재 `output_basic_qc()`는 header를 읽을 수 있고 필수 컬럼이 존재하면 PASS로 본다. row count는 기록하지만 성공 조건에 사용하지 않는다. `04_audit_prepare_batches.py`는 batch output가 읽히기만 하면 그 batch의 모든 source를 complete로 표시한다.

재구축 기준:

- 파일 존재, schema valid, row count, source coverage를 별도 필드로 저장한다.
- `SUCCESS_NONEMPTY`, `SUCCESS_EMPTY_BIOLOGICAL`, `FAILED_SCHEMA`, `FAILED_RUNTIME`, `SKIPPED_EXISTING_VALID`을 구분한다.
- header-only는 자동 성공이 아니다.
- `SUCCESS_EMPTY_BIOLOGICAL`은 입력을 정상적으로 읽었고, 각 filter 전후 row count와 0행 사유가 기록된 경우에만 허용한다.
- batch 성공은 source별 coverage를 확인한 뒤 계산한다.

### P0-3. recovery 로그의 성공과 결과가 일치하지 않는다

예: `batch_651`은 SWAP70과 TFPI raw를 처리하고 `Completed batch`를 남겼지만, original과 recovery output 모두 0행, 258 bytes이며 최신 audit은 두 source를 `prepare_complete=True`로 기록한다.

재구축 기준:

- subprocess return code 0은 실행 성공일 뿐 분석 성공이 아니다.
- stage 종료 시 output contract를 검증하고 검증 실패 시 non-zero로 종료한다.
- 각 유전자에 `raw_rows`, `standardized_rows`, `cis_rows`, `p_threshold_rows`, `strong_rows`, `clumped_rows`를 기록한다.

### P0-4. genome build가 다른데 position만으로 병합한다

- UKB-PPP exposure: GRCh38
- GIGASTROKE outcome: GRCh37
- 현재 harmonization: exposure variant ID에서 position만 추출한 뒤 outcome `base_pair_location`과 병합

재구축 기준:

- 모든 downstream 자료를 GRCh38로 통일한다.
- outcome GRCh37은 검증된 chain file로 liftover한다.
- 변환 전후 `chr`, `pos`, `ref`, `alt`, build를 모두 보존한다.
- join key는 최소 `(build, chr, pos, effect_allele, other_allele)` 검증을 포함한다.
- liftover 성공률, unmapped 수, duplicate mapping 수, reference allele 불일치 수를 보고한다.

참고: UCSC liftOver: <https://genome.ucsc.edu/cgi-bin/hgLiftOver>

### P0-5. instrument와 harmonization schema가 연결되지 않는다

`03_select_instruments.R`의 주요 출력은 `SNP`, `protein_id`, `F_stat`이다. `04_harmonize_data.R`은 `snp`, `gene_symbol`, `f_statistic` 등을 기대한다.

재구축 기준:

- 모든 단계가 하나의 canonical schema를 사용한다.
- stage마다 입력 schema를 검증하고 예상치 못한 alias를 자동 추정하지 않는다.
- alias 변환은 dataset adapter에서 한 번만 수행한다.

### P0-6. 설정에 적힌 분석과 실제 구현이 다르다

설정에는 cis-only instrument selection, LD clumping, 여러 MR 방법, sensitivity, coloc, replication이 계획되어 있으나, 현재 구현은 lead SNP 선택, 단순 allele reversal, Wald ratio/간이 IVW, EAS Wald replication 일부에 가깝다.

재구축 기준:

- `method_plan.json`은 계획 문서가 아니라 실행 가능한 stage registry로 바꾼다.
- `planned`, `implemented`, `tested`, `executed`를 구분한다.
- 최종 보고서에는 실행하지 않은 분석을 명시한다.

### P0-7. 활성 코드의 source of truth가 없다

핵심 prepare/recovery 코드가 최상위 `Backup`에 있고, `scripts`에는 일부 downstream와 audit 코드만 있다. patch 스크립트와 시점별 backup이 혼재한다.

재구축 기준:

- GitHub가 유일한 코드 source of truth다.
- Drive의 `Backup`은 읽기 전용 참고 자료로만 사용한다.
- patch script를 운영 코드로 사용하지 않는다. 필요한 변경은 본 코드와 테스트에 반영한다.
- 모든 run manifest에 Git commit SHA를 기록한다.

### P1-1. cis filter가 사실상 실행되지 않았다

대표 로그마다 `No gene coordinate table found. Cis filter will be skipped.`가 반복된다. 그럼에도 결과는 cis-pQTL pipeline처럼 취급된다.

재구축 기준:

- gene coordinate reference가 없으면 instrument stage를 실패 처리한다.
- cis window와 사용한 transcript/TSS 정의를 기록한다.
- `is_cis=NA`인 상태는 통과시키지 않는다.

### P1-2. 단일 lead SNP 고정은 계획된 MR와 맞지 않는다

현재 instrument 코드가 protein별 1개 SNP만 남기므로 대부분 Wald ratio만 가능하다.

재구축 기준:

- genome-wide significant cis variants 전체에서 LD clumping 후 독립 SNP들을 유지한다.
- 단일 lead 결과는 별도 보조 분석으로만 만든다.
- population에 맞는 LD reference를 명시한다. 참고: <https://mrcieu.github.io/ieugwasr/articles/guide.html>

### P1-3. harmonization이 불완전하다

현재 코드는 direct match와 swapped allele만 처리한다. complement strand, palindromic SNP, allele frequency 기반 판정, 불일치 drop 사유가 없다.

재구축 기준:

- effect를 동일 allele 기준으로 정렬한다.
- palindromic SNP는 EAF와 사전 정의한 임계값으로 판정한다.
- match, swap, complement, palindromic_drop, mismatch_drop 수를 보고한다.
- 가능하면 검증된 `TwoSampleMR::harmonise_data()`와 결과를 대조한다. 참고: <https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html>

### P1-4. coloc에 필요한 지역 전체 통계량이 보존되지 않는다

prepare 단계에서 p-value와 F-statistic으로 먼저 걸러진 데이터만 유지하면 coloc에 필요한 locus 주변 전체 SNP 정보를 잃는다.

재구축 기준:

- `standardized_full`과 `instrument_candidates`를 분리한다.
- coloc은 lead instrument 테이블이 아니라 동일 지역의 충분한 공통 SNP를 사용한다.
- 기본 `coloc.abf`와 prior sensitivity를 수행하고, 다중 신호가 의심되면 SuSiE 기반 방법을 고려한다. 참고: <https://chr1swallace.github.io/coloc/articles/a01_intro.html>

### P1-5. 환경이 실행 중 변경된다

로그에서 R package가 분석 실행 도중 설치된다. Python requirements와 R package/버전 lock이 분리되어 있지 않다.

재구축 기준:

- 실행 전 `environment check`를 수행한다.
- Python은 lock file, R은 `renv.lock`을 사용한다.
- 분석 중 package 설치를 금지한다.
- Colab runtime과 package version을 run manifest에 기록한다.

## 5. 재사용할 좋은 요소

- 임시 파일 작성 후 `os.replace()`하는 atomic manifest update
- source file을 정확히 지정하는 recovery filter
- Google Drive raw를 Colab local disk로 staging한 뒤 처리
- tar 전체를 풀지 않고 필요한 member를 stream으로 읽는 방식
- `current`와 timestamped snapshot을 함께 남기는 방식
- raw 파일 크기와 manifest metadata를 비교하는 QC
- 모든 실행을 dry-run으로 먼저 보여주는 방식
- 단일 gene/batch로 제한하는 CLI 옵션

단, 기존 구현을 그대로 복사하지 말고 v2 schema와 validation gate에 맞춰 테스트 후 이식한다.

## 6. 저장소 역할 분리

| 위치 | 역할 | 허용 내용 |
|---|---|---|
| GitHub | 코드의 기준 저장소 | source, tests, config template, AGENTS.md, 문서 |
| Google Drive | 데이터와 실행 산출물 | raw, reference, processed, results, run manifests |
| Colab local disk | 임시 고속 작업공간 | staged raw, decompression temp, cache |

규칙:

- Drive에서 Python/R 코드를 직접 운영하지 않는다.
- Colab은 GitHub repo를 clone하고 Drive를 mount한다.
- output에는 Git commit SHA, config hash, input checksum을 기록한다.
- raw와 결과를 GitHub에 commit하지 않는다.

## 7. 권장 `IS_Analysis_v2` 구조

```text
IS_Analysis_v2/
├── AGENTS.md
├── README.md
├── pyproject.toml
├── uv.lock 또는 requirements.lock
├── renv.lock
├── config/
│   ├── project.yaml
│   ├── datasets.yaml
│   ├── methods.yaml
│   └── schemas.yaml
├── workflow/
│   └── run_pipeline.py
├── src/is_analysis/
│   ├── inventory.py
│   ├── checkpoint.py
│   ├── schemas.py
│   ├── raw_qc.py
│   └── reporting.py
├── R/
│   ├── standardize_exposure.R
│   ├── standardize_outcome.R
│   ├── select_instruments.R
│   ├── harmonize.R
│   ├── run_mr.R
│   ├── run_sensitivity.R
│   └── run_coloc.R
├── tests/
│   ├── fixtures/test_panel_10genes/
│   ├── test_inventory.py
│   ├── test_checkpoint.py
│   ├── test_schema.py
│   └── test_panel_smoke.py
└── drive_data/
    ├── raw/
    ├── reference/
    ├── interim/
    ├── processed/
    ├── results/
    └── state/
```

`drive_data`는 실제 repo 내부 데이터가 아니라 Google Drive mount 경로를 config로 연결한다.

## 8. Canonical summary-statistics schema

모든 exposure/outcome adapter는 아래 컬럼을 생성한다.

| 컬럼 | 설명 |
|---|---|
| `dataset_id` | 데이터셋 고유 ID |
| `trait_id` | protein/phenotype ID |
| `gene_symbol` | HGNC symbol |
| `protein_id` | UniProt/Olink/aptamer ID |
| `assay_id` | OID 또는 assay ID |
| `ancestry` | EUR, EAS 등 |
| `build` | GRCh37 또는 GRCh38 |
| `chr` | `chr` prefix 없는 chromosome |
| `pos` | 해당 build의 1-based position |
| `variant_id` | 원본 variant ID |
| `rsid` | 가능한 경우 rsID |
| `effect_allele` | beta 기준 allele |
| `other_allele` | 반대 allele |
| `eaf` | effect allele frequency |
| `beta` | effect estimate |
| `se` | standard error |
| `pval` | 0–1 범위 p-value |
| `n` | sample size |
| `source_file` | 원본 파일명 |
| `source_sha256` | 원본 checksum |

추가 규칙:

- 컬럼명은 소문자 snake_case로 고정한다.
- `LOG10P`는 adapter에서 `10^(-LOG10P)`로 한 번만 변환한다.
- 원본 컬럼명과 변환 규칙을 adapter report에 기록한다.
- `variant_key`는 build와 allele normalization 이후 생성한다.

## 9. 단계별 파이프라인

### Stage 00 — inventory 및 환경 검증

입력: dataset config, Drive root, code version.  
출력: `run_manifest.json`, `inventory.tsv`, `environment_report.json`.

Gate:

- 필수 config와 reference 존재
- Python/R package lock 일치
- Drive write test 성공
- secrets가 코드/config에 없음

### Stage 01 — raw catalog 및 integrity QC

검사: ancestry-aware 경로, 파일 크기와 remote metadata, tar readability, 내부 summary-stat 파일 선택 근거, 선택적 SHA-256.

raw가 처리 후 삭제된 경우 processed output checksum과 run manifest가 유효하면 `RAW_DELETED_PROCESSED_VERIFIED`로 인정한다.

### Stage 02 — exposure 표준화

- raw full summary statistics를 canonical schema로 변환한다.
- full standardized data를 보존한다.
- p-value, allele, beta/se, chromosome/position QC를 수행한다.
- 분석 filter를 이 단계에서 적용하지 않는다.

### Stage 03 — outcome 표준화 및 build 통일

- GIGASTROKE EUR/EAS를 canonical schema로 변환한다.
- GRCh37 → GRCh38 liftover를 수행한다.
- mapping 실패와 allele/ref mismatch를 보고한다.

### Stage 04 — cis instrument selection

- gene coordinate reference를 필수로 사용한다.
- ±1 Mb cis window를 명시적으로 적용한다.
- p-value, MAF, imputation INFO, F-statistic을 순서대로 적용한다.
- EUR discovery용 LD panel로 clumping한다.
- 독립 instrument 전체와 lead-only 보조 테이블을 모두 저장한다.

### Stage 05 — harmonization

- exposure/outcome build 일치 확인
- chromosome + position + allele 기반 match
- allele swap/complement/palindromic 처리
- EAF 기반 ambiguous SNP 제거
- drop reason을 variant별로 기록

### Stage 06 — primary MR

- 1 SNP: Wald ratio
- 2개 이상: IVW를 기본 결과로 생성
- 충분한 SNP가 있을 때 weighted median, MR-Egger, weighted mode 수행
- OR, 95% CI, raw p, FDR q-value를 기록
- divide-by-zero, invalid SE, duplicated SNP를 gate에서 차단

### Stage 07 — sensitivity 및 directionality

- Cochran Q heterogeneity
- MR-Egger intercept
- leave-one-out
- Steiger directionality
- reverse MR
- phenotype/pleiotropy scan

수행할 수 없는 경우 `NOT_RUN_INSUFFICIENT_INSTRUMENTS`로 명시한다.

### Stage 08 — colocalization

- 각 protein locus의 full regional exposure/outcome statistics 사용
- 공통 SNP 수와 region coverage 보고
- `coloc.abf` PP0–PP4와 prior sensitivity 저장
- 다중 causal signal이 의심되면 conditional/SuSiE 분석 경로 사용

### Stage 09 — replication

1. exposure replication: UKB-PPP signal을 deCODE/Fenland pQTL에서 재현
2. outcome replication: EUR instrument 효과를 EAS/FinnGen/추가 stroke GWAS에서 재현
3. 방향·유의성·effect-size consistency를 별도 평가

### Stage 10 — annotation과 biological context

- Open Targets, ChEMBL, DGIdb druggability
- brain/scRNA cell-type expression context
- stroke subtype별 결과
- pleiotropic locus warning

### Stage 11 — final report

필수 표:

- dataset provenance
- stage별 input/output/QC 수
- instrument strength와 LD 기준
- harmonization flow
- MR primary/sensitivity
- coloc posterior와 prior sensitivity
- replication consistency
- 제외 항목과 미실행 분석

MR 보고는 STROBE-MR 항목을 반영한다. 참고: <https://pubmed.ncbi.nlm.nih.gov/34698778/>

## 10. Checkpoint contract

각 stage/unit는 `state/<run_id>/<stage>/<unit>.json`을 atomic하게 기록한다.

```json
{
  "run_id": "20260718_eur_eas_10gene_head100",
  "stage": "02_standardize_exposure",
  "unit": "ABO__EUR__OID30675",
  "status": "SUCCESS_NONEMPTY",
  "input": {
    "path": "...",
    "size_bytes": 550615040,
    "sha256": "..."
  },
  "output": {
    "path": "...",
    "rows": 1234,
    "schema_version": "1.0",
    "sha256": "..."
  },
  "counts": {
    "raw_rows": 0,
    "standardized_rows": 0,
    "cis_rows": 0,
    "significant_rows": 0,
    "strong_rows": 0,
    "clumped_rows": 0
  },
  "code_commit": "...",
  "config_sha256": "...",
  "started_at": "...",
  "completed_at": "..."
}
```

재시작 규칙:

- input hash, config hash, code version, output validation이 모두 같을 때만 skip한다.
- output 파일 존재만으로 skip하지 않는다.
- 실패 unit만 다시 실행할 수 있어야 한다.
- partial output은 최종 경로에 두지 않는다.

## 11. 첫 번째 smoke test: 10개 유전자 × EUR/EAS × 최대 100행

테스트 목록은 `ukb_ppp_complete_pair_genes.current.txt`에서 EUR/EAS complete pair로 기록된 유전자를 정렬 순서대로 선택한다. 선택 편향을 막기 위해 목록과 순서를 config에 고정한다.

| 순번 | gene |
|---:|---|
| 1 | A1BG |
| 2 | AAMDC |
| 3 | AARSD1 |
| 4 | ABCA2 |
| 5 | ABHD14B |
| 6 | ABL1 |
| 7 | ABO |
| 8 | ABRAXAS2 |
| 9 | ACAA1 |
| 10 | ACADM |

테스트 단위와 행 제한:

- ancestry: `EUR`, `EAS`
- 단위: `(gene_symbol, ancestry)` 20개
- source: 각 유전자에서 EUR/EAS에 공통으로 존재하는 동일 assay/source pair 1개
- 행 제한: 각 단위에서 header와 comment를 제외한 data row 최대 100행
- 전체 최대 data row: 2,000행
- sampling: 원천 파일 순서를 유지하는 deterministic head 방식
- 저장: header + 최대 100행만 압축 저장하고 full raw archive는 Drive에 새로 보존하지 않는다.

압축 tar에서는 원격 파일을 행 단위로 직접 요청할 수 없는 경우가 있다. 이때 archive member를 스트리밍으로 읽되 100번째 data row에서 reader를 중단하고, 전체 archive를 영구 저장하지 않는다. 다운로드 byte 수와 실제 추출 row 수를 checkpoint에 각각 기록한다.

실행 config는 `IS_ANALYSIS_TEST_SPEC_10GENES.yaml`을 사용한다.

실행 순서:

1. 10개 유전자의 paired EUR/EAS source manifest를 만든다.
2. 각 gene×ancestry에 source가 정확히 1개 선택됐는지 검증한다.
3. ancestry-aware raw 경로, 크기, archive member를 기록한다.
4. member를 스트리밍하고 header/comment를 분리한다.
5. data row 100개에서 읽기를 중단한다.
6. canonical schema로 변환하고 원천 순서를 보존한다.
7. unit별 checkpoint와 전체 test summary를 기록한다.
8. 같은 config로 재실행하여 hash와 row 결과가 동일한지 확인한다.

smoke acceptance criteria:

- 10개 유전자 모두 EUR와 EAS가 있어 정확히 20개 unit가 생성된다.
- 각 unit의 `0 < extracted_data_rows <= 100`이다.
- header와 comment는 100행 제한에 포함하지 않는다.
- p-value가 `(0, 1]` 범위다.
- allele은 A/C/G/T이고 effect/other가 다르다.
- ancestry, build, chromosome/position, source identity가 명시된다.
- source별 출력과 checkpoint hash가 생성된다.
- 두 번 실행한 출력 hash가 동일하다.
- full archive 또는 100행 초과 결과가 test output에 남지 않는다.

이 제한행 smoke test의 성공은 다운로드, archive reader, parser, ancestry routing, canonical schema, checkpoint 및 재시작 동작만 보장한다. 100행 표본은 cis coverage나 유효 instrument를 보장하지 않으므로 MR, sensitivity analysis, coloc의 과학적 PASS 근거로 사용하지 않는다. 과학 분석 gate는 full-data 검증 단계에서 별도로 통과해야 한다.

## 12. Codex 작업 규칙

Codex는 아래 순서를 지킨다.

1. 이 문서, `AGENTS.md`, config를 먼저 읽는다.
2. GitHub와 Drive의 대상 경로를 확정한다.
3. 기존 `IS_Analysis`는 수정하거나 삭제하지 않는다.
4. 대용량 raw는 metadata만 확인한다. TEST 실행은 config에 고정된 10개 유전자, EUR/EAS, 단위별 최대 100 data row로 제한한다.
5. 모든 small code/config/backup을 비교해 기능과 차이를 기록한다.
6. backup의 기능을 채택할 때는 테스트로 검증한다.
7. schema contract와 gate를 먼저 구현한다.
8. 10-gene EUR/EAS head-100 smoke test를 실행한다.
9. 실패하면 다음 단계로 진행하지 않고 원인과 재현 명령을 남긴다.
10. 제한행 smoke test 통과 후 full-data 과학 gate를 별도로 검증하고, 그 뒤에만 15개 후보 또는 전체 batch 확장을 제안한다.

기본 대용량 스킵 규칙:

- 50 MB 초과 파일은 자동 본문 검토에서 제외한다.
- `.tar`, `.gz`, `.bgz`, `.bgen`, `.bed`, `.vcf.gz`는 raw로 간주한다.
- 단, test config가 선택한 20개 gene×ancestry source는 스트리밍 head-100 추출에 한해 예외다.

금지 사항:

- header-only 파일을 성공으로 표시
- genome build가 다른 상태에서 position join
- ancestry를 무시한 raw lookup
- 실행 중 package install
- 원본 raw/result 삭제
- backup 파일을 검증 없이 활성 코드로 복사
- 계획만 존재하는 분석을 수행 완료로 보고

## 13. 단계별 구현 우선순위

### Milestone 1 — 신뢰 가능한 기반

- GitHub repo와 `IS_Analysis_v2` skeleton
- environment lock
- canonical schema
- ancestry-aware inventory/raw QC
- checkpoint contract
- 10-gene EUR/EAS head-100 tar adapter smoke test

### Milestone 2 — 분석 입력 완성

- UKB-PPP full standardization
- GIGASTROKE standardization
- GRCh37 → GRCh38 liftover
- gene coordinate reference
- cis + LD clumping

### Milestone 3 — MR와 검증

- harmonization
- primary MR
- sensitivity/directionality
- coloc + prior sensitivity
- EUR/EAS 및 external pQTL replication

### Milestone 4 — 보고와 확장

- 15개 후보 분석
- stroke subtype 확장
- druggability/scRNA context
- STROBE-MR 기반 final report

## 14. 완료 정의

`IS_Analysis_v2`는 다음 조건을 모두 만족할 때만 전체 실행 준비 완료로 본다.

- 10-gene EUR/EAS head-100 smoke test의 20개 unit 전체 통과
- full-data 과학 gate와 제한행 smoke 결과의 명확한 분리
- schema와 build mismatch 테스트 통과
- header-only false-success regression test 통과
- raw deleted-after-processed 상태 테스트 통과
- restart/idempotency 테스트 통과
- MR 단위 테스트와 최소 synthetic truth test 통과
- coloc 입력 coverage와 prior sensitivity 출력 확인
- run manifest에 data/code/config provenance 완성
- GitHub CI에서 small fixture 테스트 통과

## 15. 참고 방법론

- TwoSampleMR harmonization: <https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html>
- ieugwasr LD clumping: <https://mrcieu.github.io/ieugwasr/articles/guide.html>
- coloc overview: <https://chr1swallace.github.io/coloc/articles/a01_intro.html>
- coloc prior sensitivity: <https://chr1swallace.github.io/coloc/articles/a04_sensitivity.html>
- UCSC liftOver: <https://genome.ucsc.edu/cgi-bin/hgLiftOver>
- MungeSumstats paper: <https://pubmed.ncbi.nlm.nih.gov/34601555/>
- STROBE-MR: <https://pubmed.ncbi.nlm.nih.gov/34698778/>
