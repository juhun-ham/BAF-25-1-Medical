# BAF 25-1 의료팀 — 임상 발현 데이터를 활용한 임상환자의 약물 민감도 예측 및 임상분석

> 비어플(BAF) 2025년 1학기 의료팀 프로젝트 레포지토리입니다.
>
> **팀원:** 15기 이서은 · 16기 유영우 · 16기 함주헌

> 📄 **이 프로젝트는 이후 통계청 논문공모전 투고작으로 발전했습니다.**
> 논문: *「정밀 의료 및 신약 개발 정책 적용을 위한 약물 반응 예측 모델 연구 — 공공 전사체 데이터를 활용한 머신러닝 기반 환자 민감도 예측」*
> 논문 단계에서 방법론이 크게 정교화·확장되었으며, 그 차이는 아래 **[📄 논문 공모전 보완 버전](#-논문-공모전-보완-버전-프로젝트--논문)** 섹션에 정리했습니다.

---

## 📌 한 줄 요약

**세포주(CCLE) 약물반응 데이터로 "약물 민감/저항"을 학습**하고, 이 모델을 **실제 유방암(TCGA-BRCA) 환자에게 전이(transfer)** 하여 환자별 약물 반응을 예측합니다.
예측 결과는 **카이제곱 검정**(임상 변수와의 연관성)과 **생존분석**(민감군 vs 저항군 생존 차이)으로 검증합니다.

> 핵심 아이디어: *세포주에서 배운 약물 반응 패턴을 유전자 발현(특히 KEGG pathway 활성)을 매개로 환자에게 옮길 수 있는가?* → **정밀종양학(Precision Oncology) / 약물 재창출(Drug Repositioning)** 관점의 접근.

---

## 🎯 분석 목표

1. **약물 민감/저항 라벨링** — 세포주의 IC50 값으로 약물별 민감군·저항군을 정의.
2. **반응 예측 모델 구축** — 유전자 발현 기반 **pathway score**로 민감/저항을 분류.
3. **환자 단위 전이 예측** — 세포주에서 학습한 모델로 TCGA 환자의 약물 반응을 예측.
4. **임상적 검증** — 예측된 반응군이 실제 **임상 특성** 및 **생존 기간**과 유의하게 연관되는지 확인.

---

## 🗃️ 데이터 구성

본 프로젝트는 세 종류의 공개 생명정보 데이터를 결합합니다.

| 데이터 | 역할 | 설명 |
|--------|------|------|
| **CCLE** (Cancer Cell Line Encyclopedia) | 🧬 학습 피처 | 세포주 유전자 발현(RNA-seq read count) |
| **GDSC IC50** | 🏷️ 학습 라벨 | 세포주–약물 조합별 IC50 약물반응 값 (`preprocess_ic50.csv`, 약 25.7만 행) |
| **TCGA-BRCA** | 🎯 예측 대상 | 유방암 환자 유전자 발현 + 임상정보 + 생존정보 |
| **KEGG Pathway** | 🔗 도메인 지식 | 유전자–경로(pathway) 관계 (NodeInfo/EdgeInfo) → pathway score 계산용 |

- `IC50`은 약물이 세포 절반을 억제하는 농도 → 작을수록 민감. 분포가 right-skewed 하여 **`pIC50 = -log10(IC50)`** 로 변환(클수록 민감).
- TCGA 임상정보: `AJCC 병기`, `PAM50 분자아형`, `ER/PR/HER2 수용체 상태`, 생존(`OS`, `OS.time`) 등.

> ⚠️ `data/`의 일부 CSV(특히 CCLE pathway score)는 용량 문제로 샘플만 포함되어 있을 수 있습니다.

---

## 📂 레포지토리 구조

```
BAF-25-1-Medical/
├── data/preprocess/                       # 전처리 완료 데이터
│   ├── preprocess_ic50.csv                 # GDSC IC50 / pIC50 (라벨 원천)
│   ├── preprocess_CCLE_gene.csv            # CCLE 유전자 발현
│   ├── preprocess_TCGA_gene.csv            # TCGA 환자 유전자 발현
│   ├── preprocess_TCGA_clinical.csv        # TCGA 임상/생존 정보
│   ├── pathway_scores_CCLE.csv             # CCLE pathway score
│   └── pathwayscores_TCGA.csv              # TCGA pathway score
│
└── code/
    ├── preprocessing/                      # ① 전처리 파이프라인
    │   ├── EDA_Drug.ipynb                  #   - IC50 EDA & 필터링, pIC50 변환
    │   ├── EDA_TCGA.ipynb                  #   - TCGA 임상/생존 정제 (1,037명)
    │   ├── Gene_mapping.ipynb              #   - KEGG∩CCLE∩TCGA 공통 유전자 추출
    │   ├── TCGA_trans.ipynb                #   - TCGA Ensembl ID→유전자 심볼 변환
    │   ├── EdgeR_CCLE_gene.R               #   - CCLE read count edgeR 필터링
    │   └── Pathway_score.R                 #   - KEGG 기반 pathway score 계산
    │
    ├── Model/                              # ② DEG 라벨링 + 모델링
    │   ├── 설명서_READ_ME.txt               #   - 폴더/파일 구조 설명서
    │   ├── Pathway Score/                  #   - pathway score 계산(ver1 ic_label, ver2 ic50)
    │   ├── IC_50 Modeling + 카이제곱/        #   - IC50 z-score 기반 상위 10개 약물 모델
    │   │     └─ <약물명>/ ( .ipynb + train/test/valid + TCGA pathway score + 예측 )
    │   ├── IC_Label Modeling + 카이제곱/     #   - 이진분류 라벨 기반 약물 모델
    │   └── Target_sample.ipynb             #   - 예측 결과 × 임상변수 카이제곱 검정
    │
    └── Survival analysis/                  # ③ 약물별 생존분석 (Kaplan-Meier / log-rank)
        └── <약물명>_survival.ipynb
```

---

## 🔧 분석 파이프라인

### ① 전처리 (`code/preprocessing/`)

**1-1. 약물 데이터 정제 — `EDA_Drug.ipynb`**
- 세포주별 보유 약물 수를 확인 → 약물 수가 **50개 미만인 세포주는 제거**(학습 안정성 확보, rule of thumb).
- `IC50_PUBLISHED` 분포가 right-skewed → **`pIC50 = -log10(IC50)`** 로 정규화.

**1-2. CCLE 유전자 필터링 — `EdgeR_CCLE_gene.R`**
- RNA-seq read count를 정수화한 뒤, `edgeR::filterByExpr` 로 **저발현 유전자 제거**.

**1-3. 공통 유전자 매핑 — `Gene_mapping.ipynb`**
- **KEGG ∩ CCLE ∩ TCGA** 세 데이터셋에 모두 존재하는 공통 유전자만 추출.
- KEGG의 NodeInfo/EdgeInfo도 공통 유전자 기준으로 필터링(pathway 구조 정합성 유지).

**1-4. TCGA 발현 데이터 정합 — `TCGA_trans.ipynb`**
- TCGA RNA-seq(STAR-Counts)의 **Ensembl ID → 유전자 심볼** 변환(g:Profiler).
- 샘플 ID를 `TCGA-XX-XXXX-XX` 형식으로 정규화하여 임상 데이터와 매핑.

**1-5. TCGA 임상/생존 정제 — `EDA_TCGA.ipynb`**
- 결측 과다 컬럼, 생존분석 무관 변수(유전자 클러스터·행정 메타·연구 목적·중복 변수 등) 제거.
- 임상 데이터와 생존 데이터를 병합 → **공통 환자 1,037명** 확정.

### ② DEG 라벨링 & Pathway Score

**2-1. 민감/저항 라벨 정의** (`Pathway Score/pathway_score_ver2_ic50.R`)
- 약물별 `pIC50`에 **z-score 표준화** → `z > 0` 이면 **sensitive**, `z ≤ 0` 이면 **resistant**.

**2-2. DEG(차등발현유전자) 분석으로 약물 선별**
- 각 약물의 민감군 vs 저항군 간 유전자 발현 차이를 검정.
- 기준: **p-value < 0.05 & |log2 fold change| ≥ 2** 를 만족하는 유전자를 DEG로 카운트.
- **DEG가 많다 = 두 군의 생물학적 차이가 뚜렷하다 = 예측이 잘 된다**는 가설.
- DEG 개수 기준 **상위 10개 약물 선별** (`DEG_Count > 100`, `NA_Count < 200`).
- 두 트랙으로 진행: **IC_50**(z-score 연속값 기반) / **IC_Label**(이진분류 라벨 기반).

**2-3. Pathway Score 계산** (`Pathway_score.R`)
- KEGG pathway 단위로 유전자 발현을 종합한 **경로 활성 점수**를 계산(기준 유전자와의 상관 방향 × 발현 순위 기반).
- CCLE(세포주)와 TCGA(환자) 모두에 대해 동일 방식으로 산출 → 두 도메인을 잇는 **공통 피처 공간** 확보.

### ③ 모델링 (`Model/<트랙>/<약물명>/`)

약물 한 종마다 독립 노트북으로 모델을 학습합니다.

- **피처 구성**:
  1. Pathway score → `StandardScaler` 표준화 → **PCA(누적 설명력 90%)** 로 차원 축소.
  2. 해당 약물의 **DEG 상위 3개 유전자** 발현값(log2 변환 + 표준화)을 추가.
- **누수 방지**: 모든 스케일러/PCA는 **train에만 fit**, valid·test·TCGA에는 transform만.
- **모델**: Logistic Regression · Random Forest · SVM · KNN (전처리를 R로 수행해 파이프라인 묶음/CV가 어려워 **수동 튜닝**).
- **평가지표**: **AUC, F1 Score** 중심 → 최고 성능 모델 선택.
- **환자 전이 예측**: 선택 모델로 **TCGA 환자의 민감/저항 라벨 예측**.
  - `predict_proba ≥ 0.7` 인 **고신뢰 샘플만** 사용 → 임상 분석 시 노이즈 감소, 군 간 차이를 선명하게.

### ④ 임상적 연관성 검증 — 카이제곱 검정 (`Target_sample.ipynb`)

- 예측된 민감/저항 그룹과 **임상 변수**(PAM50 분자아형, ER/PR/HER2 상태, AJCC 병기)의 연관성을 **카이제곱 검정**으로 확인.

### ⑤ 생존분석 (`Survival analysis/<약물명>_survival.ipynb`)

- 예측 민감군 vs 저항군의 **전체생존(OS)** 차이를 **Kaplan-Meier 곡선 + Log-rank 검정**으로 비교.
- **층화 분석**: 전체 / **병기(AJCC Stage)별** / **분자아형(PAM50)별**.
- 정상조직·암조직을 모두 가진 환자는 제거 후 재검증.

---

## 💊 분석 대상 약물

| 트랙 | 약물 |
|------|------|
| **IC_50** (10종) | AZD7762, DASATINIB, GSK1059615, PEVONEDISTAT, REFAMETINIB, SORAFENIB, TANESPIMYCIN, TRAMETINIB, VINBLASTINE, VINORELBINE |
| **IC_Label** (10종) | AZ960, BUPARLISIB, CEDIRANIB, DASATINIB, FORETINIB, IRINOTECAN, MK-1775, PEVONEDISTAT, TOZASERTIB, TRAMETINIB |

> 각 약물 폴더 구성: ① CCLE pathway score(train/test/valid) ② 민감/저항 라벨(ic_train/test/val) ③ TCGA pathway score ④ 모델링 노트북 ⑤ (IC50 트랙) TCGA 예측 결과.

---

## 🛠️ 사용 기술 스택

| 분류 | 도구 |
|------|------|
| 데이터 처리 | Python `pandas`, `numpy` |
| 시각화 | `matplotlib`, `seaborn` |
| 머신러닝 | `scikit-learn` (PCA, LogisticRegression, RandomForest, SVM, KNN) |
| 생존분석 | `lifelines` (KaplanMeierFitter, CoxPHFitter, logrank_test) |
| 생물정보 (R) | `edgeR`, `DESeq2`, `dplyr`/`tidyr` (유전자 필터링·pathway score) |
| 도메인 DB | KEGG pathway, CCLE, GDSC, TCGA-BRCA |

### 설치 (Python)

```bash
pip install pandas numpy matplotlib seaborn scikit-learn lifelines
```

### 설치 (R)

```r
install.packages(c("dplyr", "tidyr", "stringr", "tibble", "caret"))
# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("edgeR", "DESeq2"))
```

---

## 🚀 실행 순서

```
[전처리]
EDA_Drug.ipynb → EdgeR_CCLE_gene.R → Gene_mapping.ipynb → TCGA_trans.ipynb → EDA_TCGA.ipynb

[라벨링·피처]
Pathway Score/pathway_score_ver2_ic50.R   (z-score 라벨 + DEG 상위 10약물 + pathway score)

[모델링·예측]   (약물별)
Model/IC_50 Modeling + 카이제곱/<약물>/<약물>.ipynb

[검증]
Target_sample.ipynb              (카이제곱)
Survival analysis/<약물>_survival.ipynb   (생존분석)
```

> 각 노트북/스크립트 상단의 파일 경로는 작성 당시 로컬 경로(`C:/...`, `/Users/...`)이므로 본인 환경에 맞게 수정 후 실행하세요.

---

## 📄 논문 공모전 보완 버전 (프로젝트 → 논문)

위 파이프라인은 **비어플 프로젝트 단계**의 결과물입니다. 이후 **통계청 논문공모전**에 투고하면서 방법론을 통계적으로 정교화하고, 검증 체계를 대폭 확장했습니다. 프로젝트와 논문의 차이는 다음과 같습니다.

> 📁 논문 원문: `통계청논문공모전[투고]/정밀 의료 및 신약 개발 정책 적용을 위한 약물 반응 예측 모델 연구 보완작.pdf`

### 1) 연구의 정체성 — "분석"에서 "정책 제언"으로

- 프로젝트가 *약물 반응 예측·검증*에 초점을 뒀다면, 논문은 이를 **신약 개발 R&D의 고비용·고위험 문제**와 연결.
- 제안 파이프라인을 **식품의약품안전처(MFDS) 전임상 심사 단계의 예비 스크리닝 보조도구**로 위치시키고, **규제과학(Regulatory Science)** 관점의 정책 활용 방안(통계청·식약처·제약업계 데이터 허브)을 제시.

### 2) 핵심 기여의 명확화 — DEG 기준 Pathway Score

- 기존 연구의 한계: pathway score 계산 시 **기준 유전자(reference gene)를 임의로/첫 번째 유전자로** 선택 → 약물 반응성과의 생물학적 연관성이 약함.
- **본 연구의 제안**: 각 pathway 내 **DEG를 기준 유전자로 선정**(`g_ref = G_P ∩ DEG`, DEG가 2개 이상이면 분산이 가장 큰 DEG 선택) → "이 경로가 약물 반응 결정 유전자의 발현을 얼마나 지지하는가"라는 생물학적 의미 부여.
- **기존 방식(첫 유전자 기준) vs 제안 방식(DEG 기준)을 정면 비교**하는 것이 논문의 핵심 실험 설계.

### 3) 전처리·피처 산출의 통계적 정교화

| 항목 | 프로젝트 | 논문 보완 |
|------|----------|-----------|
| **공통 유전자** | KEGG∩CCLE∩TCGA 교집합 | 각 ~20,000개 → **공통 약 16,000개** 명시 |
| **저발현 필터** | edgeR `filterByExpr` | **edgeR CPM 수식 명시**(전체 샘플 70%↑에서 CPM≥k) → 16,000 → **10,927개** → pathway 교집합 후 **약 2,700개** |
| **민감/저항 라벨** | pIC50 z-score (z>0 sensitive) | 동일하나 **수식·GDSC 기준 정합성** 명시 |
| **DEG 분석** | p<0.05 & \|logFC\|≥2 | **DESeq2 median-of-ratios 정규화 + Wald test + Benjamini-Hochberg FDR 5% 보정**, `p_adj<0.05 & \|log2FC\|≥1`, **volcano plot** 시각화 |
| **PCA 차원 축소** | `n_components=0.90` 고정 | **elbow 자동 탐지**(분산 곡선의 첫점–끝점 직선 거리 최대 지점)로 약물별 최적 주성분 수 결정 |

### 4) 배치 효과(Batch Effect) 보정 전략 추가

- 논문에서 새로 다룬 문제의식: *"세포주(CCLE)에서 도출한 결과를 환자(TCGA)에 직접 적용하는 것이 타당한가?"*
- TCGA에서는 일부 유전자의 발현 분산이 0에 가까워 **상관계수 부호(sign) 계산이 불안정**해지는 현상 관찰.
- **해법**: 기준 유전자는 **CCLE에서 고정**, 상관계수 부호 벡터는 **데이터셋별로 따로 계산** → 기준의 일관성과 데이터 특이성 반영의 균형.

### 5) 모델 구성 확장

- 프로젝트 4종(LR·RF·SVM·KNN) → 논문은 **AdaBoost를 추가한 5종**을 **소프트 보팅(예측 확률 평균) 앙상블**로 결합.
- 환자 전이 예측 시 고신뢰 필터를 비대칭으로 정교화: **민감군 proba ≥ 0.6, 저항군 proba ≥ 0.7**.

### 6) 검증 체계의 대폭 강화 (논문의 핵심 결과)

프로젝트는 약물별 개별 검증이었다면, 논문은 **총 360개 약물 전체**에 대해 두 방식을 통계적으로 비교했습니다.

- **AUC 분포 비교**: KDE + Box plot으로 360개 약물 AUC 분포를 비교 → DEG 기준 방식이 전반적으로 우상향.
- **Wilcoxon signed-rank test**: 모든 모델에서 **p < 0.001** → DEG 기준 방식의 우위가 통계적으로 유의.
- **DEG 개수 ↔ AUC 회귀분석**: DEG 수가 AUC에 **유의한 양(+)의 영향**(전 모델 p < 2e-16, R² ≈ 0.19~0.29) → "DEG가 많은(=두 군 분리가 뚜렷한) 약물일수록 예측이 잘 된다"는 가설을 실증.
- **대표 약물(DABRAFENIB)** ROC 곡선으로 기준 유전자 방식별 성능 직접 비교.

**상위 10개 약물 DEG 수 / AUC (논문 표2):**

| 약물 | DEG 수 | AUC | | 약물 | DEG 수 | AUC |
|------|:---:|:---:|---|------|:---:|:---:|
| REFAMETINIB | 217 | 0.73 | | VINBLASTINE | 150 | 0.74 |
| PEVONEDISTAT | 195 | 0.60 | | VINORELBINE | 150 | 0.70 |
| TRAMETINIB | 171 | 0.75 | | SORAFENIB | 149 | 0.76 |
| DASATINIB | 157 | 0.76 | | TANESPIMYCIN | 148 | 0.71 |
| GSK1059615 | 157 | 0.64 | | AZD7762 | 147 | 0.62 |

### 7) 임상 검증의 정교화 & 목표 환자군 도출

- **종양 vs 정상 조직 카이제곱 검정 + 방향성 분석**: 단순 연관성을 넘어 *"종양 조직에서만 선택적으로 민감한가"* 까지 판정 → **우수 약물 선정**.
  - 우수 약물: **PEVONEDISTAT · SORAFENIB · VINBLASTINE · VINORELBINE** (종양 민감군↑ & 정상 저항군↑ 방향성 충족).
- **임상 변수 5종**(암 병기·PAM50 분자아형·ER·PR·HER2)과 민감도의 독립성 카이제곱.
  - 예: **VINBLASTINE** → 분자아형/ER/PR이 유의(p<0.05) → **ER·PR 음성이며 기저형/루미널 B 환자에게 효과적**일 가능성. 이는 *기존 호르몬 치료에 반응하지 않는 환자군의 대안 치료제* 라는 임상적 함의로 연결.
- **t-SNE 군집 시각화**: VINBLASTINE·VINORELBINE에서 민감군/저항군이 뚜렷이 분리 → 세포주 학습 모델이 환자 특성 차이를 반영함을 시각적으로 입증.

### 8) 명시된 한계 (논문 결론)

- 실제 약물 투약 이력·임상 반응 정보가 포함되지 않았고, 생물학적 기전 기반 변수 설계가 제한적.
- 향후 **코호트 기반 외부 검증, 약물 표적 정보 통합, 실제 임상 반응 데이터 연계**로 고도화 필요.

---

## 💡 프로젝트 핵심 인사이트

1. **세포주 → 환자 전이가 핵심 도전** — 도메인이 다른 두 데이터를 **pathway score**라는 공통 피처로 연결.
2. **DEG가 좋은 약물이 예측도 잘 된다** — 민감/저항 군의 생물학적 분리가 뚜렷할수록 모델 신뢰도↑ → DEG 개수로 약물 선별.
3. **고신뢰 예측만 임상 검증에 사용** — `predict_proba ≥ 0.7` 로 노이즈를 줄여 군 간 차이를 선명하게.
4. **통계·생존분석으로 최종 검증** — 단순 예측 정확도를 넘어, 예측 반응군이 실제 임상·생존과 연관되는지까지 확인.
