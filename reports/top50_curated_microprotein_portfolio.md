# Curated Top 50 Human Microproteins Portfolio

**Selection Framework**: Balanced Multi-Cohort Strategy de-confounded for canonical CDS hitchhiking  
**Source Catalog**: 7,264 non-canonical smORFs from Nature 2026 (*"Expanding the human proteome with microproteins and peptideins"*, DOI: [10.1038/s41586-026-10459-x](https://doi.org/10.1038/s41586-026-10459-x))  
**Selection Metrics**: Biochemical mass-spectrometry validation, genomic autonomy (non-CDS), AlphaGenome Variant Impact (AVI) purifying selection, and multimodal functional disruption.

---

## Executive Summary of Portfolio Design

If a laboratory or drug-discovery program had to select exactly **50 microproteins** to clone, express, structurally validate (via ESMFold2/PyMOL), and interrogate phenotypically, sorting by raw evolutionary constraint would mistakenly select 40+ internal coding ORFs whose signal is merely borrowed from essential canonical proteins (*TBC1D7*, *KRR1*, *MBTD1*).

To deliver the **highest-confidence, biologically authentic candidates**, this portfolio divides the 50 spots into **four distinct strategic cohorts**:

```
                                    TOP 50 CANDIDATE PORTFOLIO
                                                │
         ┌───────────────────────────┬──────────┴────────────────┬───────────────────────────┐
         ▼                           ▼                           ▼                           ▼
┌──────────────────┐       ┌──────────────────┐       ┌──────────────────┐       ┌──────────────────┐
│ Cohort 1: (n=15) │       │ Cohort 2: (n=15) │       │ Cohort 3: (n=12) │       │ Cohort 4: (n=8)  │
│ Autonomous       │       │ Autonomous       │       │ Autonomous       │       │ Dual-Coding      │
│ Tier 1           │       │ lncRNA-ORFs      │       │ Regulatory       │       │ Tier 1           │
│ Peptideins       │       │ ("Novel Genes")  │       │ 5' UTR uORFs     │       │ Peptideins       │
└──────────────────┘       └──────────────────┘       └──────────────────┘       └──────────────────┘
```

1. **Cohort 1: Autonomous Tier 1 Mass-Spec Peptideins ($n = 15$)**: The ultimate gold standard. Directly validated by mass spectrometry (HLA/non-HLA peptides) and residing outside canonical protein-coding regions.
2. **Cohort 2: Autonomous `lncRNA-ORFs` ("Novel Gene Frontier", $n = 15$)**: Pure non-coding RNA transcripts with zero canonical CDS. Sequence constraint is 100% autonomous, characterized by high chromatin remodeling and transcriptional regulatory signatures.
3. **Cohort 3: Autonomous Regulatory 5' UTR `uORFs` ($n = 12$)**: Upstream untranslated region microproteins with extreme evolutionary constraint regulating master oncogenes, chromatin remodelers, and transcription factors.
4. **Cohort 4: Mass-Spec Validated Dual-Coding Peptideins ($n = 8$)**: Elite peptides nested inside canonical genes that possess direct mass-spec proof of translation, demonstrating genuine dual-coding function.

---

## Complete Ranked Catalog of the Top 50 Candidates

### Cohort 1: Autonomous Tier 1 Mass-Spec Peptideins ($n = 15$)
*Biochemically validated by mass spectrometry (Tier 1A/1B) + autonomous non-CDS genomic location + strong evolutionary selection.*

| Rank | ORF ID | Host Gene | Biotype | Length | Tier | Peptide Support | Median Phred | Peak Phred | Peak Hotspot Variant | Primary Mechanism | AlphaGenome Atlas Link |
| :---: | :--- | :--- | :---: | :---: | :---: | :--- | :---: | :---: | :--- | :--- | :--- |
| **1** | `c3riboseqorf106` | *ZBTB11-AS1* | lncRNA-ORF | 88 aa | 1B | HLA Peptides | **23.56** | **45.77** | `chr3:101676728:C>A` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:101676645-101676911&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **2** | `c2riboseqorf55` | *WBP1* | uORF | 36 aa | 1B | HLA Peptides | **19.61** | **33.26** | `chr2:74458516:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:74458484-74458594&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **3** | `c1riboseqorf248` | *IVNS1ABP* | uORF | 24 aa | 1B | HLA Peptides | **18.23** | **33.22** | `chr1:185316953:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:185311281-185316985&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **4** | `c7riboseqorf18` | *SNX13* | uORF | 52 aa | 1B | HLA Peptides | **17.32** | **31.59** | `chr7:17940333:G>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr7:17940325-17940483&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **5** | `c6riboseqorf128` | *BCLAF1* | uORF | 19 aa | 1B | HLA Peptides | **17.21** | **22.47** | `chr6:136289799:C>A` | Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:136289742-136289801&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **6** | `cXriboseqorf59` | *MORF4L2* | uORF | 32 aa | 1B | HLA Peptides | **16.98** | **34.03** | `chrX:103686631:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chrX:103685228-103686696&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **7** | `c8riboseqorf15` | *CNOT7* | uORF | 39 aa | 1A | Non-HLA MS | **16.78** | **30.38** | `chr8:17246687:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr8:17245236-17246782&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **8** | `c1norep182` | *IGSF3* | uORF | 35 aa | 1B | HLA Peptides | **16.65** | **22.53** | `chr1:116666677:C>A` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:116666632-116666739&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **9** | `c6norep97` | *PSMB8-AS1* | lncRNA-ORF | 63 aa | 1B | HLA Peptides | **16.04** | **45.89** | `chr6:32845662:C>A` | Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:32845552-32845743&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **10** | `c17riboseqorf78` | *IGFBP4* | uORF | 17 aa | 1B | HLA Peptides | **15.44** | **24.15** | `chr17:40443555:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr17:40443505-40443558&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **11** | `c4riboseqorf2` | *CTBP1* | uORF | 32 aa | 1B | HLA Peptides | **15.18** | **31.95** | `chr4:1241388:T>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr4:1241335-1241433&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **12** | `c21riboseqorf14` | *RUNX1* | uORF | 20 aa | 1B | HLA Peptides | **14.93** | **23.84** | `chr21:34887820:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr21:34887800-34887862&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **13** | `c19riboseqorf66` | *URI1* | uORF | 53 aa | 1B | HLA Peptides | **14.79** | **24.85** | `chr19:29942370:A>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr19:29942370-29942531&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **14** | `c3riboseqorf154` | *MBNL1* | uORF | 16 aa | 1B | HLA Peptides | **14.72** | **22.44** | `chr3:152269066:G>C` | Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:152269026-152269076&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **15** | `c16riboseqorf104` | *PSKH1* | uORF | 33 aa | 1B | HLA Peptides | **13.62** | **30.51** | `chr16:67893371:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:67893276-67908685&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |

---

### Cohort 2: Autonomous `lncRNA-ORFs` ("Novel Gene Frontier", $n = 15$)
*Pure non-coding transcripts (0 canonical CDS) + top evolutionary constraint in lncRNA catalog ($\text{Phred} \ge 15.6–22.8$) + major chromatin/transcriptional disruption.*

| Rank | ORF ID | Host Gene | Length | Tier | Peptidein? | Median Phred | Peak Phred | Peak Hotspot Variant | Primary Mechanism | AlphaGenome Atlas Link |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **16** | `c12riboseqorf13` | *ENSG00000272173* | 57 aa | 2B | No | **15.67** | **59.69** | `chr12:6944514:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr12:6944390-6944563&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **17** | `c1riboseqorf199` | *ENSG00000229953* | 43 aa | 1B | No | **19.36** | **47.98** | `chr1:156647121:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:156647059-156661350&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **18** | `c1norep185` | *WARS2-AS1* | 26 aa | 4 | No | **16.75** | **46.48** | `chr1:119140623:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:119140598-119140678&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **19** | `c16norep119` | *ZFHX3-AS1* | 154 aa | 4 | No | **20.20** | **46.24** | `chr16:72788546:C>A` | **Chromatin Opening** | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:72788372-72788836&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **20** | `c1riboseqorf33` | *EMC1-AS1* | 61 aa | 3 | No | **18.00** | **45.69** | `chr1:19240340:A>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:19240260-19240445&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **21** | `c11norep50` | *NAV2-AS2* | 103 aa | 4 | No | **19.25** | **45.55** | `chr11:20044000:C>G` | **Chromatin Opening** | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:20043884-20044195&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **22** | `c20norep33` | *ENSG00000275457* | 29 aa | 4 | No | **17.24** | **45.22** | `chr20:21303402:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr20:21303362-21303451&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **23** | `c17norep109` | *MAP3K14-AS1* | 50 aa | 4 | No | **19.13** | **44.92** | `chr17:45267160:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr17:45255181-45267227&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **24** | `c2norep270` | *ENSG00000224090* | 83 aa | 4 | No | **15.92** | **44.44** | `chr2:219014007:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:219010000-219014110&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **25** | `c2riboseqorf136` | *TTN-AS1* | 127 aa | 4 | No | **22.05** | **44.23** | `chr2:178539103:C>T` | **Chromatin Opening** | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:178538902-178539285&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **26** | `c11riboseqorf31` | *NAV2-AS2* | 31 aa | 4 | No | **22.08** | **43.27** | `chr11:20049059:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:20049028-20049123&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **27** | `c19norep130` | *TMEM147-AS1* | 26 aa | 4 | No | **21.30** | **41.83** | `chr19:35542971:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr19:35542912-35542992&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **28** | `c16norep37` | *ENSG00000260592* | 61 aa | 4 | No | **17.68** | **41.78** | `chr16:19487269:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:19487189-19487877&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **29** | `c6norep46` | *NUP153-AS1* | 96 aa | 2B | No | **17.00** | **41.61** | `chr6:17706369:C>A` | Subtle / Regulatory Shift | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:17706280-17706570&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **30** | `c2norep211` | *TTN-AS1* | 87 aa | 4 | No | **22.82** | **41.45** | `chr2:178531421:G>C` | **Chromatin Closing** | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:178531419-178537610&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |

---

### Cohort 3: Autonomous Regulatory 5' UTR `uORFs` ($n = 12$)
*Strictly located in 5' UTR (non-CDS overlapping) + extreme purifying selection ($\text{Phred} \ge 18.1–20.8$) regulating essential human genes.*

| Rank | ORF ID | Host Gene | Length | Tier | Peptidein? | Median Phred | Peak Phred | Peak Hotspot Variant | Primary Mechanism | AlphaGenome Atlas Link |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **31** | `c3riboseqorf183` | *BCL6* | 19 aa | 4 | No | **20.79** | **46.16** | `chr3:187745410:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:187734882-187745443&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **32** | `c8riboseqorf89` | *DCAF13* | 61 aa | 2B | No | **19.30** | **41.83** | `chr8:103414903:G>T` | **TF Binding Disruption** | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr8:103414754-103414939&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **33** | `c10riboseqorf118` | *SHTN1* | 21 aa | 4 | No | **19.92** | **37.97** | `chr10:117005186:C>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr10:117005181-117005246&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **34** | `c17norep67` | *EVI2A* | 35 aa | 2B | No | **18.11** | **36.43** | `chr17:31321672:A>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr17:31321566-31321673&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **35** | `c1riboseqorf209` | *DCAF8* | 42 aa | 1B | No | **18.30** | **36.27** | `chr1:160262354:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:160262315-160262443&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **36** | `cXriboseqorf74` | *DNASE1L1* | 26 aa | 4 | No | **18.78** | **35.88** | `chrX:154411928:G>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chrX:154411914-154411994&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **37** | `c14norep27` | *CIDEB* | 28 aa | 2B | No | **18.10** | **35.53** | `chr14:24311302:G>A` | Subtle / Regulatory Shift | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr14:24311264-24311350&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **38** | `c22riboseqorf6` | *DGCR8* | 16 aa | 4 | No | **20.05** | **35.06** | `chr22:20085750:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr22:20085714-20085764&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **39** | `c20riboseqorf31` | *ITCH* | 16 aa | 2B | No | **20.76** | **34.29** | `chr20:34369443:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr20:34369399-34369449&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **40** | `c3norep215` | *TBL1XR1* | 18 aa | 2B | No | **20.76** | **34.19** | `chr3:177098466:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:177064980-177098479&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **41** | `c14riboseqorf81` | *NUMB* | 32 aa | 1B | No | **18.13** | **34.18** | `chr14:73366897:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr14:73355753-73366981&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **42** | `c13norep17` | *HMGB1* | 45 aa | 2B | No | **20.12** | **34.12** | `chr13:30465842:G>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr13:30465798-30465935&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |

---

### Cohort 4: Mass-Spec Validated Dual-Coding Peptideins ($n = 8$)
*Nested in canonical CDS + direct mass spectrometry proof of translation (Tier 1B) + extreme purifying selection ($\text{Phred} \ge 22.2–28.0$).*

| Rank | ORF ID | Host Gene | Biotype | Length | Tier | Peptide Support | Median Phred | Peak Phred | Peak Hotspot Variant | Primary Mechanism | AlphaGenome Atlas Link |
| :---: | :--- | :--- | :---: | :---: | :---: | :--- | :---: | :---: | :--- | :--- | :--- |
| **43** | `c1norep256` | *PRRC2C* | intORF | 34 aa | 1B | HLA Peptides | **27.99** | **42.02** | `chr1:171513009:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:171513007-171513111&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **44** | `c1riboseqorf39` | *HNRNPR* | intORF | 45 aa | 1B | HLA Peptides | **27.50** | **47.17** | `chr1:23337831:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:23337831-23338596&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **45** | `c11norep30` | *EIF4G2* | intORF | 74 aa | 1B | HLA Peptides | **25.74** | **45.21** | `chr11:10803915:T>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:10803521-10804139&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **46** | `c11riboseqorf112` | *EMSY* | uoORF | 91 aa | 1B | HLA Peptides | **25.08** | **59.69** | `chr11:76447008:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:76445106-76453357&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **47** | `c20norep82` | *ADNP* | intORF | 65 aa | 1B | HLA Peptides | **23.65** | **42.26** | `chr20:50894419:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr20:50894233-50894430&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **48** | `c3norep241` | *P3H2* | intORF | 63 aa | 1B | HLA Peptides | **22.31** | **48.43** | `chr3:190120252:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:189995313-190120313&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **49** | `c16riboseqorf58` | *SRCAP* | uoORF | 62 aa | 1B | HLA Peptides | **22.26** | **52.90** | `chr16:30700876:C>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:30700709-30704082&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **50** | `c10norep31` | *CDC123* | intORF | 29 aa | 1B | HLA Peptides | **22.17** | **42.06** | `chr10:12198733:A>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr10:12196259-12198733&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |

---

## Data Artifacts

The curated portfolio has been exported in open scientific formats:
* **TSV Format**: [`reports/top50_curated_candidates.tsv`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/top50_curated_candidates.tsv)
* **Parquet Format**: [`reports/top50_curated_candidates.parquet`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/top50_curated_candidates.parquet)
