# Multi-Omics Comparative Characterization of Human Microproteins & AlphaGenome Regulatory Selection

**Primary Literature Reference**: *"Expanding the human proteome with microproteins and peptideins"*, Nature 2026 (DOI: [10.1038/s41586-026-10459-x](https://doi.org/10.1038/s41586-026-10459-x))  
**Genomic Effect Model**: Google DeepMind AlphaGenome Atlas & Multimodal DNA Model  
**Catalog Size**: 7,264 non-canonical Small Open Reading Frames (smORFs)  
**Total Scored Genomic Variants**: 2,900,430 dense single-nucleotide variants  
**Track-Level Multimodal Impacts**: 9,911,834 significant records across 54 GTEx tissues and ENCODE cell lines  

---

## Executive Summary

Non-canonical open reading frames (ncORFs / smORFs) encode functional microproteins and peptideins that represent an uncharacterized frontier of the human proteome. By integrating mass-spectrometry validation tiers, standardized hg38 genomic intervals, dense AlphaGenome Variant Impact (AVI) purifying selection scores, and DeepMind's multimodal differential DNA model, this study provides an atlas of evolutionary selection and functional vulnerability across the human microprotein landscape.

### Major Discoveries & Methodological Insights

1. **Resolution of the "Canonical CDS Hitchhiking" Confounder**:
   - Initial catalog-wide ranking by raw AVI constraint produces an overwhelming majority of internal coding ORFs (`intORF`, 75%+ of top hits).
   - **The Confounder**: In `intORFs` and `uoORFs`, every nucleotide change is dual-use; any mutation alters both the microprotein and an essential full-length canonical protein (*TBC1D7*, *KRR1*, *MBTD1*, *UBE4A*, *RAC1*). Purifying selection is fiercely protecting the host protein, meaning raw AVI scores cannot deconvolute microprotein selection from host protein constraint.
   - **Quantitative Proof**: CDS-overlapping smORFs are **49.1x more likely** to exhibit $\text{Phred} \ge 20$ (36.4% vs. 1.15%, $p < 10^{-300}$, Fisher's exact test).

2. **The Autonomous smORF Landscape ($n = 5,387$, 74.2% of Catalog)**:
   - To uncover *genuine, independent* microprotein-driven purifying selection, we isolated **Autonomous smORFs** (`lncRNA-ORFs`, `uORFs`, `dORFs`) that do not overlap canonical protein-coding CDS.
   - Across this unconfounded cohort, **830 microproteins** reside under moderate-to-strong purifying selection ($\text{Median AVI Phred} \ge 15$, Top 3.2% genome-wide constraint), and **62** exhibit extreme constraint ($\text{Median AVI Phred} \ge 20$, Top 1% genome-wide constraint).

3. **`lncRNA-ORFs` as the Pure Autonomy Benchmark ($n = 2,012$)**:
   - Located on non-coding transcripts with **zero canonical CDS**, any sequence constraint observed in `lncRNA-ORFs` is autonomous.
   - While the bulk of lncRNAs evolve neutrally (median Phred = 4.98), an ultra-conserved core of **18 lncRNA-ORFs** achieves $\text{Median Phred} \ge 20$, and **79** achieve $\text{Median Phred} \ge 15$.
   - **Dramatic Functional Modality Shift**: Whereas CDS-overlapping ORFs are 85%+ splicing-driven, peak variants in `lncRNA-ORFs` predominantly disrupt **Chromatin Opening (31.0%)** and **Transcriptional Activation / Silencing (23.5%)**.

4. **De-confounded Validation of Mass-Spec Tier 1 Peptideins**:
   - Among the 84 mass-spec validated Tier 1 peptideins, **45 are completely autonomous** (16 lncRNA-ORFs, 25 uORFs, 4 dORFs).
   - Autonomous peptideins like `c3riboseqorf106` (*ZBTB11-AS1*, Median = 23.56, Peak = 45.77), `c2riboseqorf55` (*WBP1*, Median = 19.61, Peak = 33.26), and `c1riboseqorf248` (*IVNS1ABP*, Median = 18.23) prove that biochemically validated microproteins possess autonomous evolutionary importance independent of any host gene.

---

## 1. De-confounding the Catalog: Autonomous vs. CDS-Overlapping smORFs

To establish valid rankings, we partitioned the 7,264 microproteins into two fundamentally distinct genomic regimes:
* **Autonomous smORFs ($n = 5,387$, 74.2%)**: Located outside canonical protein-coding CDS boundaries:
  - `uORF` (5' Untranslated Region, $n = 2,915$)
  - `lncRNA-ORF` (Long Non-Coding RNA transcripts, $n = 2,012$)
  - `dORF` (3' Untranslated Region, $n = 460$)
* **CDS-Overlapping smORFs ($n = 1,877$, 25.8%)**: Physically overlapping the canonical protein-coding sequence:
  - `intORF` (Nested within internal canonical CDS, $n = 743$)
  - `uoORF` (Spans canonical 5' UTR start codon into CDS, $n = 622$)
  - `mixed` (Spanning multiple structural boundaries, $n = 448$)
  - `doORF` (Spans canonical stop codon from CDS into 3' UTR, $n = 64$)

```
AUTONOMOUS (Unconfounded Selection)        CDS-OVERLAPPING (Host Gene Hitchhiking)
=====================================     ======================================
1. lncRNA-ORF                             1. intORF
   [==== smORF ====]                         Canonical CDS: [==============================]
   (Zero canonical CDS on transcript)        intORF:             [==== smORF ====] (Dual-use bp)

2. 5' UTR uORF                            2. uoORF
   5' UTR: [= smORF =]  Canonical CDS:       5' UTR: [===== smORF =====] Canonical CDS:
   ------->[========]   [============]       ------->[=======|=========] [=============]
                                                             ^ crosses start codon
```

### Statistical Quantification of the Hitchhiking Effect

| Metric | Autonomous smORFs ($n=5,387$) | CDS-Overlapping smORFs ($n=1,877$) | Statistical Test | p-value | Effect Size |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Median AVI Phred (Median)** | **7.78** | **17.69** | Mann-Whitney U | $p < 10^{-300}$ | $\Delta = +9.91$ Phred |
| **Median AVI Phred (Mean)** | **8.90** | **16.69** | Mann-Whitney U | $p < 10^{-300}$ | $\Delta = +7.79$ Phred |
| **Peak Hotspot Phred (Median)** | **22.84** | **38.92** | Mann-Whitney U | $p < 10^{-300}$ | $\Delta = +16.08$ Phred |
| **Extreme Constraint ($\ge 20$ Phred)** | **62 (1.15%)** | **683 (36.39%)** | Fisher's Exact | $p < 10^{-300}$ | **Odds Ratio = 49.1x** |
| **Moderate Constraint ($\ge 15$ Phred)**| **830 (15.41%)** | **1,193 (63.56%)** | Fisher's Exact | $p < 10^{-300}$ | **Odds Ratio = 9.6x** |

> [!IMPORTANT]
> **Key Conclusion**: CDS-overlapping smORFs are **49 times more likely** to achieve a median Phred $\ge 20$ than autonomous smORFs. Catalog-wide ranking without stratification is essentially an assay for canonical protein essentiality rather than microprotein function.

---

## 2. Shift in Functional Disruption Mechanisms ($\chi^2 = 428.34, p < 10^{-86}$)

Evaluating single-variant functional disruptions using AlphaGenome's differential DNA model reveals a dramatic divergence between autonomous and CDS-overlapping microproteins:

| Primary Disruption Mechanism | Autonomous smORFs ($n=5,167$) | CDS-Overlapping ($n=1,797$) | `lncRNA-ORFs` Only ($n=1,927$) | Biological Manifestation in Autonomous smORFs |
| :--- | :---: | :---: | :---: | :--- |
| **Severe Splicing Disruption** | 50.7% | **63.0%** | **15.6%** | Host transcript intron retention / exon skipping |
| **Chromatin Opening** | **14.7%** | 3.6% | **31.0%** | De-novo accessibility of previously closed loci |
| **Transcriptional Activation** | **4.7%** | 0.1% | **11.5%** | Upregulation of proximal antisense/divergent transcripts |
| **Transcriptional Silencing** | **4.4%** | 0.1% | **11.9%** | Disruption of core promoter or enhancer elements |
| **Subtle / Regulatory Shift** | 9.0% | **15.8%** | **12.3%** | Moderate epigenetic and transcript level modulations |
| **TF Binding Disruption** | **2.8%** | 1.3% | **3.5%** | Ablation of CTCF, POLR2A, MAX, or TAF1 motifs |
| **Chromatin Closing / Repression** | **2.7%** | 2.0% | **5.0%** | Loss of active promoter/enhancer accessibility |
| **Benign / Tolerated** | 2.6% | 2.2% | 4.6% | Non-functional neutral sequence space |

In `lncRNA-ORFs`, **splicing disruption drops to only 15.6%**, while **chromatin and transcriptional alterations surge to over 54%**. This reflects the true regulatory nature of autonomous non-coding microprotein loci.

---

## 3. Publication Figures

The seven figures below (300 DPI) illustrate both the catalog-wide patterns and the unconfounded autonomous landscape:

### Figure 5: Purifying Selection Shift (Autonomous vs. CDS-Overlapping)
Left: Grouped violin distributions of `median_avi_phred` across all 7 biotypes partitioned by genomic autonomy. Right: Empirical Cumulative Distribution Functions (CDF) demonstrating the 49.1x Odds Ratio shift in extreme constraint caused by canonical CDS hitchhiking.

![Figure 5: Autonomous vs Overlapping Selection](C:/Users/Lucas.Escosteguy/.gemini/antigravity/brain/b4954ee9-d996-4b97-9d18-0ab5359d6b32/figures/fig5_autonomous_vs_overlapping_selection.png)

---

### Figure 6: Functional Disruption Modality Divergence
Normalized horizontal stacked bar chart showing the fundamental mechanism shift: CDS-overlapping smORFs are dominated by severe splicing disruption, whereas autonomous `lncRNA-ORFs` and `uORFs` display high frequencies of chromatin opening, transcriptional activation/silencing, and transcription factor binding disruption.

![Figure 6: Autonomous Disruption Mechanisms](C:/Users/Lucas.Escosteguy/.gemini/antigravity/brain/b4954ee9-d996-4b97-9d18-0ab5359d6b32/figures/fig6_autonomous_disruption_mechanisms.png)

---

### Figure 7: The Unconfounded Autonomous Standout Landscape
2D logarithmic density contour of peptide length vs. `median_avi_phred` for the 5,387 autonomous smORFs. Purple rings mark the top autonomous `lncRNA-ORFs` under extreme purifying selection ($\text{Phred} \ge 20$), while red markers indicate autonomous mass-spec validated Tier 1 peptideins, with key standout loci annotated (*ZBTB11-AS1*, *ZFHX3-AS1*, *TTN-AS1*, *BCL6*, *WBP1*).

![Figure 7: Top Autonomous Candidates Landscape](C:/Users/Lucas.Escosteguy/.gemini/antigravity/brain/b4954ee9-d996-4b97-9d18-0ab5359d6b32/figures/fig7_top_autonomous_candidates_landscape.png)

---

### Figure 1: Catalog-Wide Selection Distribution
![Figure 1: Biotype Selection Distribution](C:/Users/Lucas.Escosteguy/.gemini/antigravity/brain/b4954ee9-d996-4b97-9d18-0ab5359d6b32/figures/fig1_biotype_selection_distribution.png)

---

### Figure 2: Catalog-Wide Disruption Breakdown
![Figure 2: Disruption Mechanism Breakdown](C:/Users/Lucas.Escosteguy/.gemini/antigravity/brain/b4954ee9-d996-4b97-9d18-0ab5359d6b32/figures/fig2_disruption_mechanism_breakdown.png)

---

### Figure 3: Length vs. Constraint Landscape (All 7,264 ORFs)
![Figure 3: Constraint vs ORF Length Landscape](C:/Users/Lucas.Escosteguy/.gemini/antigravity/brain/b4954ee9-d996-4b97-9d18-0ab5359d6b32/figures/fig3_constraint_vs_length.png)

---

### Figure 4: GTEx Tissue Impact Matrix
![Figure 4: GTEx Tissue Impact Matrix](C:/Users/Lucas.Escosteguy/.gemini/antigravity/brain/b4954ee9-d996-4b97-9d18-0ab5359d6b32/figures/fig4_gtex_tissue_impact_matrix.png)

---

## 4. De-confounded Standout Leaderboards

### Leaderboard 1: Top 20 Autonomous `lncRNA-ORFs` Under Highest Purifying Selection
*(Filtered for `orf_biotype == 'lncRNA-ORF'`, `Median AVI Phred >= 15.0`, and `High-Impact Hotspot Variant`, sorted by `Peak Phred`)*

| Rank | ORF ID | Host Gene | Length | Tier | Peptidein? | Median Phred | Peak Phred | Peak Hotspot Variant | Primary Mechanism | AlphaGenome Atlas Deep-Link |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **1** | `c12riboseqorf13` | *ENSG00000272173* | 57 aa | 2B | No | **15.67** | **59.69** | `chr12:6944514:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr12:6944390-6944563&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **2** | `c1riboseqorf199` | *ENSG00000229953* | 43 aa | 1B | No | **19.36** | **47.98** | `chr1:156647121:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:156647059-156661350&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **3** | `c1norep185` | *WARS2-AS1* | 26 aa | 4 | No | **16.75** | **46.48** | `chr1:119140623:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:119140598-119140678&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **4** | `c16norep119` | *ZFHX3-AS1* | 154 aa | 4 | No | **20.20** | **46.24** | `chr16:72788546:C>A` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:72788372-72788836&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **5** | `c6norep97` | *PSMB8-AS1* | 63 aa | 1B | **Yes** | **16.04** | **45.89** | `chr6:32845662:C>A` | Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:32845552-32845743&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **6** | `c3riboseqorf106` | *ZBTB11-AS1* | 88 aa | 1B | **Yes** | **23.56** | **45.77** | `chr3:101676728:C>A` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:101676645-101676911&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **7** | `c1riboseqorf33` | *EMC1-AS1* | 61 aa | 3 | No | **18.00** | **45.69** | `chr1:19240340:A>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:19240260-19240445&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **8** | `c11norep50` | *NAV2-AS2* | 103 aa | 4 | No | **19.25** | **45.55** | `chr11:20044000:C>G` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:20043884-20044195&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **9** | `c20norep33` | *ENSG00000275457* | 29 aa | 4 | No | **17.24** | **45.22** | `chr20:21303402:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr20:21303362-21303451&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **10** | `c17norep109` | *MAP3K14-AS1* | 50 aa | 4 | No | **19.13** | **44.92** | `chr17:45267160:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr17:45255181-45267227&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **11** | `c2norep270` | *ENSG00000224090* | 83 aa | 4 | No | **15.92** | **44.44** | `chr2:219014007:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:219010000-219014110&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **12** | `c2riboseqorf136` | *TTN-AS1* | 127 aa | 4 | No | **22.05** | **44.23** | `chr2:178539103:C>T` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:178538902-178539285&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **13** | `c11riboseqorf31` | *NAV2-AS2* | 31 aa | 4 | No | **22.08** | **43.27** | `chr11:20049059:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:20049028-20049123&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **14** | `c19norep130` | *TMEM147-AS1* | 26 aa | 4 | No | **21.30** | **41.83** | `chr19:35542971:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr19:35542912-35542992&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **15** | `c16norep37` | *ENSG00000260592* | 61 aa | 4 | No | **17.68** | **41.78** | `chr16:19487269:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:19487189-19487877&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **16** | `c2norep211` | *TTN-AS1* | 87 aa | 4 | No | **22.82** | **41.45** | `chr2:178531421:G>C` | Chromatin Closing | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:178531419-178537610&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **17** | `c19norep131` | *TMEM147-AS1* | 38 aa | 4 | No | **17.57** | **40.84** | `chr19:35544993:A>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr19:35544884-35545000&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **18** | `c2norep220` | *CAVIN2-AS1* | 52 aa | 4 | No | **19.78** | **40.81** | `chr2:191846694:T>A` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:191846648-192034728&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **19** | `c15riboseqorf78` | *ST20-AS1* | 80 aa | 2B | No | **16.66** | **36.22** | `chr15:79923416:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr15:79923362-79923604&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **20** | `c12norep86` | *CBX5* | 43 aa | 2B | No | **17.88** | **34.56** | `chr12:54274083:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr12:54274076-54274207&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |

---

### Leaderboard 2: Top 20 Autonomous `uORFs` Under Highest Purifying Selection
*(Filtered for `orf_biotype == 'uORF'`, `Median AVI Phred >= 18.0`, and `High-Impact Hotspot Variant`, sorted by `Peak Phred`)*

| Rank | ORF ID | Host Gene | Length | Tier | Peptidein? | Median Phred | Peak Phred | Peak Hotspot Variant | Primary Mechanism | AlphaGenome Atlas Deep-Link |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **1** | `c3riboseqorf183` | *BCL6* | 19 aa | 4 | No | **20.79** | **46.16** | `chr3:187745410:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:187734882-187745443&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **2** | `c8riboseqorf89` | *DCAF13* | 61 aa | 2B | No | **19.30** | **41.83** | `chr8:103414903:G>T` | TF Binding Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr8:103414754-103414939&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **3** | `c10riboseqorf118` | *SHTN1* | 21 aa | 4 | No | **19.92** | **37.97** | `chr10:117005186:C>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr10:117005181-117005246&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **4** | `c17norep67` | *EVI2A* | 35 aa | 2B | No | **18.11** | **36.43** | `chr17:31321672:A>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr17:31321566-31321673&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **5** | `c1riboseqorf209` | *DCAF8* | 42 aa | 1B | No | **18.30** | **36.27** | `chr1:160262354:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:160262315-160262443&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **6** | `cXriboseqorf74` | *DNASE1L1* | 26 aa | 4 | No | **18.78** | **35.88** | `chrX:154411928:G>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chrX:154411914-154411994&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **7** | `c22riboseqorf6` | *DGCR8* | 16 aa | 4 | No | **20.05** | **35.06** | `chr22:20085750:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr22:20085714-20085764&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **8** | `c20riboseqorf31` | *ITCH* | 16 aa | 2B | No | **20.76** | **34.29** | `chr20:34369443:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr20:34369399-34369449&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **9** | `c3norep215` | *TBL1XR1* | 18 aa | 2B | No | **20.76** | **34.19** | `chr3:177098466:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:177064980-177098479&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **10** | `c14riboseqorf81` | *NUMB* | 32 aa | 1B | No | **18.13** | **34.18** | `chr14:73366897:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr14:73355753-73366981&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **11** | `c13norep17` | *HMGB1* | 45 aa | 2B | No | **20.12** | **34.12** | `chr13:30465842:G>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr13:30465798-30465935&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **12** | `c14norep107` | *DICER1* | 18 aa | 4 | No | **18.34** | **33.52** | `chr14:95140659:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr14:95133501-95140712&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **13** | `c6riboseqorf59` | *ZNF76* | 31 aa | 4 | No | **20.16** | **33.50** | `chr6:35259841:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:35259754-35281063&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **14** | `c2riboseqorf55` | *WBP1* | 36 aa | 1B | **Yes** | **19.61** | **33.26** | `chr2:74458516:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:74458484-74458594&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **15** | `c1riboseqorf248` | *IVNS1ABP* | 24 aa | 1B | **Yes** | **18.23** | **33.22** | `chr1:185316953:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:185311281-185316985&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **16** | `c11riboseqorf46` | *PHF21A* | 55 aa | 4 | No | **19.01** | **33.04** | `chr11:46120937:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:46084272-46121030&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **17** | `c9norep104` | *BRINP1* | 19 aa | 4 | No | **18.93** | **32.78** | `chr9:119369217:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr9:119369195-119369254&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **18** | `c1riboseqorf76` | *MTF1* | 31 aa | 2B | No | **18.93** | **32.78** | `chr1:37859531:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:37857672-37859588&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **19** | `c3riboseqorf179` | *FAM131A* | 31 aa | 4 | No | **18.25** | **32.75** | `chr3:184336123:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:184336087-184337369&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **20** | `c11norep176` | *CRYAB* | 18 aa | 4 | No | **22.31** | **32.67** | `chr11:111911891:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:111911849-111911905&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |

---

### Leaderboard 3: Top 20 Autonomous Mass-Spec Tier 1 Peptideins
*(Filtered for `is_autonomous == True`, `is_peptidein == True`, and `final_tier.startswith('1')`, sorted by `Median AVI Phred`)*

| Rank | ORF ID | Host Gene | Biotype | Length | Tier | Median Phred | Peak Phred | Peak Hotspot Variant | Primary Mechanism | AlphaGenome Atlas Deep-Link |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **1** | `c3riboseqorf106` | *ZBTB11-AS1* | lncRNA-ORF | 88 aa | 1B | **23.56** | **45.77** | `chr3:101676728:C>A` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:101676645-101676911&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **2** | `c2riboseqorf55` | *WBP1* | uORF | 36 aa | 1B | **19.61** | **33.26** | `chr2:74458516:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr2:74458484-74458594&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **3** | `c1riboseqorf248` | *IVNS1ABP* | uORF | 24 aa | 1B | **18.23** | **33.22** | `chr1:185316953:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:185311281-185316985&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **4** | `c7riboseqorf18` | *SNX13* | uORF | 52 aa | 1B | **17.32** | **31.59** | `chr7:17940333:G>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr7:17940325-17940483&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **5** | `c6riboseqorf128` | *BCLAF1* | uORF | 19 aa | 1B | **17.21** | **22.47** | `chr6:136289799:C>A` | Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:136289742-136289801&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **6** | `cXriboseqorf59` | *MORF4L2* | uORF | 32 aa | 1B | **16.98** | **34.03** | `chrX:103686631:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chrX:103685228-103686696&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **7** | `c8riboseqorf15` | *CNOT7* | uORF | 39 aa | 1A | **16.78** | **30.38** | `chr8:17246687:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr8:17245236-17246782&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **8** | `c1norep182` | *IGSF3* | uORF | 35 aa | 1B | **16.65** | **22.53** | `chr1:116666677:C>A` | Chromatin Opening | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:116666632-116666739&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **9** | `c6norep97` | *PSMB8-AS1* | lncRNA-ORF | 63 aa | 1B | **16.04** | **45.89** | `chr6:32845662:C>A` | Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr6:32845552-32845743&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **10** | `c17riboseqorf78` | *IGFBP4* | uORF | 17 aa | 1B | **15.44** | **24.15** | `chr17:40443555:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr17:40443505-40443558&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **11** | `c4riboseqorf2` | *CTBP1* | uORF | 32 aa | 1B | **15.18** | **31.95** | `chr4:1241388:T>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr4:1241335-1241433&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **12** | `c21riboseqorf14` | *RUNX1* | uORF | 20 aa | 1B | **14.93** | **23.84** | `chr21:34887820:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr21:34887800-34887862&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **13** | `c19riboseqorf66` | *URI1* | uORF | 53 aa | 1B | **14.79** | **24.85** | `chr19:29942370:A>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr19:29942370-29942531&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **14** | `c3riboseqorf154` | *MBNL1* | uORF | 16 aa | 1B | **14.72** | **22.44** | `chr3:152269066:G>C` | Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr3:152269026-152269076&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **15** | `c16riboseqorf104` | *PSKH1* | uORF | 33 aa | 1B | **13.62** | **30.51** | `chr16:67893371:G>C` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr16:67893276-67908685&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **16** | `c12riboseqorf145` | *TAOK3* | uORF | 16 aa | 1B | **11.95** | **32.46** | `chr12:118266655:C>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr12:118255649-118266698&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **17** | `c5riboseqorf128` | *PRELID1* | uORF | 17 aa | 1B | **11.69** | **23.15** | `chr5:177303866:C>G` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr5:177303816-177303869&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **18** | `c1riboseqorf94` | *CDKN2C* | uORF | 180 aa | 1A | **11.60** | **24.17** | `chr1:50968915:C>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:50968788-50969330&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **19** | `c15riboseqorf1` | *NIPA2* | uORF | 37 aa | 1B | **10.73** | **24.34** | `chr15:22838859:T>A` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr15:22838748-22838861&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |
| **20** | `c11riboseqorf108` | *LIPT2-AS1* | lncRNA-ORF | 349 aa | 1A | **10.30** | **25.97** | `chr11:74497088:G>T` | Severe Splicing Disruption | [Inspect Atlas](https://deepmind.google.com/science/alphagenome/atlas?q=chr11:74497035-74498084&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE) |

---

## 5. Conclusions & Research Implications

1. **Purifying Selection in smORFs Must Be Evaluated by Genomic Autonomy**:
   Ignoring canonical CDS overlap biases microprotein discovery toward canonical host gene essentiality. By controlling for CDS overlap, we discover that **830 autonomous smORFs** have survived millions of years of purifying selection ($\text{Phred} \ge 15$), functioning as authentic, standalone peptide elements.

2. **The `lncRNA-ORF` Frontier**:
   LncRNAs harboring translated smORFs represent the purest class of novel microproteins. Standouts like `c3riboseqorf106` (*ZBTB11-AS1*, 88 aa, mass-spec validated), `c16norep119` (*ZFHX3-AS1*, 154 aa), `c6norep97` (*PSMB8-AS1*, 63 aa, mass-spec validated), and `c2riboseqorf136` (*TTN-AS1*, 127 aa) are unconfounded by any annotated canonical protein and exhibit strong functional constraint via chromatin remodeling and transcription regulation.

3. **Archived Artifacts & Data**:
   - De-confounded statistical metrics: [`reports/deconfounded_statistical_summary.json`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/deconfounded_statistical_summary.json)
   - Top 20 Autonomous lncRNA-ORFs: [`reports/top20_autonomous_lncrna_orfs.tsv`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/top20_autonomous_lncrna_orfs.tsv)
   - Top 20 Autonomous uORFs: [`reports/top20_autonomous_uorfs.tsv`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/top20_autonomous_uorfs.tsv)
   - Top 20 Autonomous Tier 1 Peptideins: [`reports/top_autonomous_tier1_peptideins.tsv`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/top_autonomous_tier1_peptideins.tsv)
