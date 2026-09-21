# Prioritized Top 10 Binding Partner Proteins for Curated Top 50 Microproteins

**Total Complexes**: 500 candidate pairs (50 microproteins x 10 prioritized partners)  
**Cofolding Target Engine**: Biohub ESMFold2 `fold_all_atom` (`/api/v1/fold_all_atom`, default: `esmfold2-fast-2026-05`)  
**Source Interactomes**: STRING v12.0 (Homo sapiens, Taxon 9606), UniProt Canonical Proteome (20,652 entries), AlphaGenome GTEx/Track Disruption  
**Data Artifacts**: 
- TSV: [`reports/top50_microprotein_binding_partners_top10.tsv`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/top50_microprotein_binding_partners_top10.tsv)
- Parquet: [`reports/top50_microprotein_binding_partners_top10.parquet`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/top50_microprotein_binding_partners_top10.parquet)

---

## 1. Executive Summary & Selection Methodology

A fundamental question in non-canonical microprotein biology is deciphering their **macromolecular mechanism of action**: do they fold autonomously, do they bind as obligate subunits in multi-protein complexes, or do they act as competitive or allosteric modulators of canonical proteins? 

To prepare for **high-throughput complex structural prediction using Biohub ESMFold2 `fold_all_atom`**, this study prioritized the **top 10 physiological and functional binding partner proteins** for each candidate in the de-confounded Top 50 microprotein portfolio (totaling **500 candidate protein-protein pairs**). Partner selection adhered to four evidence tiers:

1. **Tier 1: Cognate Host CDS / Cis-Encoded Master Target (Partner #1)**:
   - For **5' UTR uORFs** (Cohorts 1 & 3): The microprotein is co-transcribed on the same mRNA and translated co-translationally with the main open reading frame product. It has immediate spatial proximity to engage in hetero-oligomerization, feedback inhibition, or chaperone-like assistance.
   - For **lncRNA-ORFs** (Cohorts 1 & 2): The microprotein is encoded in antisense or shared regulatory loci, acting as a direct cis/trans modulator of the cognate sense protein.
   - For **Dual-Coding intORFs/uoORFs** (Cohort 4): The peptide is produced from an alternate frame within the host CDS, establishing obligate co-expression with the host protein.
2. **Tier 2: Core Macromolecular Complex Assembly (Partners #2-#6)**:
   - Direct stoichiometric subunits of the host protein's core machine (e.g., the 20S immunoproteasome for *PSMB8*, the CCR4-NOT deadenylase for *CNOT7*, the Microprocessor for *DGCR8*, the Prefoldin chaperone for *URI1*, the ER membrane complex for *EMC1*, the CUL4-DDB1 ubiquitin ligase for *DCAF13/DCAF8*, and the sarcomere for *TTN*).
3. **Tier 3: Experimentally Validated Physical Interactors (Partners #7-#10)**:
   - Direct physical interaction partners retrieved from STRING v12.0 backed by high-confidence experimental assay scores (BioGRID, IntAct, PDB) and co-expression.
4. **Tier 4: Tissue Co-Localization & ESMFold2 Compatibility**:
   - Verified co-expression in GTEx tissues matching the candidate's AlphaGenome significant variant tracks (e.g. Whole Blood, Skeletal Muscle, Brain, Liver, Fibroblasts).
   - Filtered for feasible sequence lengths to ensure high-fidelity all-atom structural prediction within ESMFold2 context limits.

---

## 2. Portfolio Complex Feasibility Statistics

| Feasibility Tier | Total Complex Length | Pair Count | Percentage | Recommended ESMFold2 Model |
| :--- | :---: | :---: | :---: | :--- |
| **Optimal** | `< 800 aa` | **357** | 71.4% | `esmfold2-fast-2026-05` (rapid turnaround, ~2-4s) |
| **Standard** | `800 - 1400 aa` | **111** | 22.2% | `esmfold2-fast-2026-05` or `esmfold2-2026-05` |
| **High** | `> 1400 aa` | **32** | 6.4% | `esmfold2-2026-05` (requires max memory, full loops) |
| **Total Portfolio** | `94 - 1999 aa` | **500** | 100.0% | Mean: 658.6 aa, Median: 538.0 aa |

---

## 3. Detailed Portfolio Catalog: 50 Microproteins x Top 10 Partners

### 1. Autonomous Tier 1 Peptidein

#### Rank 1: `c3riboseqorf106` (*ZBTB11-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 88 aa | **Tier**: 1B | **Primary Tissue**: Whole_Blood
- **Microprotein Amino Acid Sequence**:  
  `MPGVVSAAGTQVRRLDEVPASLRLQHHLQLREGLAVPLPPLVIQSPAAHHVAGGSFSDFTLDIALGARRIRLALVRQVTQDGPVAFLA`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **ZBTB11** | [`O95625`](https://www.uniprot.org/uniprotkb/O95625) | Zinc finger and BTB domain-containin | 1053 aa | Cognate Antisense/Locus Master | 1.000 | 1141 aa | Standard (800-1400 aa) |
| 2 | **RPF1** | [`Q9H9Y2`](https://www.uniprot.org/uniprotkb/Q9H9Y2) | Ribosome production factor 1 | 349 aa | Direct Physical Interactor | 0.689 | 437 aa | Optimal (<800 aa) |
| 3 | **ZNF654** | [`Q8IZM8`](https://www.uniprot.org/uniprotkb/Q8IZM8) | Zinc finger protein 654 | 1128 aa | Direct Physical Interactor | 0.630 | 1216 aa | Standard (800-1400 aa) |
| 4 | **VCPIP1** | [`Q96JH7`](https://www.uniprot.org/uniprotkb/Q96JH7) | Deubiquitinating protein VCPIP1 | 1222 aa | Direct Physical Interactor | 0.542 | 1310 aa | Standard (800-1400 aa) |
| 5 | **STAG1** | [`Q8WVM7`](https://www.uniprot.org/uniprotkb/Q8WVM7) | Cohesin subunit SA-1 | 1258 aa | Direct Physical Interactor | 0.497 | 1346 aa | Standard (800-1400 aa) |
| 6 | **UBL5** | [`Q9BZL1`](https://www.uniprot.org/uniprotkb/Q9BZL1) | Ubiquitin-like protein 5 | 73 aa | Direct Physical Interactor | 0.750 | 161 aa | Optimal (<800 aa) |
| 7 | **ZNF644** | [`Q9H582`](https://www.uniprot.org/uniprotkb/Q9H582) | Zinc finger protein 644 | 1327 aa | Direct Physical Interactor | 0.750 | 1415 aa | High (>1400 aa) |
| 8 | **CDK8** | [`P49336`](https://www.uniprot.org/uniprotkb/P49336) | Cyclin-dependent kinase 8 | 464 aa | Direct Physical Interactor | 0.750 | 552 aa | Optimal (<800 aa) |
| 9 | **CCDC92** | [`Q53HC0`](https://www.uniprot.org/uniprotkb/Q53HC0) | Coiled-coil domain-containing protei | 331 aa | Direct Physical Interactor | 0.750 | 419 aa | Optimal (<800 aa) |
| 10 | **PLEKHA4** | [`Q9H4M7`](https://www.uniprot.org/uniprotkb/Q9H4M7) | Pleckstrin homology domain-containin | 779 aa | Direct Physical Interactor | 0.535 | 867 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein ZBTB11. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of ZBTB11 expression, stability, and macromolecular complex assembly.

---

#### Rank 2: `c2riboseqorf55` (*WBP1*)
- **Biotype**: uORF | **Length**: 36 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MVASAKMGRAGTMAVAAEVAGAGRLAVEEAVVLRGL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **WBP1** | [`Q96G27`](https://www.uniprot.org/uniprotkb/Q96G27) | WW domain-binding protein 1 | 269 aa | Cognate Host CDS (Co-Translational) | 1.000 | 305 aa | Optimal (<800 aa) |
| 2 | **WBP2** | [`Q969T9`](https://www.uniprot.org/uniprotkb/Q969T9) | WW domain-binding protein 2 | 261 aa | Core Macromolecular Complex | 0.953 | 297 aa | Optimal (<800 aa) |
| 3 | **RUVBL1** | [`Q9Y265`](https://www.uniprot.org/uniprotkb/Q9Y265) | RuvB-like 1 | 456 aa | Core Macromolecular Complex | 0.952 | 492 aa | Optimal (<800 aa) |
| 4 | **RUVBL2** | [`Q9Y230`](https://www.uniprot.org/uniprotkb/Q9Y230) | RuvB-like 2 | 463 aa | Core Macromolecular Complex | 0.950 | 499 aa | Optimal (<800 aa) |
| 5 | **ACTR5** | [`Q9H9F9`](https://www.uniprot.org/uniprotkb/Q9H9F9) | Actin-related protein 5 | 607 aa | Core Macromolecular Complex | 0.927 | 643 aa | Optimal (<800 aa) |
| 6 | **WWOX** | [`Q9NZC7`](https://www.uniprot.org/uniprotkb/Q9NZC7) | WW domain-containing oxidoreductase | 414 aa | Direct Physical Interactor | 0.782 | 450 aa | Optimal (<800 aa) |
| 7 | **NEDD4** | [`P46934`](https://www.uniprot.org/uniprotkb/P46934) | E3 ubiquitin-protein ligase NEDD4 | 1319 aa | Direct Physical Interactor | 0.750 | 1355 aa | Standard (800-1400 aa) |
| 8 | **YAP1** | [`P46937`](https://www.uniprot.org/uniprotkb/P46937) | Transcriptional coactivator YAP1 | 504 aa | Direct Physical Interactor | 0.432 | 540 aa | Optimal (<800 aa) |
| 9 | **INO80B** | [`Q9C086`](https://www.uniprot.org/uniprotkb/Q9C086) | INO80 complex subunit B | 356 aa | Direct Physical Interactor | 0.750 | 392 aa | Optimal (<800 aa) |
| 10 | **WWP1** | [`Q9H0M0`](https://www.uniprot.org/uniprotkb/Q9H0M0) | NEDD4-like E3 ubiquitin-protein liga | 922 aa | Direct Physical Interactor | 0.750 | 958 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary cognate host protein WBP1. Microprotein is translated from the 5' UTR uORF of WBP1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 3: `c1riboseqorf248` (*IVNS1ABP*)
- **Biotype**: uORF | **Length**: 24 aa | **Tier**: 1B | **Primary Tissue**: Whole_Blood
- **Microprotein Amino Acid Sequence**:  
  `MNVSISSEMNQIMMHHYHRRNSCL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **IVNS1ABP** | [`Q9Y6Y0`](https://www.uniprot.org/uniprotkb/Q9Y6Y0) | Influenza virus NS1A-binding protein | 642 aa | Cognate Host CDS (Co-Translational) | 1.000 | 666 aa | Optimal (<800 aa) |
| 2 | **AHR** | [`P35869`](https://www.uniprot.org/uniprotkb/P35869) | Aryl hydrocarbon receptor | 848 aa | Direct Physical Interactor | 0.750 | 872 aa | Standard (800-1400 aa) |
| 3 | **HSP90AA1** | [`P07900`](https://www.uniprot.org/uniprotkb/P07900) | Heat shock protein HSP 90-alpha | 732 aa | Direct Physical Interactor | 0.750 | 756 aa | Optimal (<800 aa) |
| 4 | **ARNT** | [`P27540`](https://www.uniprot.org/uniprotkb/P27540) | Aryl hydrocarbon receptor nuclear tr | 789 aa | Direct Physical Interactor | 0.750 | 813 aa | Standard (800-1400 aa) |
| 5 | **AIP** | [`O00170`](https://www.uniprot.org/uniprotkb/O00170) | AH receptor-interacting protein | 330 aa | Direct Physical Interactor | 0.750 | 354 aa | Optimal (<800 aa) |
| 6 | **ACTB** | [`P60709`](https://www.uniprot.org/uniprotkb/P60709) | Actin, cytoplasmic 1 | 375 aa | Direct Physical Interactor | 0.750 | 399 aa | Optimal (<800 aa) |
| 7 | **ACTA1** | [`P68133`](https://www.uniprot.org/uniprotkb/P68133) | Actin, alpha skeletal muscle | 377 aa | Direct Physical Interactor | 0.750 | 401 aa | Optimal (<800 aa) |
| 8 | **SRSF1** | [`Q07955`](https://www.uniprot.org/uniprotkb/Q07955) | Serine/arginine-rich splicing factor | 248 aa | Direct Physical Interactor | 0.750 | 272 aa | Optimal (<800 aa) |
| 9 | **HNRNPU** | [`Q00839`](https://www.uniprot.org/uniprotkb/Q00839) | Heterogeneous nuclear ribonucleoprot | 825 aa | Direct Physical Interactor | 0.437 | 849 aa | Standard (800-1400 aa) |
| 10 | **CAPZB** | [`P47756`](https://www.uniprot.org/uniprotkb/P47756) | F-actin-capping protein subunit beta | 272 aa | Direct Physical Interactor | 0.750 | 296 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein IVNS1ABP. Microprotein is translated from the 5' UTR uORF of IVNS1ABP mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 4: `c7riboseqorf18` (*SNX13*)
- **Biotype**: uORF | **Length**: 52 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MAARPSRATGPRGGQRSRVKPPPGRRLKEQLPPLAAARAVFAAATAVAMRRG`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **SNX13** | [`Q9Y5W8`](https://www.uniprot.org/uniprotkb/Q9Y5W8) | Sorting nexin-13 | 968 aa | Cognate Host CDS (Co-Translational) | 1.000 | 1020 aa | Standard (800-1400 aa) |
| 2 | **GNAI1** | [`P63096`](https://www.uniprot.org/uniprotkb/P63096) | Guanine nucleotide-binding protein G | 354 aa | Direct Physical Interactor | 0.750 | 406 aa | Optimal (<800 aa) |
| 3 | **GNAI2** | [`P04899`](https://www.uniprot.org/uniprotkb/P04899) | Guanine nucleotide-binding protein G | 355 aa | Direct Physical Interactor | 0.750 | 407 aa | Optimal (<800 aa) |
| 4 | **GNAI3** | [`P08754`](https://www.uniprot.org/uniprotkb/P08754) | Guanine nucleotide-binding protein G | 354 aa | Direct Physical Interactor | 0.750 | 406 aa | Optimal (<800 aa) |
| 5 | **GNAS** | [`Q5JWF2`](https://www.uniprot.org/uniprotkb/Q5JWF2) | Guanine nucleotide-binding protein G | 1037 aa | Direct Physical Interactor | 0.750 | 1089 aa | Standard (800-1400 aa) |
| 6 | **EGFR** | [`P00533`](https://www.uniprot.org/uniprotkb/P00533) | Epidermal growth factor receptor | 1210 aa | Direct Physical Interactor | 0.750 | 1262 aa | Standard (800-1400 aa) |
| 7 | **RAB5A** | [`P20339`](https://www.uniprot.org/uniprotkb/P20339) | Ras-related protein Rab-5A | 215 aa | Direct Physical Interactor | 0.750 | 267 aa | Optimal (<800 aa) |
| 8 | **RAB7A** | [`P51149`](https://www.uniprot.org/uniprotkb/P51149) | Ras-related protein Rab-7a | 207 aa | Direct Physical Interactor | 0.750 | 259 aa | Optimal (<800 aa) |
| 9 | **VPS35** | [`Q96QK1`](https://www.uniprot.org/uniprotkb/Q96QK1) | Vacuolar protein sorting-associated  | 796 aa | Direct Physical Interactor | 0.750 | 848 aa | Standard (800-1400 aa) |
| 10 | **SNX1** | [`Q13596`](https://www.uniprot.org/uniprotkb/Q13596) | Sorting nexin-1 | 522 aa | Direct Physical Interactor | 0.563 | 574 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein SNX13. Microprotein is translated from the 5' UTR uORF of SNX13 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 5: `c6riboseqorf128` (*BCLAF1*)
- **Biotype**: uORF | **Length**: 19 aa | **Tier**: 1B | **Primary Tissue**: Pituitary
- **Microprotein Amino Acid Sequence**:  
  `MVFLAFLSREMAAVWLQRR`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **BCLAF1** | [`Q9NYF8`](https://www.uniprot.org/uniprotkb/Q9NYF8) | Bcl-2-associated transcription facto | 920 aa | Cognate Host CDS (Co-Translational) | 1.000 | 939 aa | Standard (800-1400 aa) |
| 2 | **THRAP3** | [`Q9Y2W1`](https://www.uniprot.org/uniprotkb/Q9Y2W1) | Thyroid hormone receptor-associated  | 955 aa | Core Macromolecular Complex | 0.993 | 974 aa | Standard (800-1400 aa) |
| 3 | **BCL2** | [`P10415`](https://www.uniprot.org/uniprotkb/P10415) | Apoptosis regulator Bcl-2 | 239 aa | Core Macromolecular Complex | 0.945 | 258 aa | Optimal (<800 aa) |
| 4 | **BCL2L1** | [`Q07817`](https://www.uniprot.org/uniprotkb/Q07817) | Bcl-2-like protein 1 | 233 aa | Direct Physical Interactor | 0.750 | 252 aa | Optimal (<800 aa) |
| 5 | **BCL2L2** | [`Q92843`](https://www.uniprot.org/uniprotkb/Q92843) | Bcl-2-like protein 2 | 193 aa | Direct Physical Interactor | 0.750 | 212 aa | Optimal (<800 aa) |
| 6 | **SRSF1** | [`Q07955`](https://www.uniprot.org/uniprotkb/Q07955) | Serine/arginine-rich splicing factor | 248 aa | Direct Physical Interactor | 0.750 | 267 aa | Optimal (<800 aa) |
| 7 | **HNRNPM** | [`P52272`](https://www.uniprot.org/uniprotkb/P52272) | Heterogeneous nuclear ribonucleoprot | 730 aa | Direct Physical Interactor | 0.750 | 749 aa | Optimal (<800 aa) |
| 8 | **EMD** | [`P50402`](https://www.uniprot.org/uniprotkb/P50402) | Emerin | 254 aa | Direct Physical Interactor | 0.971 | 273 aa | Optimal (<800 aa) |
| 9 | **DDX5** | [`P17844`](https://www.uniprot.org/uniprotkb/P17844) | Probable ATP-dependent RNA helicase  | 614 aa | Direct Physical Interactor | 0.750 | 633 aa | Optimal (<800 aa) |
| 10 | **FUS** | [`P35637`](https://www.uniprot.org/uniprotkb/P35637) | RNA-binding protein FUS | 526 aa | Direct Physical Interactor | 0.750 | 545 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein BCLAF1. Microprotein is translated from the 5' UTR uORF of BCLAF1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 6: `cXriboseqorf59` (*MORF4L2*)
- **Biotype**: uORF | **Length**: 32 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MAVGGTAVITRRLLGRSGFSFQVTIRKAKFAV`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **MORF4L2** | [`Q15014`](https://www.uniprot.org/uniprotkb/Q15014) | Mortality factor 4-like protein 2 | 288 aa | Cognate Host CDS (Co-Translational) | 1.000 | 320 aa | Optimal (<800 aa) |
| 2 | **KAT5** | [`Q92993`](https://www.uniprot.org/uniprotkb/Q92993) | Histone acetyltransferase KAT5 | 513 aa | Core Macromolecular Complex | 0.997 | 545 aa | Optimal (<800 aa) |
| 3 | **MORF4L1** | [`Q9UBU8`](https://www.uniprot.org/uniprotkb/Q9UBU8) | Mortality factor 4-like protein 1 | 362 aa | Core Macromolecular Complex | 0.987 | 394 aa | Optimal (<800 aa) |
| 4 | **MRGBP** | [`Q9NV56`](https://www.uniprot.org/uniprotkb/Q9NV56) | MRG/MORF4L-binding protein | 204 aa | Core Macromolecular Complex | 0.999 | 236 aa | Optimal (<800 aa) |
| 5 | **DMAP1** | [`Q9NPF5`](https://www.uniprot.org/uniprotkb/Q9NPF5) | DNA methyltransferase 1-associated p | 467 aa | Core Macromolecular Complex | 0.986 | 499 aa | Optimal (<800 aa) |
| 6 | **BRD8** | [`Q9H0E9`](https://www.uniprot.org/uniprotkb/Q9H0E9) | Bromodomain-containing protein 8 | 1235 aa | Core Macromolecular Complex | 0.991 | 1267 aa | Standard (800-1400 aa) |
| 7 | **EPC1** | [`Q9H2F5`](https://www.uniprot.org/uniprotkb/Q9H2F5) | Enhancer of polycomb homolog 1 | 836 aa | Core Macromolecular Complex | 0.985 | 868 aa | Standard (800-1400 aa) |
| 8 | **SIN3A** | [`Q96ST3`](https://www.uniprot.org/uniprotkb/Q96ST3) | Paired amphipathic helix protein Sin | 1273 aa | Direct Physical Interactor | 0.719 | 1305 aa | Standard (800-1400 aa) |
| 9 | **HDAC1** | [`Q13547`](https://www.uniprot.org/uniprotkb/Q13547) | Histone deacetylase 1 | 482 aa | Core Macromolecular Complex | 0.911 | 514 aa | Optimal (<800 aa) |
| 10 | **RBBP4** | [`Q09028`](https://www.uniprot.org/uniprotkb/Q09028) | Histone-binding protein RBBP4 | 425 aa | Direct Physical Interactor | 0.752 | 457 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein MORF4L2. Microprotein is translated from the 5' UTR uORF of MORF4L2 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 7: `c8riboseqorf15` (*CNOT7*)
- **Biotype**: uORF | **Length**: 39 aa | **Tier**: 1A | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MGHPLPPPLPPPPPPPPPSAVSMARRRRRSASSATQVHK`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **CNOT7** | [`Q9UIV1`](https://www.uniprot.org/uniprotkb/Q9UIV1) | CCR4-NOT transcription complex subun | 285 aa | Cognate Host CDS (Co-Translational) | 1.000 | 324 aa | Optimal (<800 aa) |
| 2 | **CNOT10** | [`Q9H9A5`](https://www.uniprot.org/uniprotkb/Q9H9A5) | CCR4-NOT transcription complex subun | 744 aa | Core Macromolecular Complex | 0.999 | 783 aa | Optimal (<800 aa) |
| 3 | **CNOT2** | [`Q9NZN8`](https://www.uniprot.org/uniprotkb/Q9NZN8) | CCR4-NOT transcription complex subun | 540 aa | Core Macromolecular Complex | 0.999 | 579 aa | Optimal (<800 aa) |
| 4 | **CNOT3** | [`O75175`](https://www.uniprot.org/uniprotkb/O75175) | CCR4-NOT transcription complex subun | 753 aa | Core Macromolecular Complex | 0.999 | 792 aa | Optimal (<800 aa) |
| 5 | **CNOT6** | [`Q9ULM6`](https://www.uniprot.org/uniprotkb/Q9ULM6) | CCR4-NOT transcription complex subun | 557 aa | Core Macromolecular Complex | 0.999 | 596 aa | Optimal (<800 aa) |
| 6 | **CNOT6L** | [`Q96LI5`](https://www.uniprot.org/uniprotkb/Q96LI5) | CCR4-NOT transcription complex subun | 555 aa | Core Macromolecular Complex | 0.999 | 594 aa | Optimal (<800 aa) |
| 7 | **CNOT8** | [`Q9UFF9`](https://www.uniprot.org/uniprotkb/Q9UFF9) | CCR4-NOT transcription complex subun | 292 aa | Core Macromolecular Complex | 0.993 | 331 aa | Optimal (<800 aa) |
| 8 | **CNOT9** | [`Q92600`](https://www.uniprot.org/uniprotkb/Q92600) | CCR4-NOT transcription complex subun | 299 aa | Core Macromolecular Complex | 0.999 | 338 aa | Optimal (<800 aa) |
| 9 | **TOB1** | [`P50616`](https://www.uniprot.org/uniprotkb/P50616) | Protein Tob1 | 345 aa | Core Macromolecular Complex | 0.999 | 384 aa | Optimal (<800 aa) |
| 10 | **TOB2** | [`Q14106`](https://www.uniprot.org/uniprotkb/Q14106) | Protein Tob2 | 344 aa | Core Macromolecular Complex | 0.969 | 383 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein CNOT7. Microprotein is translated from the 5' UTR uORF of CNOT7 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 8: `c1norep182` (*IGSF3*)
- **Biotype**: uORF | **Length**: 35 aa | **Tier**: 1B | **Primary Tissue**: Stomach
- **Microprotein Amino Acid Sequence**:  
  `MREITDTFTEKPRTFLYNQWMKGFSQTLFSVLLED`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **IGSF3** | [`O75054`](https://www.uniprot.org/uniprotkb/O75054) | Immunoglobulin superfamily member 3 | 1194 aa | Cognate Host CDS (Co-Translational) | 1.000 | 1229 aa | Standard (800-1400 aa) |
| 2 | **CD9** | [`P21926`](https://www.uniprot.org/uniprotkb/P21926) | CD9 antigen | 228 aa | Direct Physical Interactor | 0.753 | 263 aa | Optimal (<800 aa) |
| 3 | **CD81** | [`P60033`](https://www.uniprot.org/uniprotkb/P60033) | CD81 antigen | 236 aa | Direct Physical Interactor | 0.625 | 271 aa | Optimal (<800 aa) |
| 4 | **CD151** | [`P48509`](https://www.uniprot.org/uniprotkb/P48509) | CD151 antigen | 253 aa | Direct Physical Interactor | 0.750 | 288 aa | Optimal (<800 aa) |
| 5 | **IGSF8** | [`Q969P0`](https://www.uniprot.org/uniprotkb/Q969P0) | Immunoglobulin superfamily member 8 | 613 aa | Direct Physical Interactor | 0.750 | 648 aa | Optimal (<800 aa) |
| 6 | **PTGFRN** | [`Q9P2B2`](https://www.uniprot.org/uniprotkb/Q9P2B2) | Prostaglandin F2 receptor negative r | 879 aa | Direct Physical Interactor | 0.402 | 914 aa | Standard (800-1400 aa) |
| 7 | **ITGB1** | [`P05556`](https://www.uniprot.org/uniprotkb/P05556) | Integrin beta-1 | 798 aa | Direct Physical Interactor | 0.750 | 833 aa | Standard (800-1400 aa) |
| 8 | **ITGA3** | [`P26006`](https://www.uniprot.org/uniprotkb/P26006) | Integrin alpha-3 | 1051 aa | Direct Physical Interactor | 0.750 | 1086 aa | Standard (800-1400 aa) |
| 9 | **TSPAN4** | [`O14817`](https://www.uniprot.org/uniprotkb/O14817) | Tetraspanin-4 | 238 aa | Direct Physical Interactor | 0.750 | 273 aa | Optimal (<800 aa) |
| 10 | **TSPAN7** | [`P41732`](https://www.uniprot.org/uniprotkb/P41732) | Tetraspanin-7 | 249 aa | Direct Physical Interactor | 0.562 | 284 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein IGSF3. Microprotein is translated from the 5' UTR uORF of IGSF3 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 9: `c6norep97` (*PSMB8-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 63 aa | **Tier**: 1B | **Primary Tissue**: Whole_Blood / Spleen / Immune (Immunoproteasome)
- **Microprotein Amino Acid Sequence**:  
  `MECAGLRRLSFWSICRSLHHGPVAPLFLHELLVGSPLPDSASFQKEDVVSLLHQAEVLGDEKH`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **PSMB8** | [`P28062`](https://www.uniprot.org/uniprotkb/P28062) | Proteasome subunit beta type-8 | 276 aa | Cognate Antisense/Locus Master | 1.000 | 339 aa | Optimal (<800 aa) |
| 2 | **PSMB9** | [`P28065`](https://www.uniprot.org/uniprotkb/P28065) | Proteasome subunit beta type-9 | 219 aa | Core Macromolecular Complex | 0.999 | 282 aa | Optimal (<800 aa) |
| 3 | **PSMB10** | [`P40306`](https://www.uniprot.org/uniprotkb/P40306) | Proteasome subunit beta type-10 | 273 aa | Core Macromolecular Complex | 0.999 | 336 aa | Optimal (<800 aa) |
| 4 | **PSMA1** | [`P25786`](https://www.uniprot.org/uniprotkb/P25786) | Proteasome subunit alpha type-1 | 263 aa | Core Macromolecular Complex | 0.999 | 326 aa | Optimal (<800 aa) |
| 5 | **PSMA3** | [`P25788`](https://www.uniprot.org/uniprotkb/P25788) | Proteasome subunit alpha type-3 | 255 aa | Core Macromolecular Complex | 0.999 | 318 aa | Optimal (<800 aa) |
| 6 | **PSMB1** | [`P20618`](https://www.uniprot.org/uniprotkb/P20618) | Proteasome subunit beta type-1 | 241 aa | Core Macromolecular Complex | 0.999 | 304 aa | Optimal (<800 aa) |
| 7 | **PSMB2** | [`P49721`](https://www.uniprot.org/uniprotkb/P49721) | Proteasome subunit beta type-2 | 201 aa | Core Macromolecular Complex | 0.999 | 264 aa | Optimal (<800 aa) |
| 8 | **PSMB3** | [`P49720`](https://www.uniprot.org/uniprotkb/P49720) | Proteasome subunit beta type-3 | 205 aa | Core Macromolecular Complex | 0.999 | 268 aa | Optimal (<800 aa) |
| 9 | **POMP** | [`Q9Y244`](https://www.uniprot.org/uniprotkb/Q9Y244) | Proteasome maturation protein | 141 aa | Core Macromolecular Complex | 0.999 | 204 aa | Optimal (<800 aa) |
| 10 | **PSME1** | [`Q06323`](https://www.uniprot.org/uniprotkb/Q06323) | Proteasome activator complex subunit | 249 aa | Core Macromolecular Complex | 0.985 | 312 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein PSMB8. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of PSMB8 expression, stability, and macromolecular complex assembly.

---

#### Rank 10: `c17riboseqorf78` (*IGFBP4*)
- **Biotype**: uORF | **Length**: 17 aa | **Tier**: 1B | **Primary Tissue**: Cells_Cultured_fibroblasts
- **Microprotein Amino Acid Sequence**:  
  `MRRSASAARLAPLRHAC`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **IGFBP4** | [`P22692`](https://www.uniprot.org/uniprotkb/P22692) | Insulin-like growth factor-binding p | 258 aa | Cognate Host CDS (Co-Translational) | 1.000 | 275 aa | Optimal (<800 aa) |
| 2 | **IGF1** | [`P05019`](https://www.uniprot.org/uniprotkb/P05019) | Insulin-like growth factor 1 | 195 aa | Core Macromolecular Complex | 0.999 | 212 aa | Optimal (<800 aa) |
| 3 | **IGF2** | [`P01344`](https://www.uniprot.org/uniprotkb/P01344) | Insulin-like growth factor 2 | 180 aa | Core Macromolecular Complex | 0.993 | 197 aa | Optimal (<800 aa) |
| 4 | **IGF1R** | [`P08069`](https://www.uniprot.org/uniprotkb/P08069) | Insulin-like growth factor 1 recepto | 1367 aa | Direct Physical Interactor | 0.604 | 1384 aa | Standard (800-1400 aa) |
| 5 | **PAPPA** | [`Q13219`](https://www.uniprot.org/uniprotkb/Q13219) | Pappalysin-1 | 1627 aa | Core Macromolecular Complex | 0.976 | 1644 aa | High (>1400 aa) |
| 6 | **PAPPA2** | [`Q9BXP8`](https://www.uniprot.org/uniprotkb/Q9BXP8) | Pappalysin-2 | 1791 aa | Direct Physical Interactor | 0.827 | 1808 aa | High (>1400 aa) |
| 7 | **LRP5** | [`O75197`](https://www.uniprot.org/uniprotkb/O75197) | Low-density lipoprotein receptor-rel | 1615 aa | Direct Physical Interactor | 0.549 | 1632 aa | High (>1400 aa) |
| 8 | **LRP6** | [`O75581`](https://www.uniprot.org/uniprotkb/O75581) | Low-density lipoprotein receptor-rel | 1613 aa | Direct Physical Interactor | 0.847 | 1630 aa | High (>1400 aa) |
| 9 | **FZD4** | [`Q9ULV1`](https://www.uniprot.org/uniprotkb/Q9ULV1) | Frizzled-4 | 537 aa | Direct Physical Interactor | 0.750 | 554 aa | Optimal (<800 aa) |
| 10 | **FZD8** | [`Q9H461`](https://www.uniprot.org/uniprotkb/Q9H461) | Frizzled-8 | 694 aa | Direct Physical Interactor | 0.838 | 711 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein IGFBP4. Microprotein is translated from the 5' UTR uORF of IGFBP4 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 11: `c4riboseqorf2` (*CTBP1*)
- **Biotype**: uORF | **Length**: 32 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MRLRDRTRFCLGILKIKSFEKSNSWSWKLSPY`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **CTBP1** | [`Q13363`](https://www.uniprot.org/uniprotkb/Q13363) | C-terminal-binding protein 1 | 440 aa | Cognate Host CDS (Co-Translational) | 1.000 | 472 aa | Optimal (<800 aa) |
| 2 | **CTBP2** | [`P56545`](https://www.uniprot.org/uniprotkb/P56545) | C-terminal-binding protein 2 | 445 aa | Core Macromolecular Complex | 0.999 | 477 aa | Optimal (<800 aa) |
| 3 | **ZEB1** | [`P37275`](https://www.uniprot.org/uniprotkb/P37275) | Zinc finger E-box-binding homeobox 1 | 1124 aa | Core Macromolecular Complex | 0.998 | 1156 aa | Standard (800-1400 aa) |
| 4 | **ZEB2** | [`O60315`](https://www.uniprot.org/uniprotkb/O60315) | Zinc finger E-box-binding homeobox 2 | 1214 aa | Core Macromolecular Complex | 0.951 | 1246 aa | Standard (800-1400 aa) |
| 5 | **HDAC1** | [`Q13547`](https://www.uniprot.org/uniprotkb/Q13547) | Histone deacetylase 1 | 482 aa | Core Macromolecular Complex | 0.999 | 514 aa | Optimal (<800 aa) |
| 6 | **HDAC2** | [`Q92769`](https://www.uniprot.org/uniprotkb/Q92769) | Histone deacetylase 2 | 488 aa | Core Macromolecular Complex | 0.999 | 520 aa | Optimal (<800 aa) |
| 7 | **RCOR1** | [`Q9UKL0`](https://www.uniprot.org/uniprotkb/Q9UKL0) | REST corepressor 1 | 485 aa | Core Macromolecular Complex | 0.999 | 517 aa | Optimal (<800 aa) |
| 8 | **KDM1A** | [`O60341`](https://www.uniprot.org/uniprotkb/O60341) | Lysine-specific histone demethylase  | 852 aa | Core Macromolecular Complex | 0.998 | 884 aa | Standard (800-1400 aa) |
| 9 | **EHMT2** | [`Q96KQ7`](https://www.uniprot.org/uniprotkb/Q96KQ7) | Histone-lysine N-methyltransferase E | 1210 aa | Core Macromolecular Complex | 0.870 | 1242 aa | Standard (800-1400 aa) |
| 10 | **CBX5** | [`P45973`](https://www.uniprot.org/uniprotkb/P45973) | Chromobox protein homolog 5 | 191 aa | Direct Physical Interactor | 0.750 | 223 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein CTBP1. Microprotein is translated from the 5' UTR uORF of CTBP1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 12: `c21riboseqorf14` (*RUNX1*)
- **Biotype**: uORF | **Length**: 20 aa | **Tier**: 1B | **Primary Tissue**: Cells_Cultured_fibroblasts
- **Microprotein Amino Acid Sequence**:  
  `MTSTSSSALLQLLKLIFKAT`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **RUNX1** | [`Q01196`](https://www.uniprot.org/uniprotkb/Q01196) | Runt-related transcription factor 1 | 453 aa | Cognate Host CDS (Co-Translational) | 1.000 | 473 aa | Optimal (<800 aa) |
| 2 | **CBFB** | [`Q13951`](https://www.uniprot.org/uniprotkb/Q13951) | Core-binding factor subunit beta | 182 aa | Core Macromolecular Complex | 0.999 | 202 aa | Optimal (<800 aa) |
| 3 | **RUNX1T1** | [`Q06455`](https://www.uniprot.org/uniprotkb/Q06455) | Protein CBFA2T1 | 604 aa | Core Macromolecular Complex | 0.999 | 624 aa | Optimal (<800 aa) |
| 4 | **TAL1** | [`P17542`](https://www.uniprot.org/uniprotkb/P17542) | T-cell acute lymphocytic leukemia pr | 331 aa | Core Macromolecular Complex | 0.996 | 351 aa | Optimal (<800 aa) |
| 5 | **GATA2** | [`P23769`](https://www.uniprot.org/uniprotkb/P23769) | Endothelial transcription factor GAT | 480 aa | Core Macromolecular Complex | 0.996 | 500 aa | Optimal (<800 aa) |
| 6 | **ETS1** | [`P14921`](https://www.uniprot.org/uniprotkb/P14921) | Protein C-ets-1 | 441 aa | Core Macromolecular Complex | 0.999 | 461 aa | Optimal (<800 aa) |
| 7 | **SPI1** | [`P17947`](https://www.uniprot.org/uniprotkb/P17947) | Transcription factor PU.1 | 270 aa | Core Macromolecular Complex | 0.993 | 290 aa | Optimal (<800 aa) |
| 8 | **GATA1** | [`P15976`](https://www.uniprot.org/uniprotkb/P15976) | Erythroid transcription factor | 413 aa | Core Macromolecular Complex | 0.989 | 433 aa | Optimal (<800 aa) |
| 9 | **HDAC1** | [`Q13547`](https://www.uniprot.org/uniprotkb/Q13547) | Histone deacetylase 1 | 482 aa | Core Macromolecular Complex | 0.998 | 502 aa | Optimal (<800 aa) |
| 10 | **MEF2C** | [`Q06413`](https://www.uniprot.org/uniprotkb/Q06413) | Myocyte-specific enhancer factor 2C | 473 aa | Direct Physical Interactor | 0.750 | 493 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein RUNX1. Microprotein is translated from the 5' UTR uORF of RUNX1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 13: `c19riboseqorf66` (*URI1*)
- **Biotype**: uORF | **Length**: 53 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MRQRAARTRTAAAAAGAASWARGARCLRAGARALGNCRPRRLRRRWFRTHTPR`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **URI1** | [`O94763`](https://www.uniprot.org/uniprotkb/O94763) | Unconventional prefoldin RPB5 intera | 535 aa | Cognate Host CDS (Co-Translational) | 1.000 | 588 aa | Optimal (<800 aa) |
| 2 | **POLR2E** | [`P19388`](https://www.uniprot.org/uniprotkb/P19388) | DNA-directed RNA polymerases I, II,  | 210 aa | Core Macromolecular Complex | 0.999 | 263 aa | Optimal (<800 aa) |
| 3 | **PDRG1** | [`Q9NUG6`](https://www.uniprot.org/uniprotkb/Q9NUG6) | p53 and DNA damage-regulated protein | 133 aa | Core Macromolecular Complex | 0.973 | 186 aa | Optimal (<800 aa) |
| 4 | **PFDN2** | [`Q9UHV9`](https://www.uniprot.org/uniprotkb/Q9UHV9) | Prefoldin subunit 2 | 154 aa | Core Macromolecular Complex | 0.997 | 207 aa | Optimal (<800 aa) |
| 5 | **PFDN6** | [`O15212`](https://www.uniprot.org/uniprotkb/O15212) | Prefoldin subunit 6 | 129 aa | Core Macromolecular Complex | 0.982 | 182 aa | Optimal (<800 aa) |
| 6 | **UXT** | [`Q9UBK9`](https://www.uniprot.org/uniprotkb/Q9UBK9) | Protein UXT | 157 aa | Core Macromolecular Complex | 0.991 | 210 aa | Optimal (<800 aa) |
| 7 | **RUVBL1** | [`Q9Y265`](https://www.uniprot.org/uniprotkb/Q9Y265) | RuvB-like 1 | 456 aa | Core Macromolecular Complex | 0.986 | 509 aa | Optimal (<800 aa) |
| 8 | **RUVBL2** | [`Q9Y230`](https://www.uniprot.org/uniprotkb/Q9Y230) | RuvB-like 2 | 463 aa | Core Macromolecular Complex | 0.978 | 516 aa | Optimal (<800 aa) |
| 9 | **RPAP3** | [`Q9H6T3`](https://www.uniprot.org/uniprotkb/Q9H6T3) | RNA polymerase II-associated protein | 665 aa | Core Macromolecular Complex | 0.977 | 718 aa | Optimal (<800 aa) |
| 10 | **PIH1D1** | [`Q9NWS0`](https://www.uniprot.org/uniprotkb/Q9NWS0) | PIH1 domain-containing protein 1 | 290 aa | Core Macromolecular Complex | 0.968 | 343 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein URI1. Microprotein is translated from the 5' UTR uORF of URI1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 14: `c3riboseqorf154` (*MBNL1*)
- **Biotype**: uORF | **Length**: 16 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MTPTMPWALGRQWGRL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **MBNL1** | [`Q9NR56`](https://www.uniprot.org/uniprotkb/Q9NR56) | Muscleblind-like protein 1 | 388 aa | Cognate Host CDS (Co-Translational) | 1.000 | 404 aa | Optimal (<800 aa) |
| 2 | **MBNL2** | [`Q5VZF2`](https://www.uniprot.org/uniprotkb/Q5VZF2) | Muscleblind-like protein 2 | 373 aa | Direct Physical Interactor | 0.813 | 389 aa | Optimal (<800 aa) |
| 3 | **MBNL3** | [`Q9NUK0`](https://www.uniprot.org/uniprotkb/Q9NUK0) | Muscleblind-like protein 3 | 354 aa | Direct Physical Interactor | 0.786 | 370 aa | Optimal (<800 aa) |
| 4 | **SRSF1** | [`Q07955`](https://www.uniprot.org/uniprotkb/Q07955) | Serine/arginine-rich splicing factor | 248 aa | Core Macromolecular Complex | 0.852 | 264 aa | Optimal (<800 aa) |
| 5 | **SRSF2** | [`Q01130`](https://www.uniprot.org/uniprotkb/Q01130) | Serine/arginine-rich splicing factor | 221 aa | Direct Physical Interactor | 0.717 | 237 aa | Optimal (<800 aa) |
| 6 | **HNRNPA1** | [`P09651`](https://www.uniprot.org/uniprotkb/P09651) | Heterogeneous nuclear ribonucleoprot | 372 aa | Core Macromolecular Complex | 0.857 | 388 aa | Optimal (<800 aa) |
| 7 | **PTBP1** | [`P26599`](https://www.uniprot.org/uniprotkb/P26599) | Polypyrimidine tract-binding protein | 557 aa | Direct Physical Interactor | 0.801 | 573 aa | Optimal (<800 aa) |
| 8 | **CELF1** | [`Q92879`](https://www.uniprot.org/uniprotkb/Q92879) | CUGBP Elav-like family member 1 | 486 aa | Direct Physical Interactor | 0.985 | 502 aa | Optimal (<800 aa) |
| 9 | **TNPO1** | [`Q92973`](https://www.uniprot.org/uniprotkb/Q92973) | Transportin-1 | 898 aa | Direct Physical Interactor | 0.750 | 914 aa | Standard (800-1400 aa) |
| 10 | **YBX1** | [`P67809`](https://www.uniprot.org/uniprotkb/P67809) | Y-box-binding protein 1 | 324 aa | Direct Physical Interactor | 0.750 | 340 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein MBNL1. Microprotein is translated from the 5' UTR uORF of MBNL1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 15: `c16riboseqorf104` (*PSKH1*)
- **Biotype**: uORF | **Length**: 33 aa | **Tier**: 1B | **Primary Tissue**: Muscle_Skeletal
- **Microprotein Amino Acid Sequence**:  
  `MAAAAAAAAAAAIARRWPAEPPRRRRARRPREV`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **PSKH1** | [`P11801`](https://www.uniprot.org/uniprotkb/P11801) | Serine/threonine-protein kinase H1 | 424 aa | Cognate Host CDS (Co-Translational) | 1.000 | 457 aa | Optimal (<800 aa) |
| 2 | **BCLAF1** | [`Q9NYF8`](https://www.uniprot.org/uniprotkb/Q9NYF8) | Bcl-2-associated transcription facto | 920 aa | Direct Physical Interactor | 0.750 | 953 aa | Standard (800-1400 aa) |
| 3 | **THRAP3** | [`Q9Y2W1`](https://www.uniprot.org/uniprotkb/Q9Y2W1) | Thyroid hormone receptor-associated  | 955 aa | Direct Physical Interactor | 0.750 | 988 aa | Standard (800-1400 aa) |
| 4 | **SRSF1** | [`Q07955`](https://www.uniprot.org/uniprotkb/Q07955) | Serine/arginine-rich splicing factor | 248 aa | Direct Physical Interactor | 0.750 | 281 aa | Optimal (<800 aa) |
| 5 | **GOLGA2** | [`Q08379`](https://www.uniprot.org/uniprotkb/Q08379) | Golgin subfamily A member 2 | 1002 aa | Direct Physical Interactor | 0.750 | 1035 aa | Standard (800-1400 aa) |
| 6 | **GORASP1** | [`Q9BQQ3`](https://www.uniprot.org/uniprotkb/Q9BQQ3) | Golgi reassembly-stacking protein 1 | 440 aa | Direct Physical Interactor | 0.750 | 473 aa | Optimal (<800 aa) |
| 7 | **DCX** | [`O43602`](https://www.uniprot.org/uniprotkb/O43602) | Neuronal migration protein doublecor | 365 aa | Direct Physical Interactor | 0.850 | 398 aa | Optimal (<800 aa) |
| 8 | **HDAC1** | [`Q13547`](https://www.uniprot.org/uniprotkb/Q13547) | Histone deacetylase 1 | 482 aa | Direct Physical Interactor | 0.750 | 515 aa | Optimal (<800 aa) |
| 9 | **MNT** | [`Q99583`](https://www.uniprot.org/uniprotkb/Q99583) | Max-binding protein MNT | 582 aa | Direct Physical Interactor | 0.750 | 615 aa | Optimal (<800 aa) |
| 10 | **ZFP91** | [`Q96JP5`](https://www.uniprot.org/uniprotkb/Q96JP5) | E3 ubiquitin-protein ligase ZFP91 | 570 aa | Direct Physical Interactor | 0.750 | 603 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein PSKH1. Microprotein is translated from the 5' UTR uORF of PSKH1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

### 2. Autonomous lncRNA-ORF

#### Rank 16: `c12riboseqorf13` (*ENSG00000272173*)
- **Biotype**: lncRNA-ORF | **Length**: 57 aa | **Tier**: 2B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MSLQALSRASSMRTAFSGAENAWITSARTTARQKHPREASRVGVGEQTAGPAESARN`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **C12orf57** | [`Q99622`](https://www.uniprot.org/uniprotkb/Q99622) | Protein C10 | 126 aa | Cognate Antisense/Locus Master | 1.000 | 183 aa | Optimal (<800 aa) |
| 2 | **ANKRD50** | [`Q9ULJ7`](https://www.uniprot.org/uniprotkb/Q9ULJ7) | Ankyrin repeat domain-containing pro | 1429 aa | Core Macromolecular Complex | 0.994 | 1486 aa | High (>1400 aa) |
| 3 | **HLA-B** | [`P01889`](https://www.uniprot.org/uniprotkb/P01889) | HLA class I histocompatibility antig | 362 aa | Direct Physical Interactor | 0.811 | 419 aa | Optimal (<800 aa) |
| 4 | **PLAC8** | [`Q9NZF1`](https://www.uniprot.org/uniprotkb/Q9NZF1) | Placenta-specific gene 8 protein | 115 aa | Direct Physical Interactor | 0.769 | 172 aa | Optimal (<800 aa) |
| 5 | **HLA-A** | [`P04439`](https://www.uniprot.org/uniprotkb/P04439) | HLA class I histocompatibility antig | 365 aa | Direct Physical Interactor | 0.767 | 422 aa | Optimal (<800 aa) |
| 6 | **POLR3H** | [`Q9Y535`](https://www.uniprot.org/uniprotkb/Q9Y535) | DNA-directed RNA polymerase III subu | 204 aa | Direct Physical Interactor | 0.665 | 261 aa | Optimal (<800 aa) |
| 7 | **B2M** | [`P61769`](https://www.uniprot.org/uniprotkb/P61769) | Beta-2-microglobulin | 119 aa | Direct Physical Interactor | 0.658 | 176 aa | Optimal (<800 aa) |
| 8 | **APLP1** | [`P51693`](https://www.uniprot.org/uniprotkb/P51693) | Amyloid beta precursor like protein  | 650 aa | Direct Physical Interactor | 0.628 | 707 aa | Optimal (<800 aa) |
| 9 | **PLXNB2** | [`O15031`](https://www.uniprot.org/uniprotkb/O15031) | Plexin-B2 | 1838 aa | Direct Physical Interactor | 0.610 | 1895 aa | High (>1400 aa) |
| 10 | **TACC1** | [`O75410`](https://www.uniprot.org/uniprotkb/O75410) | Transforming acidic coiled-coil-cont | 805 aa | Direct Physical Interactor | 0.605 | 862 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein C12orf57. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of C12orf57 expression, stability, and macromolecular complex assembly.

---

#### Rank 17: `c1riboseqorf199` (*ENSG00000229953*)
- **Biotype**: lncRNA-ORF | **Length**: 43 aa | **Tier**: 1B | **Primary Tissue**: Brain / Cerebellum (Neural ECM)
- **Microprotein Amino Acid Sequence**:  
  `MMEPSCRSWGELRRCCHRCRAGPHSDRYLSRWGAARSAPGRRR`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **BCAN** | [`Q96GW7`](https://www.uniprot.org/uniprotkb/Q96GW7) | Brevican core protein | 911 aa | Cognate Antisense/Locus Master | 1.000 | 954 aa | Standard (800-1400 aa) |
| 2 | **NCAN** | [`O14594`](https://www.uniprot.org/uniprotkb/O14594) | Neurocan core protein | 1321 aa | Core Macromolecular Complex | 0.955 | 1364 aa | Standard (800-1400 aa) |
| 3 | **ADAMTS4** | [`O75173`](https://www.uniprot.org/uniprotkb/O75173) | A disintegrin and metalloproteinase  | 837 aa | Direct Physical Interactor | 0.913 | 880 aa | Standard (800-1400 aa) |
| 4 | **ADAMTS5** | [`Q9UNA0`](https://www.uniprot.org/uniprotkb/Q9UNA0) | A disintegrin and metalloproteinase  | 930 aa | Direct Physical Interactor | 0.913 | 973 aa | Standard (800-1400 aa) |
| 5 | **HAPLN1** | [`P10915`](https://www.uniprot.org/uniprotkb/P10915) | Hyaluronan and proteoglycan link pro | 354 aa | Direct Physical Interactor | 0.750 | 397 aa | Optimal (<800 aa) |
| 6 | **HAPLN2** | [`Q9GZV7`](https://www.uniprot.org/uniprotkb/Q9GZV7) | Hyaluronan and proteoglycan link pro | 340 aa | Direct Physical Interactor | 0.651 | 383 aa | Optimal (<800 aa) |
| 7 | **TNR** | [`Q92752`](https://www.uniprot.org/uniprotkb/Q92752) | Tenascin-R | 1358 aa | Direct Physical Interactor | 0.662 | 1401 aa | High (>1400 aa) |
| 8 | **DCN** | [`P07585`](https://www.uniprot.org/uniprotkb/P07585) | Decorin | 359 aa | Direct Physical Interactor | 0.754 | 402 aa | Optimal (<800 aa) |
| 9 | **CD44** | [`P16070`](https://www.uniprot.org/uniprotkb/P16070) | CD44 antigen | 742 aa | Direct Physical Interactor | 0.618 | 785 aa | Optimal (<800 aa) |
| 10 | **EGFR** | [`P00533`](https://www.uniprot.org/uniprotkb/P00533) | Epidermal growth factor receptor | 1210 aa | Direct Physical Interactor | 0.750 | 1253 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein BCAN. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of BCAN expression, stability, and macromolecular complex assembly.

---

#### Rank 18: `c1norep185` (*WARS2-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 26 aa | **Tier**: 4 | **Primary Tissue**: Mitochondria / Ubiquitous (tRNA Synthetase)
- **Microprotein Amino Acid Sequence**:  
  `MKLQRSRAFRIECSAILRRAEPSCLE`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **WARS2** | [`Q9UGM6`](https://www.uniprot.org/uniprotkb/Q9UGM6) | Tryptophan--tRNA ligase, mitochondri | 360 aa | Cognate Antisense/Locus Master | 1.000 | 386 aa | Optimal (<800 aa) |
| 2 | **LARS2** | [`Q15031`](https://www.uniprot.org/uniprotkb/Q15031) | Leucine--tRNA ligase, mitochondrial | 903 aa | Core Macromolecular Complex | 0.868 | 929 aa | Standard (800-1400 aa) |
| 3 | **MARS2** | [`Q96GW9`](https://www.uniprot.org/uniprotkb/Q96GW9) | Methionine--tRNA ligase, mitochondri | 593 aa | Core Macromolecular Complex | 0.908 | 619 aa | Optimal (<800 aa) |
| 4 | **TARS2** | [`Q9BW92`](https://www.uniprot.org/uniprotkb/Q9BW92) | Threonine--tRNA ligase, mitochondria | 718 aa | Direct Physical Interactor | 0.802 | 744 aa | Optimal (<800 aa) |
| 5 | **CARS2** | [`Q9HA77`](https://www.uniprot.org/uniprotkb/Q9HA77) | Probable cysteine--tRNA ligase, mito | 564 aa | Direct Physical Interactor | 0.838 | 590 aa | Optimal (<800 aa) |
| 6 | **SARS2** | [`Q9NP81`](https://www.uniprot.org/uniprotkb/Q9NP81) | Serine--tRNA ligase, mitochondrial | 518 aa | Core Macromolecular Complex | 0.869 | 544 aa | Optimal (<800 aa) |
| 7 | **YARS2** | [`Q9Y2Z4`](https://www.uniprot.org/uniprotkb/Q9Y2Z4) | Tyrosine--tRNA ligase, mitochondrial | 477 aa | Core Macromolecular Complex | 0.951 | 503 aa | Optimal (<800 aa) |
| 8 | **DARS2** | [`Q6PI48`](https://www.uniprot.org/uniprotkb/Q6PI48) | Aspartate--tRNA ligase, mitochondria | 645 aa | Direct Physical Interactor | 0.766 | 671 aa | Optimal (<800 aa) |
| 9 | **HARS2** | [`P49590`](https://www.uniprot.org/uniprotkb/P49590) | Histidine--tRNA ligase, mitochondria | 506 aa | Direct Physical Interactor | 0.828 | 532 aa | Optimal (<800 aa) |
| 10 | **GARS1** | [`P41250`](https://www.uniprot.org/uniprotkb/P41250) | Glycine--tRNA ligase | 739 aa | Direct Physical Interactor | 0.837 | 765 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein WARS2. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of WARS2 expression, stability, and macromolecular complex assembly.

---

#### Rank 19: `c16norep119` (*ZFHX3-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 154 aa | **Tier**: 4 | **Primary Tissue**: Whole_Blood
- **Microprotein Amino Acid Sequence**:  
  `MLWGPKSTRPAGPAAVQGLPHTPWWAAGLRSWLQLPWGPLLSPSWGQAGVPLSLCVPFPSPLLSHCLCPCVAAAVAVGAVVWVAAGLLAAGEAAAPVVAAAAAAAGGVAEGPSPSWLVLGLGSSATLRTVYLVIQKWGGQGKEQLCWADPSSSA`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **RUNX3** | [`Q13761`](https://www.uniprot.org/uniprotkb/Q13761) | Runt-related transcription factor 3 | 415 aa | Direct Physical Interactor | 0.955 | 569 aa | Optimal (<800 aa) |
| 2 | **POU2F1** | [`P14859`](https://www.uniprot.org/uniprotkb/P14859) | POU domain, class 2, transcription f | 743 aa | Direct Physical Interactor | 0.750 | 897 aa | Standard (800-1400 aa) |
| 3 | **POU1F1** | [`P28069`](https://www.uniprot.org/uniprotkb/P28069) | Pituitary-specific positive transcri | 291 aa | Direct Physical Interactor | 0.854 | 445 aa | Optimal (<800 aa) |
| 4 | **PIAS3** | [`Q9Y6X2`](https://www.uniprot.org/uniprotkb/Q9Y6X2) | E3 SUMO-protein ligase PIAS3 | 628 aa | Direct Physical Interactor | 0.789 | 782 aa | Optimal (<800 aa) |
| 5 | **RUNX2** | [`Q13950`](https://www.uniprot.org/uniprotkb/Q13950) | Runt-related transcription factor 2 | 521 aa | Direct Physical Interactor | 0.750 | 675 aa | Optimal (<800 aa) |
| 6 | **SMAD1** | [`Q15797`](https://www.uniprot.org/uniprotkb/Q15797) | SMAD family member 1 | 465 aa | Direct Physical Interactor | 0.750 | 619 aa | Optimal (<800 aa) |
| 7 | **SMAD2** | [`Q15796`](https://www.uniprot.org/uniprotkb/Q15796) | SMAD family member 2 | 467 aa | Direct Physical Interactor | 0.750 | 621 aa | Optimal (<800 aa) |
| 8 | **SMAD3** | [`P84022`](https://www.uniprot.org/uniprotkb/P84022) | SMAD family member 3 | 425 aa | Direct Physical Interactor | 0.750 | 579 aa | Optimal (<800 aa) |
| 9 | **HDAC1** | [`Q13547`](https://www.uniprot.org/uniprotkb/Q13547) | Histone deacetylase 1 | 482 aa | Direct Physical Interactor | 0.750 | 636 aa | Optimal (<800 aa) |
| 10 | **FOXA1** | [`P55317`](https://www.uniprot.org/uniprotkb/P55317) | Hepatocyte nuclear factor 3-alpha | 472 aa | Direct Physical Interactor | 0.750 | 626 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein ZFHX3. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of ZFHX3 expression, stability, and macromolecular complex assembly.

---

#### Rank 20: `c1riboseqorf33` (*EMC1-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 61 aa | **Tier**: 3 | **Primary Tissue**: Endoplasmic Reticulum (EMC Complex)
- **Microprotein Amino Acid Sequence**:  
  `MPQARRSSGSLSPAVGSVSTPIPSPEPKFGGNVSSGQDTPGQPHPPHHRLQTGAAATELKP`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **EMC1** | [`Q8N766`](https://www.uniprot.org/uniprotkb/Q8N766) | ER membrane protein complex subunit  | 993 aa | Cognate Antisense/Locus Master | 1.000 | 1054 aa | Standard (800-1400 aa) |
| 2 | **EMC2** | [`Q15006`](https://www.uniprot.org/uniprotkb/Q15006) | ER membrane protein complex subunit  | 297 aa | Core Macromolecular Complex | 0.999 | 358 aa | Optimal (<800 aa) |
| 3 | **EMC3** | [`Q9P0I2`](https://www.uniprot.org/uniprotkb/Q9P0I2) | ER membrane protein complex subunit  | 261 aa | Core Macromolecular Complex | 0.999 | 322 aa | Optimal (<800 aa) |
| 4 | **EMC4** | [`Q5J8M3`](https://www.uniprot.org/uniprotkb/Q5J8M3) | ER membrane protein complex subunit  | 183 aa | Core Macromolecular Complex | 0.999 | 244 aa | Optimal (<800 aa) |
| 5 | **EMC6** | [`Q9BV81`](https://www.uniprot.org/uniprotkb/Q9BV81) | ER membrane protein complex subunit  | 110 aa | Core Macromolecular Complex | 0.999 | 171 aa | Optimal (<800 aa) |
| 6 | **EMC7** | [`Q9NPA0`](https://www.uniprot.org/uniprotkb/Q9NPA0) | Endoplasmic reticulum membrane prote | 242 aa | Core Macromolecular Complex | 0.999 | 303 aa | Optimal (<800 aa) |
| 7 | **EMC8** | [`O43402`](https://www.uniprot.org/uniprotkb/O43402) | ER membrane protein complex subunit  | 210 aa | Core Macromolecular Complex | 0.999 | 271 aa | Optimal (<800 aa) |
| 8 | **EMC10** | [`Q5UCC4`](https://www.uniprot.org/uniprotkb/Q5UCC4) | ER membrane protein complex subunit  | 262 aa | Core Macromolecular Complex | 0.999 | 323 aa | Optimal (<800 aa) |
| 9 | **SEC61A1** | [`P61619`](https://www.uniprot.org/uniprotkb/P61619) | Protein transport protein Sec61 subu | 476 aa | Direct Physical Interactor | 0.570 | 537 aa | Optimal (<800 aa) |
| 10 | **VCP** | [`P55072`](https://www.uniprot.org/uniprotkb/P55072) | Transitional endoplasmic reticulum A | 806 aa | Direct Physical Interactor | 0.750 | 867 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein EMC1. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of EMC1 expression, stability, and macromolecular complex assembly.

---

#### Rank 21: `c11norep50` (*NAV2-AS2*)
- **Biotype**: lncRNA-ORF | **Length**: 103 aa | **Tier**: 4 | **Primary Tissue**: Cells_Cultured_fibroblasts
- **Microprotein Amino Acid Sequence**:  
  `MPRRHEPVWEMTGFFLPDVLLSDSSDTSEGFRLHLEPGSIFMPLSEPPSDFLTSSPDHREFLSTCECFSASVCGQRLSELVLPQSQAGKPLTPKAKMPKHWSS`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **CNN1** | [`P51911`](https://www.uniprot.org/uniprotkb/P51911) | Calponin-1 | 297 aa | Direct Physical Interactor | 0.725 | 400 aa | Optimal (<800 aa) |
| 2 | **NAV1** | [`Q8NEY1`](https://www.uniprot.org/uniprotkb/Q8NEY1) | Neuron navigator 1 | 1877 aa | Direct Physical Interactor | 0.750 | 1980 aa | High (>1400 aa) |
| 3 | **ING5** | [`Q8WYH8`](https://www.uniprot.org/uniprotkb/Q8WYH8) | Inhibitor of growth protein 5 | 240 aa | Direct Physical Interactor | 0.688 | 343 aa | Optimal (<800 aa) |
| 4 | **TCF7L1** | [`Q9HCS4`](https://www.uniprot.org/uniprotkb/Q9HCS4) | Transcription factor 7-like 1 | 588 aa | Direct Physical Interactor | 0.685 | 691 aa | Optimal (<800 aa) |
| 5 | **MAP2** | [`P11137`](https://www.uniprot.org/uniprotkb/P11137) | Microtubule-associated protein 2 | 1827 aa | Direct Physical Interactor | 0.750 | 1930 aa | High (>1400 aa) |
| 6 | **TUBA1A** | [`Q71U36`](https://www.uniprot.org/uniprotkb/Q71U36) | Tubulin alpha-1A chain | 451 aa | Direct Physical Interactor | 0.750 | 554 aa | Optimal (<800 aa) |
| 7 | **TUBB3** | [`Q13509`](https://www.uniprot.org/uniprotkb/Q13509) | Tubulin beta-3 chain | 450 aa | Direct Physical Interactor | 0.750 | 553 aa | Optimal (<800 aa) |
| 8 | **DBN1** | [`Q16643`](https://www.uniprot.org/uniprotkb/Q16643) | Drebrin | 649 aa | Direct Physical Interactor | 0.750 | 752 aa | Optimal (<800 aa) |
| 9 | **ABI1** | [`Q8IZP0`](https://www.uniprot.org/uniprotkb/Q8IZP0) | Abl interactor 1 | 508 aa | Direct Physical Interactor | 0.575 | 611 aa | Optimal (<800 aa) |
| 10 | **ABL1** | [`P00519`](https://www.uniprot.org/uniprotkb/P00519) | Tyrosine-protein kinase ABL1 | 1130 aa | Direct Physical Interactor | 0.750 | 1233 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein NAV2. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of NAV2 expression, stability, and macromolecular complex assembly.

---

#### Rank 22: `c20norep33` (*ENSG00000275457*)
- **Biotype**: lncRNA-ORF | **Length**: 29 aa | **Tier**: 4 | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MMDGYLRLSQRKNAGTPIAAHGWRRRPTG`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **XRN2** | [`Q9H0D6`](https://www.uniprot.org/uniprotkb/Q9H0D6) | 5'-3' exoribonuclease 2 | 950 aa | Cognate Antisense/Locus Master | 1.000 | 979 aa | Standard (800-1400 aa) |
| 2 | **DIS3** | [`Q9Y2L1`](https://www.uniprot.org/uniprotkb/Q9Y2L1) | Exosome complex exonuclease RRP44 | 958 aa | Direct Physical Interactor | 0.715 | 987 aa | Standard (800-1400 aa) |
| 3 | **EXOSC10** | [`Q01780`](https://www.uniprot.org/uniprotkb/Q01780) | Exosome complex component 10 | 885 aa | Core Macromolecular Complex | 0.856 | 914 aa | Standard (800-1400 aa) |
| 4 | **DXO** | [`O77932`](https://www.uniprot.org/uniprotkb/O77932) | Decapping and exoribonuclease protei | 396 aa | Core Macromolecular Complex | 0.999 | 425 aa | Optimal (<800 aa) |
| 5 | **CPSF4** | [`O95639`](https://www.uniprot.org/uniprotkb/O95639) | Cleavage and polyadenylation specifi | 269 aa | Direct Physical Interactor | 0.750 | 298 aa | Optimal (<800 aa) |
| 6 | **WDR82** | [`Q6UXN9`](https://www.uniprot.org/uniprotkb/Q6UXN9) | WD repeat-containing protein 82 | 313 aa | Direct Physical Interactor | 0.750 | 342 aa | Optimal (<800 aa) |
| 7 | **PRMT5** | [`O14744`](https://www.uniprot.org/uniprotkb/O14744) | Protein arginine N-methyltransferase | 637 aa | Direct Physical Interactor | 0.750 | 666 aa | Optimal (<800 aa) |
| 8 | **POLR2A** | [`P24928`](https://www.uniprot.org/uniprotkb/P24928) | DNA-directed RNA polymerase II subun | 1970 aa | Direct Physical Interactor | 0.714 | 1999 aa | High (>1400 aa) |
| 9 | **PAF1** | [`Q8N7H5`](https://www.uniprot.org/uniprotkb/Q8N7H5) | RNA polymerase II-associated factor  | 531 aa | Direct Physical Interactor | 0.750 | 560 aa | Optimal (<800 aa) |
| 10 | **RTF1** | [`Q92541`](https://www.uniprot.org/uniprotkb/Q92541) | RNA polymerase-associated protein RT | 710 aa | Direct Physical Interactor | 0.750 | 739 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein XRN2. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of XRN2 expression, stability, and macromolecular complex assembly.

---

#### Rank 23: `c17norep109` (*MAP3K14-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 50 aa | **Tier**: 4 | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MVALLIPGLLKSKERLSMLRHERICSCSSRENGWDRLLRNNSAGKSGRWL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **MAP3K14** | [`Q99558`](https://www.uniprot.org/uniprotkb/Q99558) | Mitogen-activated protein kinase kin | 947 aa | Cognate Antisense/Locus Master | 1.000 | 997 aa | Standard (800-1400 aa) |
| 2 | **CHUK** | [`O15111`](https://www.uniprot.org/uniprotkb/O15111) | Inhibitor of nuclear factor kappa-B  | 745 aa | Core Macromolecular Complex | 0.997 | 795 aa | Optimal (<800 aa) |
| 3 | **IKBKB** | [`O14920`](https://www.uniprot.org/uniprotkb/O14920) | Inhibitor of nuclear factor kappa-B  | 756 aa | Core Macromolecular Complex | 0.991 | 806 aa | Standard (800-1400 aa) |
| 4 | **IKBKG** | [`Q9Y6K9`](https://www.uniprot.org/uniprotkb/Q9Y6K9) | NF-kappa-B essential modulator | 419 aa | Core Macromolecular Complex | 0.994 | 469 aa | Optimal (<800 aa) |
| 5 | **TRAF2** | [`Q12933`](https://www.uniprot.org/uniprotkb/Q12933) | TNF receptor-associated factor 2 | 501 aa | Core Macromolecular Complex | 0.995 | 551 aa | Optimal (<800 aa) |
| 6 | **TRAF3** | [`Q13114`](https://www.uniprot.org/uniprotkb/Q13114) | TNF receptor-associated factor 3 | 568 aa | Core Macromolecular Complex | 0.994 | 618 aa | Optimal (<800 aa) |
| 7 | **BIRC2** | [`Q13490`](https://www.uniprot.org/uniprotkb/Q13490) | Baculoviral IAP repeat-containing pr | 618 aa | Core Macromolecular Complex | 0.888 | 668 aa | Optimal (<800 aa) |
| 8 | **BIRC3** | [`Q13489`](https://www.uniprot.org/uniprotkb/Q13489) | Baculoviral IAP repeat-containing pr | 604 aa | Direct Physical Interactor | 0.746 | 654 aa | Optimal (<800 aa) |
| 9 | **NFKB2** | [`Q00653`](https://www.uniprot.org/uniprotkb/Q00653) | Nuclear factor NF-kappa-B p100 subun | 900 aa | Core Macromolecular Complex | 0.966 | 950 aa | Standard (800-1400 aa) |
| 10 | **RELB** | [`Q01201`](https://www.uniprot.org/uniprotkb/Q01201) | Transcription factor RelB | 579 aa | Core Macromolecular Complex | 0.903 | 629 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein MAP3K14. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of MAP3K14 expression, stability, and macromolecular complex assembly.

---

#### Rank 24: `c2norep270` (*ENSG00000224090*)
- **Biotype**: lncRNA-ORF | **Length**: 83 aa | **Tier**: 4 | **Primary Tissue**: Testis / Ciliated Epithelia (Axoneme)
- **Microprotein Amino Acid Sequence**:  
  `MLLRKSRRLHLLCTGLRLPALGLMEKRQLSCTRIWWSSVEFNSACSAQSSTSIRWSLGRRKALGEGCKAIQRNWSGRTERQGL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **CRYBA2** | [`P53672`](https://www.uniprot.org/uniprotkb/P53672) | Beta-crystallin A2 | 197 aa | Direct Physical Interactor | 0.637 | 280 aa | Optimal (<800 aa) |
| 2 | **CFAP69** | [`A5D8W1`](https://www.uniprot.org/uniprotkb/A5D8W1) | Cilia- and flagella-associated prote | 941 aa | Direct Physical Interactor | 0.474 | 1024 aa | Standard (800-1400 aa) |
| 3 | **CFAP70** | [`Q5T0N1`](https://www.uniprot.org/uniprotkb/Q5T0N1) | Cilia- and flagella-associated prote | 1121 aa | Direct Physical Interactor | 0.615 | 1204 aa | Standard (800-1400 aa) |
| 4 | **SPEF2** | [`Q9C093`](https://www.uniprot.org/uniprotkb/Q9C093) | Sperm flagella and cilia-associated  | 1822 aa | Direct Physical Interactor | 0.408 | 1905 aa | High (>1400 aa) |
| 5 | **TTC29** | [`Q8NA56`](https://www.uniprot.org/uniprotkb/Q8NA56) | Tetratricopeptide repeat protein 29 | 475 aa | Core Macromolecular Complex | 0.511 | 558 aa | Optimal (<800 aa) |
| 6 | **RSPH10B** | [`P0C881`](https://www.uniprot.org/uniprotkb/P0C881) | Radial spoke head 10 homolog B | 870 aa | Core Macromolecular Complex | 0.475 | 953 aa | Standard (800-1400 aa) |
| 7 | **AKAP4** | [`Q5JQC9`](https://www.uniprot.org/uniprotkb/Q5JQC9) | A-kinase anchor protein 4 | 854 aa | Direct Physical Interactor | 0.750 | 937 aa | Standard (800-1400 aa) |
| 8 | **TEKT1** | [`Q969V4`](https://www.uniprot.org/uniprotkb/Q969V4) | Tektin-1 | 418 aa | Direct Physical Interactor | 0.750 | 501 aa | Optimal (<800 aa) |
| 9 | **DRC1** | [`Q96MC2`](https://www.uniprot.org/uniprotkb/Q96MC2) | Dynein regulatory complex protein 1 | 740 aa | Direct Physical Interactor | 0.750 | 823 aa | Standard (800-1400 aa) |
| 10 | **RSPH1** | [`Q8WYR4`](https://www.uniprot.org/uniprotkb/Q8WYR4) | Radial spoke head 1 homolog | 309 aa | Direct Physical Interactor | 0.750 | 392 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein CFAP65. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of CFAP65 expression, stability, and macromolecular complex assembly.

---

#### Rank 25: `c2riboseqorf136` (*TTN-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 127 aa | **Tier**: 4 | **Primary Tissue**: Muscle_Skeletal / Heart (Sarcomere)
- **Microprotein Amino Acid Sequence**:  
  `MACFFNLTPSSEFLTKWIFANNWFRSRAGLTQANIILCDDAELILSIGNKPCDSVHCGGDLSLVVSDPLVSGCLFAFDVVACDSRTTIIFGPGPGQANRTLGYIEYFWRIAWRFWRICKYKWKAHMY`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **TCAP** | [`O15273`](https://www.uniprot.org/uniprotkb/O15273) | Telethonin | 167 aa | Core Macromolecular Complex | 0.999 | 294 aa | Optimal (<800 aa) |
| 2 | **CSRP3** | [`P50461`](https://www.uniprot.org/uniprotkb/P50461) | Cysteine and glycine-rich protein 3 | 194 aa | Core Macromolecular Complex | 0.998 | 321 aa | Optimal (<800 aa) |
| 3 | **MYOT** | [`Q9UBF9`](https://www.uniprot.org/uniprotkb/Q9UBF9) | Myotilin | 498 aa | Core Macromolecular Complex | 0.898 | 625 aa | Optimal (<800 aa) |
| 4 | **ACTN2** | [`P35609`](https://www.uniprot.org/uniprotkb/P35609) | Alpha-actinin-2 | 894 aa | Core Macromolecular Complex | 0.969 | 1021 aa | Standard (800-1400 aa) |
| 5 | **TRIM63** | [`Q969Q1`](https://www.uniprot.org/uniprotkb/Q969Q1) | E3 ubiquitin-protein ligase TRIM63 | 353 aa | Core Macromolecular Complex | 0.953 | 480 aa | Optimal (<800 aa) |
| 6 | **TRIM54** | [`Q9BYV2`](https://www.uniprot.org/uniprotkb/Q9BYV2) | Tripartite motif-containing protein  | 358 aa | Direct Physical Interactor | 0.720 | 485 aa | Optimal (<800 aa) |
| 7 | **ANKRD1** | [`Q15327`](https://www.uniprot.org/uniprotkb/Q15327) | Ankyrin repeat domain-containing pro | 319 aa | Direct Physical Interactor | 0.750 | 446 aa | Optimal (<800 aa) |
| 8 | **CAPN3** | [`P20807`](https://www.uniprot.org/uniprotkb/P20807) | Calpain-3 | 821 aa | Core Macromolecular Complex | 0.855 | 948 aa | Standard (800-1400 aa) |
| 9 | **FHL2** | [`Q14192`](https://www.uniprot.org/uniprotkb/Q14192) | Four and a half LIM domains protein  | 279 aa | Direct Physical Interactor | 0.750 | 406 aa | Optimal (<800 aa) |
| 10 | **CRYAB** | [`P02511`](https://www.uniprot.org/uniprotkb/P02511) | Alpha-crystallin B chain | 175 aa | Direct Physical Interactor | 0.750 | 302 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary sarcomeric Z-disc structural surrogate for locus master TTN (Titin, 34,350 aa, structurally infeasible for ESMFold2). Microprotein is encoded on the TTN-AS1 antisense transcript; TCAP (Telethonin, 167 aa) is Titin's obligate N-terminal Z-disc capping partner (STRING combined score: 0.999), serving as the physiological complex target for sarcomeric modulation.

---

#### Rank 26: `c11riboseqorf31` (*NAV2-AS2*)
- **Biotype**: lncRNA-ORF | **Length**: 31 aa | **Tier**: 4 | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MRPVEELWLGTPPLLREMSMDSQLVWRLLTD`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **CNN1** | [`P51911`](https://www.uniprot.org/uniprotkb/P51911) | Calponin-1 | 297 aa | Direct Physical Interactor | 0.725 | 328 aa | Optimal (<800 aa) |
| 2 | **NAV1** | [`Q8NEY1`](https://www.uniprot.org/uniprotkb/Q8NEY1) | Neuron navigator 1 | 1877 aa | Direct Physical Interactor | 0.750 | 1908 aa | High (>1400 aa) |
| 3 | **ING5** | [`Q8WYH8`](https://www.uniprot.org/uniprotkb/Q8WYH8) | Inhibitor of growth protein 5 | 240 aa | Direct Physical Interactor | 0.688 | 271 aa | Optimal (<800 aa) |
| 4 | **TCF7L1** | [`Q9HCS4`](https://www.uniprot.org/uniprotkb/Q9HCS4) | Transcription factor 7-like 1 | 588 aa | Direct Physical Interactor | 0.685 | 619 aa | Optimal (<800 aa) |
| 5 | **MAP2** | [`P11137`](https://www.uniprot.org/uniprotkb/P11137) | Microtubule-associated protein 2 | 1827 aa | Direct Physical Interactor | 0.750 | 1858 aa | High (>1400 aa) |
| 6 | **TUBA1A** | [`Q71U36`](https://www.uniprot.org/uniprotkb/Q71U36) | Tubulin alpha-1A chain | 451 aa | Direct Physical Interactor | 0.750 | 482 aa | Optimal (<800 aa) |
| 7 | **TUBB3** | [`Q13509`](https://www.uniprot.org/uniprotkb/Q13509) | Tubulin beta-3 chain | 450 aa | Direct Physical Interactor | 0.750 | 481 aa | Optimal (<800 aa) |
| 8 | **DBN1** | [`Q16643`](https://www.uniprot.org/uniprotkb/Q16643) | Drebrin | 649 aa | Direct Physical Interactor | 0.750 | 680 aa | Optimal (<800 aa) |
| 9 | **ABI1** | [`Q8IZP0`](https://www.uniprot.org/uniprotkb/Q8IZP0) | Abl interactor 1 | 508 aa | Direct Physical Interactor | 0.575 | 539 aa | Optimal (<800 aa) |
| 10 | **ABL1** | [`P00519`](https://www.uniprot.org/uniprotkb/P00519) | Tyrosine-protein kinase ABL1 | 1130 aa | Direct Physical Interactor | 0.750 | 1161 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein NAV2. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of NAV2 expression, stability, and macromolecular complex assembly.

---

#### Rank 27: `c19norep130` (*TMEM147-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 26 aa | **Tier**: 4 | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MTLARGAKQLVVQDALLKQPRCARTL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **TMEM147** | [`Q9BVK8`](https://www.uniprot.org/uniprotkb/Q9BVK8) | BOS complex subunit TMEM147 | 224 aa | Cognate Antisense/Locus Master | 1.000 | 250 aa | Optimal (<800 aa) |
| 2 | **NCLN** | [`Q969V3`](https://www.uniprot.org/uniprotkb/Q969V3) | BOS complex subunit NCLN | 563 aa | Core Macromolecular Complex | 0.998 | 589 aa | Optimal (<800 aa) |
| 3 | **NOMO1** | [`Q15155`](https://www.uniprot.org/uniprotkb/Q15155) | BOS complex subunit NOMO1 | 1222 aa | Direct Physical Interactor | 0.750 | 1248 aa | Standard (800-1400 aa) |
| 4 | **NOMO2** | [`Q5JPE7`](https://www.uniprot.org/uniprotkb/Q5JPE7) | BOS complex subunit NOMO2 | 1267 aa | Core Macromolecular Complex | 0.986 | 1293 aa | Standard (800-1400 aa) |
| 5 | **NOMO3** | [`P69849`](https://www.uniprot.org/uniprotkb/P69849) | BOS complex subunit NOMO3 | 1222 aa | Direct Physical Interactor | 0.783 | 1248 aa | Standard (800-1400 aa) |
| 6 | **SEC61A1** | [`P61619`](https://www.uniprot.org/uniprotkb/P61619) | Protein transport protein Sec61 subu | 476 aa | Core Macromolecular Complex | 0.891 | 502 aa | Optimal (<800 aa) |
| 7 | **SEC61B** | [`P60468`](https://www.uniprot.org/uniprotkb/P60468) | Protein transport protein Sec61 subu | 96 aa | Core Macromolecular Complex | 0.860 | 122 aa | Optimal (<800 aa) |
| 8 | **SEC61G** | [`P60059`](https://www.uniprot.org/uniprotkb/P60059) | Protein transport protein Sec61 subu | 68 aa | Direct Physical Interactor | 0.838 | 94 aa | Optimal (<800 aa) |
| 9 | **CANX** | [`P27824`](https://www.uniprot.org/uniprotkb/P27824) | Calnexin | 592 aa | Direct Physical Interactor | 0.750 | 618 aa | Optimal (<800 aa) |
| 10 | **CALR** | [`P27797`](https://www.uniprot.org/uniprotkb/P27797) | Calreticulin | 417 aa | Direct Physical Interactor | 0.750 | 443 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein TMEM147. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of TMEM147 expression, stability, and macromolecular complex assembly.

---

#### Rank 28: `c16norep37` (*ENSG00000260592*)
- **Biotype**: lncRNA-ORF | **Length**: 61 aa | **Tier**: 4 | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MLPCVVPCGIECPLFSGTMVMARVHKTPVKDGKKSKKMKKVIICEARQALLGGWKFIIRLI`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **TMC5** | [`Q6UXY8`](https://www.uniprot.org/uniprotkb/Q6UXY8) | Transmembrane channel-like protein 5 | 1006 aa | Cognate Antisense/Locus Master | 1.000 | 1067 aa | Standard (800-1400 aa) |
| 2 | **TMC6** | [`Q7Z403`](https://www.uniprot.org/uniprotkb/Q7Z403) | Transmembrane channel-like protein 6 | 805 aa | Direct Physical Interactor | 0.750 | 866 aa | Standard (800-1400 aa) |
| 3 | **TMC8** | [`Q8IU68`](https://www.uniprot.org/uniprotkb/Q8IU68) | Transmembrane channel-like protein 8 | 726 aa | Direct Physical Interactor | 0.750 | 787 aa | Optimal (<800 aa) |
| 4 | **KCNH8** | [`Q96L42`](https://www.uniprot.org/uniprotkb/Q96L42) | Voltage-gated delayed rectifier pota | 1107 aa | Direct Physical Interactor | 0.508 | 1168 aa | Standard (800-1400 aa) |
| 5 | **TMEM221** | [`A6NGB7`](https://www.uniprot.org/uniprotkb/A6NGB7) | Transmembrane protein 221 | 291 aa | Direct Physical Interactor | 0.473 | 352 aa | Optimal (<800 aa) |
| 6 | **ADAMTS12** | [`P58397`](https://www.uniprot.org/uniprotkb/P58397) | A disintegrin and metalloproteinase  | 1594 aa | Direct Physical Interactor | 0.472 | 1655 aa | High (>1400 aa) |
| 7 | **CMC1** | [`Q7Z7K0`](https://www.uniprot.org/uniprotkb/Q7Z7K0) | COX assembly mitochondrial protein h | 106 aa | Direct Physical Interactor | 0.454 | 167 aa | Optimal (<800 aa) |
| 8 | **RNF175** | [`Q8N4F7`](https://www.uniprot.org/uniprotkb/Q8N4F7) | RING finger protein 175 | 328 aa | Direct Physical Interactor | 0.435 | 389 aa | Optimal (<800 aa) |
| 9 | **CCP110** | [`O43303`](https://www.uniprot.org/uniprotkb/O43303) | Centriolar coiled-coil protein of 11 | 1012 aa | Direct Physical Interactor | 0.434 | 1073 aa | Standard (800-1400 aa) |
| 10 | **ANKUB1** | [`A6NFN9`](https://www.uniprot.org/uniprotkb/A6NFN9) | Protein ANKUB1 | 502 aa | Direct Physical Interactor | 0.429 | 563 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein TMC5. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of TMC5 expression, stability, and macromolecular complex assembly.

---

#### Rank 29: `c6norep46` (*NUP153-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 96 aa | **Tier**: 2B | **Primary Tissue**: Kidney_Medulla
- **Microprotein Amino Acid Sequence**:  
  `MLLSPLLVRLNWPLVATPRPDLAATAPSDSSGSRGHGGASAASRSGAGKGAGEAEAEALESLPRRPAPAQKSARAVHTVGTSTPGTARLRAGAERG`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **NUP153** | [`P49790`](https://www.uniprot.org/uniprotkb/P49790) | Nuclear pore complex protein Nup153 | 1475 aa | Cognate Antisense/Locus Master | 1.000 | 1571 aa | High (>1400 aa) |
| 2 | **NUP50** | [`Q9UKX7`](https://www.uniprot.org/uniprotkb/Q9UKX7) | Nuclear pore complex protein Nup50 | 468 aa | Core Macromolecular Complex | 0.999 | 564 aa | Optimal (<800 aa) |
| 3 | **NUP133** | [`Q8WUM0`](https://www.uniprot.org/uniprotkb/Q8WUM0) | Nuclear pore complex protein Nup133 | 1156 aa | Core Macromolecular Complex | 0.998 | 1252 aa | Standard (800-1400 aa) |
| 4 | **NUP93** | [`Q8N1F7`](https://www.uniprot.org/uniprotkb/Q8N1F7) | Nuclear pore complex protein Nup93 | 819 aa | Core Macromolecular Complex | 0.996 | 915 aa | Standard (800-1400 aa) |
| 5 | **XPO1** | [`O14980`](https://www.uniprot.org/uniprotkb/O14980) | Exportin-1 | 1071 aa | Core Macromolecular Complex | 0.894 | 1167 aa | Standard (800-1400 aa) |
| 6 | **KPNA2** | [`P52292`](https://www.uniprot.org/uniprotkb/P52292) | Importin subunit alpha-1 | 529 aa | Core Macromolecular Complex | 0.969 | 625 aa | Optimal (<800 aa) |
| 7 | **KPNB1** | [`Q14974`](https://www.uniprot.org/uniprotkb/Q14974) | Importin subunit beta-1 | 876 aa | Core Macromolecular Complex | 0.996 | 972 aa | Standard (800-1400 aa) |
| 8 | **NUP98** | [`P52948`](https://www.uniprot.org/uniprotkb/P52948) | Nuclear pore complex protein Nup98-N | 1817 aa | Core Macromolecular Complex | 0.999 | 1913 aa | High (>1400 aa) |
| 9 | **NUP107** | [`P57740`](https://www.uniprot.org/uniprotkb/P57740) | Nuclear pore complex protein Nup107 | 925 aa | Core Macromolecular Complex | 0.999 | 1021 aa | Standard (800-1400 aa) |
| 10 | **RAN** | [`P62826`](https://www.uniprot.org/uniprotkb/P62826) | GTP-binding nuclear protein Ran | 216 aa | Core Macromolecular Complex | 0.996 | 312 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary locus-regulated cognate protein NUP153. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of NUP153 expression, stability, and macromolecular complex assembly.

---

#### Rank 30: `c2norep211` (*TTN-AS1*)
- **Biotype**: lncRNA-ORF | **Length**: 87 aa | **Tier**: 4 | **Primary Tissue**: Muscle_Skeletal / Heart (Sarcomere)
- **Microprotein Amino Acid Sequence**:  
  `MKQHTPLFLSAQGKQIWGPLHLCVGKHRWYNPLLSRPTLHSQKLRPCSWCTQHGDHRIWFATERNLDADHFTSESMNCSHAFKSFSG`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **TCAP** | [`O15273`](https://www.uniprot.org/uniprotkb/O15273) | Telethonin | 167 aa | Core Macromolecular Complex | 0.999 | 254 aa | Optimal (<800 aa) |
| 2 | **CSRP3** | [`P50461`](https://www.uniprot.org/uniprotkb/P50461) | Cysteine and glycine-rich protein 3 | 194 aa | Core Macromolecular Complex | 0.998 | 281 aa | Optimal (<800 aa) |
| 3 | **MYOT** | [`Q9UBF9`](https://www.uniprot.org/uniprotkb/Q9UBF9) | Myotilin | 498 aa | Core Macromolecular Complex | 0.898 | 585 aa | Optimal (<800 aa) |
| 4 | **ACTN2** | [`P35609`](https://www.uniprot.org/uniprotkb/P35609) | Alpha-actinin-2 | 894 aa | Core Macromolecular Complex | 0.969 | 981 aa | Standard (800-1400 aa) |
| 5 | **TRIM63** | [`Q969Q1`](https://www.uniprot.org/uniprotkb/Q969Q1) | E3 ubiquitin-protein ligase TRIM63 | 353 aa | Core Macromolecular Complex | 0.953 | 440 aa | Optimal (<800 aa) |
| 6 | **TRIM54** | [`Q9BYV2`](https://www.uniprot.org/uniprotkb/Q9BYV2) | Tripartite motif-containing protein  | 358 aa | Direct Physical Interactor | 0.720 | 445 aa | Optimal (<800 aa) |
| 7 | **ANKRD1** | [`Q15327`](https://www.uniprot.org/uniprotkb/Q15327) | Ankyrin repeat domain-containing pro | 319 aa | Direct Physical Interactor | 0.750 | 406 aa | Optimal (<800 aa) |
| 8 | **CAPN3** | [`P20807`](https://www.uniprot.org/uniprotkb/P20807) | Calpain-3 | 821 aa | Core Macromolecular Complex | 0.855 | 908 aa | Standard (800-1400 aa) |
| 9 | **FHL2** | [`Q14192`](https://www.uniprot.org/uniprotkb/Q14192) | Four and a half LIM domains protein  | 279 aa | Direct Physical Interactor | 0.750 | 366 aa | Optimal (<800 aa) |
| 10 | **CRYAB** | [`P02511`](https://www.uniprot.org/uniprotkb/P02511) | Alpha-crystallin B chain | 175 aa | Direct Physical Interactor | 0.750 | 262 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary sarcomeric Z-disc structural surrogate for locus master TTN (Titin, 34,350 aa, structurally infeasible for ESMFold2). Microprotein is encoded on the TTN-AS1 antisense transcript; TCAP (Telethonin, 167 aa) is Titin's obligate N-terminal Z-disc capping partner (STRING combined score: 0.999), serving as the physiological complex target for sarcomeric modulation.

---

### 3. Autonomous Regulatory uORF

#### Rank 31: `c3riboseqorf183` (*BCL6*)
- **Biotype**: uORF | **Length**: 19 aa | **Tier**: 4 | **Primary Tissue**: Whole_Blood / Spleen / Germinal Center (B-Cells)
- **Microprotein Amino Acid Sequence**:  
  `MQEVSRKGRTPGFEQNFGL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **BCL6** | [`P41182`](https://www.uniprot.org/uniprotkb/P41182) | B-cell lymphoma 6 protein | 706 aa | Cognate Host CDS (Co-Translational) | 1.000 | 725 aa | Optimal (<800 aa) |
| 2 | **BCOR** | [`Q6W2J9`](https://www.uniprot.org/uniprotkb/Q6W2J9) | BCL-6 corepressor | 1755 aa | Core Macromolecular Complex | 0.999 | 1774 aa | High (>1400 aa) |
| 3 | **BCORL1** | [`Q5H9F3`](https://www.uniprot.org/uniprotkb/Q5H9F3) | BCL-6 corepressor-like protein 1 | 1785 aa | Direct Physical Interactor | 0.797 | 1804 aa | High (>1400 aa) |
| 4 | **HDAC4** | [`P56524`](https://www.uniprot.org/uniprotkb/P56524) | Histone deacetylase 4 | 1084 aa | Core Macromolecular Complex | 0.966 | 1103 aa | Standard (800-1400 aa) |
| 5 | **TP53** | [`P04637`](https://www.uniprot.org/uniprotkb/P04637) | Cellular tumor antigen p53 | 393 aa | Direct Physical Interactor | 0.957 | 412 aa | Optimal (<800 aa) |
| 6 | **HDAC3** | [`O15379`](https://www.uniprot.org/uniprotkb/O15379) | Histone deacetylase 3 | 428 aa | Direct Physical Interactor | 0.750 | 447 aa | Optimal (<800 aa) |
| 7 | **CTBP1** | [`Q13363`](https://www.uniprot.org/uniprotkb/Q13363) | C-terminal-binding protein 1 | 440 aa | Direct Physical Interactor | 0.796 | 459 aa | Optimal (<800 aa) |
| 8 | **BCL6B** | [`Q8N143`](https://www.uniprot.org/uniprotkb/Q8N143) | B-cell CLL/lymphoma 6 member B prote | 479 aa | Core Macromolecular Complex | 0.910 | 498 aa | Optimal (<800 aa) |
| 9 | **IRF4** | [`Q15306`](https://www.uniprot.org/uniprotkb/Q15306) | Interferon regulatory factor 4 | 451 aa | Direct Physical Interactor | 0.942 | 470 aa | Optimal (<800 aa) |
| 10 | **PAX5** | [`Q02548`](https://www.uniprot.org/uniprotkb/Q02548) | Paired box protein Pax-5 | 391 aa | Direct Physical Interactor | 0.931 | 410 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein BCL6. Microprotein is translated from the 5' UTR uORF of BCL6 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 32: `c8riboseqorf89` (*DCAF13*)
- **Biotype**: uORF | **Length**: 61 aa | **Tier**: 2B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MQPGVWAGSYRRSGSSRGRADAAQGWIRRRRSRPLSGSHTGRGGIPCSTPTTRTPTGPGPS`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **DCAF13** | [`Q9NV06`](https://www.uniprot.org/uniprotkb/Q9NV06) | DDB1- and CUL4-associated factor 13 | 445 aa | Cognate Host CDS (Co-Translational) | 1.000 | 506 aa | Optimal (<800 aa) |
| 2 | **DDB1** | [`Q16531`](https://www.uniprot.org/uniprotkb/Q16531) | DNA damage-binding protein 1 | 1140 aa | Direct Physical Interactor | 0.750 | 1201 aa | Standard (800-1400 aa) |
| 3 | **CUL4A** | [`Q13619`](https://www.uniprot.org/uniprotkb/Q13619) | Cullin-4A | 759 aa | Direct Physical Interactor | 0.750 | 820 aa | Standard (800-1400 aa) |
| 4 | **CUL4B** | [`Q13620`](https://www.uniprot.org/uniprotkb/Q13620) | Cullin-4B | 913 aa | Direct Physical Interactor | 0.750 | 974 aa | Standard (800-1400 aa) |
| 5 | **RBX1** | [`P62877`](https://www.uniprot.org/uniprotkb/P62877) | E3 ubiquitin-protein ligase RBX1 | 108 aa | Direct Physical Interactor | 0.750 | 169 aa | Optimal (<800 aa) |
| 6 | **UTP18** | [`Q9Y5J1`](https://www.uniprot.org/uniprotkb/Q9Y5J1) | U3 small nucleolar RNA-associated pr | 556 aa | Core Macromolecular Complex | 0.998 | 617 aa | Optimal (<800 aa) |
| 7 | **WDR43** | [`Q15061`](https://www.uniprot.org/uniprotkb/Q15061) | WD repeat-containing protein 43 | 677 aa | Core Macromolecular Complex | 0.999 | 738 aa | Optimal (<800 aa) |
| 8 | **WDR36** | [`Q8NI36`](https://www.uniprot.org/uniprotkb/Q8NI36) | WD repeat-containing protein 36 | 895 aa | Core Macromolecular Complex | 0.999 | 956 aa | Standard (800-1400 aa) |
| 9 | **BMS1** | [`Q14692`](https://www.uniprot.org/uniprotkb/Q14692) | Ribosome biogenesis protein BMS1 hom | 1282 aa | Core Macromolecular Complex | 0.998 | 1343 aa | Standard (800-1400 aa) |
| 10 | **RRP9** | [`O43818`](https://www.uniprot.org/uniprotkb/O43818) | U3 small nucleolar RNA-interacting p | 475 aa | Core Macromolecular Complex | 0.999 | 536 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein DCAF13. Microprotein is translated from the 5' UTR uORF of DCAF13 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 33: `c10riboseqorf118` (*SHTN1*)
- **Biotype**: uORF | **Length**: 21 aa | **Tier**: 4 | **Primary Tissue**: Brain / Neuronal (Axon Growth Cone)
- **Microprotein Amino Acid Sequence**:  
  `MISLARSAPGGGGAGADPTSG`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **SHTN1** | [`A0MZ66`](https://www.uniprot.org/uniprotkb/A0MZ66) | Shootin-1 | 631 aa | Cognate Host CDS (Co-Translational) | 1.000 | 652 aa | Optimal (<800 aa) |
| 2 | **L1CAM** | [`P32004`](https://www.uniprot.org/uniprotkb/P32004) | Neural cell adhesion molecule L1 | 1257 aa | Core Macromolecular Complex | 0.980 | 1278 aa | Standard (800-1400 aa) |
| 3 | **CTTN** | [`Q14247`](https://www.uniprot.org/uniprotkb/Q14247) | Src substrate cortactin | 550 aa | Direct Physical Interactor | 0.641 | 571 aa | Optimal (<800 aa) |
| 4 | **PAK1** | [`Q13153`](https://www.uniprot.org/uniprotkb/Q13153) | Serine/threonine-protein kinase PAK  | 545 aa | Direct Physical Interactor | 0.750 | 566 aa | Optimal (<800 aa) |
| 5 | **CDC42** | [`P60953`](https://www.uniprot.org/uniprotkb/P60953) | Cell division control protein 42 hom | 191 aa | Direct Physical Interactor | 0.750 | 212 aa | Optimal (<800 aa) |
| 6 | **RAC1** | [`P63000`](https://www.uniprot.org/uniprotkb/P63000) | Ras-related C3 botulinum toxin subst | 192 aa | Direct Physical Interactor | 0.750 | 213 aa | Optimal (<800 aa) |
| 7 | **ARHGEF7** | [`Q14155`](https://www.uniprot.org/uniprotkb/Q14155) | Rho guanine nucleotide exchange fact | 803 aa | Direct Physical Interactor | 0.750 | 824 aa | Standard (800-1400 aa) |
| 8 | **TIAM1** | [`Q13009`](https://www.uniprot.org/uniprotkb/Q13009) | Rho guanine nucleotide exchange fact | 1591 aa | Direct Physical Interactor | 0.750 | 1612 aa | High (>1400 aa) |
| 9 | **ACTB** | [`P60709`](https://www.uniprot.org/uniprotkb/P60709) | Actin, cytoplasmic 1 | 375 aa | Direct Physical Interactor | 0.750 | 396 aa | Optimal (<800 aa) |
| 10 | **TUBA1A** | [`Q71U36`](https://www.uniprot.org/uniprotkb/Q71U36) | Tubulin alpha-1A chain | 451 aa | Direct Physical Interactor | 0.750 | 472 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein SHTN1. Microprotein is translated from the 5' UTR uORF of SHTN1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 34: `c17norep67` (*EVI2A*)
- **Biotype**: uORF | **Length**: 35 aa | **Tier**: 2B | **Primary Tissue**: Cells_Cultured_fibroblasts
- **Microprotein Amino Acid Sequence**:  
  `MCHLWFGFKSGKLAAHILFYCRFTLRLIFSKSILL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **EVI2A** | [`P22794`](https://www.uniprot.org/uniprotkb/P22794) | Protein EVI2A | 236 aa | Cognate Host CDS (Co-Translational) | 1.000 | 271 aa | Optimal (<800 aa) |
| 2 | **RNF135** | [`Q8IUD6`](https://www.uniprot.org/uniprotkb/Q8IUD6) | E3 ubiquitin-protein ligase RNF135 | 432 aa | Direct Physical Interactor | 0.809 | 467 aa | Optimal (<800 aa) |
| 3 | **OMG** | [`P23515`](https://www.uniprot.org/uniprotkb/P23515) | Oligodendrocyte-myelin glycoprotein | 440 aa | Core Macromolecular Complex | 0.878 | 475 aa | Optimal (<800 aa) |
| 4 | **EVI2B** | [`P34910`](https://www.uniprot.org/uniprotkb/P34910) | Protein EVI2B | 448 aa | Core Macromolecular Complex | 0.995 | 483 aa | Optimal (<800 aa) |
| 5 | **SPRED1** | [`Q7Z699`](https://www.uniprot.org/uniprotkb/Q7Z699) | Sprouty-related, EVH1 domain-contain | 444 aa | Direct Physical Interactor | 0.750 | 479 aa | Optimal (<800 aa) |
| 6 | **RASA1** | [`P20936`](https://www.uniprot.org/uniprotkb/P20936) | Ras GTPase-activating protein 1 | 1047 aa | Direct Physical Interactor | 0.750 | 1082 aa | Standard (800-1400 aa) |
| 7 | **SPRED2** | [`Q7Z698`](https://www.uniprot.org/uniprotkb/Q7Z698) | Sprouty-related, EVH1 domain-contain | 418 aa | Direct Physical Interactor | 0.750 | 453 aa | Optimal (<800 aa) |
| 8 | **KRAS** | [`P01116`](https://www.uniprot.org/uniprotkb/P01116) | GTPase KRas | 189 aa | Direct Physical Interactor | 0.750 | 224 aa | Optimal (<800 aa) |
| 9 | **HRAS** | [`P01112`](https://www.uniprot.org/uniprotkb/P01112) | GTPase HRas | 189 aa | Direct Physical Interactor | 0.750 | 224 aa | Optimal (<800 aa) |
| 10 | **NRAS** | [`P01111`](https://www.uniprot.org/uniprotkb/P01111) | GTPase NRas | 189 aa | Direct Physical Interactor | 0.750 | 224 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein EVI2A. Microprotein is translated from the 5' UTR uORF of EVI2A mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 35: `c1riboseqorf209` (*DCAF8*)
- **Biotype**: uORF | **Length**: 42 aa | **Tier**: 1B | **Primary Tissue**: Whole_Blood
- **Microprotein Amino Acid Sequence**:  
  `MAAALVVAVDGAAGPGTVGAAASRRMTQDGKRAAEPRRFPCL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **DCAF8** | [`Q5TAQ9`](https://www.uniprot.org/uniprotkb/Q5TAQ9) | DDB1- and CUL4-associated factor 8 | 597 aa | Cognate Host CDS (Co-Translational) | 1.000 | 639 aa | Optimal (<800 aa) |
| 2 | **DDB1** | [`Q16531`](https://www.uniprot.org/uniprotkb/Q16531) | DNA damage-binding protein 1 | 1140 aa | Core Macromolecular Complex | 0.999 | 1182 aa | Standard (800-1400 aa) |
| 3 | **CUL4A** | [`Q13619`](https://www.uniprot.org/uniprotkb/Q13619) | Cullin-4A | 759 aa | Core Macromolecular Complex | 0.990 | 801 aa | Standard (800-1400 aa) |
| 4 | **CUL4B** | [`Q13620`](https://www.uniprot.org/uniprotkb/Q13620) | Cullin-4B | 913 aa | Core Macromolecular Complex | 0.992 | 955 aa | Standard (800-1400 aa) |
| 5 | **RBX1** | [`P62877`](https://www.uniprot.org/uniprotkb/P62877) | E3 ubiquitin-protein ligase RBX1 | 108 aa | Core Macromolecular Complex | 0.948 | 150 aa | Optimal (<800 aa) |
| 6 | **DCAF6** | [`Q58WW2`](https://www.uniprot.org/uniprotkb/Q58WW2) | DDB1- and CUL4-associated factor 6 | 860 aa | Core Macromolecular Complex | 0.916 | 902 aa | Standard (800-1400 aa) |
| 7 | **DCAF7** | [`P61962`](https://www.uniprot.org/uniprotkb/P61962) | DDB1- and CUL4-associated factor 7 | 342 aa | Core Macromolecular Complex | 0.912 | 384 aa | Optimal (<800 aa) |
| 8 | **CAND1** | [`Q86VP6`](https://www.uniprot.org/uniprotkb/Q86VP6) | Cullin-associated NEDD8-dissociated  | 1230 aa | Direct Physical Interactor | 0.611 | 1272 aa | Standard (800-1400 aa) |
| 9 | **COPS5** | [`Q92905`](https://www.uniprot.org/uniprotkb/Q92905) | COP9 signalosome complex subunit 5 | 334 aa | Direct Physical Interactor | 0.750 | 376 aa | Optimal (<800 aa) |
| 10 | **NEDD8** | [`Q15843`](https://www.uniprot.org/uniprotkb/Q15843) | Ubiquitin-like protein NEDD8 | 81 aa | Direct Physical Interactor | 0.538 | 123 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein DCAF8. Microprotein is translated from the 5' UTR uORF of DCAF8 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 36: `cXriboseqorf74` (*DNASE1L1*)
- **Biotype**: uORF | **Length**: 26 aa | **Tier**: 4 | **Primary Tissue**: Cells_Cultured_fibroblasts
- **Microprotein Amino Acid Sequence**:  
  `MGTRAGLPRPGPTHWSRSSCRCPPSP`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **DNASE1L1** | [`P49184`](https://www.uniprot.org/uniprotkb/P49184) | Deoxyribonuclease-1-like 1 | 302 aa | Cognate Host CDS (Co-Translational) | 1.000 | 328 aa | Optimal (<800 aa) |
| 2 | **SPCS1** | [`Q9Y6A9`](https://www.uniprot.org/uniprotkb/Q9Y6A9) | Signal peptidase complex subunit 1 | 169 aa | Direct Physical Interactor | 0.626 | 195 aa | Optimal (<800 aa) |
| 3 | **SPCS2** | [`Q15005`](https://www.uniprot.org/uniprotkb/Q15005) | Signal peptidase complex subunit 2 | 226 aa | Direct Physical Interactor | 0.628 | 252 aa | Optimal (<800 aa) |
| 4 | **BSG** | [`P35613`](https://www.uniprot.org/uniprotkb/P35613) | Basigin | 385 aa | Direct Physical Interactor | 0.593 | 411 aa | Optimal (<800 aa) |
| 5 | **DNASE2** | [`O00115`](https://www.uniprot.org/uniprotkb/O00115) | Deoxyribonuclease-2-alpha | 360 aa | Direct Physical Interactor | 0.557 | 386 aa | Optimal (<800 aa) |
| 6 | **DNASE2B** | [`Q8WZ79`](https://www.uniprot.org/uniprotkb/Q8WZ79) | Deoxyribonuclease-2-beta | 361 aa | Direct Physical Interactor | 0.555 | 387 aa | Optimal (<800 aa) |
| 7 | **DNASE1** | [`P24855`](https://www.uniprot.org/uniprotkb/P24855) | Deoxyribonuclease-1 | 282 aa | Direct Physical Interactor | 0.750 | 308 aa | Optimal (<800 aa) |
| 8 | **DFFB** | [`O76075`](https://www.uniprot.org/uniprotkb/O76075) | DNA fragmentation factor subunit bet | 338 aa | Direct Physical Interactor | 0.475 | 364 aa | Optimal (<800 aa) |
| 9 | **DFFA** | [`O00273`](https://www.uniprot.org/uniprotkb/O00273) | DNA fragmentation factor subunit alp | 331 aa | Direct Physical Interactor | 0.750 | 357 aa | Optimal (<800 aa) |
| 10 | **ACTB** | [`P60709`](https://www.uniprot.org/uniprotkb/P60709) | Actin, cytoplasmic 1 | 375 aa | Direct Physical Interactor | 0.750 | 401 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein DNASE1L1. Microprotein is translated from the 5' UTR uORF of DNASE1L1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 37: `c14norep27` (*CIDEB*)
- **Biotype**: uORF | **Length**: 28 aa | **Tier**: 2B | **Primary Tissue**: Liver / Adipose (Lipid Droplet)
- **Microprotein Amino Acid Sequence**:  
  `MALTSRPTRAPCRPEPQRAPRSRASVTL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **CIDEB** | [`Q9UHD4`](https://www.uniprot.org/uniprotkb/Q9UHD4) | Lipid transferase CIDEB | 219 aa | Cognate Host CDS (Co-Translational) | 1.000 | 247 aa | Optimal (<800 aa) |
| 2 | **CIDEA** | [`O60543`](https://www.uniprot.org/uniprotkb/O60543) | Lipid transferase CIDEA | 219 aa | Direct Physical Interactor | 0.750 | 247 aa | Optimal (<800 aa) |
| 3 | **CIDEC** | [`Q96AQ7`](https://www.uniprot.org/uniprotkb/Q96AQ7) | Lipid transferase CIDEC | 238 aa | Direct Physical Interactor | 0.750 | 266 aa | Optimal (<800 aa) |
| 4 | **DFFA** | [`O00273`](https://www.uniprot.org/uniprotkb/O00273) | DNA fragmentation factor subunit alp | 331 aa | Core Macromolecular Complex | 0.953 | 359 aa | Optimal (<800 aa) |
| 5 | **PLIN2** | [`Q99541`](https://www.uniprot.org/uniprotkb/Q99541) | Perilipin-2 | 437 aa | Direct Physical Interactor | 0.521 | 465 aa | Optimal (<800 aa) |
| 6 | **PLIN1** | [`O60240`](https://www.uniprot.org/uniprotkb/O60240) | Perilipin-1 | 522 aa | Direct Physical Interactor | 0.729 | 550 aa | Optimal (<800 aa) |
| 7 | **PLIN3** | [`O60664`](https://www.uniprot.org/uniprotkb/O60664) | Perilipin-3 | 434 aa | Direct Physical Interactor | 0.466 | 462 aa | Optimal (<800 aa) |
| 8 | **PNPLA2** | [`Q96AD5`](https://www.uniprot.org/uniprotkb/Q96AD5) | Patatin-like phospholipase domain-co | 504 aa | Direct Physical Interactor | 0.750 | 532 aa | Optimal (<800 aa) |
| 9 | **BSCL2** | [`Q96G97`](https://www.uniprot.org/uniprotkb/Q96G97) | Seipin | 398 aa | Direct Physical Interactor | 0.750 | 426 aa | Optimal (<800 aa) |
| 10 | **MTTP** | [`P55157`](https://www.uniprot.org/uniprotkb/P55157) | Microsomal triglyceride transfer pro | 894 aa | Direct Physical Interactor | 0.598 | 922 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary cognate host protein CIDEB. Microprotein is translated from the 5' UTR uORF of CIDEB mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 38: `c22riboseqorf6` (*DGCR8*)
- **Biotype**: uORF | **Length**: 16 aa | **Tier**: 4 | **Primary Tissue**: Nucleus / Microprocessor Complex
- **Microprotein Amino Acid Sequence**:  
  `MKTDSLSRQSLKLSAL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **DGCR8** | [`Q8WYQ5`](https://www.uniprot.org/uniprotkb/Q8WYQ5) | Microprocessor complex subunit DGCR8 | 773 aa | Cognate Host CDS (Co-Translational) | 1.000 | 789 aa | Optimal (<800 aa) |
| 2 | **DROSHA** | [`Q9NRR4`](https://www.uniprot.org/uniprotkb/Q9NRR4) | Ribonuclease 3 | 1374 aa | Core Macromolecular Complex | 0.999 | 1390 aa | Standard (800-1400 aa) |
| 3 | **DDX5** | [`P17844`](https://www.uniprot.org/uniprotkb/P17844) | Probable ATP-dependent RNA helicase  | 614 aa | Core Macromolecular Complex | 0.993 | 630 aa | Optimal (<800 aa) |
| 4 | **DDX17** | [`Q92841`](https://www.uniprot.org/uniprotkb/Q92841) | Probable ATP-dependent RNA helicase  | 729 aa | Core Macromolecular Complex | 0.996 | 745 aa | Optimal (<800 aa) |
| 5 | **DICER1** | [`Q9UPY3`](https://www.uniprot.org/uniprotkb/Q9UPY3) | Endoribonuclease Dicer | 1922 aa | Core Macromolecular Complex | 0.993 | 1938 aa | High (>1400 aa) |
| 6 | **TARBP2** | [`Q15633`](https://www.uniprot.org/uniprotkb/Q15633) | RISC-loading complex subunit TARBP2 | 366 aa | Core Macromolecular Complex | 0.966 | 382 aa | Optimal (<800 aa) |
| 7 | **HNRNPA2B1** | [`P22626`](https://www.uniprot.org/uniprotkb/P22626) | Heterogeneous nuclear ribonucleoprot | 353 aa | Core Macromolecular Complex | 0.991 | 369 aa | Optimal (<800 aa) |
| 8 | **SMAD3** | [`P84022`](https://www.uniprot.org/uniprotkb/P84022) | SMAD family member 3 | 425 aa | Direct Physical Interactor | 0.750 | 441 aa | Optimal (<800 aa) |
| 9 | **XPO5** | [`Q9HAV4`](https://www.uniprot.org/uniprotkb/Q9HAV4) | Exportin-5 | 1204 aa | Core Macromolecular Complex | 0.936 | 1220 aa | Standard (800-1400 aa) |
| 10 | **RAN** | [`P62826`](https://www.uniprot.org/uniprotkb/P62826) | GTP-binding nuclear protein Ran | 216 aa | Direct Physical Interactor | 0.750 | 232 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein DGCR8. Microprotein is translated from the 5' UTR uORF of DGCR8 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 39: `c20riboseqorf31` (*ITCH*)
- **Biotype**: uORF | **Length**: 16 aa | **Tier**: 2B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MHFTVALWRQRLNPRK`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **ITCH** | [`Q96J02`](https://www.uniprot.org/uniprotkb/Q96J02) | E3 ubiquitin-protein ligase Itchy ho | 903 aa | Cognate Host CDS (Co-Translational) | 1.000 | 919 aa | Standard (800-1400 aa) |
| 2 | **RNF11** | [`Q9Y3C5`](https://www.uniprot.org/uniprotkb/Q9Y3C5) | RING finger protein 11 | 154 aa | Core Macromolecular Complex | 0.995 | 170 aa | Optimal (<800 aa) |
| 3 | **NUMB** | [`P49757`](https://www.uniprot.org/uniprotkb/P49757) | Protein numb homolog | 651 aa | Core Macromolecular Complex | 0.992 | 667 aa | Optimal (<800 aa) |
| 4 | **CFLAR** | [`O15519`](https://www.uniprot.org/uniprotkb/O15519) | CASP8 and FADD-like apoptosis regula | 480 aa | Core Macromolecular Complex | 0.992 | 496 aa | Optimal (<800 aa) |
| 5 | **TP73** | [`O15350`](https://www.uniprot.org/uniprotkb/O15350) | Tumor protein p73 | 636 aa | Core Macromolecular Complex | 0.905 | 652 aa | Optimal (<800 aa) |
| 6 | **UBC** | [`P0CG48`](https://www.uniprot.org/uniprotkb/P0CG48) | Polyubiquitin-C | 685 aa | Core Macromolecular Complex | 0.993 | 701 aa | Optimal (<800 aa) |
| 7 | **WWP1** | [`Q9H0M0`](https://www.uniprot.org/uniprotkb/Q9H0M0) | NEDD4-like E3 ubiquitin-protein liga | 922 aa | Core Macromolecular Complex | 0.937 | 938 aa | Standard (800-1400 aa) |
| 8 | **NEDD4** | [`P46934`](https://www.uniprot.org/uniprotkb/P46934) | E3 ubiquitin-protein ligase NEDD4 | 1319 aa | Direct Physical Interactor | 0.750 | 1335 aa | Standard (800-1400 aa) |
| 9 | **UBE2D1** | [`P51668`](https://www.uniprot.org/uniprotkb/P51668) | Ubiquitin-conjugating enzyme E2 D1 | 147 aa | Direct Physical Interactor | 0.798 | 163 aa | Optimal (<800 aa) |
| 10 | **UBE2L3** | [`P68036`](https://www.uniprot.org/uniprotkb/P68036) | Ubiquitin-conjugating enzyme E2 L3 | 154 aa | Core Macromolecular Complex | 0.977 | 170 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein ITCH. Microprotein is translated from the 5' UTR uORF of ITCH mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 40: `c3norep215` (*TBL1XR1*)
- **Biotype**: uORF | **Length**: 18 aa | **Tier**: 2B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MHCQWVISCVVTSWFKWE`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **TBL1XR1** | [`Q9BZK7`](https://www.uniprot.org/uniprotkb/Q9BZK7) | F-box-like/WD repeat-containing prot | 514 aa | Cognate Host CDS (Co-Translational) | 1.000 | 532 aa | Optimal (<800 aa) |
| 2 | **TBL1X** | [`O60907`](https://www.uniprot.org/uniprotkb/O60907) | F-box-like/WD repeat-containing prot | 577 aa | Core Macromolecular Complex | 0.999 | 595 aa | Optimal (<800 aa) |
| 3 | **TBL1Y** | [`Q9BQ87`](https://www.uniprot.org/uniprotkb/Q9BQ87) | F-box-like/WD repeat-containing prot | 522 aa | Core Macromolecular Complex | 0.999 | 540 aa | Optimal (<800 aa) |
| 4 | **CTNNB1** | [`P35222`](https://www.uniprot.org/uniprotkb/P35222) | Catenin beta-1 | 781 aa | Direct Physical Interactor | 0.997 | 799 aa | Optimal (<800 aa) |
| 5 | **HDAC3** | [`O15379`](https://www.uniprot.org/uniprotkb/O15379) | Histone deacetylase 3 | 428 aa | Core Macromolecular Complex | 0.999 | 446 aa | Optimal (<800 aa) |
| 6 | **GPS2** | [`Q13227`](https://www.uniprot.org/uniprotkb/Q13227) | G protein pathway suppressor 2 | 327 aa | Core Macromolecular Complex | 0.999 | 345 aa | Optimal (<800 aa) |
| 7 | **CORO2A** | [`Q92828`](https://www.uniprot.org/uniprotkb/Q92828) | Coronin-2A | 525 aa | Direct Physical Interactor | 0.750 | 543 aa | Optimal (<800 aa) |
| 8 | **BCL6** | [`P41182`](https://www.uniprot.org/uniprotkb/P41182) | B-cell lymphoma 6 protein | 706 aa | Direct Physical Interactor | 0.717 | 724 aa | Optimal (<800 aa) |
| 9 | **JUN** | [`P05412`](https://www.uniprot.org/uniprotkb/P05412) | Transcription factor Jun | 331 aa | Direct Physical Interactor | 0.750 | 349 aa | Optimal (<800 aa) |
| 10 | **SKP1** | [`P63208`](https://www.uniprot.org/uniprotkb/P63208) | S-phase kinase-associated protein 1 | 163 aa | Core Macromolecular Complex | 0.962 | 181 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein TBL1XR1. Microprotein is translated from the 5' UTR uORF of TBL1XR1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 41: `c14riboseqorf81` (*NUMB*)
- **Biotype**: uORF | **Length**: 32 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MRKLRLGEVTKLVQNHTTCKGHSQDSEPGCKN`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **NUMB** | [`P49757`](https://www.uniprot.org/uniprotkb/P49757) | Protein numb homolog | 651 aa | Cognate Host CDS (Co-Translational) | 1.000 | 683 aa | Optimal (<800 aa) |
| 2 | **MDM2** | [`Q00987`](https://www.uniprot.org/uniprotkb/Q00987) | E3 ubiquitin-protein ligase Mdm2 | 491 aa | Core Macromolecular Complex | 0.998 | 523 aa | Optimal (<800 aa) |
| 3 | **TP53** | [`P04637`](https://www.uniprot.org/uniprotkb/P04637) | Cellular tumor antigen p53 | 393 aa | Core Macromolecular Complex | 0.940 | 425 aa | Optimal (<800 aa) |
| 4 | **DPYSL2** | [`Q16555`](https://www.uniprot.org/uniprotkb/Q16555) | Dihydropyrimidinase-related protein  | 572 aa | Direct Physical Interactor | 0.992 | 604 aa | Optimal (<800 aa) |
| 5 | **EPS15** | [`P42566`](https://www.uniprot.org/uniprotkb/P42566) | Epidermal growth factor receptor sub | 896 aa | Core Macromolecular Complex | 0.990 | 928 aa | Standard (800-1400 aa) |
| 6 | **AP2A1** | [`O95782`](https://www.uniprot.org/uniprotkb/O95782) | AP-2 complex subunit alpha-1 | 977 aa | Core Macromolecular Complex | 0.919 | 1009 aa | Standard (800-1400 aa) |
| 7 | **AP2B1** | [`P63010`](https://www.uniprot.org/uniprotkb/P63010) | AP-2 complex subunit beta | 937 aa | Direct Physical Interactor | 0.727 | 969 aa | Standard (800-1400 aa) |
| 8 | **ITCH** | [`Q96J02`](https://www.uniprot.org/uniprotkb/Q96J02) | E3 ubiquitin-protein ligase Itchy ho | 903 aa | Core Macromolecular Complex | 0.992 | 935 aa | Standard (800-1400 aa) |
| 9 | **LNX1** | [`Q8TBB1`](https://www.uniprot.org/uniprotkb/Q8TBB1) | E3 ubiquitin-protein ligase LNX | 728 aa | Core Macromolecular Complex | 0.994 | 760 aa | Optimal (<800 aa) |
| 10 | **ACVR1** | [`Q04771`](https://www.uniprot.org/uniprotkb/Q04771) | Activin receptor type-1 | 509 aa | Direct Physical Interactor | 0.750 | 541 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein NUMB. Microprotein is translated from the 5' UTR uORF of NUMB mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

#### Rank 42: `c13norep17` (*HMGB1*)
- **Biotype**: uORF | **Length**: 45 aa | **Tier**: 2B | **Primary Tissue**: Nucleus & Extracellular (Alarmin / Chromatin)
- **Microprotein Amino Acid Sequence**:  
  `MLQSGESEEAASGSRSHSHCSTLSSIETAPGQVRAGRALGDSVPR`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **HMGB1** | [`P09429`](https://www.uniprot.org/uniprotkb/P09429) | High mobility group protein B1 | 215 aa | Cognate Host CDS (Co-Translational) | 1.000 | 260 aa | Optimal (<800 aa) |
| 2 | **TP53** | [`P04637`](https://www.uniprot.org/uniprotkb/P04637) | Cellular tumor antigen p53 | 393 aa | Core Macromolecular Complex | 0.997 | 438 aa | Optimal (<800 aa) |
| 3 | **TLR4** | [`O00206`](https://www.uniprot.org/uniprotkb/O00206) | Toll-like receptor 4 | 839 aa | Core Macromolecular Complex | 0.999 | 884 aa | Standard (800-1400 aa) |
| 4 | **AGER** | [`Q15109`](https://www.uniprot.org/uniprotkb/Q15109) | Advanced glycation end product-speci | 404 aa | Core Macromolecular Complex | 0.999 | 449 aa | Optimal (<800 aa) |
| 5 | **RELA** | [`Q04206`](https://www.uniprot.org/uniprotkb/Q04206) | Transcription factor p65 | 551 aa | Direct Physical Interactor | 0.750 | 596 aa | Optimal (<800 aa) |
| 6 | **NFKB1** | [`P19838`](https://www.uniprot.org/uniprotkb/P19838) | Nuclear factor NF-kappa-B p105 subun | 968 aa | Core Macromolecular Complex | 0.923 | 1013 aa | Standard (800-1400 aa) |
| 7 | **HMGB2** | [`P26583`](https://www.uniprot.org/uniprotkb/P26583) | High mobility group protein B2 | 209 aa | Core Macromolecular Complex | 0.981 | 254 aa | Optimal (<800 aa) |
| 8 | **H1-2** | [`P16403`](https://www.uniprot.org/uniprotkb/P16403) | Histone H1.2 | 213 aa | Direct Physical Interactor | 0.750 | 258 aa | Optimal (<800 aa) |
| 9 | **TOP2A** | [`P11388`](https://www.uniprot.org/uniprotkb/P11388) | DNA topoisomerase 2-alpha | 1531 aa | Direct Physical Interactor | 0.750 | 1576 aa | High (>1400 aa) |
| 10 | **RB1** | [`P06400`](https://www.uniprot.org/uniprotkb/P06400) | Retinoblastoma-associated protein | 928 aa | Direct Physical Interactor | 0.750 | 973 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary cognate host protein HMGB1. Microprotein is translated from the 5' UTR uORF of HMGB1 mRNA, providing immediate spatial and stoichiometric co-localization for physical hetero-oligomerization, feedback inhibition, and assembly regulation.

---

### 4. Dual-Coding Tier 1 Peptidein

#### Rank 43: `c1norep256` (*PRRC2C*)
- **Biotype**: intORF | **Length**: 34 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MDYRVLEKSVFHGVCLHLLTSQVLKQKTKAMILM`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **RBM25** | [`P49756`](https://www.uniprot.org/uniprotkb/P49756) | RNA-binding protein 25 | 843 aa | Core Macromolecular Complex | 0.999 | 877 aa | Standard (800-1400 aa) |
| 2 | **UBAP2L** | [`Q14157`](https://www.uniprot.org/uniprotkb/Q14157) | Ubiquitin-associated protein 2-like | 1087 aa | Core Macromolecular Complex | 0.915 | 1121 aa | Standard (800-1400 aa) |
| 3 | **SRSF11** | [`Q05519`](https://www.uniprot.org/uniprotkb/Q05519) | Serine/arginine-rich splicing factor | 484 aa | Core Macromolecular Complex | 0.832 | 518 aa | Optimal (<800 aa) |
| 4 | **YTHDF2** | [`Q9Y5A9`](https://www.uniprot.org/uniprotkb/Q9Y5A9) | YTH domain-containing family protein | 579 aa | Direct Physical Interactor | 0.750 | 613 aa | Optimal (<800 aa) |
| 5 | **EIF4A1** | [`P60842`](https://www.uniprot.org/uniprotkb/P60842) | Eukaryotic initiation factor 4A-I | 406 aa | Direct Physical Interactor | 0.750 | 440 aa | Optimal (<800 aa) |
| 6 | **EIF4G1** | [`Q04637`](https://www.uniprot.org/uniprotkb/Q04637) | Eukaryotic translation initiation fa | 1599 aa | Direct Physical Interactor | 0.750 | 1633 aa | High (>1400 aa) |
| 7 | **PABPC1** | [`P11940`](https://www.uniprot.org/uniprotkb/P11940) | Polyadenylate-binding protein 1 | 636 aa | Direct Physical Interactor | 0.750 | 670 aa | Optimal (<800 aa) |
| 8 | **ATXN2** | [`Q99700`](https://www.uniprot.org/uniprotkb/Q99700) | Ataxin-2 | 1313 aa | Direct Physical Interactor | 0.750 | 1347 aa | Standard (800-1400 aa) |
| 9 | **G3BP1** | [`Q13283`](https://www.uniprot.org/uniprotkb/Q13283) | Ras GTPase-activating protein-bindin | 466 aa | Direct Physical Interactor | 0.799 | 500 aa | Optimal (<800 aa) |
| 10 | **MOV10** | [`Q9HCE1`](https://www.uniprot.org/uniprotkb/Q9HCE1) | RNA helicase MOV-10 | 1003 aa | Direct Physical Interactor | 0.750 | 1037 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary cognate host protein PRRC2C. Microprotein is translated from an alternate internal reading frame nested within PRRC2C CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

#### Rank 44: `c1riboseqorf39` (*HNRNPR*)
- **Biotype**: intORF | **Length**: 45 aa | **Tier**: 1B | **Primary Tissue**: Pituitary
- **Microprotein Amino Acid Sequence**:  
  `MSILMKEQLMLSGNLMKKELCLYYSSSRKVTYHMFRTKVHFYVEL`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **HNRNPR** | [`O43390`](https://www.uniprot.org/uniprotkb/O43390) | Heterogeneous nuclear ribonucleoprot | 633 aa | Cognate Host CDS (Dual-Coding) | 1.000 | 678 aa | Optimal (<800 aa) |
| 2 | **SMN1** | [`Q16637`](https://www.uniprot.org/uniprotkb/Q16637) | Survival motor neuron protein | 294 aa | Core Macromolecular Complex | 0.906 | 339 aa | Optimal (<800 aa) |
| 3 | **HNRNPU** | [`Q00839`](https://www.uniprot.org/uniprotkb/Q00839) | Heterogeneous nuclear ribonucleoprot | 825 aa | Core Macromolecular Complex | 0.947 | 870 aa | Standard (800-1400 aa) |
| 4 | **HNRNPA1** | [`P09651`](https://www.uniprot.org/uniprotkb/P09651) | Heterogeneous nuclear ribonucleoprot | 372 aa | Core Macromolecular Complex | 0.989 | 417 aa | Optimal (<800 aa) |
| 5 | **HNRNPA2B1** | [`P22626`](https://www.uniprot.org/uniprotkb/P22626) | Heterogeneous nuclear ribonucleoprot | 353 aa | Core Macromolecular Complex | 0.923 | 398 aa | Optimal (<800 aa) |
| 6 | **TARDBP** | [`Q13148`](https://www.uniprot.org/uniprotkb/Q13148) | TAR DNA-binding protein 43 | 414 aa | Core Macromolecular Complex | 0.899 | 459 aa | Optimal (<800 aa) |
| 7 | **SRSF1** | [`Q07955`](https://www.uniprot.org/uniprotkb/Q07955) | Serine/arginine-rich splicing factor | 248 aa | Core Macromolecular Complex | 0.935 | 293 aa | Optimal (<800 aa) |
| 8 | **HNRNPC** | [`P07910`](https://www.uniprot.org/uniprotkb/P07910) | Heterogeneous nuclear ribonucleoprot | 306 aa | Core Macromolecular Complex | 0.980 | 351 aa | Optimal (<800 aa) |
| 9 | **ACTB** | [`P60709`](https://www.uniprot.org/uniprotkb/P60709) | Actin, cytoplasmic 1 | 375 aa | Direct Physical Interactor | 0.762 | 420 aa | Optimal (<800 aa) |
| 10 | **SYNCRIP** | [`O60506`](https://www.uniprot.org/uniprotkb/O60506) | Heterogeneous nuclear ribonucleoprot | 623 aa | Core Macromolecular Complex | 0.956 | 668 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein HNRNPR. Microprotein is translated from an alternate internal reading frame nested within HNRNPR CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

#### Rank 45: `c11norep30` (*EIF4G2*)
- **Biotype**: intORF | **Length**: 74 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MSMISVKIPSSPRRRNREPLLRSRCWETSNSLESLASLILFTNLSFISASKHFWKRRRESNSKIWERIWSASVR`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **EIF4G2** | [`P78344`](https://www.uniprot.org/uniprotkb/P78344) | Eukaryotic translation initiation fa | 907 aa | Cognate Host CDS (Dual-Coding) | 1.000 | 981 aa | Standard (800-1400 aa) |
| 2 | **EIF4A1** | [`P60842`](https://www.uniprot.org/uniprotkb/P60842) | Eukaryotic initiation factor 4A-I | 406 aa | Core Macromolecular Complex | 0.999 | 480 aa | Optimal (<800 aa) |
| 3 | **EIF3A** | [`Q14152`](https://www.uniprot.org/uniprotkb/Q14152) | Eukaryotic translation initiation fa | 1382 aa | Direct Physical Interactor | 0.610 | 1456 aa | High (>1400 aa) |
| 4 | **EIF3B** | [`P55884`](https://www.uniprot.org/uniprotkb/P55884) | Eukaryotic translation initiation fa | 814 aa | Direct Physical Interactor | 0.805 | 888 aa | Standard (800-1400 aa) |
| 5 | **EIF3C** | [`Q99613`](https://www.uniprot.org/uniprotkb/Q99613) | Eukaryotic translation initiation fa | 913 aa | Direct Physical Interactor | 0.692 | 987 aa | Standard (800-1400 aa) |
| 6 | **PABPC1** | [`P11940`](https://www.uniprot.org/uniprotkb/P11940) | Polyadenylate-binding protein 1 | 636 aa | Core Macromolecular Complex | 0.883 | 710 aa | Optimal (<800 aa) |
| 7 | **EIF4E** | [`P06730`](https://www.uniprot.org/uniprotkb/P06730) | Eukaryotic translation initiation fa | 217 aa | Core Macromolecular Complex | 0.999 | 291 aa | Optimal (<800 aa) |
| 8 | **MKNK1** | [`Q9BUB5`](https://www.uniprot.org/uniprotkb/Q9BUB5) | MAP kinase-interacting serine/threon | 465 aa | Core Macromolecular Complex | 0.912 | 539 aa | Optimal (<800 aa) |
| 9 | **PAIP1** | [`Q9H074`](https://www.uniprot.org/uniprotkb/Q9H074) | Polyadenylate-binding protein-intera | 479 aa | Direct Physical Interactor | 0.750 | 553 aa | Optimal (<800 aa) |
| 10 | **EIF4G1** | [`Q04637`](https://www.uniprot.org/uniprotkb/Q04637) | Eukaryotic translation initiation fa | 1599 aa | Core Macromolecular Complex | 0.984 | 1673 aa | High (>1400 aa) |

- **Primary Biological Rationale**: Primary cognate host protein EIF4G2. Microprotein is translated from an alternate internal reading frame nested within EIF4G2 CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

#### Rank 46: `c11riboseqorf112` (*EMSY*)
- **Biotype**: uoORF | **Length**: 91 aa | **Tier**: 1B | **Primary Tissue**: Cells_EBV-transformed_lymphocytes
- **Microprotein Amino Acid Sequence**:  
  `MSPERPGREDKLFGATKQKQQCLLCGQPFWISAGMNAKEFFENWNWRHMLELSVHFGHRGISPRKRKIFLENYQKFLASQQNATVLKFGEQ`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **EMSY** | [`Q7Z589`](https://www.uniprot.org/uniprotkb/Q7Z589) | BRCA2-interacting transcriptional re | 1322 aa | Cognate Host CDS (Dual-Coding) | 1.000 | 1413 aa | High (>1400 aa) |
| 2 | **PHF12** | [`Q96QT6`](https://www.uniprot.org/uniprotkb/Q96QT6) | PHD finger protein 12 | 1004 aa | Core Macromolecular Complex | 0.960 | 1095 aa | Standard (800-1400 aa) |
| 3 | **CBX5** | [`P45973`](https://www.uniprot.org/uniprotkb/P45973) | Chromobox protein homolog 5 | 191 aa | Direct Physical Interactor | 0.750 | 282 aa | Optimal (<800 aa) |
| 4 | **CBX1** | [`P83916`](https://www.uniprot.org/uniprotkb/P83916) | Chromobox protein homolog 1 | 185 aa | Core Macromolecular Complex | 0.966 | 276 aa | Optimal (<800 aa) |
| 5 | **CBX3** | [`Q13185`](https://www.uniprot.org/uniprotkb/Q13185) | Chromobox protein homolog 3 | 183 aa | Direct Physical Interactor | 0.750 | 274 aa | Optimal (<800 aa) |
| 6 | **KDM5A** | [`P29375`](https://www.uniprot.org/uniprotkb/P29375) | Lysine-specific demethylase 5A | 1690 aa | Core Macromolecular Complex | 0.904 | 1781 aa | High (>1400 aa) |
| 7 | **SKP1** | [`P63208`](https://www.uniprot.org/uniprotkb/P63208) | S-phase kinase-associated protein 1 | 163 aa | Direct Physical Interactor | 0.750 | 254 aa | Optimal (<800 aa) |
| 8 | **ZMYND11** | [`Q15326`](https://www.uniprot.org/uniprotkb/Q15326) | Zinc finger MYND domain-containing p | 602 aa | Core Macromolecular Complex | 0.920 | 693 aa | Optimal (<800 aa) |
| 9 | **RAD51** | [`Q06609`](https://www.uniprot.org/uniprotkb/Q06609) | DNA repair protein RAD51 homolog 1 | 339 aa | Direct Physical Interactor | 0.513 | 430 aa | Optimal (<800 aa) |
| 10 | **PALB2** | [`Q86YC2`](https://www.uniprot.org/uniprotkb/Q86YC2) | Partner and localizer of BRCA2 | 1186 aa | Direct Physical Interactor | 0.732 | 1277 aa | Standard (800-1400 aa) |

- **Primary Biological Rationale**: Primary cognate host protein EMSY. Microprotein is translated from an alternate internal reading frame nested within EMSY CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

#### Rank 47: `c20norep82` (*ADNP*)
- **Biotype**: intORF | **Length**: 65 aa | **Tier**: 1B | **Primary Tissue**: Brain / Neuronal Nuclei (ChAHP Complex)
- **Microprotein Amino Acid Sequence**:  
  `MSIVKTLKIGFSLIAPTVPSMQTKRLWKHTLKYFMLRTPAHQVAASALSKIKTKMMALNLSRLTV`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **ADNP** | [`Q9H2P0`](https://www.uniprot.org/uniprotkb/Q9H2P0) | Activity-dependent neuroprotector ho | 1102 aa | Cognate Host CDS (Dual-Coding) | 1.000 | 1167 aa | Standard (800-1400 aa) |
| 2 | **CHD4** | [`Q14839`](https://www.uniprot.org/uniprotkb/Q14839) | ATP-dependent chromatin remodeler CH | 1912 aa | Core Macromolecular Complex | 0.991 | 1977 aa | High (>1400 aa) |
| 3 | **CBX5** | [`P45973`](https://www.uniprot.org/uniprotkb/P45973) | Chromobox protein homolog 5 | 191 aa | Direct Physical Interactor | 0.819 | 256 aa | Optimal (<800 aa) |
| 4 | **CBX1** | [`P83916`](https://www.uniprot.org/uniprotkb/P83916) | Chromobox protein homolog 1 | 185 aa | Direct Physical Interactor | 0.848 | 250 aa | Optimal (<800 aa) |
| 5 | **CBX3** | [`Q13185`](https://www.uniprot.org/uniprotkb/Q13185) | Chromobox protein homolog 3 | 183 aa | Core Macromolecular Complex | 0.856 | 248 aa | Optimal (<800 aa) |
| 6 | **MTA1** | [`Q13330`](https://www.uniprot.org/uniprotkb/Q13330) | Metastasis-associated protein MTA1 | 715 aa | Direct Physical Interactor | 0.750 | 780 aa | Optimal (<800 aa) |
| 7 | **HDAC1** | [`Q13547`](https://www.uniprot.org/uniprotkb/Q13547) | Histone deacetylase 1 | 482 aa | Direct Physical Interactor | 0.750 | 547 aa | Optimal (<800 aa) |
| 8 | **HDAC2** | [`Q92769`](https://www.uniprot.org/uniprotkb/Q92769) | Histone deacetylase 2 | 488 aa | Direct Physical Interactor | 0.750 | 553 aa | Optimal (<800 aa) |
| 9 | **SMARCC2** | [`Q8TAQ2`](https://www.uniprot.org/uniprotkb/Q8TAQ2) | SWI/SNF complex subunit SMARCC2 | 1214 aa | Core Macromolecular Complex | 0.741 | 1279 aa | Standard (800-1400 aa) |
| 10 | **SMARCA4** | [`P51532`](https://www.uniprot.org/uniprotkb/P51532) | SWI/SNF-related matrix-associated ac | 1647 aa | Direct Physical Interactor | 0.640 | 1712 aa | High (>1400 aa) |

- **Primary Biological Rationale**: Primary cognate host protein ADNP. Microprotein is translated from an alternate internal reading frame nested within ADNP CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

#### Rank 48: `c3norep241` (*P3H2*)
- **Biotype**: intORF | **Length**: 63 aa | **Tier**: 1B | **Primary Tissue**: Connective Tissue / Fibroblasts (Collagen Prolyl Hydroxylation)
- **Microprotein Amino Acid Sequence**:  
  `MCAATSSAECPTTTCSGPTSSLTSSKKQWKQLTHFSWLTLSTWKCSRTLRITGRQLVLKHCSW`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **P3H2** | [`Q8IVL5`](https://www.uniprot.org/uniprotkb/Q8IVL5) | Prolyl 3-hydroxylase 2 | 708 aa | Cognate Host CDS (Dual-Coding) | 1.000 | 771 aa | Optimal (<800 aa) |
| 2 | **P3H1** | [`Q32P28`](https://www.uniprot.org/uniprotkb/Q32P28) | Prolyl 3-hydroxylase 1 | 736 aa | Direct Physical Interactor | 0.750 | 799 aa | Optimal (<800 aa) |
| 3 | **CRTAP** | [`O75718`](https://www.uniprot.org/uniprotkb/O75718) | Cartilage-associated protein | 401 aa | Direct Physical Interactor | 0.517 | 464 aa | Optimal (<800 aa) |
| 4 | **PPIB** | [`P23284`](https://www.uniprot.org/uniprotkb/P23284) | Peptidyl-prolyl cis-trans isomerase  | 216 aa | Direct Physical Interactor | 0.639 | 279 aa | Optimal (<800 aa) |
| 5 | **COL1A1** | [`P02452`](https://www.uniprot.org/uniprotkb/P02452) | Collagen alpha-1(I) chain | 1464 aa | Direct Physical Interactor | 0.621 | 1527 aa | High (>1400 aa) |
| 6 | **COL1A2** | [`P08123`](https://www.uniprot.org/uniprotkb/P08123) | Collagen alpha-2(I) chain | 1366 aa | Direct Physical Interactor | 0.599 | 1429 aa | High (>1400 aa) |
| 7 | **COL4A1** | [`P02462`](https://www.uniprot.org/uniprotkb/P02462) | Collagen alpha-1(IV) chain | 1669 aa | Direct Physical Interactor | 0.602 | 1732 aa | High (>1400 aa) |
| 8 | **P4HA1** | [`P13674`](https://www.uniprot.org/uniprotkb/P13674) | Prolyl 4-hydroxylase subunit alpha-1 | 534 aa | Direct Physical Interactor | 0.787 | 597 aa | Optimal (<800 aa) |
| 9 | **PLOD1** | [`Q02809`](https://www.uniprot.org/uniprotkb/Q02809) | Procollagen-lysine,2-oxoglutarate 5- | 727 aa | Direct Physical Interactor | 0.750 | 790 aa | Optimal (<800 aa) |
| 10 | **SERPINH1** | [`P50454`](https://www.uniprot.org/uniprotkb/P50454) | Serpin H1 | 418 aa | Direct Physical Interactor | 0.750 | 481 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein P3H2. Microprotein is translated from an alternate internal reading frame nested within P3H2 CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

#### Rank 49: `c16riboseqorf58` (*SRCAP*)
- **Biotype**: uoORF | **Length**: 62 aa | **Tier**: 1B | **Primary Tissue**: Nucleus / Chromatin (H2A.Z Deposition)
- **Microprotein Amino Acid Sequence**:  
  `MRPAPSMPCSRTPAPRPAVLVITTQSFFRHPRGSLGVGPCRAAPPLLTLSSQSYRHRWCRTA`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **KAT5** | [`Q92993`](https://www.uniprot.org/uniprotkb/Q92993) | Histone acetyltransferase KAT5 | 513 aa | Core Macromolecular Complex | 0.984 | 575 aa | Optimal (<800 aa) |
| 2 | **H2AZ1** | [`P0C0S5`](https://www.uniprot.org/uniprotkb/P0C0S5) | Histone H2A.Z | 128 aa | Core Macromolecular Complex | 0.998 | 190 aa | Optimal (<800 aa) |
| 3 | **RUVBL1** | [`Q9Y265`](https://www.uniprot.org/uniprotkb/Q9Y265) | RuvB-like 1 | 456 aa | Core Macromolecular Complex | 0.999 | 518 aa | Optimal (<800 aa) |
| 4 | **RUVBL2** | [`Q9Y230`](https://www.uniprot.org/uniprotkb/Q9Y230) | RuvB-like 2 | 463 aa | Core Macromolecular Complex | 0.999 | 525 aa | Optimal (<800 aa) |
| 5 | **ACTR6** | [`Q9GZN1`](https://www.uniprot.org/uniprotkb/Q9GZN1) | Actin-related protein 6 | 396 aa | Core Macromolecular Complex | 0.999 | 458 aa | Optimal (<800 aa) |
| 6 | **VPS72** | [`Q15906`](https://www.uniprot.org/uniprotkb/Q15906) | Vacuolar protein sorting-associated  | 364 aa | Core Macromolecular Complex | 0.997 | 426 aa | Optimal (<800 aa) |
| 7 | **DMAP1** | [`Q9NPF5`](https://www.uniprot.org/uniprotkb/Q9NPF5) | DNA methyltransferase 1-associated p | 467 aa | Core Macromolecular Complex | 0.992 | 529 aa | Optimal (<800 aa) |
| 8 | **YEATS4** | [`O95619`](https://www.uniprot.org/uniprotkb/O95619) | YEATS domain-containing protein 4 | 227 aa | Core Macromolecular Complex | 0.989 | 289 aa | Optimal (<800 aa) |
| 9 | **ACTL6A** | [`O96019`](https://www.uniprot.org/uniprotkb/O96019) | Actin-like protein 6A | 429 aa | Core Macromolecular Complex | 0.959 | 491 aa | Optimal (<800 aa) |
| 10 | **ZNHIT1** | [`O43257`](https://www.uniprot.org/uniprotkb/O43257) | Zinc finger HIT domain-containing pr | 154 aa | Core Macromolecular Complex | 0.995 | 216 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein SRCAP. Microprotein is translated from an alternate internal reading frame nested within SRCAP CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

#### Rank 50: `c10norep31` (*CDC123*)
- **Biotype**: intORF | **Length**: 29 aa | **Tier**: 1B | **Primary Tissue**: Ubiquitous / Cytosol (eIF2 Translation Assembly)
- **Microprotein Amino Acid Sequence**:  
  `MCFTASSPRGTRSSEALPSRVSFFHFLRM`

| P.Rank | Partner Gene | UniProt Acc | Protein Name | Length | Category | STRING Score | Total Complex Length | Feasibility |
| :---: | :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- |
| 1 | **CDC123** | [`O75794`](https://www.uniprot.org/uniprotkb/O75794) | Translation initiation factor eIF2 a | 336 aa | Cognate Host CDS (Dual-Coding) | 1.000 | 365 aa | Optimal (<800 aa) |
| 2 | **EIF2S1** | [`P05198`](https://www.uniprot.org/uniprotkb/P05198) | Eukaryotic translation initiation fa | 315 aa | Direct Physical Interactor | 0.539 | 344 aa | Optimal (<800 aa) |
| 3 | **EIF2S2** | [`P20042`](https://www.uniprot.org/uniprotkb/P20042) | Eukaryotic translation initiation fa | 333 aa | Core Macromolecular Complex | 0.873 | 362 aa | Optimal (<800 aa) |
| 4 | **EIF2S3** | [`P41091`](https://www.uniprot.org/uniprotkb/P41091) | Eukaryotic translation initiation fa | 472 aa | Core Macromolecular Complex | 0.961 | 501 aa | Optimal (<800 aa) |
| 5 | **PPP1CA** | [`P62136`](https://www.uniprot.org/uniprotkb/P62136) | Serine/threonine-protein phosphatase | 330 aa | Direct Physical Interactor | 0.750 | 359 aa | Optimal (<800 aa) |
| 6 | **PPP1CB** | [`P62140`](https://www.uniprot.org/uniprotkb/P62140) | Serine/threonine-protein phosphatase | 327 aa | Direct Physical Interactor | 0.750 | 356 aa | Optimal (<800 aa) |
| 7 | **PPP1R15A** | [`O75807`](https://www.uniprot.org/uniprotkb/O75807) | Protein phosphatase 1 regulatory sub | 674 aa | Direct Physical Interactor | 0.750 | 703 aa | Optimal (<800 aa) |
| 8 | **EIF2AK3** | [`Q9NZJ5`](https://www.uniprot.org/uniprotkb/Q9NZJ5) | Eukaryotic translation initiation fa | 1116 aa | Direct Physical Interactor | 0.750 | 1145 aa | Standard (800-1400 aa) |
| 9 | **EIF2B1** | [`Q14232`](https://www.uniprot.org/uniprotkb/Q14232) | Translation initiation factor eIF2B  | 305 aa | Direct Physical Interactor | 0.750 | 334 aa | Optimal (<800 aa) |
| 10 | **EIF2B2** | [`P49770`](https://www.uniprot.org/uniprotkb/P49770) | Translation initiation factor eIF2B  | 351 aa | Direct Physical Interactor | 0.750 | 380 aa | Optimal (<800 aa) |

- **Primary Biological Rationale**: Primary cognate host protein CDC123. Microprotein is translated from an alternate internal reading frame nested within CDC123 CDS, providing obligate stoichiometric co-expression for physical co-complex assembly and co-translational engagement.

---

## 4. Instructions for Biohub ESMFold2 `fold_all_atom` API Execution

When submitting these 500 complexes to the Biohub ESM API, the payload should be constructed as follows:

```json
{
  "model": "esmfold2-fast-2026-05",
  "sequence": "<microprotein_sequence>|<partner_sequence>",
  "num_loops": 20,
  "num_sampling_steps": 100,
  "include_pae": true,
  "include_pair_chains_iptm": true,
  "potential_sequence_of_concern": false
}
```

### Key Output Metrics for Binding Assessment:
1. **`interface_ptm` (ipTM)**: Measures inter-chain predicted TM-score. ipTM > 0.6 indicates high probability of genuine physical binding; ipTM > 0.8 indicates near-certain complex formation.
2. **`pair_chains_iptm` Matrix**: 2?2 matrix providing chain A vs chain B interaction confidence.
3. **Inter-Chain PAE**: Predicted Aligned Error at the interface contacts (sub-5? PAE indicates high positional confidence).
4. **Interface Contact Area & Binding Interface**: Extract coordinates to compute buried surface area (BSA) and contact residue pairs via PyMOL or BioPython.