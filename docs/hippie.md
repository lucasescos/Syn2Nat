# HIPPIE: Human Integrated Protein-Protein Interaction rEference

> **Official Web Server**: [https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/)  
> **Information & Documentation**: [https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/information.php](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/information.php)  
> **Local Datasets**: [`data/hippie/`](../data/hippie/)  
> **Python Client**: [`scripts/hippie_client.py`](../scripts/hippie_client.py)  
> **Related Framework**: [MAPPIE (`docs/mappie.md`)](mappie.md) (uses high-confidence HIPPIE PPIs $s \ge 0.65$ as ground truth)

---

## 1. Overview & System Architecture

**HIPPIE** (Human Integrated Protein-Protein Interaction rEference) is a curated, integrated repository of human protein-protein interactions (PPIs) developed and maintained by the **Computational Biology & Data Mining (CBDM)** group at Johannes Gutenberg University (JGU) Mainz.

### Key Highlights
* **Current Release**: **v2.4** (April 2026 update, incorporating 275,000+ new interactions).
* **Multi-Source Integration**: Integrates physical interactions from 10 major primary PPI databases (BioGRID, IntAct, MINT, HPRD, DIP, BIND, I2D, STRING, etc.) and individual high-throughput screens.
* **Confidence Scoring**: Every interaction is evaluated and assigned an experiment-based confidence score between **0.00** and **1.00**, reflecting experimental reliability, reproducibility across studies, and cross-species conservation.
* **Context Annotations**: Enriched with GTEx tissue expression profiles (RNA-seq RPKM $\ge 1$), Gene Ontology / MeSH terms, interaction types (direct physical vs. complex association), edge directionality, and predicted activating/inhibitory effects.

```
Primary PPI Databases (BioGRID, IntAct, MINT, HPRD, DIP, etc.)
                          │
                          ▼
            [ HIPPIE Integration Pipeline ]
     - ID Harmonization (UniProt, Entrez, Gene Symbols)
     - Quality Weighting of Experimental Techniques
     - Multi-Study & Orthology Verification
                          │
                          ▼
         [ Confidence Scored Interactome ]
                 Score: 0.00 – 1.00
    ├── Low Confidence:     0.00 <= score < 0.63
    ├── Medium Confidence:  0.63 <= score < 0.72 (Q2 quartile)
    └── High Confidence:    0.72 <= score <= 1.00 (Q3 quartile)
                          │
                          ▼
         [ Access Points & Downstream Tools ]
    ├── REST API (`fast_query_tissue.php` / `queryHIPPIE.php`)
    ├── Flat Files (`hippie_current.txt`, `HIPPIE-current.mitab.txt`)
    └── MAPPIE Deep Learning Embeddings (`esm2_650M` + Autoencoder)
```

---

## 2. Confidence Scoring Scheme

A central feature of HIPPIE is its confidence scoring formula, parameterized via joint optimization by domain experts and computer algorithms.

### Mathematical Formulation
The overall score $S \in [0, 1]$ is a weighted sum of three subscores:

$$S = w_s \cdot S_{\text{study}} + w_t \cdot S_{\text{technique}} + w_o \cdot S_{\text{orthology}}$$

Where the optimized default weights are:
* **Study Weight ($w_s$)**: $0.60$ (number of independent PubMed studies reporting the interaction)
* **Technique Weight ($w_t$)**: $0.30$ (reliability scores of the experimental methods used)
* **Orthology Weight ($w_o$)**: $0.10$ (reproduction of the interaction in non-human organisms)

$$\sum w = w_s + w_t + w_o = 0.60 + 0.30 + 0.10 = 1.00$$

### Subscores & Saturation Parameters
Each subscore saturates asymptotically as evidence accumulates, governed by parameter $a$:

$$S_x = 1 - e^{-a_x \cdot x}$$

* **Study Saturation ($a_s$)**: $2.3$ (reflecting rapid confidence increase from multiple independent publications)
* **Technique Saturation ($a_t$)**: $0.2$ (summing the quality scores assigned to each experimental technique)
* **Orthology Saturation ($a_o$)**: $1.6$ (evaluating conservation in model organisms like *M. musculus*, *D. melanogaster*, *C. elegans*, *S. cerevisiae*)

### Experimental Technique Quality Scores
HIPPIE assigns a quality score (from 1 to 10) to each experimental technique based on its false-positive rate and biochemical rigor. The complete list of 126 technique scores is preserved locally at [`data/hippie/experimental_scores.tsv`](../data/hippie/experimental_scores.tsv).

Selected examples:
| Technique | PSI-MI ID | Quality Score |
| :--- | :--- | :--- |
| X-ray Crystallography | `MI:0114` | 9.0 |
| Nuclear Magnetic Resonance (NMR) | `MI:0077` | 8.0 |
| Surface Plasmon Resonance (SPR) | `MI:0107` | 8.0 |
| Isothermal Titration Calorimetry (ITC) | `MI:0065` | 8.0 |
| Atomic Force Microscopy | `MI:0872` | 9.0 |
| Two-hybrid (Y2H) | `MI:0018` | 5.0 |
| Affinity Chromatography / Pull Down | `MI:0004` / `MI:0096` | 5.0 |
| Tandem Affinity Purification (TAP) | `MI:0676` | 5.0 |
| BioID / Proximity Labeling | `MI:1313` | 5.0 |
| Cross-linking Study | `MI:0030` | 5.0 |
| Array Technology | `MI:0008` | 3.0 |

### Confidence Thresholds
* **Medium Confidence**: $\ge 0.63$ (second quartile / median of HIPPIE distribution).
* **High Confidence**: $\ge 0.72$ (or $0.73$, third quartile of HIPPIE distribution).
* **MAPPIE Benchmark Cutoff**: $s \ge 0.65$ was selected by Cihan et al. (2026) to train the MAPPIE latent autoencoder.

---

## 3. REST API Reference

HIPPIE provides programmatic query access via HTTP endpoints hosted at `cbdm-01.zdv.uni-mainz.de`.

### Endpoints

#### 1. Direct Data Endpoint (Recommended for Pipelines)
* **URL**: `https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/fast_query_tissue.php`
* **Method**: `POST`
* **Behavior**: Returns the raw interaction data immediately in the requested format (`conc_file` TSV or `mitab` PSI-MI TAB).

#### 2. Query Wrapper Endpoint
* **URL**: `https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/queryHIPPIE.php`
* **Method**: `GET`
* **Behavior**: Designed for web browsers; returns an HTML page containing a JavaScript auto-submitting POST form targeting `fast_query_tissue.php`.

### Query Parameters

| Parameter | Type | Required | Values / Format | Description |
| :--- | :--- | :--- | :--- | :--- |
| `query_genes` (or `proteins`) | string | **Yes** | Multiline string (`\n`), comma, semicolon, or pipe separated | Identifier(s) of interest: UniProt accession (e.g. `P04637`), UniProt mnemonic (`P53_HUMAN`), Gene symbol (`TP53`), or Entrez ID (`7157`). Max 100 entries per batch. |
| `layers` | integer | No | `0` or `1` (default: `1`) | `0`: interactions strictly *within* the input set.<br>`1`: interactions between the input set and all HIPPIE proteins. |
| `conf_thres` | float | No | `0.0` – `1.0` (default: `0.0`) | Minimum confidence score threshold. Only interactions with $s \ge \text{conf\_thres}$ are returned. |
| `out_type` | string | No | `conc_file`, `mitab`, `browser`, `viz` | Output format (default: `conc_file`). |
| `n_links` | integer | No | $\ge 1$ (default: `1`) | Minimum number of interactions to the query set. Setting to $\ge 2$ filters for common interactors. |
| `direction_type` | integer | No | `0`, `1`, `2`, `3` (default: `0`) | `0`: None.<br>`1`: Unweighted shortest paths.<br>`2`: Confidence-weighted shortest paths.<br>`3`: KEGG direction. |
| `effect_type` | integer | No | `0`, `1`, `2` (default: `0`) | `0`: None.<br>`1`: Predicted effect (activating vs inhibitory).<br>`2`: KEGG effect. |
| `assoc` | flag | No | present / absent | Physical or functional association. |
| `phys_assoc` | flag | No | present / absent | Physical association. |
| `dir_int` | flag | No | present / absent | Direct interaction only. |
| `coloc` | flag | No | present / absent | Colocalization. |

### Output Formats

#### 1. `conc_file` (HIPPIE Concise TAB)
Tab-separated table with 7 columns:
```tsv
uniprot id 1	entrez gene id 1	gene name 1	uniprot id 2	entrez gene id 2	gene name 2	score
P53_HUMAN	7157	TP53	MDM2_HUMAN	4193	MDM2	0.97
P53_HUMAN	7157	TP53	DAXX_HUMAN	1616	DAXX	0.99
```

#### 2. `mitab` (PSI-MI TAB 2.5)
Standardized 15-column format specified by the Proteomics Standards Initiative (PSI), detailing:
1. Identifier interactor A (e.g. `entrez gene:7157`)
2. Identifier interactor B (e.g. `entrez gene:4193`)
3. Alternative identifier interactor A (e.g. `uniprotkb:P53_HUMAN`)
4. Alternative identifier interactor B (e.g. `uniprotkb:MDM2_HUMAN`)
5. Aliases interactor A
6. Aliases interactor B
7. Interaction detection methods (PSI-MI terms and free text, e.g. `MI:0018(Two-hybrid)`)
8. First author / publication
9. Publication identifiers (PubMed IDs, e.g. `pubmed:12667443|pubmed:16151013`)
10. Taxid interactor A (`taxid:9606(Homo sapiens)`)
11. Taxid interactor B (`taxid:9606(Homo sapiens)`)
12. Interaction types
13. Source databases (`BioGRID|HPRD|IntAct|MINT|STRING`)
14. Interaction identifiers
15. Confidence score (e.g. `0.97`)

---

## 4. Bulk Datasets & Offline Tools

For large-scale offline analyses (such as whole-interactome graph learning or batch scoring), HIPPIE provides flat files and Java tools:

### Available Downloads
* **HIPPIE TAB Flat File**: [`https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt`](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt) (~138 MB uncompressed)
  * Format: `UniProt_A \t Entrez_A \t UniProt_B \t Entrez_B \t Score \t Evidence`
* **HIPPIE PSI-MI TAB 2.5**: [`https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/HIPPIE-current.mitab.txt`](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/HIPPIE-current.mitab.txt)
* **Experimental Scores Table**: Stored locally at [`data/hippie/experimental_scores.tsv`](../data/hippie/experimental_scores.tsv).
* **Network Construction Tool**: [`HIPPIE_NC.jar`](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/NC/HIPPIE_NC.jar) (Documentation: [`data/hippie/README_NC.txt`](../data/hippie/README_NC.txt)).
* **Rescoring Tool**: [`HIPPIE_RS.jar`](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/RS/HIPPIE_RS.jar) (Documentation: [`data/hippie/README_RS.txt`](../data/hippie/README_RS.txt)).

---

## 5. Python Client & Usage Examples

A complete, production-ready Python client is provided in [`scripts/hippie_client.py`](../scripts/hippie_client.py).

### CLI Usage Examples

```powershell
# Query interactors for a single protein (confidence >= 0.72)
python scripts/hippie_client.py --proteins TP53 --conf 0.72

# Check interactions strictly WITHIN a query set (Layer 0)
python scripts/hippie_client.py --proteins TP53,MDM2,CDKN1A,ATM --layer 0

# Query by UniProt accession
python scripts/hippie_client.py --proteins P04637 --conf 0.85

# Export to CSV or JSON
python scripts/hippie_client.py --proteins BRCA1,BARD1 --format mitab --output ppi.tsv
```

### Python Programmatic Examples

```python
from scripts.hippie_client import HippieClient

client = HippieClient()

# 1. Query interactors of TP53 above high-confidence threshold
df_p53 = client.query(
    proteins=["TP53"],
    conf_thres=0.72,
    layer=1
)
print(f"Found {len(df_p53)} high-confidence interactors for TP53")
print(df_p53[["gene_1", "gene_2", "score"]].head())

# 2. Check if two specific proteins interact and get their score
interaction = client.query_pair("TP53", "MDM2")
if interaction is not None:
    print(f"Interaction Score: {interaction['score']}")

# 3. Construct an internal subnetwork (Layer 0) from candidate microprotein targets
candidates = ["TP53", "MDM2", "CDKN1A", "BAX", "PML", "EP300"]
subnetwork = client.query_subnetwork(candidates, conf_thres=0.63)
print(subnetwork)
```

---

## 6. Tissue Specificity (GTEx Integration)

HIPPIE integrates RNA-Seq expression data from the **Genotype-Tissue Expression (GTEx)** project across 53 non-diseased human tissues. A gene is defined as expressed in a given tissue when $\text{RPKM} \ge 1.0$.

When the tissue filter is activated, **both** interacting partners must be expressed in at least one of the selected tissues for the interaction edge to be preserved in the constructed network.

Available tissue filters include:
* Adipose Tissue (Subcutaneous, Visceral)
* Adrenal Gland
* Brain (13 specific subregions: Cortex, Cerebellum, Hippocampus, Striatum, Hypothalamus, etc.)
* Breast (Mammary Tissue)
* Colon (Sigmoid, Transverse)
* Esophagus (Mucosa, Muscularis, Gastroesophageal Junction)
* Heart (Left Ventricle, Atrial Appendage)
* Kidney, Liver, Lung
* Muscle (Skeletal)
* Nerve (Tibial)
* Ovary, Pancreas, Pituitary, Prostate
* Skin (Suprapubic, Lower leg)
* Spleen, Stomach, Testis, Thyroid, Uterus, Vagina, Whole Blood.

---

## 7. Related Manuscripts & Citations

1. **HIPPIE v1.0 (Scoring & Expert System)**:
   > Schaefer MH, Fontaine J-F, Vinayagam A, Porras P, Wanker EE, Andrade-Navarro MA (2012). *HIPPIE: Integrating Protein Interaction Networks with Experiment Based Quality Scores*. **PLoS ONE**, 7(2): e31826.  
   > [DOI: 10.1371/journal.pone.0031826](https://doi.org/10.1371/journal.pone.0031826)

2. **HIPPIE Context & Tissue Filtering**:
   > Schaefer MH, Lopes TJS, Mah N, Shoemaker JE, Matsuoka Y, Fontaine JF, et al. (2013). *Adding Protein Context to the Human Protein-Protein Interaction Network to Reveal Meaningful Interactions*. **PLoS Computational Biology**, 9(1): e1002860.  
   > [DOI: 10.1371/journal.pcbi.1002860](https://doi.org/10.1371/journal.pcbi.1002860)

3. **Phenotypic Effects & Directionality**:
   > Suratanee A, Schaefer MH, Betts M, Soons Z, Mannsperger H, Harder N, et al. (2014). *Characterizing Protein Interactions Employing a Genome-Wide siRNA Cellular Phenotyping Screen*. **PLoS Computational Biology**, 10(9): e1003814.  
   > [DOI: 10.1371/journal.pcbi.1003814](https://doi.org/10.1371/journal.pcbi.1003814)

4. **HIPPIE v2.0**:
   > Alanis-Lobato G, Andrade-Navarro MA, Schaefer MH (2016). *HIPPIE v2.0: enhancing meaningfulness and reliability of protein-protein interaction networks*. **Nucleic Acids Research**, 45(D1): D408–D414.  
   > [DOI: 10.1093/nar/gkw985](https://doi.org/10.1093/nar/gkw985)

5. **MAPPIE (Embedding HIPPIE into ESM-2 Latent Space)**:
   > Cihan M, Distler U, Andrade-Navarro MA (2026). *A map of human protein-protein interaction embeddings for functional discovery*. **bioRxiv**, 2026.08.07.743440.  
   > [DOI: 10.64898/2026.08.07.743440](https://doi.org/10.64898/2026.08.07.743440)
