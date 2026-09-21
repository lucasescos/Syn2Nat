# Strategic Roadmap: Large-Scale Microprotein Representation Engine

## 1. Executive Directive & Strategic Pivot

> **Key Context**: The prior genomic selection analysis (Google DeepMind AlphaGenome AVI Phred purifying selection) and CDS-overlapping hitchhiking de-confounding analyses were **exploratory studies**. 

### The Core Objective
The primary mission of this project is:
1. **Curate and Scale**: Acquire a massive, comprehensive catalog of human microproteins (expanding from 7,264 to **617,462 candidates** via the Human Microprotein Atlas / HMPA).
2. **Representation Extraction**: Systematically extract deep protein language model embeddings (penultimate layer representations, e.g., Biohub `esmc-6b-2024-12` layer 79 `[2560]` vectors and `[L, 2560]` per-residue matrices).
3. **Quota-Managed Daily Cadence**: Because Biohub ESM / cloud model inference operates under daily API quotas, extraction is executed via automated, checkpointed daily batches.
4. **Decoupled Architecture Design**: Decouple the data collection and embedding extraction phase from final downstream architecture decisions. By banking a rich, high-dimensional representation store first, downstream neural architectures (binding affinity predictors, functional classifiers, structural decoders) can be rapidly tested and iterated without data bottlenecks.

---

## 2. Data Source: The Human Microprotein Atlas (HMPA)

- **Publication**: *"Comprehensive annotation and analysis of human microproteins by human microprotein atlas platform"*, *Communications Chemistry* (2026) 9:188.
- **Scale**: 617,462 non-redundant human smORFs / microproteins.
- **Repository**: Zenodo (DOI: [10.5281/zenodo.19547725](https://zenodo.org/records/19547725)) & Web Portal ([cuilab.cn/microaf](http://www.cuilab.cn/microaf)).
- **Key Master Files**:
  - `hmpa_structure_alphafold_scores.csv` (139.4 MB): Complete 617,462 entries with sequences, genomic coordinates, AlphaFold2 pLDDT/pTM/PAE, and nearest gene.
  - `smORFs.pro.fa` (26.8 MB): Clean FASTA protein sequences.
  - `smORFs.merge.hg38.bed` (25 MB): hg38 genomic coordinates.
  - `microprotein_subcellular_location.csv` (52.8 MB): LocPro 10-compartment probabilities.
  - `hmpa_pathogenicity_esm2_scores.csv` (13.4 MB): Residue/sequence mutational vulnerability.
  - `hmpa_structure_classification.csv` (43.1 MB): CATH fold categories.

---

## 3. Quota-Aware Daily Representation Workflow

The existing engine in `bindpred/extract_representations.py` is configured to extract:
1. **Sequence-level vectors**: `[2560]` float32 representations saved in Apache Parquet.
2. **Per-residue matrices**: `[L, 2560]` float16 representations saved in compressed NPZ files.

### Batch Execution Plan
```
  [ 617,462 Microproteins ]
              │
              ▼
   [ Daily Quota Batcher ]  ───►  Daily Limit: N sequences / day
              │
              ▼
   [ Biohub ESMC 6B API ]   ───►  Layer 79 representations
              │
              ▼
  [ Checkpointed Store ]   ───►  Parquet + NPZ batch archives
                                 in bindpred/representations/
```

- Each daily run processes a fixed chunk (e.g. batch size of 200–500), automatically saving state in `bindpred/representations/` so no tokens or API quotas are ever wasted.
- As the representation database grows, clustering and exploratory probing can proceed in parallel.
