# Variant scoring

A genomic variant is a difference between an individual's genome and the
reference. Most have no appreciable effect; scoring is how you pick out the ones
that do. AlphaGenome treats a variant as a pair of sequences — reference (`REF`)
and alternate (`ALT`) — and estimates the effect by comparing predictions for the
two across modalities.

## How a score is computed

### 1. Predict REF and ALT for one modality

Generate predictions for both alleles, restricted to a modality of interest
(`RNA_SEQ`, `ATAC`, …). The inputs are REF and ALT sequences whose interval
contains the variant.

### 2. (Indels only) Align ALT to REF coordinates

The ALT profile is aligned back to the REF coordinate space: inserted bases are
summarised by the **maximum** over the inserted segment; deleted bases are
treated as **zero signal** in the ALT context. This keeps positional comparisons
consistent. (Same strategy as SpliceAI.)

### 3. Apply a spatial mask

The mask defines the region of interest inside the interval. It is either
**centred on the variant** (a window of fixed width) or **gene-shaped** — gene
body, exons, or TSS, from a GTF. Values outside the mask are discarded.

### 4. Aggregate spatially and take ALT − REF

- Reduce along the spatial axis (`mean`, `sum`, …).
- Optionally scale (log or L2 transform).
- Take the difference `ALT − REF`.

The result is **one scalar per track**.

Aggregation logic lives in `variant_scorers.AggregationType`. Names read
right-to-left in order of application: `DIFF_SUM_LOG2` applies a log2 transform,
then a sum, then the ALT − REF difference. Some options apply the same steps in a
different order. Every aggregation type returns a single scalar per track.

Available: `DIFF_MEAN`, `DIFF_SUM`, `DIFF_SUM_LOG2`, `DIFF_LOG2_SUM`, `L2_DIFF`,
`L2_DIFF_LOG1P`, `ACTIVE_MEAN`, `ACTIVE_SUM`.

### 5. (Optional) Aggregate across tracks

Afterwards you can take a mean/max/sum across all tracks or a subset, or just
pick the single track for the sample you care about.

## Recommended scorers

Accessible as a dict, or via `get_recommended_scorers(organism)`:

```python
from alphagenome.models import variant_scorers

variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())   # all of them
variant_scorers.get_recommended_scorers(
    dna_client.Organism.HOMO_SAPIENS.to_proto())   # takes the PROTO enum
```

> **Note on `get_recommended_scorers`:** This function requires the protobuf enum
> (`.to_proto()`). Passing `dna_client.Organism.HOMO_SAPIENS` directly will
> silently return an empty list `[]` without raising an error.


### Differential scorers

| Key | Scorer | Comparison · mask · aggregation |
| :--- | :--- | :--- |
| `RNA_SEQ` | `GeneMaskLFCScorer(RNA_SEQ)` | RNA coverage; exons of a gene; log-fold change `log(mean(ALT)+0.001) − log(mean(REF)+0.001)` |
| `CAGE` | `CenterMaskScorer(CAGE, 501, DIFF_LOG2_SUM)` | CAGE coverage; 501 bp centred on variant; `log2[(sum(ALT)+1)/(sum(REF)+1)]` |
| `PROCAP` | `CenterMaskScorer(PROCAP, 501, DIFF_LOG2_SUM)` | as above, PRO-cap |
| `ATAC` | `CenterMaskScorer(ATAC, 501, DIFF_LOG2_SUM)` | ATAC coverage; 501 bp window; log2 ratio of sums |
| `DNASE` | `CenterMaskScorer(DNASE, 501, DIFF_LOG2_SUM)` | DNase coverage; 501 bp window; log2 ratio of sums |
| `CHIP_TF` | `CenterMaskScorer(CHIP_TF, 501, DIFF_LOG2_SUM)` | TF binding intensity; 501 bp window; log2 ratio of sums |
| `CHIP_HISTONE` | `CenterMaskScorer(CHIP_HISTONE, 2001, DIFF_LOG2_SUM)` | histone marks; **2001 bp** window; log2 ratio of sums |
| `SPLICE_SITES` | `GeneMaskSplicingScorer(SPLICE_SITES)` | splice-site class probabilities; gene body; `max(|ALT − REF|)` |
| `SPLICE_SITE_USAGE` | `GeneMaskSplicingScorer(SPLICE_SITE_USAGE)` | splice site usage; gene body; `max(|ALT − REF|)` |
| `SPLICE_JUNCTIONS` | `SpliceJunctionScorer()` | paired junction counts; top-k splice sites for the gene (annotated and predicted); `max(|log(ALT) − log(REF)|)` |
| `POLYADENYLATION` | `PolyadenylationScorer()` | RNA coverage; 400 bp windows around 3' cleavage junctions; max absolute log-fold change of distal/proximal PAS usage ratios over all splits. Follows Borzoi's paQTL method. **Human only.** |
| `CONTACT_MAPS` | `ContactMapScorer()` | contact frequencies; 1 Mb window centred on the variant; mean absolute difference over all interactions involving the variant's bin |

### Active-allele scorers

These capture the **absolute activity level** of one allele rather than the
change between alleles: the maximum of the aggregated REF and ALT signals over
the masked window or gene region.

| Key | Definition |
| :--- | :--- |
| `RNA_SEQ_ACTIVE` | `max(mean(ALT), mean(REF))` across exons of a gene |
| `CAGE_ACTIVE`, `PROCAP_ACTIVE` | `max(sum(ALT), sum(REF))` in a 501 bp window |
| `ATAC_ACTIVE`, `DNASE_ACTIVE` | `max(sum(ALT), sum(REF))` in a 501 bp window |
| `CHIP_TF_ACTIVE` | `max(sum(ALT), sum(REF))` in a 501 bp window |
| `CHIP_HISTONE_ACTIVE` | `max(sum(ALT), sum(REF))` in a 2001 bp window |

## Scorer classes

Seven base scorer types (`BaseVariantScorer`): `CENTER_MASK`, `CONTACT_MAP`,
`GENE_MASK_LFC`, `GENE_MASK_ACTIVE`, `GENE_MASK_SPLICING`, `PA_QTL`,
`SPLICE_JUNCTION`. `BaseVariantScorer` is itself an enum; the seven **scorer
classes** below are frozen dataclasses, each exposing `.name`, `.is_signed`,
`.requested_output`, and `.base_variant_scorer`.

```python
variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.DNASE,
    width=501,
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
)
variant_scorers.GeneMaskLFCScorer(requested_output=OutputType.RNA_SEQ)
variant_scorers.GeneMaskActiveScorer(requested_output=OutputType.RNA_SEQ)
variant_scorers.GeneMaskSplicingScorer(requested_output=OutputType.SPLICE_SITES, width=None)
variant_scorers.SpliceJunctionScorer()
variant_scorers.PolyadenylationScorer()
variant_scorers.ContactMapScorer()
```

### Constraints

| Base scorer | Supported output types | Widths |
| :--- | :--- | :--- |
| `CENTER_MASK` | ATAC, CAGE, DNASE, PROCAP, RNA_SEQ, CHIP_HISTONE, CHIP_TF, SPLICE_SITES, SPLICE_SITE_USAGE | `None`, 501, 2001, 10001, 100001, 200001 |
| `CONTACT_MAP` | CONTACT_MAPS | — |
| `GENE_MASK_LFC` | ATAC, CAGE, DNASE, PROCAP, RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE | — |
| `GENE_MASK_ACTIVE` | same as `GENE_MASK_LFC` | — |
| `GENE_MASK_SPLICING` | SPLICE_SITES, SPLICE_SITE_USAGE | `None`, 101, 1001, 10001 |
| `PA_QTL` | — | human only |
| `SPLICE_JUNCTION` | — | — |

`CENTER_MASK` accepts every `AggregationType`. Check organism support with
`variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]`.

**Which scorer for which modality?** In practice most strategies work for any
modality; the table above is the recommendation based on DeepMind's evaluations.
Custom scoring strategies are not supported by the API, but since scoring is
just aggregation of REF and ALT track predictions, you can write your own
aggregation over `predict_variant` output.

## Running a scorer

```python
scores = dna_model.score_variant(
    interval=variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
    variant=variant,
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
    organism=dna_client.Organism.HOMO_SAPIENS,
)
# list[AnnData], one per scorer
```

Note: `score_variant` computes scores for all tracks associated with the scorer
and does **not** accept `ontology_terms`. Filter the output with `tidy_scores`
by `ontology_curie` or `biosample_name`.


Many variants at once:

```python
scores = dna_model.score_variants(
    intervals=intervals, variants=variants,
    variant_scorers=selected_scorers,
    max_workers=2,
)
# list[list[AnnData]]: outer = variants, inner = scorers
```

Maximum 20 scorers per request (`MAX_VARIANT_SCORERS_PER_REQUEST`). Use 1 Mb
sequence length for gene-level scoring — it gives the best predictions.

## Tidying scores

`tidy_scores` flattens `AnnData` into a long DataFrame, one row per
(variant, gene, scorer, ontology). This is the recommended way to work with
scores.

```python
df = variant_scorers.tidy_scores(scores, match_gene_strand=True)
df.to_csv('variant_scores.csv', index=False)
```

Accepts either a flat sequence of `AnnData` (from `score_variant`) or a nested
sequence (from `score_variants`); results from multiple scorers and variants are
concatenated, taking the union of applicable columns.

- `match_gene_strand=True` drops rows where the track strand doesn't match the
  gene's strand, for gene-centric scorers.
- `include_extended_metadata=False` keeps only `track_name` and `track_strand`
  instead of the full biosample metadata.
- `tidy_anndata(adata, ...)` does the same for a single `AnnData`.

### Columns

`variant_id`, `scored_interval`, `gene_id` (ENSEMBL, no version), `gene_name`
(HGNC symbol), `gene_type`, `gene_strand`, `output_type`, `variant_scorer` (or
`interval_scorer`), `track_name`, `track_strand`, `ontology_curie`,
`gtex_tissue`, `Assay title`, `biosample_name`, `biosample_type`,
`transcription_factor`, `histone_mark`, `raw_score`, `quantile_score`.

Gene columns are `None` for non-gene-centric scorers.

## `raw_score` vs `quantile_score`

**`raw_score`** is the scorer's direct output. Different tracks and modalities
produce scores on different scales — the splice-site-usage scorer returns values
in [0, 1], while the RNA-seq scorer is unbounded and signed — so raw scores are
not comparable across scorers.

**`quantile_score`** places a raw score in an empirical background distribution.
The background is built from scores for **common variants** (MAF > 0.01 in any
gnomAD v3 population), per scorer and per track. A quantile score of 0.99 means
the raw score sits at the 99th percentile of common variants.

- For **signed** scorers the [0, 1] quantile is linearly rescaled to [−1, 1] so
  direction is preserved: 0th percentile (most negative) → −1, 50th → 0,
  100th → +1.
- The magnitude never exceeds 0.999990 (or −0.999990), because ~300K variants
  were used to compute the quantiles.
- Quantile scores are only available for the **recommended** scorers.

**Recommended use:** treat the quantile score as the indicator of whether a raw
score is unusually large, and the raw score as the measure of effect magnitude
for that scorer and track.

> **Caution — The "High Quantile + Low Raw Score" Trap:**
> In genes with low or near-zero baseline expression in a biosample, empirical
> variance across common variants is minimal. Even negligible numerical noise can
> produce an extreme quantile score (e.g. `|quantile_score| > 0.999`) while the
> actual biological effect is nonexistent (`|raw_score| < 0.1`).
>
> General rules of thumb for RNA-seq log-fold change raw scores:
> - `|raw_score| < 0.1`: Background noise / No significant molecular effect.
> - `0.1 <= |raw_score| < 0.5`: Subtle or weak change (1.07x to 1.4x fold change).
> - `0.5 <= |raw_score| < 1.0`: Moderate effect (1.4x to 2x fold change).
> - `|raw_score| >= 1.0`: Strong molecular effect (>2x fold change).
>
> Always cross-reference high quantile scores against the active allele score
> (`RNA_SEQ_ACTIVE`) or track prediction to verify that the target gene is
> genuinely transcribed in that tissue.


## Interval scorers

A parallel, less-used API scores an interval without a variant:

```python
from alphagenome.models import interval_scorers

interval_scorers.RECOMMENDED_INTERVAL_SCORERS['RNA_SEQ']
# GeneMaskScorer(RNA_SEQ, width=200001, aggregation_type=IntervalAggregationType.MEAN)

dna_model.score_interval(interval=interval, interval_scorers=[scorer])
dna_model.score_intervals(intervals=intervals, interval_scorers=[scorer])
```

`tidy_scores` handles interval-scorer output too, producing an
`interval_scorer` column instead of `variant_scorer`.

## See also

- Splicing-specific scoring and the merged splicing score →
  [08-splicing.md](08-splicing.md)
- Scoring a VCF in bulk → [11-recipes.md](11-recipes.md#batch-score-a-vcf)
- ISM (scoring every possible SNV) →
  [10-interpretation-ism.md](10-interpretation-ism.md)
