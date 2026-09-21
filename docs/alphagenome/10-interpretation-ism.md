# Interpretation: in silico mutagenesis

ISM highlights which bases in a region drive a prediction, by systematically
scoring **every possible single-nucleotide variant** in a window and reading off
how much each one changes the output.

## The shape of an ISM run

Two intervals:

- **`sequence_interval`** — the contextual sequence the model sees when making
  each prediction. Must be a supported length.
- **`ism_interval`** — the (much smaller) region whose bases get mutated.

```python
from alphagenome.data import genome
from alphagenome.interpretation import ism
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

# 16 kb of context.
sequence_interval = genome.Interval('chr20', 3_753_000, 3_753_400)
sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)

# Mutate the central 256 bp.
ism_interval = sequence_interval.resize(256)
```

Pick a scorer — what "effect" means for this run:

```python
dnase_variant_scorer = variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.DNASE,
    width=501,
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
)
```

Run it:

```python
variant_scores = dna_model.score_ism_variants(
    interval=sequence_interval,
    ism_interval=ism_interval,
    variant_scorers=[dnase_variant_scorer],
)
len(variant_scores)          # 768 = 256 positions x 3 alternative bases
variant_scores[0][0].X.shape # (1, 305) — non-gene-centric scorer, 305 DNASE tracks
```

The return is `list[list[AnnData]]`: outer index = variant, inner = scorer.

> **This is expensive.** Every base means three scored variants (one per
> alternative allele); the client batches them into RPCs of
> `MAX_ISM_INTERVAL_WIDTH` = 10 bp. Use the shortest
> context sequence that still makes sense and the narrowest ISM window you can.

Extra arguments: `interval_variant` (apply a background variant to the whole
context first), `progress_bar`, `max_workers`, `merge_stranded_gene_tracks`.

## Reducing to one value per variant

A sequence logo needs a single scalar per variant, so pick a track or average
across several:

```python
def extract_k562(adata):
    values = adata.X[:, adata.var['ontology_curie'] == 'EFO:0002067']
    assert values.size == 1
    return values.flatten()[0]

ism_result = ism.ism_matrix(
    [extract_k562(x[0]) for x in variant_scores],
    variants=[v[0].uns['variant'] for v in variant_scores],
)
ism_result.shape   # (256, 4) — one score per position per base
```

`ism.ism_matrix(variant_scores, variants, interval=None, multiply_by_sequence=True,
vocabulary='ACGT', require_fully_filled=True)`

With `multiply_by_sequence=True` (the default), the output array is non-zero
only at the bases present in the reference sequence — the standard convention for
contribution-score logos.

`ism.ism_variants(interval, sequence, vocabulary='ACGT', skip_n=False)` builds the
list of all single-nucleotide variants for an interval, if you want to drive the
scoring loop yourself.

## Plotting

```python
plot_components.plot(
    [
        plot_components.SeqLogo(
            scores=ism_result,
            scores_interval=ism_interval,
            ylabel='ISM K562 DNase',
        )
    ],
    interval=ism_interval,
    fig_width=35,
)
```

Tall letters mark positions where changing the base most changes the prediction.
Positive letter heights represent activating effects (signal increases when the
reference base is present), while negative heights represent repressive effects.
In the worked example, positions ~225–240 have the strongest effect on predicted
nearby DNase signal in K562 cells.

**Reverse complements:** Because transcription factors can bind in either
orientation on double-stranded DNA, always check the reverse complement sequence
when matching sequence motifs (e.g. TATA box `TATAAA` <-> `TTTATA`; PolyA signal
`AATAAA` <-> `TTTATT`).

A wide `fig_width` is usually necessary for the letters to be readable. Note that
`SeqLogo` draws **nothing at all** — a silent blank panel — if the plotted
interval is wider than `max_width` (default 1000); raise `max_width` or narrow
the plot interval.

## Downstream: finding motifs

Contribution scores can be used to systematically discover motifs that matter for
different modalities and cell types, identify the transcription factors binding
them, and map motif instances across the genome:

- [tfmodisco-lite](https://github.com/jmschrei/tfmodisco-lite/) — motif discovery
  from contribution scores
- [tangermeme](https://github.com/jmschrei/tangermeme) — motif analysis toolkit
- [tomtom](https://meme-suite.org/meme/tools/tomtom) — match discovered motifs
  against known databases

## A cheaper alternative

If you only want to know whether a *specific* set of candidate sequences matters,
scoring those variants directly with `score_variants` is far cheaper than a full
ISM sweep. See the shuffled-background approach in
[11-recipes.md](11-recipes.md#3-compare-real-variants-against-shuffled-backgrounds),
which asks "does this exact sequence do something a random sequence of the same
length wouldn't?" using a few dozen calls instead of hundreds.
