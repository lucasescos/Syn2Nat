# Splicing

AlphaGenome predicts splicing through three output types, and the recommended
way to ask "does this variant cause aberrant splicing?" is to combine three
scorers into a single **merged splicing score**.

## The three splicing outputs

| Output type | What it is |
| :--- | :--- |
| `SPLICE_SITES` | Per-base probability that a position is a donor or acceptor splice site, per strand. Tissue-agnostic — 4 tracks total. |
| `SPLICE_SITE_USAGE` | Fraction of transcripts spanning a site that actually use it. Tissue-specific. |
| `SPLICE_JUNCTIONS` | Predicted spliced read counts for donor/acceptor pairs. Returned as `JunctionData`, not `TrackData`. |

The model predicts acceptor sites near the beginnings of exons and donor sites
near the ends; splice-site *usage* has peaks on both sides of exons. `RNA_SEQ`
also carries splicing information, so plotting it alongside is useful.

## The merged splicing score

This is the approach used in the AlphaGenome paper to score ClinVar variants for
missplicing, and the recommended method for assessing aberrant splicing.

Three scorers, each aggregated by taking the **maximum absolute score across all
tissues and genes**:

- **Splice sites** — changes in splice site class assignment probabilities
  (donor, acceptor) between ALT and REF.
- **Splice site usage** — changes in the relative usage of splice sites.
- **Splice junctions** — log-fold changes in predicted splice junction counts.

Combined as:

```
alphagenome_splicing = max(splice_sites)
                     + max(splice_site_usage)
                     + max(splice_junctions) / 5
```

Full weight goes to changes in splice site identity and usage; junction-level
changes get a weight of 0.2 because their magnitude is larger (the junction
scorer computes a log fold change).

### Interpretation

In theory the score is unbounded from 0 to infinity, but empirically most
variants fall in **[0, 6]**:
- `< 0.2`: Baseline / background alternative splicing (typically neutral).
- `0.2 - 1.0`: Moderate or ambiguous splicing alteration.
- `> 1.0`: Strong predicted aberrant splicing disruption (e.g. canonical site
  loss, robust exon skipping, or pseudo-exon activation).


## Scoring a variant for splicing

```python
from alphagenome import colab_utils
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
import pandas as pd

dna_model = dna_client.create(colab_utils.get_api_key())

splicing_scorers = [
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITES'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_SITE_USAGE'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_JUNCTIONS'],
]
for scorer in splicing_scorers:
    print(f'{scorer.name} (signed={scorer.is_signed})')
```

Under the hood these are `GeneMaskSplicingScorer(SPLICE_SITES)`,
`GeneMaskSplicingScorer(SPLICE_SITE_USAGE)` and `SpliceJunctionScorer()`.

```python
variant = genome.Variant(
    chromosome='chr13', position=32316462,
    reference_bases='T', alternate_bases='G',   # near a BRCA2 splice site
)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=splicing_scorers,
    organism=dna_client.Organism.HOMO_SAPIENS,
)
df_scores = variant_scorers.tidy_scores([scores])
```

Use 1 Mb — gene-level scoring is best at full context length.

### Computing the merged score

```python
def compute_merged_splicing_score(df: pd.DataFrame) -> pd.DataFrame:
    """max raw_score per variant per output type, combined with paper weights."""
    df['variant_id'] = df['variant_id'].map(str)
    max_scores = (
        df.groupby(['variant_id', 'output_type'])['raw_score']
        .max()
        .reset_index()
        .pivot(index='variant_id', columns='output_type', values='raw_score')
        .fillna(0.0)
    )
    max_scores['alphagenome_splicing'] = (
        max_scores.get('SPLICE_SITES', 0.0)
        + max_scores.get('SPLICE_SITE_USAGE', 0.0)
        + max_scores.get('SPLICE_JUNCTIONS', 0.0) / 5.0
    )
    return max_scores.reset_index()

merged_scores = compute_merged_splicing_score(df_scores)
```

### A batch of variants

```python
variants = [
    genome.Variant('chr13', 32316462, 'T', 'G'),
    genome.Variant('chr17', 43092919, 'G', 'A'),
    genome.Variant('chr7',  117559593, 'T', 'A'),
]

all_scores = []
for variant in variants:
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    all_scores.append(dna_model.score_variant(
        interval=interval, variant=variant,
        variant_scorers=splicing_scorers,
        organism=dna_client.Organism.HOMO_SAPIENS,
    ))

merged_scores = compute_merged_splicing_score(
    variant_scorers.tidy_scores(all_scores)
)
```

## Deriving PSI values

AlphaGenome's splicing scorers predict splice site usage and junction counts but
do **not** output Percent Spliced In directly. You can derive PSI3/PSI5 from
predicted junction counts.

- **PSI5** — usage of a 5' splice site (donor) relative to alternative donors:
  normalise the donor × acceptor count matrix **row-wise**.
- **PSI3** — usage of a 3' splice site (acceptor) relative to alternative
  acceptors: normalise **column-wise**.

```python
output = dna_model.predict_variant(
    interval=variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
    variant=variant,
    requested_outputs=[dna_client.OutputType.SPLICE_JUNCTIONS],
    ontology_terms=['CL:0000084'],   # T cell
)

output.reference.splice_junctions.junctions.shape   # (num_junctions,)
output.reference.splice_junctions.values.shape      # (num_junctions, num_tracks)
output.reference.splice_junctions.metadata
```

Two tracks come back for one ontology term because junction counts derive from
RNA-seq and two assay types are available:

- **polyA RNA-seq** — enriches for fully processed transcripts, better for
  mature mRNA splicing.
- **Total RNA-seq** — captures pre-mRNA and nascent transcripts, so it can look
  "messier" because splicing is still in progress.

```python
def compute_psi_from_raw(junction_data, track_idx):
    data = []
    for i, j in enumerate(junction_data.junctions):
        # Uses built-in strand-aware donor and acceptor properties on genome.Junction
        data.append({'donor': j.donor, 'acceptor': j.acceptor,
                     'count': junction_data.values[i, track_idx]})

    df = pd.DataFrame(data)
    matrix = df.pivot_table(index='donor', columns='acceptor',
                            values='count', aggfunc='sum', fill_value=0)

    psi5 = matrix.div(matrix.sum(axis=1), axis=0).fillna(0)  # row-wise (donor PSI)
    psi3 = matrix.div(matrix.sum(axis=0), axis=1).fillna(0)  # column-wise (acceptor PSI)
    return psi5, psi3
```

Pick a track by assay title:

```python
titles = output.reference.splice_junctions.metadata['Assay title'].astype(str).tolist()
polya_idx = next((i for i, t in enumerate(titles) if 'polya' in t.lower()), None)
total_idx = next((i for i, t in enumerate(titles) if 'total' in t.lower()), None)

psi5, psi3 = compute_psi_from_raw(output.reference.splice_junctions, polya_idx)
```

Reading the results: **PSI3** is the proportion of total splicing activity at a
given acceptor contributed by a given donor; **PSI5** is the proportion of
activity at a given donor directed to a given acceptor.

## Visualising splicing

Sashimi plots draw junctions as arcs whose thickness scales with the predicted
counts. The numbers annotated above the arcs are **predicted normalised read
counts** per donor/acceptor pair, not PSI values.

> **Note on `Sashimi`:** `plot_components.Sashimi` does **not** accept a `strand`
> argument directly. Filter the `JunctionData` input by strand beforehand (e.g.
> `output.reference.splice_junctions.filter_to_strand('-')`).

```python
from alphagenome.visualization import plot_components

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(transcript_extractor.extract(interval)),
        plot_components.Sashimi(
            output.reference.splice_junctions.filter_to_strand('-'),
            ylabel_template='Reference {biosample_name} ({strand})\n{name}',
            color='skyblue',
        ),
        plot_components.Sashimi(
            output.alternate.splice_junctions.filter_to_strand('-'),
            ylabel_template='Alternate {biosample_name} ({strand})\n{name}',
            color='red',
        ),
    ],
    annotations=[plot_components.VariantAnnotation([variant])],
    interval=interval.resize(40000),
)
```

For splicing work, consider extracting more than just MANE/protein-coding
transcripts — though including everything makes the annotation crowded:

```python
transcript_extractor_all     = transcript_utils.TranscriptExtractor(gtf)
transcript_extractor_protein = transcript_utils.TranscriptExtractor(
    gene_annotation.filter_protein_coding(gtf))
transcript_extractor_mane    = transcript_utils.TranscriptExtractor(
    gene_annotation.filter_to_mane_select_transcript(gtf))
```

Filter junctions before plotting to keep things readable:

```python
output.reference.splice_junctions.filter_to_strand('-').filter_by_tissue('Colon_Transverse')
```

A full REF/ALT splicing panel — sashimi, RNA-seq, splice sites and usage —
appears in [09-visualization.md](09-visualization.md#splicing).
