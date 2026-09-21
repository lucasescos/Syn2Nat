# Quick start

The shortest path from an API key to predictions, plots, and variant scores.
Condensed from the official Quick Start Colab.

## Setup

```python
from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.interpretation import ism
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import pandas as pd

dna_model = dna_client.create(colab_utils.get_api_key())
```

Available output types:

```python
[output.name for output in dna_client.OutputType]
# ATAC, CAGE, DNASE, RNA_SEQ, CHIP_HISTONE, CHIP_TF,
# SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, CONTACT_MAPS, PROCAP
```

## 1. Predict from a raw sequence

Request only the output types and tissues you need — predictions for subsets are
much cheaper than the full set.

```python
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),  # pad to a valid length
    requested_outputs=[dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048'],  # Lung.
)

dnase = output.dnase          # a TrackData object
dnase.values.shape            # (sequence_length, num_tracks)
dnase.metadata                # one row per track
```

Multiple assays and tissues at once:

```python
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    requested_outputs=[dna_client.OutputType.CAGE, dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048', 'UBERON:0000955'],  # Lung, Brain.
)
output.dnase.values.shape  # 2 tracks
output.cage.values.shape   # 4 tracks: 2 tissues x 2 strands (CAGE is stranded)
```

Ontology terms come from UBERON (anatomy), CL (cell ontology), EFO, CLO and NTR.
See [11-recipes.md](11-recipes.md#find-the-ontology-term-for-a-tissue) to look
one up.

## 2. Predict for a reference-genome interval

Load GENCODE annotations, find a gene, resize to a supported length:

```python
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

# MANE Select = one curated transcript per locus.
gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
gtf_transcripts = gene_annotation.filter_to_mane_select_transcript(gtf_transcripts)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)

interval = gene_annotation.get_gene_interval(gtf, gene_symbol='CYP2B6')
interval = interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

output = dna_model.predict_interval(
    interval=interval,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001114'],  # Right liver lobe.
)
```

`.resize()` expands or contracts around the interval's centre using **real
genomic sequence**, not padding. Supported lengths:

```python
dna_client.SUPPORTED_SEQUENCE_LENGTHS.keys()
# SEQUENCE_LENGTH_16KB, SEQUENCE_LENGTH_100KB, SEQUENCE_LENGTH_500KB, SEQUENCE_LENGTH_1MB
```

## 3. Plot it

```python
transcripts = transcript_extractor.extract(interval)

plot_components.plot(
    components=[
        plot_components.TranscriptAnnotation(transcripts),
        plot_components.Tracks(output.rna_seq),
    ],
    interval=output.rna_seq.interval,
)
plt.show()
```

Zoom by resizing the plotting interval only — the data is unchanged:

```python
plot_components.plot(
    components=[
        plot_components.TranscriptAnnotation(transcripts, fig_height=0.1),
        plot_components.Tracks(output.rna_seq),
    ],
    interval=output.rna_seq.interval.resize(2**15),
)
```

Predicted RNA-seq aligns with exons, and stranded tracks show signal on the
strand the gene is transcribed from.

## 4. Predict a variant's effect

```python
variant = genome.Variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',   # may differ from the true reference base
    alternate_bases='C',
)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

variant_output = dna_model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157'],  # Colon - Transverse.
)
# -> VariantOutput with .reference and .alternate, each an Output
```

Overlay REF and ALT:

```python
plot_components.plot(
    [
        plot_components.TranscriptAnnotation(transcript_extractor.extract(interval)),
        plot_components.OverlaidTracks(
            tdata={'REF': variant_output.reference.rna_seq,
                   'ALT': variant_output.alternate.rna_seq},
            colors={'REF': 'dimgrey', 'ALT': 'red'},
        ),
    ],
    interval=variant_output.reference.rna_seq.interval.resize(2**15),
    annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
)
plt.show()
```

For this variant the ALT allele is associated with lower expression and an exon
skipping event in *APOL4* (on the negative strand).

## 5. Score the variant

Scoring makes REF and ALT predictions and aggregates the difference into one
scalar per track.

```python
variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']

variant_scores = dna_model.score_variant(
    interval=interval, variant=variant, variant_scorers=[variant_scorer]
)
scores = variant_scores[0]   # one AnnData per scorer
scores.X.shape               # (num_genes, num_tracks), e.g. (37, 371)
scores.obs.head()            # gene metadata (None for non-gene-centric scorers)
scores.var                   # track metadata
scores.uns['interval'], scores.uns['variant'], scores.uns['variant_scorer']
```

Flatten to a tidy DataFrame — one row per (variant, gene, scorer, ontology):

```python
df = variant_scorers.tidy_scores(variant_scores, match_gene_strand=True)
```

Unlike `predict_variant`, `score_variant` computes scores across all tracks for
the chosen scorer and does not accept `ontology_terms` directly. Filter the
resulting DataFrame instead (e.g. `df[df['ontology_curie'] == 'UBERON:0001157']`).

`raw_score` is the scorer's own output; `quantile_score` ranks it against a
background distribution of common variants so that different scorers are
comparable. See [07-variant-scoring.md](07-variant-scoring.md).

For a non-gene-centric scorer (e.g. `CenterMaskScorer`), `X` has shape
`(1, num_tracks)` and there is no gene metadata.

## 6. In silico mutagenesis

Score every possible SNV in a small window to find the bases that matter.

```python
sequence_interval = genome.Interval('chr20', 3_753_000, 3_753_400)
sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
ism_interval = sequence_interval.resize(256)   # mutate the central 256 bp

dnase_variant_scorer = variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.DNASE,
    width=501,
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
)

variant_scores = dna_model.score_ism_variants(
    interval=sequence_interval,
    ism_interval=ism_interval,
    variant_scorers=[dnase_variant_scorer],
)
len(variant_scores)  # 768 = 256 positions x 3 alternative bases
```

This is expensive — keep the context interval and the ISM window small. Collapse
to one value per variant and plot a sequence logo:

```python
def extract_k562(adata):
  values = adata.X[:, adata.var['ontology_curie'] == 'EFO:0002067']
  assert values.size == 1
  return values.flatten()[0]

ism_result = ism.ism_matrix(
    [extract_k562(x[0]) for x in variant_scores],
    variants=[v[0].uns['variant'] for v in variant_scores],
)  # shape (256, 4)

plot_components.plot(
    [plot_components.SeqLogo(scores=ism_result,
                             scores_interval=ism_interval,
                             ylabel='ISM K562 DNase')],
    interval=ism_interval,
    fig_width=35,
)
plt.show()
```

See [10-interpretation-ism.md](10-interpretation-ism.md).

## 7. Mouse predictions

Pass `organism=dna_client.Organism.MUS_MUSCULUS`. Supported ontology terms
differ between species, and `PROCAP` is not available for mouse.

```python
output = dna_model.predict_interval(
    interval=genome.Interval('chr1', 3_000_000, 3_000_001).resize(
        dna_client.SEQUENCE_LENGTH_1MB),
    organism=dna_client.Organism.MUS_MUSCULUS,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0002048'],
)
```

## Where to go next

- Objects and their methods → [04-core-objects.md](04-core-objects.md)
- All prediction entry points → [05-making-predictions.md](05-making-predictions.md)
- Scoring in depth → [07-variant-scoring.md](07-variant-scoring.md)
- Plotting every modality → [09-visualization.md](09-visualization.md)
