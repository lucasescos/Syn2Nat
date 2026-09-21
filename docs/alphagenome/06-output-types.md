# Output types and metadata

AlphaGenome returns predictions for **11 output types**. Counts below are for the
**human** model. For precise definitions, units and normalisation, see the
Methods section of the AlphaGenome paper.

## The 11 output types

| `OutputType` | Description | Units | Resolution | Biosamples | Tracks |
| :--- | :--- | :--- | :--- | ---: | ---: |
| `RNA_SEQ` | RNA expression by RNA-seq. Mixture of PolyA+ and Total RNA assays; some tracks stranded. | Normalized read signal | 1 bp | 285 | 667 |
| `CAGE` | Expression at TSSs by Cap Analysis Gene Expression. | Normalized read signal | 1 bp | 264 | 546 |
| `PROCAP` | Expression at TSSs by Precision Run-On sequencing and capping. | Normalized read signal | 1 bp | 6 | 12 |
| `DNASE` | Chromatin accessibility by DNase I hypersensitive sites sequencing. | Normalized insertion signal | 1 bp | 305 | 305 |
| `ATAC` | Chromatin accessibility by ATAC-seq. | Normalized insertion signal | 1 bp | 167 | 167 |
| `CHIP_HISTONE` | Histone modification abundance by ChIP-seq, [24 markers](https://www.encodeproject.org/chip-seq/histone/) (e.g. H3K27ac). | Fold-change over control, summed per 128 bp bin | 128 bp | 219 | 1116 |
| `CHIP_TF` | DNA-bound transcription factor abundance by ChIP-seq, [43 proteins](https://www.encodeproject.org/chip-seq/transcription_factor/). | Fold-change over control, summed per 128 bp bin | 128 bp | 163 | 1617 |
| `SPLICE_SITES` | Donor/acceptor splice site probability, both strands. | Predicted probability | 1 bp | NA | 4 |
| `SPLICE_JUNCTIONS` | Spliced read counts per junction, for all pairings of at most 512 donors and 512 acceptors per strand in the interval. | Normalized junction signal | 1 bp | 282 | 734 |
| `SPLICE_SITE_USAGE` | Fraction of transcripts using a splice site, among reads spanning it. | Fraction | 1 bp | 282 | 734 |
| `CONTACT_MAPS` | Relative frequency of physical contact between pairwise positions (symmetric), from Micro-C and Hi-C. Coarse-grained and normalized by removing off-diagonal power-law decay. | Log-fold over distance-based expectation | 2048 bp | 12 | 28 |

Notes:

- `SPLICE_SITES` is **tissue-agnostic** — 4 tracks (donor/acceptor × strand), and
  the `ontology_terms` filter does not apply to it.
- `SPLICE_JUNCTIONS` metadata has **half** the rows in the table, because strand
  is a property of a junction rather than of a track.
- `CHIP_HISTONE` and `CHIP_TF` come back at 128 bp resolution. `CONTACT_MAPS`
  comes back at 2048 bp resolution as a pairwise interaction matrix per track,
  stored in `TrackData.values` with shape `(num_bins, num_bins, num_tracks)`.
  Plot intervals are adjusted to a compatible multiple automatically.
- `PROCAP` has no mouse tracks, and polyadenylation scoring (`PA_QTL`) is human only.
- `SPLICE_SITE_USAGE` predicts the tissue-specific fraction of spanning transcripts
  using each splice site, whereas `SPLICE_SITES` predicts intrinsic donor/acceptor
  sequence probabilities independent of tissue.

```python
[o.name for o in dna_client.OutputType]
```

## Track metadata

```python
output_metadata = dna_model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
output_metadata.rna_seq       # a pandas DataFrame, one row per track
```

Core columns, present for every output type:

- **`name`** — track name, e.g. `CL:0000047 polyA plus RNA-seq`
- **`strand`** — `+`, `-`, or `.` (unstranded)
- **`ontology_curie`** — biosample ontology ID, e.g. `CL:0000100`
- **`biosample_name`** — plain-text biosample, e.g. `motor neuron`

Additional columns on some types:

- **`gtex_tissue`** — on `rna_seq` and `splice_sites`, populated for tracks
  matching tissues sampled in the [GTEx project](https://gtexportal.org/home/samplingSitePage).
  One exception: for 'Brain - Cerebellar hemisphere' AlphaGenome uses
  `UBERON:0002245` (Uberon's ID for cerebellar hemisphere) rather than the
  `UBERON:0002037` given in the GTEx documentation.
- **`transcription_factor`** — on `chip_tf`, e.g. `CTCF`
- **`histone_mark`** — on `chip_histone`, e.g. `H3K4ME3`
- **`biosample_type`** — e.g. `tissue`, `primary cell`
- **`Assay title`** — assay subtype, e.g. `total RNA-seq`, `polyA plus RNA-seq`

### Browsing metadata

```python
all_metadata = dna_model.output_metadata(
    dna_client.Organism.HOMO_SAPIENS
).concatenate()

# All colon-related ChIP-histone tracks
output_metadata.chip_histone[
    output_metadata.chip_histone['biosample_name'].str.contains('colon')
]
```

Track counts per output type, human vs mouse:

```python
human = (dna_model.output_metadata(dna_client.Organism.HOMO_SAPIENS)
         .concatenate().groupby('output_type').size().rename('# Human tracks'))
mouse = (dna_model.output_metadata(dna_client.Organism.MUS_MUSCULUS)
         .concatenate().groupby('output_type').size().rename('# Mouse tracks'))
pd.concat([human, mouse], axis=1).astype(pd.Int64Dtype())
```

## Strandedness

DNA is double-stranded; by convention one molecule is the forward/positive
strand (5'→3') and the other the reverse/negative strand (3'→5').

- **Unstranded assays** do not distinguish which strand a measurement came from —
  ATAC-seq, for example, gives unstranded accessibility. One track per biosample,
  marked `.`.
- **Stranded (strand-specific) assays** annotate each measurement with its
  strand. Two tracks per biosample, `+` and `-`. Important for transcriptional
  assays, e.g. to separate two transcripts sharing a TSS on opposite strands.

Not all RNA-seq is stranded — older experiments often are not, and GTEx RNA-seq
is unstranded. Background reading:
[strand-specific protocols](https://www.ecseq.com/support/ngs/how-do-strand-specific-sequencing-protocols-work).

Convenience filters on `TrackData`: `filter_to_positive_strand()`,
`filter_to_negative_strand()`, `filter_to_unstranded()`,
`filter_to_nonpositive_strand()`, `filter_to_nonnegative_strand()`,
`filter_to_stranded()`.

## How many tracks, and what do they mean?

Track counts range from 4 (`SPLICE_SITES`) to 1617 (`CHIP_TF`) per output type —
see the table above. Each track corresponds to a particular cell type or tissue,
plus other properties such as strand or, for `CHIP_TF`, a specific transcription
factor. To map a tissue name to a CURIE, see
[11-recipes.md](11-recipes.md#find-the-ontology-term-for-a-tissue).
