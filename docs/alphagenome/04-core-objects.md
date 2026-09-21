# Core objects

The data types you pass in and get back. Condensed from the Essential Commands
Colab plus the `alphagenome.data` API.

---

## `genome.Interval` — a genomic region

```python
from alphagenome.data import genome

interval = genome.Interval(chromosome='chr1', start=1_000, end=1_010)
```

Fields: `chromosome`, `start`, `end`, `strand` (`'+'`, `'-'`, `'.'`, default
`'.'`), `name`, `info` (free-form dict).

**Indexing is 0-based and half-open.** The interval includes the base at `start`
up to the base at `end - 1`. So `genome.Interval('chr1', 0, 1)` is the first base
of chr1, with width 1, and
`genome.Interval('chr1', 0, 1).overlaps(genome.Interval('chr1', 1, 2))` is
`False`. Default assembly is human hg38.

### Properties and geometry

```python
interval.width          # 10
interval.center()       # centre coordinate (strand-aware by default)
interval.negative_strand
interval.copy()
```

### Resizing and shifting

```python
interval.resize(100)                      # new interval, centred on the original centre
interval.resize_inplace(100)              # mutate in place
interval.shift(500)                       # move along the chromosome
interval.pad(start_pad, end_pad)          # asymmetric expansion
interval.boundary_shift(start_offset=, end_offset=)
interval.truncate(reference_length)       # clip to the reference
interval.within_reference(reference_length)
```

`resize` is how you snap an interval to a model-compatible length. It expands
using **real genomic sequence**, not padding:

```python
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
```

### Comparing intervals

```python
a = genome.Interval('chr1', 1_000, 1_010)
b = genome.Interval('chr1', 1_005, 1_015)

a.overlaps(b)    # True
a.contains(b)    # False
a.intersect(b)   # Interval('chr1', 1005, 1010) or None
```

### Conversion and coverage

```python
genome.Interval.from_str('chr1:100-200:+')
interval.to_interval_dict() / genome.Interval.from_interval_dict(d)
interval.to_pyranges_dict() / genome.Interval.from_pyranges_dict(row)
interval.swap_strand(); interval.as_unstranded()
interval.coverage(intervals, bin_size=1)          # coverage array over this interval
interval.binary_mask(intervals, bin_size=1)       # bool mask
interval.coverage_stranded(...), interval.binary_mask_stranded(...)
```

Module functions: `genome.intersect_intervals`, `genome.union_intervals`,
`genome.merge_overlapping_intervals`.

---

## `genome.Variant` — a genetic variant

```python
variant = genome.Variant(
    chromosome='chr3', position=10_000,
    reference_bases='A', alternate_bases='C',
    name='rs123',            # optional
)
```

**`position` is 1-based**, to match VCF, dbSNP and other public formats.
Internally it is converted to 0-based: `variant.start == variant.position - 1`,
and `variant.end` is 0-based too.

`reference_bases` need **not** match the actual hg38 base at that position. The
model inserts whatever REF and ALT you give it and predicts on those sequences;
`predict_variant` is agnostic to the true reference allele.

### Indels

Indels are left-aligned and share the prefix anchor base (consistent with VCF):

```python
# Insertion: 1-based pos 10,000, anchor 'T' followed by inserted bases
genome.Variant('chr3', 10_000, 'T', 'TCGTCAA')
# Deletion: 1-based pos 10,000, anchor 'A' followed by deleted bases
genome.Variant('chr3', 10_000, 'AGGGATC', 'A')
```

Note: If the first base differs between REF and ALT (e.g. `'AGGGATC'` -> `'C'`),
the variant represents a simultaneous substitution and deletion (delins).


For scoring, AlphaGenome adopts SpliceAI's indel alignment: inserted bases are
summarised by the **maximum** over the inserted segment, and deleted bases are
treated as **zero signal** in the ALT context, so positions stay comparable.

### Properties

```python
variant.start, variant.end              # 0-based
variant.reference_interval              # Interval spanning the REF bases
variant.is_snv, variant.is_insertion, variant.is_deletion
variant.is_indel, variant.is_frameshift, variant.is_structural
variant.reference_overlaps(interval)
variant.alternate_overlaps(interval)
variant.as_truncated_str(max_length=50)
variant.split(anchor)                   # split into two at a position
```

### Parsing from strings

```python
genome.Variant.from_str('chr22:36201698:A>C')                       # DEFAULT
genome.Variant.from_str('chr22_36201698_A_C_b38', genome.VariantFormat.GTEX)
# Also: OPEN_TARGETS ('22_36201698_A_C'), OPEN_TARGETS_BIGQUERY ('22:36201698:A:C'),
#       GNOMAD ('22-36201698-A-C')
```

Standard interval-from-variant idiom:

```python
input_interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
```

---

## `genome.Junction` and `Strand`

`Junction` subclasses `Interval` and adds `k` (the number of reads supporting the
junction — the value `normalize_values` and `k_threshold` act on) plus `.donor` /
`.acceptor` accessors — for `'+'` strand junctions the donor is `start` and the
acceptor is `end`; for `'-'` they swap.

`genome.Strand` is an `IntEnum` with `POSITIVE`, `NEGATIVE`, `UNSTRANDED`, plus
`Strand.from_str('+')`. Constants: `STRAND_POSITIVE = '+'`,
`STRAND_NEGATIVE = '-'`, `STRAND_UNSTRANDED = '.'`.

---

## `track_data.TrackData` — model predictions

Returned for almost every output type. Three parts:

- `tdata.values` — `numpy.ndarray`, shape `(sequence_length / resolution, num_tracks)`
  for 1D tracks (`RNA_SEQ`, `DNASE`, `ATAC`, etc.); or 3D shape
  `(num_bins, num_bins, num_tracks)` for pairwise `CONTACT_MAPS` (2048 bp bins).
- `tdata.metadata` — `pandas.DataFrame`, one row per track
- `tdata.uns` — dict of extra unstructured metadata

Plus `tdata.resolution` (bp per bin) and `tdata.interval`.

### Building one yourself

Metadata must contain at least `name` and `strand`:

```python
from alphagenome.data import track_data
import numpy as np, pandas as pd

values = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]).astype(np.float32)
metadata = pd.DataFrame({'name': ['track1', 'track1', 'track2'],
                         'strand': ['+', '-', '.']})

tdata = track_data.TrackData(values=values, metadata=metadata)

# With genomic context — len(values) must match interval.width / resolution
interval = genome.Interval('chr1', 1_000, 1_004)
tdata = track_data.TrackData(values=values, metadata=metadata,
                             resolution=1, interval=interval)

# 128bp-resolution version of the same 4 values spans 512 bp
interval = genome.Interval('chr1', 1_000, 1_512)
tdata = track_data.TrackData(values=values, metadata=metadata,
                             resolution=128, interval=interval)
```

### Changing resolution

Downsampling sums adjacent values; upsampling repeats while preserving the sum.

```python
tdata = tdata.change_resolution(resolution=2)   # 4 positions -> 2
tdata = tdata.change_resolution(resolution=1)   # back to 4
# Explicit: tdata.downsample(...) / tdata.upsample(...),
# aggregation_type=track_data.AggregationType.SUM (default) or MAX
```

### Filtering by strand

```python
tdata.filter_to_positive_strand()
tdata.filter_to_negative_strand()
tdata.filter_to_unstranded()
tdata.filter_to_nonpositive_strand()   # negative + unstranded
tdata.filter_to_nonnegative_strand()   # positive + unstranded
tdata.filter_to_stranded()             # drop unstranded
```

Use `filter_to_nonpositive_strand()` when looking at a gene on the `-` strand —
it keeps the negative and unstranded tracks and drops the irrelevant `+` ones.

### Resizing, slicing, subsetting

```python
tdata.resize(width=2)        # crop, fixed centre
tdata.resize(width=8)        # pad with zeros
tdata.pad(start_pad, end_pad)

tdata.slice_by_positions(start=2, end=4)
tdata.slice_by_interval(genome.Interval('chr1', 1_002, 1_004), match_resolution=True)

tdata.select_tracks_by_name(names='track1')      # metadata is filtered too
tdata.select_tracks_by_index(idx)
tdata.filter_tracks(mask)                        # boolean mask over tracks
tdata.groupby('biosample_name')                  # -> dict[str, TrackData]
```

`[]` access generalises position and track slicing, numpy-style:

```python
tdata[genome.Interval('chr1', 1_002, 1_004), (tdata.metadata.name == 'track1').values]
tdata[2:4, 0:2]
tdata[2:4, ['track1']]
```

### Other helpers

```python
tdata.num_tracks, tdata.width, tdata.names, tdata.strands, tdata.ontology_terms
tdata.with_name_suffix(' (REF)')
tdata.with_metadata_column('group', 'colon')
tdata.bin_index(relative_position)
tdata.reverse_complement()   # strand-aware; needs a stranded or None interval
tdata.copy()

track_data.concat([t1, t2])                       # along the track axis
track_data.interleave([t1, t2], name_prefixes=['REF', 'ALT'])
```

`TrackData` supports arithmetic, so `alt.rna_seq - ref.rna_seq` gives a
difference track you can plot directly.

---

## `junction_data.JunctionData` — splice junctions

Returned as `output.splice_junctions`. Sparse, 2D-ish: predictions for pairs of
donor/acceptor positions rather than a value per base.

- `.junctions` — array of `genome.Junction`, shape `(num_junctions,)`
- `.values` — shape `(num_junctions, num_tracks)`
- `.metadata`, `.interval`, `.uns`

```python
jd.filter_to_strand('-')          # also filter_to_positive_strand/negative_strand
jd.filter_by_tissue('Colon_Transverse')
jd.filter_by_ontology('CL:0000084')
jd.filter_by_name(name)
jd.filter_tracks(mask)
jd.normalize_values(total_k=10.0)
jd.intersect_with_interval(interval)
jd.num_tracks, jd.names, jd.strands, jd.possible_strands, jd.ontology_terms

junction_data.get_junctions_to_plot(predictions=jd, name=..., strand='-', k_threshold=0.0)
```

Note: strand is a property of a *junction*, not of a track, so
`SPLICE_JUNCTIONS` metadata has half as many rows as the total track count.

---

## `ontology.OntologyTerm` — tissues and cell types

CURIEs (Compact Uniform Resource Identifiers) are short standardised codes such
as `UBERON:0001114` (liver) or `EFO:0002067` (K562). They are sourced from the
IDs in the training data, restricted to **UBERON, CL, CLO, EFO and NTR**
following ENCODE practice. Use EBI's
[Ontology Lookup Service](https://www.ebi.ac.uk/ols4) to explore relationships
between terms.

```python
from alphagenome.data import ontology

term = ontology.from_curie('UBERON:0001114')
terms = ontology.from_curies(['UBERON:0001114', 'EFO:0002067'])
term.ontology_curie   # 'UBERON:0001114'
term.type             # OntologyType.UBERON
```

Anywhere the API takes `ontology_terms`, plain CURIE strings work too.

---

## `anndata.AnnData` — variant scoring output

`score_variant` returns a list of `AnnData`, one per scorer. The format pairs a
score matrix with row and column annotations.

- `.X` — `numpy.ndarray` of shape `(num_genes, num_tracks)` for gene-centric
  scorers, or `(1, num_tracks)` otherwise
- `.obs` — gene metadata `DataFrame` (`None` for non-gene-centric scorers)
- `.var` — track metadata `DataFrame` (same frame as the output metadata)
- `.uns` — provenance: `['variant']`, `['interval']`, `['variant_scorer']`

```python
scores = dna_model.score_variant(...)[0]
scores.X.shape
scores.obs.head()
scores.var
scores.uns['variant']
```

Prefer `variant_scorers.tidy_scores(...)` over reading `AnnData` directly — see
[07-variant-scoring.md](07-variant-scoring.md#tidying-scores).

Constructing one by hand (rarely needed):

```python
import anndata
variant_scores = anndata.AnnData(
    X=np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]]),
    obs=pd.DataFrame({'gene_id': ['ENSG0001', 'ENSG0002', 'ENSG0003']}),
    var=pd.DataFrame({'name': ['track1', 'track2'], 'strand': ['+', '-']}),
)
```

---

## `transcript.Transcript` and `TranscriptExtractor`

```python
from alphagenome.data import gene_annotation, transcript as transcript_utils

gtf = pd.read_feather(HG38_GTF_FEATHER)
gtf_t = gene_annotation.filter_protein_coding(gtf)
gtf_t = gene_annotation.filter_to_mane_select_transcript(gtf_t)

extractor = transcript_utils.TranscriptExtractor(gtf_t)
extractor.cache_transcripts()               # optional speed-up
transcripts = extractor.extract(interval)   # list[Transcript]
```

`Transcript` exposes `exons`, `cds`, `start_codon`, `stop_codon`,
`transcript_id`, `gene_id`, `protein_id`, `uniprot_id`, `info`, plus derived
properties: `chromosome`, `strand`, `is_positive_strand`, `transcript_interval`,
`introns`, `is_coding`, `utr5`, `utr3`, `cds_including_stop_codon`,
`splice_regions`, `splice_donors`, `splice_acceptors`,
`splice_donor_sites`, `splice_acceptor_sites`, `selenocysteines`,
`offset_in_cds(genome_position)`.

Gene name is in `t.info['gene_name']`:

```python
[t.info['gene_name'] for t in transcripts if t.strand == '+']
```

### GTF filters (`gene_annotation`)

```python
gene_annotation.filter_protein_coding(gtf, include_gene_entries=False)
gene_annotation.filter_to_longest_transcript(gtf)
gene_annotation.filter_to_mane_select_transcript(gtf)
gene_annotation.filter_transcript_support_level(gtf, ['1'])
gene_annotation.filter_transcript_type(gtf, (TranscriptType.PROTEIN_CODING,))
gene_annotation.extract_tss(gtf, feature='transcript')
gene_annotation.get_gene_interval(gtf, gene_symbol='APOL4')      # or gene_id=
gene_annotation.get_gene_intervals(gtf, gene_symbols=[...])
gene_annotation.upgrade_annotation_ids(old_ids, new_ids)
```

Transcript annotations come from GENCODE GTFs: hg38 release 46 (human) and mm10
release M23 (mouse).

> **GTF Feather Column Names:** The pre-compiled `.feather` files use **PascalCase**
> for genomic coordinates and structural features (`Chromosome`, `Start`, `End`,
> `Strand`, `Feature`), while metadata attributes remain lowercase (`gene_name`,
> `gene_id`, `gene_type`). Use `gtf[gtf['Feature'] == 'exon']` rather than
> lowercase `feature`. Note that `gene_annotation.get_gene_interval` takes
> `gene_symbol=`, whereas `tidy_scores()` and `Transcript.info` use `gene_name`.


---

## `io.fasta.FastaExtractor`

Read reference sequence directly — used for haplotype construction.

```python
from alphagenome.io import fasta

extractor = fasta.FastaExtractor(fasta_path)
seq = extractor.extract(genome.Interval('chr5', 1_295_112, 1_295_135)).upper()
extractor.sequence_names
extractor.get_length_for_sequence_name('chr1')
fasta.reverse_complement('ACGT')   # 'ACGT' -> 'ACGT' reverse complemented
```
