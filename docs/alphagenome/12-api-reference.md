# API reference

Complete public surface of the `alphagenome` package, condensed from the
auto-generated reference and the source. Signatures show defaults; `*` marks the
start of keyword-only arguments.

- [`alphagenome.models`](#alphagenomemodels)
- [`alphagenome.data`](#alphagenomedata)
- [`alphagenome.visualization`](#alphagenomevisualization)
- [`alphagenome.interpretation`](#alphagenomeinterpretation)
- [`alphagenome.atlas`](#alphagenomeatlas)
- [`alphagenome.io`](#alphagenomeio)
- [`alphagenome.colab_utils`](#alphagenomecolab_utils)

---

## `alphagenome.models`

### `models.dna_client`

**Constants**

```python
SEQUENCE_LENGTH_16KB  = 2**14
SEQUENCE_LENGTH_100KB = 2**17
SEQUENCE_LENGTH_500KB = 2**19
SEQUENCE_LENGTH_1MB   = 2**20
SUPPORTED_SEQUENCE_LENGTHS          # dict name -> value
MAX_VARIANT_SCORERS_PER_REQUEST = 20
MAX_ISM_INTERVAL_WIDTH = 10
```

**Functions**

```python
create(api_key: str, *, model_version=None, timeout=None, address=None) -> DnaClient
validate_sequence_length(length: int)          # raises if unsupported
construct_output_metadata(responses) -> OutputMetadata
retry_rpc(function, *, max_attempts=5, initial_backoff=1.25, backoff_multiplier=1.5,
          retry_status_codes=frozenset({RESOURCE_EXHAUSTED, UNAVAILABLE}), jitter=0.2)
```

**`class DnaClient`** (implements `dna_model.DnaModel`)

```python
predict_sequence(sequence, *, organism=Organism.HOMO_SAPIENS, requested_outputs,
                 ontology_terms, interval=None) -> Output
predict_interval(interval, *, organism=..., requested_outputs, ontology_terms) -> Output
predict_variant(interval, variant, *, organism=..., requested_outputs,
                ontology_terms) -> VariantOutput
score_interval(interval, interval_scorers=(), *, organism=...,
               merge_stranded_gene_tracks=True) -> list[AnnData]
score_variant(interval, variant, variant_scorers=(), *, organism=...,
              merge_stranded_gene_tracks=True) -> list[AnnData]
score_ism_variants(interval, ism_interval, variant_scorers=(), *, organism=...,
                   interval_variant=None, progress_bar=True, max_workers=5,
                   merge_stranded_gene_tracks=True) -> list[list[AnnData]]
output_metadata(organism=Organism.HOMO_SAPIENS) -> OutputMetadata
```

Batch methods inherited from `dna_model.DnaModel` (each takes
`progress_bar=True`, `max_workers=5`):

```python
predict_sequences(sequences, *, organism=..., requested_outputs, ontology_terms,
                  intervals=None) -> list[Output]
predict_intervals(intervals, *, ...) -> list[Output]
predict_variants(intervals, variants, *, ...) -> list[VariantOutput]
score_intervals(intervals, interval_scorers=(), *, ...) -> list[list[AnnData]]
score_variants(intervals, variants, variant_scorers=(), *, ...) -> list[list[AnnData]]
```

`intervals` on `predict_variants` / `score_variants` may be a single `Interval`
shared by all variants, or one per variant. `DEFAULT_MAX_WORKERS = 5`.

**`enum ModelVersion`** — `ALL_FOLDS`, `FOLD_0`, `FOLD_1`, `FOLD_2`, `FOLD_3`.

**`enum Organism`** — `HOMO_SAPIENS`, `MUS_MUSCULUS`. Method: `.to_proto()`.

### `models.dna_output`

**`enum OutputType`** — `ATAC`, `CAGE`, `DNASE`, `RNA_SEQ`, `CHIP_HISTONE`,
`CHIP_TF`, `SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`,
`CONTACT_MAPS`, `PROCAP`. Method: `.to_proto()`.

**`class Output`** (frozen dataclass) — one field per output type; `None` if not
requested.

```python
.atac .cage .dnase .rna_seq .chip_histone .chip_tf
.splice_sites .splice_site_usage .contact_maps .procap   # TrackData | None
.splice_junctions                                        # JunctionData | None

.get(output_type) -> TrackData | JunctionData | None
.filter_to_strand(strand: '+' | '-' | '.') -> Output
.filter_ontology_terms(ontology_terms) -> Output
.filter_output_type(output_types) -> Output
.resize(width) -> Output
.map_track_data(fn) -> Output
```

**`class OutputMetadata`** (frozen dataclass) — same field names, each a
`TrackMetadata` (`pandas.DataFrame`) or `None`.

```python
.get(output: OutputType) -> TrackMetadata | None     # note: 'output', not 'output_type'
.concatenate() -> TrackMetadata                 # all output types in one frame
OutputMetadata.from_outputs(mapping)            # classmethod
```

**`class VariantOutput`** (frozen dataclass) — `.reference: Output`,
`.alternate: Output`.

### `models.variant_scorers`

**`enum AggregationType`** — `DIFF_MEAN`, `DIFF_SUM`, `DIFF_SUM_LOG2`,
`DIFF_LOG2_SUM`, `L2_DIFF`, `L2_DIFF_LOG1P`, `ACTIVE_MEAN`, `ACTIVE_SUM`. Names
read right-to-left in order of application.

**`enum BaseVariantScorer`** — `CENTER_MASK`, `CONTACT_MAP`, `GENE_MASK_LFC`,
`GENE_MASK_ACTIVE`, `GENE_MASK_SPLICING`, `PA_QTL`, `SPLICE_JUNCTION`.

**Scorer classes** (all frozen dataclasses, all with `.base_variant_scorer`,
`.name`, `.is_signed`, `.requested_output`, `.to_proto()`):

```python
CenterMaskScorer(requested_output: OutputType, width: int | None,
                 aggregation_type: AggregationType)
ContactMapScorer()
GeneMaskLFCScorer(requested_output: OutputType)
GeneMaskActiveScorer(requested_output: OutputType)
GeneMaskSplicingScorer(requested_output: OutputType, width: int | None)
PolyadenylationScorer()
SpliceJunctionScorer()
```

**Constants**

```python
RECOMMENDED_VARIANT_SCORERS   # dict: key -> scorer instance (see 07-variant-scoring.md)
SUPPORTED_ORGANISMS           # BaseVariantScorer -> organisms
SUPPORTED_OUTPUT_TYPES        # BaseVariantScorer -> [OutputType]
SUPPORTED_WIDTHS              # CENTER_MASK: [None,501,2001,10001,100001,200001]
                              # GENE_MASK_SPLICING: [None,101,1001,10001]
SUPPORTED_AGGREGATIONS        # CENTER_MASK: all AggregationTypes
```

**Functions**

```python
get_recommended_scorers(organism: dna_model_pb2.Organism) -> list[VariantScorerTypes]
    # NOTE: takes the proto enum, i.e. Organism.HOMO_SAPIENS.to_proto().
    # Passing dna_client.Organism directly returns an empty list, silently.
tidy_scores(scores, match_gene_strand=True, include_extended_metadata=True) -> pd.DataFrame | None
tidy_anndata(adata, match_gene_strand=True, include_extended_metadata=True) -> pd.DataFrame
```

`tidy_scores` accepts a flat or nested sequence of `AnnData` and concatenates
across scorers and variants. Columns: `variant_id`, `scored_interval`, `gene_id`,
`gene_name`, `gene_type`, `gene_strand`, `output_type`, `variant_scorer` /
`interval_scorer`, `track_name`, `track_strand`, `ontology_curie`, `gtex_tissue`,
`Assay title`, `biosample_name`, `biosample_type`, `transcription_factor`,
`histone_mark`, `raw_score`, `quantile_score`.

### `models.interval_scorers`

```python
enum IntervalAggregationType
enum BaseIntervalScorer                          # GENE_MASK
GeneMaskScorer(requested_output, width, aggregation_type)
RECOMMENDED_INTERVAL_SCORERS = {
    'RNA_SEQ': GeneMaskScorer(RNA_SEQ, width=200001, aggregation_type=MEAN)}
SUPPORTED_OUTPUT_TYPES, SUPPORTED_WIDTHS
```

---

## `alphagenome.data`

### `data.genome`

**Constants** — `STRAND_POSITIVE='+'`, `STRAND_NEGATIVE='-'`,
`STRAND_UNSTRANDED='.'`, `STRAND_OPTIONS`, `VALID_VARIANT_BASES=frozenset('ACGTN')`,
`PYRANGES_INTERVAL_COLUMNS`.

**`enum Strand`** (`IntEnum`) — `POSITIVE`, `NEGATIVE`, `UNSTRANDED`;
`.from_str(s)`, `.to_proto()`, `.from_proto(p)`.

**`class Interval`** (ordered dataclass)

```python
Interval(chromosome: str, start: int, end: int, strand: str = '.',
         name: str = '', info: dict = {})

# properties
.width  .negative_strand
# geometry
.center(use_strand=True) -> int
.resize(width, use_strand=True) -> Interval        .resize_inplace(...)
.shift(offset, use_strand=True)                    .boundary_shift(start_offset=0, end_offset=0)
.pad(start_pad, end_pad, *, use_strand=True)       .pad_inplace(...)
.truncate(reference_length=sys.maxsize)            .within_reference(reference_length)
.swap_strand()  .as_unstranded()  .copy()
# relations
.overlaps(interval) -> bool   .contains(interval) -> bool
.intersect(interval) -> Interval | None
.coverage(intervals, *, bin_size=1) -> np.ndarray
.coverage_stranded(intervals, *, bin_size=1)
.binary_mask(intervals, bin_size=1)  .binary_mask_stranded(...)
# conversion
Interval.from_str('chr1:100-200:+')      .to_proto() / Interval.from_proto(p)
.to_interval_dict() / Interval.from_interval_dict(d)
.to_pyranges_dict() / Interval.from_pyranges_dict(row, ignore_info=False)
```

**`class Variant`** (dataclass)

```python
Variant(chromosome: str, position: int, reference_bases: str,
        alternate_bases: str, name: str = '', info: dict = {})
# position is 1-BASED; .start / .end are 0-based

.start  .end  .reference_interval
.is_snv  .is_insertion  .is_deletion  .is_indel  .is_frameshift  .is_structural
.reference_overlaps(interval) -> bool     .alternate_overlaps(interval) -> bool
.as_truncated_str(max_length=50)          .split(anchor) -> tuple[Variant|None, Variant|None]
.copy()  .to_dict() / Variant.from_dict(d)  .to_proto() / Variant.from_proto(p)
Variant.from_str(string, variant_format=VariantFormat.DEFAULT)
```

**`enum VariantFormat`** — `DEFAULT` (`chr1:100:A>C`), `GTEX` (`chr1_100_A_C_b38`),
`OPEN_TARGETS` (`1_100_A_C`), `OPEN_TARGETS_BIGQUERY` (`1:100:A:C`),
`GNOMAD` (`1-100-A-C`). Method `.to_regex()`.

**`class Junction(Interval)`** — adds `.k`, `.donor`, `.acceptor` (strand-aware).

**Module functions** — `normalize_variant`, `intersect_intervals`,
`union_intervals`, `merge_overlapping_intervals`.

### `data.track_data`

**`enum AggregationType`** — `SUM`, `MAX` (for resolution changes).

**`class TrackData`** (frozen dataclass)

```python
TrackData(values: np.ndarray,            # (positions, num_tracks)
          metadata: pd.DataFrame,        # needs at least 'name' and 'strand'
          resolution: int = 1,
          interval: Interval | None = None,
          uns: dict | None = None)

# properties
.num_tracks  .width  .names  .strands  .ontology_terms  .positional_axes
# resolution
.change_resolution(resolution, aggregation_type=SUM)
.upsample(resolution, ...)   .downsample(resolution, ...)
# shape
.resize(width)   .pad(start_pad, end_pad)
.slice_by_positions(start, end)
.slice_by_interval(interval, match_resolution=False)
.bin_index(relative_position)
# track selection
.filter_tracks(mask)  .select_tracks_by_index(idx)  .select_tracks_by_name(names)
.filter_to_positive_strand()   .filter_to_negative_strand()   .filter_to_unstranded()
.filter_to_nonpositive_strand()  .filter_to_nonnegative_strand()  .filter_to_stranded()
.groupby(column) -> dict[str, TrackData]
# misc
.with_name_suffix(suffix)  .with_metadata_column(name, value)
.reverse_complement()  .copy()
# numpy-style indexing
tdata[interval_or_slice, track_mask_or_names_or_slice]
```

**Functions** — `concat(track_datas, extra_metadata_name_and_keys=None)`,
`interleave(track_datas, name_prefixes)`.

### `data.junction_data`

**`class JunctionData`** (frozen dataclass)

```python
JunctionData(junctions: np.ndarray,     # of genome.Junction, (num_junctions,)
             values: np.ndarray,        # (num_junctions, num_tracks)
             metadata: pd.DataFrame, interval=None, uns=None)

.num_tracks  .names  .strands  .possible_strands  .ontology_terms
.filter_tracks(mask)         .filter_to_strand(strand)
.filter_to_positive_strand() .filter_to_negative_strand()
.filter_by_tissue(tissue)    .filter_by_ontology(ontology_curie)  .filter_by_name(name)
.normalize_values(total_k=10.0)
.intersect_with_interval(interval)
```

**Function** —
`get_junctions_to_plot(*, predictions, name, strand, k_threshold=0.0) -> list[Junction]`

### `data.ontology`

**`enum OntologyType`** — `CLO`, `UBERON`, `CL`, `EFO`, `NTR`.

**`class OntologyTerm`** (frozen dataclass) — `.type`, `.id`, `.ontology_curie`,
`.to_proto()`.

**Functions** — `from_curie(curie)`, `from_curies(curies)`, `from_proto(proto)`.

### `data.gene_annotation`

**`enum TranscriptType`** — the full GENCODE set: `PROTEIN_CODING`, `LNCRNA`,
`MIRNA`, `SNRNA`, `SNORNA`, `RRNA`, `MISC_RNA`, `RETAINED_INTRON`,
`PROCESSED_TRANSCRIPT`, `NONSENSE_MEDIATED_DECAY`, `NON_STOP_DECAY`,
`PROCESSED_PSEUDOGENE`, `UNPROCESSED_PSEUDOGENE`, `TRANSCRIBED_*_PSEUDOGENE`,
`IG_*`, `TR_*`, `MT_RRNA`, `MT_TRNA`, `ARTIFACT`, `TEC`, and others.

```python
extract_tss(gtf, feature='transcript') -> pd.DataFrame
filter_transcript_type(gtf, transcript_types=None)
filter_protein_coding(gtf, include_gene_entries=False)
filter_to_longest_transcript(gtf)
filter_to_mane_select_transcript(gtf)
filter_transcript_support_level(gtf, transcript_support_levels)     # e.g. ['1']
upgrade_annotation_ids(old_ids, new_ids, patchless=False)
get_gene_interval(gtf, gene_symbol=None, gene_id=None) -> Interval
get_gene_intervals(gtf, gene_symbols=None, gene_ids=None) -> list[Interval]
```

### `data.transcript`

**`class Transcript`** (frozen dataclass)

```python
Transcript(exons, cds=None, start_codon=None, stop_codon=None,
           transcript_id=None, gene_id=None, protein_id=None,
           uniprot_id=None, info={})

.chromosome  .strand  .strand_int  .is_positive_strand  .is_negative_strand
.is_mitochondrial  .is_coding  .transcript_interval
.introns  .utr5  .utr3  .cds_including_stop_codon
.splice_regions  .splice_donors  .splice_acceptors
.splice_donor_sites  .splice_acceptor_sites
.selenocysteines  .selenocysteine_pos_in_protein
.offset_in_cds(genome_position) -> int | None
Transcript.from_gtf_df(transcript_df, ignore_info=True, fix_truncation=False)
Transcript.fix_truncation(transcript)
```

Gene name lives in `transcript.info['gene_name']`.

**`class TranscriptExtractor`**

```python
TranscriptExtractor(gtf_df: pd.DataFrame)
.cache_transcripts() -> None            # speeds up repeated extract() calls
.extract(interval) -> list[Transcript]
```

`MITOCHONDRIAL_CHROMS = ['M', 'chrM', 'MT']`

### `data.fold_intervals`

Training/validation/test splits used by the model folds.

```python
enum Subset            # TRAIN, VALID, TEST
get_all_folds() -> list[str]
get_fold_names(model_version, subset) -> list[str]
get_fold_intervals(model_version, organism, subset,
                   example_regions_path=None) -> pd.DataFrame
```

---

## `alphagenome.visualization`

### `visualization.plot_components`

```python
plot(components, interval, fig_width=20, fig_height_scale=1.0, title=None,
     despine=True, despine_keep_bottom=False, annotations=None,
     annotation_offset_range=(0.1, 0.6), hspace=0.3,
     xlabel=None) -> matplotlib.figure.Figure
```

**Components** — all subclass `AbstractComponent` (`.plot_ax(ax, axis_index,
interval)`, `.num_axes`, `.get_ax_height(axis_index)`, `.total_height`).

```python
Tracks(tdata, cmap='viridis', truncate_cmap=True, track_height=1.0, filled=False,
       ylabel_template='{name}:{strand}', ylabel_horizontal=True,
       shared_y_scale=False, global_ylims=None, max_num_tracks=50,
       track_colors=None, **kwargs)

OverlaidTracks(tdata: Mapping[str, TrackData], colors=None, cmap='viridis',
               track_height=1.0, ylabel_template='{name}:{strand}',
               ylabel_horizontal=True, shared_y_scale=False, global_ylims=None,
               yticks=None, yticklabels=None, alpha=0.8,
               order_tdata_by_mean=True, max_num_tracks=50,
               legend_loc='upper right', **kwargs)

ContactMaps(tdata, track_height=10.0, vmin=-1.0, vmax=2.0, norm=None,
            ylabel_horizontal=True, ylabel_template='{name}', cmap=None,
            max_num_tracks=10, **kwargs)

ContactMapsDiff(tdata, track_height=10.0, vmin=-1.0, vmax=1.0,
                ylabel_horizontal=True, ylabel_template='{name}',
                cmap='RdBu_r', max_num_tracks=10, **kwargs)

TranscriptAnnotation(transcripts, adaptive_fig_height=True, fig_height=1.0,
                     transcript_style=TranscriptStylePreset.MINIMAL.value,
                     plot_labels_once=True, label_name='gene_name', **kwargs)

SeqLogo(scores, scores_interval, fig_height=1.0, alphabet='ACGT',
        max_width=1000, ylabel='', ylabel_horizontal=True, ylim=None, **kwargs)

Sashimi(junction_track, fig_height=1.0, filter_threshold=None,
        ylabel_template='{name}', ylabel_horizontal=True, annotate_counts=True,
        normalize_values=True, interval_contained=True, rng=None, color=None)

EmptyComponent(fig_height=1.0)
```

**Annotations** — subclass `AbstractAnnotation` (`.plot_ax`, `.plot_labels`,
`.is_variant`, `.has_labels`, `.add_label`).

```python
IntervalAnnotation(intervals, colors='darkgray', alpha=0.2, labels=None,
                   use_default_labels=True, label_angle=15)
VariantAnnotation(variants, colors='orange', alpha=0.8, labels=None,
                  use_default_labels=True, label_angle=15, label_position='left')
```

### `visualization.plot`

```python
seqlogo(letter_heights, alphabet='ACGT', one_based=True, start=0, ax=None,
        letter_colors='default')                      # or 'factorbook'
plot_contact_map(contact_map, vmin=None, vmax=None, square=True,
                 cbar_shrink=0.4, ax=None, **kwargs)
plot_track(arr, ax, x=None, legend=False, ylim=None, color=None, filled=False)
plot_tracks(tracks, x=None, title=None, legend=False, fig_width=20,
            fig_track_height=1.5, ylim='auto', yticks_min_max_only=False,
            ylab=True, color=None, horizontal_ylab=True, filled_tracks=None,
            despine=True, despine_keep_bottom=False, plot_track_fn=plot_track)
sashimi_plot(junctions, ax, interval=None, filter_threshold=0.01,
             annotate_counts=True, rng=None, color=None)
pad_track(track, new_len, value=0) -> np.ndarray
```

### `visualization.plot_transcripts`

```python
class TranscriptStyle:      # dataclass
    cds_height, utr_height, cds_color, utr5_color, utr3_color,
    first_noncoding_exon_color, label_color, xlim_pad

enum TranscriptStylePreset  # STANDARD, MINIMAL

plot_transcripts(ax, transcripts, interval, zero_origin=False, label_name=None,
                 transcript_style=TranscriptStylePreset.STANDARD.value,
                 plot_labels_once=False, **kwargs)
draw_transcript(ax, transcript, interval, y, cds_height=0.7, utr_height=0.35,
                cds_color='#7f7f7f', utr5_color='#ff7f0e', utr3_color='#1f77b4',
                first_noncoding_exon_color='#2ca02c', shift=0, label=None,
                label_color='#7f7f7f', num_transcripts=1, **kwargs)
draw_strand_arrows(ax, transcript, interval, y, color, *, cds_height=0.22,
                   num_transcripts=1, max_arrows_per_intron=5)
draw_interval(ax, interval, y, label=None, height=0.5, shift=0,
              label_color='#7f7f7f', **kwargs)
```

---

## `alphagenome.interpretation`

### `interpretation.ism`

```python
ism_variants(interval, sequence, vocabulary='ACGT',
             skip_n=False) -> list[genome.Variant]
ism_matrix(variant_scores, variants, interval=None, multiply_by_sequence=True,
           vocabulary='ACGT', require_fully_filled=True) -> np.ndarray
```

`ism_matrix` returns shape `(interval.width, 4)` — the second axis is hard-coded
to 4, so a non-4-letter `vocabulary` is not really supported. With
`multiply_by_sequence=True` only the reference base at each position is non-zero.

---

## `alphagenome.atlas`

### `atlas.atlas`

```python
create(api_key, *, timeout=None, address=None) -> AtlasClient

class ScorerMetadata:       # frozen dataclass
    name: str
    is_signed: bool
    track_metadata: pd.DataFrame

class AtlasClient:
    query_variant(variant, *, requested_scorers, ontology_terms=None,
                  gene_ids=None, gene_names=None) -> Mapping[str, AnnData]
    query_variants(variants, *, requested_scorers, ontology_terms=None,
                   gene_ids=None, gene_names=None, progress_bar=True,
                   max_workers=10) -> Mapping[str, AnnData]
    query_interval(interval, *, requested_scorers, ontology_terms=None,
                   gene_ids=None, gene_names=None, progress_bar=True,
                   max_workers=10) -> Mapping[str, AnnData]
    scorer_metadata() -> Mapping[str, ScorerMetadata]
```

Results are keyed by scorer name. `DEFAULT_MAX_WORKERS = 10`.

---

## `alphagenome.io`

### `io.fasta`

```python
class FastaExtractor:
    FastaExtractor(fasta_path)
    .extract(interval) -> str
    .sequence_names -> list[str]
    .get_length_for_sequence_name(name) -> int

reverse_complement(sequence: str) -> str
```

---

## `alphagenome.colab_utils`

```python
get_api_key(secret='ALPHA_GENOME_API_KEY') -> str
```

Checks the environment variable first, then Colab secrets. Raises `ValueError`
with setup instructions if neither is available.
