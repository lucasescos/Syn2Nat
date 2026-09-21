# Visualization

`alphagenome.visualization` turns model outputs into matplotlib figures. Three
concepts: **plot**, **components**, **annotations**.

## The model

```python
from alphagenome.visualization import plot_components

fig = plot_components.plot(
    components=[...],           # vertically stacked panels
    interval=interval,          # shared x-axis
    annotations=[...],          # overlaid across all panels
    title='...',
)
```

- **`plot`** takes a list of components and returns a `matplotlib.figure.Figure`.
- A **component** is a light wrapper around one model output plus plot
  aesthetics. Each maps to one vertically stacked subplot. Each has its own
  y-axis but they **share an x-axis** — the DNA interval, in base pairs.
- An **annotation** is a figure element tied to the interval but outside any
  component — variant positions, promoter regions — drawn across all panels.

### `plot` arguments

`components`, `interval`, `fig_width=20`, `fig_height_scale=1.0`, `title=None`,
`despine=True`, `despine_keep_bottom=False`, `annotations=None`,
`annotation_offset_range=(0.1, 0.6)`, `hspace=0.3`, `xlabel=None`.

> **Caution:** Always pass `VariantAnnotation` and `IntervalAnnotation` to the
> `annotations=[...]` argument of `plot()`, **never** in `components=[...]`.
> Placing annotations in `components` raises:
> `AttributeError: 'VariantAnnotation' object has no attribute 'num_axes'`.

**Zooming** is done by passing a smaller `interval` to `plot` — the data is
untouched:


```python
interval=output.rna_seq.interval.resize(2**15)
interval=interval.shift(offset).resize(width)
```

## Components

| Component | Draws | Data shape | Use for | Good for variants? |
| :--- | :--- | :--- | :--- | :--- |
| `Tracks` | Line plot, one value per position | 1D | All except `SPLICE_JUNCTIONS`, `CONTACT_MAPS` | No |
| `OverlaidTracks` | Two-plus lines on one axis in different colours | 1D × N | Same, e.g. REF vs ALT | Yes |
| `Sashimi` | Arcs between position pairs, thickness ∝ value | 2D sparse | `SPLICE_JUNCTIONS` | Yes |
| `SeqLogo` | Letters with heights ∝ per-position score | 1D + sequence | ISM contribution scores | Yes |
| `ContactMaps` | Heatmap of a position × position matrix | 2D | `CONTACT_MAPS` | No |
| `ContactMapsDiff` | Same, diverging colormap centred on zero (white) | 2D | `CONTACT_MAPS` differences (ALT − REF) | Yes |
| `TranscriptAnnotation` | Horizontal transcript models — exons, introns, UTRs, direction | Interval(s) | Gene context | No |
| `VariantAnnotation` | Semi-transparent rectangle/vertical line spanning all panels | Variant(s) | Marking variants | Yes |
| `IntervalAnnotation` | Semi-transparent rectangles spanning all panels | Interval(s) | Promoters, pA sites, TSSs | — |
| `EmptyComponent` | A blank panel | — | Spacing | — |
| `AbstractComponent` | Base class to subclass | — | Custom components | — |

### Key component arguments

```python
plot_components.Tracks(
    tdata,                          # TrackData
    ylabel_template='{name}:{strand}',
    track_height=1.0,
    filled=False,                   # fill under the line — good for ChIP
    shared_y_scale=False,
    global_ylims=None,
    track_colors=None,              # sequence of colors, or one color
    cmap='viridis',
    max_num_tracks=50,
)

plot_components.OverlaidTracks(
    tdata={'REF': ..., 'ALT': ...},   # mapping label -> TrackData
    colors={'REF': 'dimgrey', 'ALT': 'red'},
    alpha=0.8,
    legend_loc='upper right',        # None to hide
    order_tdata_by_mean=True,
)

plot_components.TranscriptAnnotation(
    transcripts,
    fig_height=1.0,
    adaptive_fig_height=True,
    label_name='gene_name',
    transcript_style=plot_transcripts.TranscriptStylePreset.MINIMAL.value,
)

plot_components.Sashimi(
    junction_track,                  # JunctionData
    ylabel_template='{name}',
    color=None,
    filter_threshold=None,
    annotate_counts=True,
    normalize_values=True,
)

plot_components.ContactMaps(tdata, vmin=-1.0, vmax=2.0, cmap=None, track_height=10.0)
plot_components.ContactMapsDiff(tdata, vmin=-1.0, vmax=1.0, cmap='RdBu_r')

plot_components.SeqLogo(scores, scores_interval, ylabel='', alphabet='ACGT',
                        max_width=1000, ylim=None)
```

`ylabel_template` interpolates any metadata column: `{name}`, `{strand}`,
`{biosample_name}`, `{transcription_factor}`, `{histone_mark}`.

### Annotations

```python
plot_components.VariantAnnotation(
    variants, colors='orange', alpha=0.8,
    labels=None, use_default_labels=True, label_angle=15, label_position='left',
)
plot_components.IntervalAnnotation(
    intervals, colors='darkgray', alpha=0.2,
    labels=None, use_default_labels=True, label_angle=15,
)
```

### Custom components

Extend `AbstractComponent` (workhorse method: `plot_ax(ax, axis_index, interval)`,
plus `num_axes` and `get_ax_height`) or `AbstractAnnotation`. Any data of your
own can go through the library as long as it is in the expected format — e.g. a
`TrackData` for `Tracks`. Figures are plain matplotlib, so they are extendable,
and you can always work with the raw arrays and plot them yourself.

## Recipes by modality

Setup shared by all of them:

```python
gtf = pd.read_feather(HG38_GTF_FEATHER)
gtf_transcript = gene_annotation.filter_transcript_support_level(
    gene_annotation.filter_protein_coding(gtf), ['1'])
transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
longest_transcript_extractor = transcript.TranscriptExtractor(
    gene_annotation.filter_to_longest_transcript(gtf_transcript))

interval = genome.Interval('chr22', 36_150_498, 36_252_898).resize(
    dna_client.SEQUENCE_LENGTH_1MB)
longest_transcripts = longest_transcript_extractor.extract(interval)
ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}
```

### Gene expression (RNA_SEQ, CAGE)

```python
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.RNA_SEQ, dna_client.OutputType.CAGE},
    ontology_terms=['UBERON:0001159', 'UBERON:0001155'],  # sigmoid, transverse colon
)

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.Tracks(output.rna_seq,
            ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}'),
        plot_components.Tracks(output.cage,
            ylabel_template='CAGE: {biosample_name} ({strand})\n{name}'),
    ],
    interval=interval,
    title='Predicted RNA expression for colon tissue',
)
```

Positive- and negative-strand tracks have clearly different signals; check which
genes are on which strand with
`[t.info['gene_name'] for t in longest_transcripts if t.strand == '+']`.

**PROCAP (Precision Run-On capped sequencing):** Like CAGE, `PROCAP` captures
transcription start sites at 1 bp resolution on positive and negative strands.
Plot it identically using `plot_components.Tracks(output.procap)` or
`OverlaidTracks`. Note that PROCAP is available for human only (6 biosamples,
12 tracks) and has no mouse tracks.


### A variant's effect on expression

For a gene on the negative strand, drop the positive-strand tracks and zoom to
the gene:

```python
variant = genome.Variant.from_str('chr22:36201698:A>C')
output = dna_model.predict_variant(
    interval=interval, variant=variant,
    requested_outputs={dna_client.OutputType.RNA_SEQ, dna_client.OutputType.CAGE},
    ontology_terms=ontology_terms,
)

apol4_interval = gene_annotation.get_gene_interval(gtf, gene_symbol='APOL4')
apol4_interval.resize_inplace(apol4_interval.width + 1000)   # 1 kb of flank

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.OverlaidTracks(
            tdata={'REF': output.reference.rna_seq.filter_to_nonpositive_strand(),
                   'ALT': output.alternate.rna_seq.filter_to_nonpositive_strand()},
            colors=ref_alt_colors,
            ylabel_template='{biosample_name} ({strand})\n{name}'),
        plot_components.OverlaidTracks(
            tdata={'REF': output.reference.cage.filter_to_nonpositive_strand(),
                   'ALT': output.alternate.cage.filter_to_nonpositive_strand()},
            colors=ref_alt_colors,
            ylabel_template='{biosample_name} ({strand})\n{name}'),
    ],
    annotations=[plot_components.VariantAnnotation([variant])],
    interval=apol4_interval,
    title='Effect of variant on predicted RNA expression in colon tissue',
)
```

ALT below REF across the gene means reduced expression; an exon where ALT is zero
but REF has a peak suggests an exon skipping event.

### Chromatin accessibility (DNASE, ATAC), with custom annotation

```python
output = dna_model.predict_interval(
    interval,
    requested_outputs={dna_client.OutputType.DNASE, dna_client.OutputType.ATAC},
    ontology_terms=['UBERON:0000317', 'UBERON:0001155', 'UBERON:0001157',
                    'UBERON:0001159', 'UBERON:0004992', 'UBERON:0008971'],
)

promoter_intervals = [
    genome.Interval('chr22', 36_201_799, 36_202_681,
                    name='Ensembl_promoter:ENSR00001367790'),
    genome.Interval('chr22', 36_204_705, 36_205_330,
                    name='Ensembl_promoter:ENSR00001367792'),
]

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.Tracks(output.dnase,
            ylabel_template='DNASE: {biosample_name} ({strand})\n{name}'),
        plot_components.Tracks(output.atac,
            ylabel_template='ATAC: {biosample_name} ({strand})\n{name}'),
    ],
    interval=variant.reference_interval.resize(8000),
    annotations=[plot_components.VariantAnnotation([variant]),
                 plot_components.IntervalAnnotation(promoter_intervals)],
)
```

DNASE and ATAC signals are similar but not identical — the assays differ.

Any 1 bp feature can be annotated the same way; e.g. polyadenylation sites from
[PolyADB v3](https://exon.apps.wistar.org/polya_db/v3/misc/download.php):

```python
apol4_pAs = [genome.Interval('chr22', 36_189_128, 36_189_129, '-'),
             genome.Interval('chr22', 36_190_089, 36_190_090, '-'),
             genome.Interval('chr22', 36_190_144, 36_190_145, '-')]

annotations=[plot_components.IntervalAnnotation(
    apol4_pAs, alpha=1, labels=['pA_3', 'pA_2', 'pA_1'], label_angle=90)]
```

Close annotations may visually overlap depending on zoom level.

### Splicing

`SPLICE_SITES` is tissue-agnostic, so `ontology_terms` is not applied to it.

```python
output = dna_model.predict_variant(
    interval=interval, variant=variant,
    requested_outputs={dna_client.OutputType.RNA_SEQ,
                       dna_client.OutputType.SPLICE_SITES,
                       dna_client.OutputType.SPLICE_SITE_USAGE,
                       dna_client.OutputType.SPLICE_JUNCTIONS},
    ontology_terms=['UBERON:0001157', 'UBERON:0001159'],
)
ref_output, alt_output = output.reference, output.alternate

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(transcript_extractor.extract(interval)),
        plot_components.Sashimi(
            ref_output.splice_junctions.filter_to_strand('-')
                      .filter_by_tissue('Colon_Transverse'),
            ylabel_template='Reference {biosample_name} ({strand})\n{name}'),
        plot_components.Sashimi(
            alt_output.splice_junctions.filter_to_strand('-')
                      .filter_by_tissue('Colon_Transverse'),
            ylabel_template='Alternate {biosample_name} ({strand})\n{name}'),
        plot_components.OverlaidTracks(
            tdata={'REF': ref_output.rna_seq.filter_to_nonpositive_strand(),
                   'ALT': alt_output.rna_seq.filter_to_nonpositive_strand()},
            colors=ref_alt_colors,
            ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}'),
        plot_components.OverlaidTracks(
            tdata={'REF': ref_output.splice_sites.filter_to_nonpositive_strand(),
                   'ALT': alt_output.splice_sites.filter_to_nonpositive_strand()},
            colors=ref_alt_colors,
            ylabel_template='SPLICE SITES: {name} ({strand})'),
        plot_components.OverlaidTracks(
            tdata={'REF': ref_output.splice_site_usage.filter_to_nonpositive_strand(),
                   'ALT': alt_output.splice_site_usage.filter_to_nonpositive_strand()},
            colors=ref_alt_colors,
            ylabel_template='SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}'),
    ],
    interval=apol4_interval,
    annotations=[plot_components.VariantAnnotation([variant])],
)
```

### Histone marks (CHIP_HISTONE)

Predictions come back at 128 bp resolution, so the plotted interval width must be
a multiple of 128 — the library handles this automatically. Five major marks
(H3K4me3, H3K4me1, H3K27ac, H3K27me3, H3K36me3) are available for ~40% of
biosamples; check coverage first:

```python
output_metadata.chip_histone[
    output_metadata.chip_histone['biosample_name'].str.contains('colon')]
```

Group tracks by mark and colour them:

```python
reordered = output.chip_histone.select_tracks_by_index(
    output.chip_histone.metadata.sort_values('histone_mark').index)

histone_to_color = {'H3K27AC': '#e41a1c', 'H3K36ME3': '#ff7f00',
                    'H3K4ME1': '#377eb8', 'H3K4ME3': '#984ea3',
                    'H3K9AC': '#4daf4a', 'H3K27ME3': '#ffc0cb'}
track_colors = (reordered.metadata['histone_mark']
                .map(lambda x: histone_to_color.get(x.upper(), '#000000')).values)

# TSS positions from GENCODE, as annotation.
gtf_tss = gene_annotation.extract_tss(gtf_longest_transcript)
tss_as_intervals = [
    genome.Interval(chromosome=row.Chromosome, start=row.Start,
                    end=row.End + 1000,      # widen so TSSs are visible
                    name=row.gene_name)
    for _, row in gtf_tss.iterrows()
]

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.Tracks(reordered, filled=True, track_colors=track_colors,
            ylabel_template='CHIP HISTONE: {biosample_name} ({strand})\n{histone_mark}'),
    ],
    interval=interval,
    annotations=[plot_components.IntervalAnnotation(tss_as_intervals,
                                                    alpha=0.5, colors='blue')],
    despine_keep_bottom=True,
)
```

Promoter-associated marks (H3K4me3, H3K27ac, H3K9ac) peak at TSSs; H3K36me3 and
H3K4me1 track gene bodies and enhancers instead.

### TF binding (CHIP_TF)

Most biosamples only cover CTCF and POLR2A; HepG2 and K562 cover many more (501
and 269 respectively), so include those cell lines for good TF coverage.

```python
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.CHIP_TF},
    ontology_terms=['UBERON:0001159', 'UBERON:0001157',
                    'EFO:0002067',   # K562
                    'EFO:0001187'],  # HepG2
)

# Keep only the cell lines.
output_chip_tf = output.chip_tf.filter_tracks(
    output.chip_tf.metadata['ontology_curie'].isin(['EFO:0002067', 'EFO:0001187']).values)

# Keep only tracks with a strong peak.
output_filtered = output_chip_tf.filter_tracks(output_chip_tf.values.max(axis=0) > 8000)

# Or: the 10 strongest tracks within a gene.
max_predictions = output.chip_tf.slice_by_interval(
    apol4_interval, match_resolution=True).values.max(axis=0)
output_filtered = output.chip_tf.filter_tracks(
    max_predictions >= np.sort(max_predictions)[-10])
```

Averaging a TF's signal across tissues means building a new `TrackData`:

```python
mean_ctcf = output_filtered.values[
    :, output_filtered.metadata['transcription_factor'] == 'CTCF'].mean(axis=1)

tdata_mean_ctcf = track_data.TrackData(
    values=mean_ctcf[:, None],
    metadata=pd.DataFrame({'transcription_factor': ['CTCF'],
                           'name': ['mean'], 'strand': ['.']}),
    interval=output_filtered.interval,
    resolution=output_filtered.resolution,
)
```

### Contact maps

2048 bp resolution, so only intervals much wider than that are meaningful.

```python
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.CONTACT_MAPS},
    ontology_terms=['EFO:0002824'],   # HCT116 colon carcinoma
)

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.ContactMaps(output.contact_maps,
            ylabel_template='{biosample_name}\n{name}', cmap='autumn_r', vmax=1.0),
    ],
    interval=interval,
)
```

Blocks along the diagonal are candidate topologically-associated domains (TADs).
For a variant's effect, use `ContactMapsDiff` on the difference:

```python
plot_components.ContactMapsDiff(tdata=alt.contact_maps - ref.contact_maps)
```

## Lower-level plotting functions

`visualization.plot` holds the primitives the components are built on, for use
with your own axes:

```python
plot.seqlogo(letter_heights, alphabet='ACGT', ax=None,
             letter_colors='default')          # or 'factorbook'
plot.plot_contact_map(contact_map, vmin=None, vmax=None, square=True, ax=None)
plot.plot_track(arr, ax, color=None, filled=False, ylim=None)
plot.plot_tracks(tracks, fig_width=20, fig_track_height=1.5, ylim='auto')
plot.sashimi_plot(junctions, ax, interval=None, filter_threshold=0.01)
plot.pad_track(track, new_len, value=0)
```

`visualization.plot_transcripts` draws transcript models:

```python
plot_transcripts.plot_transcripts(ax, transcripts, interval, label_name='gene_name',
                                  transcript_style=TranscriptStylePreset.STANDARD.value)
plot_transcripts.draw_transcript(ax, transcript, interval, y, ...)
plot_transcripts.draw_interval(ax, interval, y, label=None, height=0.5)
```

`TranscriptStyle` fields: `cds_height`, `utr_height`, `cds_color`, `utr5_color`,
`utr3_color`, `first_noncoding_exon_color`, `label_color`, `xlim_pad`. Presets:
`STANDARD`, `MINIMAL`.

## Practical notes

- Transcript annotations come from GENCODE: hg38 release 46 (human), mm10
  release M23 (mouse).
- You are **not** limited to protein-coding genes or the longest transcript —
  just drop the `filter_protein_coding` / `filter_to_longest_transcript` calls.
  More transcripts makes for a busier plot; raise `fig_height` on
  `TranscriptAnnotation` to keep it legible.
- Figures are matplotlib, so `plt.show()`, `fig.savefig(...)` and any further
  customisation work as usual.
