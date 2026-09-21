# FAQ

Condensed from the official FAQ. Cross-references point at the relevant file here.

## Model inputs

**How do I make predictions for a specific genomic region?**
Define any region in the human or mouse genome and call `predict_interval`. See
[03-quickstart.md](03-quickstart.md).

**How do I specify a genomic region?**
With `genome.Interval(chromosome, start, end)`. It uses **0-based indexing**,
consistent with Python: the interval includes the base at `start` up to the base
at `end - 1`. So `genome.Interval('chr1', 0, 1)` is the first base of chr1, width
1. For overlaps, remember the end is exclusive:
`genome.Interval('chr1', 0, 1).overlaps(genome.Interval('chr1', 1, 2))` is
`False`.

**Which reference genome versions?**
Human hg38 (`GRCh38.p13.genome.fa`) and mouse mm10 (`GRCm38.p6.genome.fa`). For
other builds such as hg19, use
[LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert to hg38 first.

**Can I predict for any arbitrary DNA sequence?**
Yes, as long as it's within a supported length. But predictions have only been
evaluated on sequences that differ from the reference by a relatively small
amount (SNPs and indels). Large departures — structural variants, heavy padding,
synthetic sequences, artificial constructs — may be unreliable.

**Can I predict for other species?**
Yes, with the caveat that the model was trained only on mouse and human DNA.
Quality likely degrades with evolutionary distance, though this has not been
formally benchmarked.

**What is the longest input sequence?**
1 Mb, precisely 2^20 bp. Also supported: ~16 kb, ~100 kb, ~500 kb. Use 1 Mb where
possible for the best results.

**What if my sequence length isn't supported?**
Use `genome.Interval.resize` to crop or expand to the nearest supported length.
`.resize` expands using actual surrounding genomic data, not padding.

## Model outputs

**How many tracks per output type, and what do they represent?**
From 4 (`SPLICE_SITES`) to 1617 (`CHIP_TF`). Each track corresponds to a cell
type or tissue, plus properties such as strand or a specific transcription factor
(for `CHIP_TF`).
Full table in [06-output-types.md](06-output-types.md).

**How do I find out what tissue or cell type a track refers to?**
Look at the output metadata, where biosample names and ontology CURIEs are
listed. See [11-recipes.md](11-recipes.md#find-the-ontology-term-for-a-tissue).

**What is an ontology CURIE?**
A Compact Uniform Resource Identifier — a standardised abbreviated code (e.g.
`UBERON:0001114` for liver) that uniquely identifies an ontology term.

**Where do the CURIEs come from?**
From the IDs in the source training data, restricted to UBERON, CL, CLO, EFO and
NTR following ENCODE practice. Use EBI's
[Ontology Lookup Service](https://www.ebi.ac.uk/ols4) to understand
relationships between terms.

**What is strandedness?**
DNA is double-stranded; by convention one molecule is the forward/positive strand
(5'→3') and the other the reverse/negative strand (3'→5'). *Unstranded* assays
don't distinguish which strand a measurement came from (ATAC-seq, for example).
*Stranded* assays annotate each measurement with its strand — important for
transcriptional assays, e.g. distinguishing two transcripts sharing a TSS on
opposite strands. Not all RNA-seq is stranded; GTEx RNA-seq is unstranded.

**How is strandedness handled in the output metadata?**
`+` positive, `-` negative, `.` unstranded. Stranded assays get two tracks per
biosample; unstranded assays get one. Convenience methods on `TrackData`:
`filter_to_negative_strand()`, `filter_to_nonpositive_strand()`, and friends.

**How do I save model outputs?**
*Variant scores*: convert with `tidy_scores` to a pandas DataFrame and export to
CSV. *Track predictions*: `TrackData.values` is a NumPy array — use `numpy.save`
(`.npy`) or `numpy.savez_compressed` (`.npz`).

**What are the model's limitations?**
See [01-overview.md](01-overview.md#limitations) — tissue-specificity and
long-range interactions remain challenging, the species scope is human and mouse,
personal genomes are unbenchmarked, the scope is molecular rather than
phenotypic, and the model is not diploid-aware (trained on unphased data, so
variant effects don't inherently model heterozygous states).

## Visualization

**How do I visualize predictions?**
Any tool works on the numerical output, but `alphagenome.visualization` produces
matplotlib figures directly from API outputs. See
[09-visualization.md](09-visualization.md).

**Can I design my own visualizations?**
Yes. The figures are matplotlib, so they're extendable, and you can always work
with the raw output arrays. To add new component types, extend
`plot_components.AbstractComponent` / `AbstractAnnotation`.

**Where do the transcript annotations come from?**
GENCODE GTF files: hg38 release 46 (human) and mm10 release M23 (mouse).

**Am I limited to protein-coding genes and the longest transcript?**
No. Remove the `gene_annotation.filter_protein_coding(gtf)` and
`gene_annotation.filter_to_longest_transcript(gtf)` calls to include everything.
More transcripts makes the plot busier — increase `fig_height` on
`TranscriptAnnotation` to keep it legible.

## Variant scoring

**How do I score splicing variants?**
Combine the three splicing scorers into the merged splicing score described in
the paper. See [08-splicing.md](08-splicing.md).

**How do I define a variant?**
Create a `genome.Variant`. Note the indexing asymmetry: `Variant.start` and
`Variant.end` are 0-based, but most public databases (dbSNP and others) are
1-based, so `Variant` is *initialised* with a **1-based `position`** and converts
internally — `start` returns `position - 1`.

**Are there tools to help define variants and run inference?**
Yes — the single-variant scoring flow, the batch VCF flow, and the splicing flow
are all in [11-recipes.md](11-recipes.md) and [08-splicing.md](08-splicing.md).

**Do `reference_bases` have to match the reference genome?**
No. You can pass any sequence. `predict_variant` is agnostic to the true
reference allele and uses the REF/ALT you specify.

**Are indels supported?**
Yes. Indels are specified with left-alignment. For scoring, AlphaGenome adopts
SpliceAI's indel alignment strategy: inserted bases are summarised by the maximum
over the inserted segment, deleted bases are treated as zero signal in the ALT
context, enabling consistent positional comparisons.

**Which scorer should I use for a given modality?**
In practice most scoring strategies work for any modality. The recommended
pairings, based on DeepMind's evaluations, are in
[07-variant-scoring.md](07-variant-scoring.md#recommended-scorers).

**Can I write my own variant scoring strategy?**
Not through the API. But since scoring is just aggregation of REF and ALT track
predictions, you can write your own aggregation over `predict_variant` output.

**What's the difference between `quantile_score` and `raw_score`?**
`raw_score` is the scorer's direct output, but different tracks and modalities
produce different scales — splice site usage returns values in [0, 1] while gene
expression is unbounded and signed. To compare across them, AlphaGenome estimates
a background distribution per scorer and track from common variants (MAF > 0.01
in any gnomAD v3 population) and converts a raw score to its rank within that
background. A quantile score of 0.99 sits at the 99th percentile of common
variants. For signed scorers, the [0, 1] quantile is linearly rescaled to
[−1, 1] to preserve direction (0th percentile → −1, 50th → 0, 100th → +1). The
magnitude never exceeds 0.999990, because ~300K variants were used. Recommended
use: **quantile score** to judge whether a raw score is unusually large,
**raw score** as the magnitude of effect for that scorer and track. Quantile
scores are only available for the recommended scorers.

**How do I convert splicing scores to PSI values?**
The splicing scorers predict splice site usage and junction counts, not PSI
directly, but PSI3/PSI5 can be derived from the predicted junction counts. See
[08-splicing.md](08-splicing.md#deriving-psi-values).

**Can I predict the combined effect of multiple variants in a haplotype?**
Yes, with a workaround. The API takes one variant at a time, but you can build an
alternative sequence containing all variants of interest and compare its
prediction to those of the individual variants. See
[11-recipes.md](11-recipes.md#haplotype-analysis-combined-effect-of-several-variants).

## Common errors & troubleshooting

**Why does `score_variant` raise `TypeError: unexpected keyword argument 'ontology_terms'`?**
`score_variant` and `score_variants` evaluate all tracks associated with the
scorer and do **not** accept `ontology_terms` (unlike `predict_variant`). Filter
the resulting DataFrame afterwards: `df[df['ontology_curie'] == 'UBERON:...']`.

**Why does a low-expression gene have a quantile score > 0.999 but no visual effect on tracks?**
This is a statistical artifact caused by variance stabilization. In genes with
minimal baseline expression in a tissue, common-variant variance is near zero,
causing tiny numerical fluctuations to appear in the extreme quantile tail.
Always check that `|raw_score| >= 0.1` and verify baseline expression with
`RNA_SEQ_ACTIVE`.

**Why did `get_recommended_scorers` return an empty list?**
`get_recommended_scorers` expects the protobuf enum: pass
`dna_client.Organism.HOMO_SAPIENS.to_proto()`. Passing `dna_client.Organism`
directly silently returns an empty list `[]`.

**Why did `gtf[gtf['feature'] == 'exon']` raise `KeyError: 'feature'`?**
The pre-compiled GENCODE `.feather` file uses PascalCase for coordinate columns:
`Chromosome`, `Start`, `End`, `Strand`, and `Feature`. Metadata columns remain
lowercase (`gene_name`, `gene_id`, `gene_type`). Use `gtf['Feature']`.

**Why did `plot_components.plot` raise `AttributeError: 'VariantAnnotation' object has no attribute 'num_axes'`?**
`VariantAnnotation` and `IntervalAnnotation` must be passed to the `annotations`
argument of `plot()`, not `components`. E.g.
`plot(components=[...], annotations=[VariantAnnotation([variant])])`.


## Other

**What terms of use apply to outputs?**
The API is for **non-commercial use only**, subject to the
[Terms of Service](https://deepmind.google.com/science/alphagenome/terms).
Outputs must not be used to train other machine-learning models.

**How should I cite AlphaGenome?**
See the BibTeX entry in [01-overview.md](01-overview.md#citing).

**Who do I contact with issues, enquiries or feedback?**
Bugs and code issues on
[GitHub](https://github.com/google-deepmind/alphagenome). For general feedback,
usage questions and feature requests, the
[community forum](https://www.alphagenomecommunity.com) is actively monitored by
the team and usually the fastest route. Otherwise <alphagenome@google.com> —
responses may be delayed given volume.
