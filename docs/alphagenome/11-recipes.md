# Recipes

Task-shaped extracts from the tutorial Colabs.

---

## Find the ontology term for a tissue

Every prediction call takes `ontology_terms`; this is how you get the CURIE.

```python
output_metadata = dna_model.output_metadata(
    dna_client.Organism.HOMO_SAPIENS
).concatenate()

# Search for a tissue by name.
output_metadata[output_metadata['biosample_name'].str.contains('brain', case=False)]
```

In Colab, `data_table.enable_dataframe_formatter()` makes the frame interactive —
click **Filter** at the top right and type a tissue name into "Search by all
fields" to find the matching `ontology_curie` and output type.

Track counts per output type, human vs mouse:

```python
human_tracks = (dna_model.output_metadata(dna_client.Organism.HOMO_SAPIENS)
                .concatenate().groupby('output_type').size().rename('# Human tracks'))
mouse_tracks = (dna_model.output_metadata(dna_client.Organism.MUS_MUSCULUS)
                .concatenate().groupby('output_type').size().rename('# Mouse tracks'))
pd.concat([human_tracks, mouse_tracks], axis=1).astype(pd.Int64Dtype())
```

`PROCAP` has no mouse tracks. Use EBI's
[Ontology Lookup Service](https://www.ebi.ac.uk/ols4) to explore relationships
between terms (e.g. to find a parent tissue when the exact one is absent).

---

## Score a single variant across every modality

```python
variant = genome.Variant('chr22', 36201698, 'A', 'C')
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

variant_scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
)
df_scores = variant_scorers.tidy_scores(variant_scores)
df_scores.to_csv(f'{variant}_scores.csv', index=False)
```

`RECOMMENDED_VARIANT_SCORERS` has 19 entries, which fits under the
20-scorer-per-request limit (`MAX_VARIANT_SCORERS_PER_REQUEST`), so passing the
whole dict works. Add your own scorers on top and you may need to subset.

Sequence lengths for the interval come from
`dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{name}']`, where `name`
is one of `16KB`, `100KB`, `500KB`, `1MB`.

---

## Batch score a VCF

Prepare a tab-separated file with these columns:

| Column | Meaning |
| :--- | :--- |
| `variant_id` | Unique identifier |
| `CHROM` | Chromosome as a string starting with `chr` (e.g. `chr1`) |
| `POS` | 1-based base-pair position, hg38 (human) or mm10 (mouse) |
| `REF` | Reference nucleotide sequence at that position |
| `ALT` | Alternate nucleotide sequence |

```python
from io import StringIO
from tqdm import tqdm

vcf_file = """variant_id\tCHROM\tPOS\tREF\tALT
chr3_58394738_A_T_b38\tchr3\t58394738\tA\tT
chr8_28520_G_C_b38\tchr8\t28520\tG\tC
chr16_636337_G_A_b38\tchr16\t636337\tG\tA
chr16_1135446_G_T_b38\tchr16\t1135446\tG\tT
"""
vcf = pd.read_csv(StringIO(vcf_file), sep='\t')

for column in ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']:
    if column not in vcf.columns:
        raise ValueError(f'VCF file is missing required column: {column}.')

organism = dna_client.Organism.HOMO_SAPIENS
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
selected_scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
```

Drop scorers that don't support your organism — otherwise the request fails:

```python
unsupported = [
    s for s in selected_scorers
    if (organism.value not in variant_scorers.SUPPORTED_ORGANISMS[s.base_variant_scorer])
    | ((s.requested_output == dna_client.OutputType.PROCAP)
       & (organism == dna_client.Organism.MUS_MUSCULUS))
]
for s in unsupported:
    selected_scorers.remove(s)
```

(Polyadenylation is human-only; `PROCAP` has no mouse tracks.)

```python
results = []
for _, vcf_row in tqdm(vcf.iterrows(), total=len(vcf)):
    variant = genome.Variant(
        chromosome=str(vcf_row.CHROM),
        position=int(vcf_row.POS),
        reference_bases=vcf_row.REF,
        alternate_bases=vcf_row.ALT,
        name=vcf_row.variant_id,
    )
    interval = variant.reference_interval.resize(sequence_length)
    results.append(dna_model.score_variant(
        interval=interval, variant=variant,
        variant_scorers=selected_scorers, organism=organism,
    ))

df_scores = variant_scorers.tidy_scores(results)
df_scores.to_csv('variant_scores.csv', index=False)
```

The result can be very large with many scorers. Filter with plain pandas:

```python
# Just the effects on T cells.
columns = [c for c in df_scores.columns if c != 'ontology_curie']
df_scores[df_scores['ontology_curie'] == 'CL:0000084'][columns]
```

For higher throughput, `dna_model.score_variants(intervals=..., variants=...,
max_workers=2)` runs the batch concurrently instead of looping.

---

## Haplotype analysis: combined effect of several variants

The API takes one variant at a time, but you can build a single synthetic
"variant" spanning all of them: fetch the reference stretch, apply every
mutation, and treat the whole span as one REF→ALT pair. Comparing that to the
individual variants reveals cooperative effects.

Constraints: all mutations on the same chromosome, span at most 1 Mb.

```python
from alphagenome.io import fasta

def create_haplotype(mutations: list[genome.Variant]) -> dict:
    """Returns {'combined_variant': Variant, 'individual_variants': {pos: Variant}}."""
    chroms = {m.chromosome for m in mutations}
    if len(chroms) > 1:
        raise ValueError(f'Multiple chromosomes detected: {chroms}.')
    target_chr = list(chroms)[0]

    positions = [m.position for m in mutations]
    min_pos, max_pos = min(positions), max(positions)
    if (max_pos - min_pos) > 1_000_000:
        raise ValueError(f'Range too large: {max_pos - min_pos} bp. Max 1MB allowed.')

    fasta_path = ('https://storage.googleapis.com/alphagenome/reference/'
                  'gencode/hg38/GRCh38.p13.genome.fa')
    extractor = fasta.FastaExtractor(fasta_path)

    # FastaExtractor takes a 0-indexed, half-open Interval.
    raw = extractor.extract(
        genome.Interval(chromosome=target_chr, start=min_pos - 1, end=max_pos)).upper()

    ref_list, alt_list = list(raw), list(raw)
    individual_variants = {}
    for mut in mutations:
        rel_index = mut.position - min_pos
        ref_list[rel_index] = mut.reference_bases.upper()
        alt_list[rel_index] = mut.alternate_bases.upper()
        individual_variants[mut.position] = mut

    return {
        'combined_variant': genome.Variant(
            chromosome=target_chr, position=min_pos,
            reference_bases=''.join(ref_list),
            alternate_bases=''.join(alt_list)),
        'individual_variants': individual_variants,
    }
```

Predict the haplotype and each variant separately with identical parameters:

```python
mutation_data = [
    genome.Variant('chr5', 1295113, 'G', 'A'),
    genome.Variant('chr5', 1295135, 'G', 'A'),
]
haplotype = create_haplotype(mutation_data)

interval = haplotype['combined_variant'].reference_interval.resize(
    dna_client.SEQUENCE_LENGTH_1MB)

prediction_params = {
    'interval': interval,
    'requested_outputs': [dna_client.OutputType.RNA_SEQ, dna_client.OutputType.ATAC],
    'ontology_terms': ['UBERON:0002107'],   # liver
}

all_predictions = {'combined': dna_model.predict_variant(
    variant=haplotype['combined_variant'], **prediction_params)}
for pos, var_obj in haplotype['individual_variants'].items():
    all_predictions[pos] = dna_model.predict_variant(variant=var_obj, **prediction_params)
```

Overlay REF, each single variant, and the combination on one axis:

```python
palette = ['#e377c2', '#17becf', '#1f77b4', '#2ca02c', '#9467bd', '#8c564b']
tdata_rna, tdata_atac, colors, annotations = {}, {}, {}, []

tdata_rna['REF'] = all_predictions['combined'].reference.rna_seq.filter_to_nonpositive_strand()
tdata_atac['REF'] = all_predictions['combined'].reference.atac
colors['REF'] = 'dimgrey'

tdata_rna['ALT_Combined'] = all_predictions['combined'].alternate.rna_seq.filter_to_nonpositive_strand()
tdata_atac['ALT_Combined'] = all_predictions['combined'].alternate.atac
colors['ALT_Combined'] = 'red'

for i, (pos, var_obj) in enumerate(haplotype['individual_variants'].items()):
    key = f'ALT_{pos}'
    tdata_rna[key] = all_predictions[pos].alternate.rna_seq.filter_to_nonpositive_strand()
    tdata_atac[key] = all_predictions[pos].alternate.atac
    colors[key] = palette[i % len(palette)]
    annotations.append(plot_components.VariantAnnotation([var_obj], alpha=0.8))

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(transcript_extractor.extract(interval)),
        plot_components.OverlaidTracks(tdata=tdata_rna, colors=colors, legend_loc=None,
            ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}'),
        plot_components.OverlaidTracks(tdata=tdata_atac, colors=colors,
            ylabel_template='ATAC: {biosample_name}\n{name}'),
    ],
    interval=interval.resize(5000),
    annotations=annotations,
)
```

In the worked example (TERT promoter mutations rs1242535815 and rs1242535816),
each mutation individually creates a de novo binding site for ETS transcription
factors (GABPA). The two are almost always mutually exclusive in cancer, since
one is normally enough to grant immortality — and AlphaGenome predicts that
having both would have a *greater* effect, consistent with cooperative GABPA
binding or a wider window of increased accessibility.

Caveat: the model is not diploid-aware and was trained on unphased data, so this
is a workaround for phased analysis rather than native support.

---

## A full worked analysis: the TAL1 locus

T-cell acute lymphoblastic leukemia (T-ALL) often arises from aberrant
upregulation of the *TAL1* oncogene, driven by non-coding variants several
kilobases away. The pattern below generalises to any "are these variants doing
something real?" question.

### 1. Place the variants in context

```python
tal1_interval = genome.Interval('chr1', 47209255, 47242023, strand='-')

plot_components.plot(
    [plot_components.TranscriptAnnotation(transcript_extractor.extract(tal1_interval))],
    annotations=[plot_components.VariantAnnotation(
        [genome.Variant('chr1', x, 'N', 'N') for x in unique_positions],
        labels=labels, use_default_labels=False)],
    interval=tal1_interval,
    title='Positions of variants near TAL1',
)
```

Three positional groups emerge: a 5' cluster upstream of *TAL1* (the MuTE,
"mutation of the TAL1 enhancer", site), an intronic variant, and a 3' cluster
downstream. All three were shown in the original studies to converge on the same
mechanism — upregulation of *TAL1*.

### 2. Predict one variant across three modalities

Choose an ontology term matching the reported tissue of origin — here
`CL:0001059` ("common myeloid progenitor, CD34-positive"), a close match to the
purified CD34+ hematopoietic stem cells used in the source study.

```python
output = dna_model.predict_variant(
    interval=tal1_interval.resize(2**20),
    variant=variant,
    requested_outputs={dna_client.OutputType.RNA_SEQ,
                       dna_client.OutputType.CHIP_HISTONE,
                       dna_client.OutputType.DNASE},
    ontology_terms=['CL:0001059'],
)
```

Plot ALT − REF directly, so the panels show *change* rather than absolute level:

```python
plot_components.plot(
    [
        plot_components.TranscriptAnnotation(transcript_extractor.extract(tal1_interval)),
        plot_components.Tracks(
            tdata=output.alternate.rna_seq.filter_to_nonpositive_strand()
                - output.reference.rna_seq.filter_to_nonpositive_strand(),
            ylabel_template='{biosample_name} ({strand})\n{name}', filled=True),
        plot_components.Tracks(
            tdata=output.alternate.dnase.filter_to_nonpositive_strand()
                - output.reference.dnase.filter_to_nonpositive_strand(),
            ylabel_template='{biosample_name} ({strand})\n{name}', filled=True),
        plot_components.Tracks(
            tdata=output.alternate.chip_histone.filter_to_nonpositive_strand()
                - output.reference.chip_histone.filter_to_nonpositive_strand(),
            ylabel_template='{biosample_name} ({strand})\n{name}', filled=True),
    ],
    annotations=[plot_components.VariantAnnotation([variant])],
    interval=tal1_interval,
)
```

Reading the result: RNA-seq coverage above the horizontal means *TAL1* is
increased with the ALT allele; DNase shows changed accessibility at the variant
and increased accessibility at the *TAL1* TSS; H3K27ac and H3K4me1 increase
directly on the variant (active enhancer marks); H3K4me3 increases near the TSS
(active promoter); H3K36me3 increases over the gene body (active transcription);
H3K27me3 and H3K9me3 decrease (silencing marks). Together: the variant creates a
de novo active enhancer that activates the *TAL1* promoter.

### 3. Compare real variants against shuffled backgrounds

The key control. For each oncogenic variant, generate synthetic variants with the
**same position and the same ALT length** but shuffled sequence, then ask whether
the real variant's effect stands out from that background.

```python
def generate_background_variants(variant: genome.Variant,
                                 max_number: int = 100) -> pd.DataFrame:
    """Random ALT alleles of the same length at the same position."""
    nucleotides = np.array(list('ACGT'), dtype='<U1')
    n = len(variant.alternate_bases)

    if 4 ** n < max_number:
        permutations = [''.join(p) for p in itertools.product(nucleotides, repeat=n)]
    else:
        rng = np.random.default_rng(42)
        generated = set()
        while len(generated) < max_number:
            s = ''.join(nucleotides[rng.integers(0, 4, size=n)])
            if s != variant.alternate_bases:
                generated.add(s)
        permutations = list(generated)

    df = pd.DataFrame({
        'ID': ['mut_' + str(variant.position) + '_' + x for x in permutations],
        'CHROM': variant.chromosome, 'POS': variant.position,
        'REF': variant.reference_bases, 'ALT': permutations,
        'output': 0.0, 'original_variant': variant.name,
    })
    return df[df['REF'] != df['ALT']]
```

Score everything in one batch, then pull out the gene and cell type of interest:

```python
scores = dna_model.score_variants(
    intervals=eval_df['interval'].to_list(),
    variants=eval_df['variant'].to_list(),
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
    max_workers=2,
)

gene_index = scores[0][0].obs.query('gene_name == "TAL1"').index[0]
cell_type_index = scores[0][0].var.query('ontology_curie == "CL:0001059"').index[0]

eval_df['tal1_diff_in_cd34'] = [
    x[0][gene_index, cell_type_index].X[0, 0] for x in scores
]
```

Group variants by position and ALT length so each real variant is compared
against length-matched shuffles, then plot the real scores (highlighted) against
the background density. In the published analysis the cancer-associated variants
sit in the right tail of their matched background distributions — the specific
inserted sequence matters, not merely the insertion length.

---

## Save predictions to disk

```python
# Variant scores -> CSV.
variant_scorers.tidy_scores(results).to_csv('variant_scores.csv', index=False)

# Track predictions -> NumPy.
np.savez_compressed('rna_seq.npz', values=output.rna_seq.values)
output.rna_seq.metadata.to_csv('rna_seq_metadata.csv', index=False)
```

---

## Cache predictions during interactive work

Model calls are the slow, rate-limited part of any notebook. A dict keyed on the
call's arguments avoids repeats when you're only changing plot settings:

```python
_prediction_cache = {}

def predict_variant_cached(interval, variant, organism, requested_outputs, ontology_terms):
    key = (str(interval), str(variant), str(organism),
           tuple(requested_outputs), tuple(ontology_terms))
    if key not in _prediction_cache:
        _prediction_cache[key] = dna_model.predict_variant(
            interval=interval, variant=variant, organism=organism,
            requested_outputs=requested_outputs, ontology_terms=ontology_terms)
    return _prediction_cache[key]
```

The same trick is worth applying to `TranscriptExtractor` construction, which
involves parsing a large GTF.

---

## Filter tracks to specific transcription factors

```python
def filter_to_tfs(tdata, transcription_factors: list[str]):
    if transcription_factors is None:
        return tdata
    tf_rows = tdata.metadata.index[
        tdata.metadata['transcription_factor'].isin(transcription_factors)]
    if not tf_rows.any():
        print('No tracks found for the specified TFs and ontology terms.')
        return None
    missing = set(transcription_factors) - set(tdata.metadata['transcription_factor'])
    if missing:
        print(f'Could not find tracks for: {missing}')
    return tdata.select_tracks_by_index(tf_rows)
```

---

## Query the Atlas

Pre-computed genome-wide variant effect scores, with a higher query rate than
live prediction. Available for **human hg38 SNVs**.

```python
from alphagenome.atlas import atlas

client = atlas.create(api_key)
client.scorer_metadata()     # -> {scorer_name: ScorerMetadata(name, is_signed, track_metadata)}

# Primary composite impact score and biological feature attributions:
atlas_scores = client.query_variant(
    variant,
    requested_scorers=['AVI_SCORE', 'AVI_SCORE_FEATURE_IMPORTANCE'],
)

# Region scan (fast saturation mutagenesis across up to 1,000 bp):
region_scores = client.query_interval(
    interval,
    requested_scorers=['AVI_SCORE'],
)
```

All query methods return `Mapping[str, anndata.AnnData]`, keyed by scorer name.
Individual assay scorers (e.g. `RNA_SEQ`, `SPLICE_JUNCTIONS`, `DNASE`) can also
be queried precomputed.

**Interpreting AVI Scores:**
The Atlas calibrates tail quantiles to a Phred impact scale:
`Phred = -10 * log10(1.0 - quantile)`
- `Phred >= 40`: Top 0.01% predicted impact genome-wide.
- `Phred >= 30`: Top 0.10% predicted impact genome-wide.
- `Phred >= 20`: Top 1.00% predicted impact genome-wide.
- `Phred < 10`: Bottom 90% of variants (typically benign/neutral).

For indels or mouse mm10, use live prediction via `dna_client.create()`.
More on the dataset: <https://deepmind.google.com/science/alphagenome/atlas>.

