# Making predictions

## Creating a client

```python
from alphagenome.models import dna_client

dna_model = dna_client.create(api_key)
# optional: model_version=..., timeout=..., address=...
```

The default API address is `dns:///gdmscience.googleapis.com:443` (communicating
via gRPC over TLS on port 443). If your environment requires an HTTP/gRPC proxy,
configure standard `grpc_proxy` / `https_proxy` environment variables.

In Colab: `dna_client.create(colab_utils.get_api_key())`.

`ModelVersion` selects which training fold to use: `ALL_FOLDS` (default),
`FOLD_0`, `FOLD_1`, `FOLD_2`, `FOLD_3`.

## Sequence-length constraints

```python
dna_client.SEQUENCE_LENGTH_16KB    # 2**14
dna_client.SEQUENCE_LENGTH_100KB   # 2**17
dna_client.SEQUENCE_LENGTH_500KB   # 2**19
dna_client.SEQUENCE_LENGTH_1MB     # 2**20  (the maximum, and the recommended default)
dna_client.SUPPORTED_SEQUENCE_LENGTHS   # dict of name -> value
dna_client.validate_sequence_length(length)
```

Anything else raises. Snap to a supported length with `Interval.resize()`, which
expands using real genomic sequence rather than padding.

## The five prediction methods

| Method | Input | Returns |
| :--- | :--- | :--- |
| `predict_sequence` | a DNA string | `Output` |
| `predict_interval` | a `genome.Interval` (reference genome) | `Output` |
| `predict_variant` | interval + `genome.Variant` | `VariantOutput` (`.reference`, `.alternate`) |
| `score_variant` | interval + variant + scorers | `list[AnnData]`, one per scorer |
| `score_variants` | many intervals + variants + scorers | `list[list[AnnData]]` |

Plural forms — `predict_sequences`, `predict_intervals`, `predict_variants`,
`score_intervals` — run batches concurrently with a progress bar.

### Common arguments

- `requested_outputs` — iterable of `dna_client.OutputType` (for `predict_*` methods).
  Request only what you need; compute cost scales with it.
- `ontology_terms` — iterable of CURIE strings or `OntologyTerm`s, or `None` for
  all tracks (for `predict_*` methods). `SPLICE_SITES` is tissue-agnostic and
  ignores this filter. **Note:** `score_variant` and `score_variants` do **not**
  accept `ontology_terms` — they score all available tracks for the requested
  scorers; filter the returned `AnnData` or tidy DataFrame afterwards.
- `organism` — `Organism.HOMO_SAPIENS` (default) or `Organism.MUS_MUSCULUS`.
- `progress_bar`, `max_workers` (default 5) — on the batch methods.
- `merge_stranded_gene_tracks` (default `True`) — on `score_variant`,
  `score_interval` and `score_ism_variants` only. The batch forms
  `score_variants` / `score_intervals` do not accept it.


### `predict_sequence`

```python
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    requested_outputs=[dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048'],
    interval=None,      # optional: attach genomic coordinates to the result
)
```

Valid characters are `ACGTN`. You can predict on any sequence within a supported
length, but predictions have only been evaluated on sequences close to the
reference genome — heavy padding, synthetic sequences and structural variants
are unreliable.

### `predict_interval`

```python
output = dna_model.predict_interval(
    interval=interval,                       # must already be a supported width
    requested_outputs={dna_client.OutputType.RNA_SEQ, dna_client.OutputType.CAGE},
    ontology_terms=['UBERON:0001159', 'UBERON:0001155'],
)
```

### `predict_variant`

```python
variant_output = dna_model.predict_variant(
    interval=variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157'],
)
ref, alt = variant_output.reference, variant_output.alternate
```

Both are `Output` objects with the same structure, so you can diff them:
`alt.rna_seq - ref.rna_seq`.

### Batch prediction

```python
outputs = dna_model.predict_variants(
    intervals=intervals,        # one Interval, or one per variant
    variants=variants,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['CL:0000084'],
    max_workers=2,
)
```

## The `Output` object

One optional `TrackData` field per output type, plus `splice_junctions` which is
a `JunctionData`:

```python
output.atac, output.cage, output.dnase, output.rna_seq
output.chip_histone, output.chip_tf
output.splice_sites, output.splice_site_usage, output.splice_junctions
output.contact_maps, output.procap
```

Fields you did not request are `None`. Helpers:

```python
output.get(dna_client.OutputType.RNA_SEQ)
output.filter_to_strand('-')                  # '+', '-' or '.'
output.filter_ontology_terms(terms)
output.filter_output_type([OutputType.ATAC])
output.resize(width)
output.map_track_data(fn)                     # apply fn to every TrackData
```

## Output metadata

```python
output_metadata = dna_model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)

output_metadata.rna_seq          # DataFrame, one row per RNA-seq track
output_metadata.get(dna_client.OutputType.ATAC)
output_metadata.concatenate()    # every output type in one DataFrame
```

See [06-output-types.md](06-output-types.md) for the metadata columns.

## Mouse

```python
output = dna_model.predict_sequence(
    sequence=...,
    organism=dna_client.Organism.MUS_MUSCULUS,
    requested_outputs=[dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048'],
)
```

Caveats: supported ontology terms differ between species, `PROCAP` has no mouse
tracks, and the polyadenylation scorer is human-only. Coordinates are mm10.

## Saving results

- **Variant scores** — convert to a DataFrame with
  `variant_scorers.tidy_scores(...)`, then `.to_csv(...)`.
- **Track predictions** — `TrackData.values` is a plain NumPy array; use
  `numpy.save` (`.npy`) or `numpy.savez_compressed` (`.npz`) and keep
  `tdata.metadata` alongside as CSV.

## Retries and errors

The client retries RPCs that fail with `RESOURCE_EXHAUSTED` or `UNAVAILABLE`
(5 attempts, exponential backoff with jitter) via `dna_client.retry_rpc`. Other
gRPC errors surface as Python exceptions. If you are hitting rate limits, reduce
`max_workers`, request fewer output types, or use shorter sequence lengths.

Other relevant limits: `MAX_VARIANT_SCORERS_PER_REQUEST = 20`, and
`MAX_ISM_INTERVAL_WIDTH = 10` — the client automatically splits a wider ISM
interval into chunks of this width.
