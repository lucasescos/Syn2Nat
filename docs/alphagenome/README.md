# AlphaGenome — condensed documentation

A compressed version of the official AlphaGenome API docs (<https://www.alphagenomedocs.com>),
the `google-deepmind/alphagenome` repository, and the ten tutorial Colabs.
Everything you need to use the API is here; the Colabs are summarised as
runnable recipes rather than reproduced cell by cell.

**AlphaGenome** predicts functional genomic outputs — gene expression, splicing,
chromatin accessibility, TF and histone binding, and 3D contact maps — directly
from DNA sequence, for inputs up to 1 Mb, at up to single-base-pair resolution.

## Files

| File | What's in it |
| :--- | :--- |
| [01-overview.md](01-overview.md) | What the model does, access and terms, limitations, how to cite |
| [02-installation.md](02-installation.md) | pip, Colab API-key setup, local install, updating |
| [03-quickstart.md](03-quickstart.md) | Shortest end-to-end path: predict, visualize, score |
| [04-core-objects.md](04-core-objects.md) | `Interval`, `Variant`, `TrackData`, `JunctionData`, `AnnData`, ontology terms |
| [05-making-predictions.md](05-making-predictions.md) | `predict_sequence` / `predict_interval` / `predict_variant`, batching, mouse |
| [06-output-types.md](06-output-types.md) | The 11 output types, units, resolutions, track counts, metadata |
| [07-variant-scoring.md](07-variant-scoring.md) | How scoring works, all recommended scorers, `tidy_scores`, quantile scores |
| [08-splicing.md](08-splicing.md) | Merged splicing score, splicing scorers, deriving PSI |
| [09-visualization.md](09-visualization.md) | `plot_components`, per-modality plotting recipes |
| [10-interpretation-ism.md](10-interpretation-ism.md) | In silico mutagenesis and sequence logos |
| [11-recipes.md](11-recipes.md) | Batch VCF scoring, haplotypes, ontology lookup, a full worked analysis |
| [12-api-reference.md](12-api-reference.md) | Complete API surface: every public class, method, and constant |
| [13-faq.md](13-faq.md) | Condensed FAQ |

## 60-second orientation

```bash
pip install -U alphagenome
```

```python
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

dna_model = dna_client.create(API_KEY)

variant = genome.Variant('chr22', 36201698, 'A', 'C')
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

# Tracks for REF and ALT.
output = dna_model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157'],  # Colon - Transverse.
)

# A scalar effect size per gene per track.
scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
)
df = variant_scorers.tidy_scores(scores)
```

Three verbs cover almost every use: **`predict_*`** returns full tracks,
**`score_*`** collapses REF-vs-ALT into a number per track, and
**`plot_components.plot`** draws either.

## Key facts

- **Genome builds**: human hg38 (GRCh38.p13), mouse mm10 (GRCm38.p6). Use
  [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) for anything else.
- **Supported input lengths**: 2^14 (~16 kb), 2^17 (~100 kb), 2^19 (~500 kb),
  2^20 (~1 Mb). Prefer 1 Mb. Use `Interval.resize()` to snap to one.
- **Indexing**: `Interval` is 0-based half-open; `Variant.position` is 1-based
  (to match VCF/dbSNP) but `.start`/`.end` are 0-based.
- **Python**: >= 3.10.
- **Terms**: free, non-commercial use only; outputs may not be used to train
  other ML models. See <https://deepmind.google.com/science/alphagenome/terms>.

## Sources

- [AlphaGenome documentation](https://www.alphagenomedocs.com/)
- [google-deepmind/alphagenome on GitHub](https://github.com/google-deepmind/alphagenome)
- Avsec et al., *Advancing regulatory variant effect prediction with AlphaGenome*,
  Nature 649:1206–1218 (2026), [doi:10.1038/s41586-025-10014-0](https://doi.org/10.1038/s41586-025-10014-0)
- [Get an API key](https://deepmind.google.com/science/alphagenome) ·
  [Community forum](https://www.alphagenomecommunity.com)
