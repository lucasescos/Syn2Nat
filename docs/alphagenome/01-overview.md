# Overview

AlphaGenome is Google DeepMind's unifying model for deciphering the regulatory
code within DNA sequences. It makes **multimodal predictions** from raw DNA:
gene expression, splicing patterns, chromatin features, TF binding, and 3D
contact maps. It accepts sequences up to **1 million base pairs** and delivers
predictions at **single base-pair resolution** for most outputs, achieving
state-of-the-art performance across a range of genomic prediction benchmarks,
including many variant effect prediction tasks.

The API also provides access to the **AlphaGenome Atlas**: pre-computed variant
effect scores, AlphaGenome Variant Impact (AVI) scores, and feature importances
across the entire human genome
(<https://deepmind.google.com/science/alphagenome/atlas>).

## What you can ask it

| Question | Tool |
| :--- | :--- |
| What does this region look like functionally in tissue X? | `predict_interval` / `predict_sequence` |
| What does this variant change? | `predict_variant` (tracks) |
| How big is that change, as a number? | `score_variant` / `score_variants` |
| Which bases in this region matter? | `score_ism_variants` + `ism.ism_matrix` |
| Has this variant already been scored genome-wide? | Atlas (`atlas.create`) |

## Access and rate limits

- Free service, **non-commercial use only**, subject to the
  [Terms of Service](https://deepmind.google.com/science/alphagenome/terms).
- Outputs must **not** be used to train other machine-learning models.
- Query rates vary with demand. Pre-computed Atlas predictions typically allow a
  higher query rate.
- Suited to small- and medium-scale analyses — a limited number of regions or
  variants requiring thousands of predictions. **Not** suited to analyses needing
  more than ~1 million predictions.
- [Get an API key](https://deepmind.google.com/science/alphagenome).

## Limitations

- **Tissue-specificity and long-range interactions.** Improved over previous
  models, but accurately capturing tissue-specific effects and long-range
  interactions remains hard for deep learning in genomics.
- **Species scope.** Trained and evaluated on human and mouse only. Performance
  on other species is undetermined and degrades with evolutionary distance.
- **Personal genomes.** Not yet benchmarked for predicting individual (personal)
  human genomes.
- **Molecular scope.** Predicts molecular consequences of genetic variation.
  Direct applicability to complex traits is limited, since those also involve
  gene function, development, and environment.
- **Unphased training, single-sequence input.** The model sees one sequence at a
  time and is not diploid-aware. It was trained on unphased data and cannot
  distinguish maternal from paternal alleles, so variant effect predictions do
  not inherently model heterozygous states.
- **Sequence realism.** Predictions have only been evaluated on sequences close
  to the reference genome (SNPs and indels). Structural variants, heavy padding,
  synthetic sequences, and artificial constructs may give unreliable results.

## Citing

```bibtex
@article{alphagenome,
  title={Advancing regulatory variant effect prediction with {AlphaGenome}},
  author={Avsec, {\v Z}iga and Latysheva, Natasha and Cheng, Jun and Novati, Guido
    and Taylor, Kyle R. and Ward, Tom and Bycroft, Clare and Nicolaisen, Lauren
    and Arvaniti, Eirini and Pan, Joshua and Thomas, Raina and Dutordoir, Vincent
    and Perino, Matteo and De, Soham and Karollus, Alexander and Gayoso, Adam
    and Sargeant, Toby and Mottram, Anne and Wong, Lai Hong and Drot{\'a}r, Pavol
    and Kosiorek, Adam and Senior, Andrew and Tanburn, Richard
    and Applebaum, Taylor and Basu, Souradeep and Hassabis, Demis and Kohli, Pushmeet},
  journal={Nature},
  volume={649}, number={8099}, pages={1206--1218}, year={2026},
  doi={10.1038/s41586-025-10014-0},
  publisher={Nature Publishing Group UK London}
}
```

## Getting help

- Bugs and code issues: <https://github.com/google-deepmind/alphagenome/issues>
- Usage questions, feedback, feature requests: <https://www.alphagenomecommunity.com>
  (actively monitored by the team — usually the fastest route)
- Direct contact: <alphagenome@google.com>
- Video intro: [AlphaGenome 101](https://youtu.be/Xbvloe13nak)

## Other works referenced by the docs

- **Borzoi** — Linder et al., *Predicting RNA-seq coverage from DNA sequence as a
  unifying model of gene regulation*, Nature Genetics (2025). Basis of the
  polyadenylation (paQTL) scoring method.
- **SpliceAI** — Jaganathan et al., Cell (2019). Basis of the indel alignment
  strategy used in scoring.
- **GTEx Consortium**, Science (2020). Source of the `gtex_tissue` metadata column.
- **Zhou (2022)**, Nature Genetics. Contact-map normalization approach.
