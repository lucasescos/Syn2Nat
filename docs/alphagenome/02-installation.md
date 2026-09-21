# Installation

Requires **Python >= 3.10**.

## From PyPI (recommended)

```bash
pip install -U alphagenome
```

Optionally create a [virtual environment](https://docs.python.org/3/tutorial/venv.html)
first to avoid conflicts with your system Python.

## Google Colab

Every tutorial notebook opens with an install cell:

```python
from IPython.display import clear_output
!pip install alphagenome
clear_output()
```

### Adding your API key to Colab secrets

To make model requests from a notebook, store the key in Colab secrets:

1. Open the notebook and click the 🔑 **Secrets** tab in the left panel.
2. Create a secret named `ALPHA_GENOME_API_KEY`.
3. Paste your API key into its **Value** box.
4. Toggle the switch on the left to allow notebook access to the secret.

Then:

```python
from alphagenome import colab_utils
from alphagenome.models import dna_client

dna_model = dna_client.create(colab_utils.get_api_key())
```

`colab_utils.get_api_key(secret='ALPHA_GENOME_API_KEY')` checks the environment
variable of that name first, then Colab secrets. Outside Colab, set the
environment variable or pass the key string directly to `dna_client.create()`.
Note that external tooling and `.env` files frequently use `ALPHAGENOME_API_KEY`
(without the second underscore); ensure you use the name expected by your code
or pass `api_key=os.environ.get('ALPHAGENOME_API_KEY')`.

## Local install from source

```bash
rm -rf ./alphagenome
git clone https://github.com/google-deepmind/alphagenome.git
pip install -e ./alphagenome
```

A virtual environment manager such as
[miniconda](https://docs.anaconda.com/miniconda/) or
[uv](https://docs.astral.sh/uv/pip/environments/) is strongly recommended. With
miniconda:

```bash
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env
pip install -e ./alphagenome
```

## Updating

With the relevant environment already activated:

```bash
cd ./alphagenome
git pull
pip install --upgrade .
```

Or, for the PyPI install: `pip install -U alphagenome`.

## Dependencies

Installed automatically: `absl-py`, `anndata`, `fsspec`, `grpcio>=1.67.1`,
`immutabledict`, `jaxtyping`, `matplotlib`, `ml_dtypes`, `numpy`, `pandas`,
`protobuf>=5.28.3`, `pyarrow`, `pyfaidx`, `scipy`, `seaborn`, `tqdm`,
`typeguard`, `typing_extensions`, `zstandard`.

So reading the GENCODE `.feather` files (`pyarrow`) and the reference FASTA
(`pyfaidx`) works out of the box. One tutorial — the TAL1 analysis — additionally
needs `plotnine`, which is not a dependency.

## Reference annotation files

The tutorials pull GENCODE annotations from public Google Cloud Storage:

```python
HG38_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)
MM10_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'mm10/gencode.vM23.annotation.gtf.gz.feather'
)
gtf = pd.read_feather(HG38_GTF_FEATHER)
```

A reference FASTA (used by `io.fasta.FastaExtractor`) lives at
`https://storage.googleapis.com/alphagenome/reference/gencode/hg38/GRCh38.p13.genome.fa`.

> **Memory advisory:** The GENCODE feather file is ~318 MB compressed, but expands
> to **~4.4 GB in RAM** once loaded into a DataFrame by pandas/pyarrow. On systems
> with limited memory, filter the DataFrame immediately (e.g. using
> `gene_annotation.filter_to_mane_select_transcript(gtf)`) or cache the filtered
> subset to disk to avoid Out-Of-Memory kernel terminations. Loading from remote
> `https://` URLs requires `fsspec` installed.

