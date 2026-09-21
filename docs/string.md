# STRING Database API: Architecture & Programmatic Access Guide

> **Official Portal**: [https://string-db.org](https://string-db.org)  
> **API Documentation**: [https://string-db.org/help/api/](https://string-db.org/help/api/)  
> **MCP Server**: [https://mcp.string-db.org/](https://mcp.string-db.org/) (GitHub: [meringlab/string-mcp](https://github.com/meringlab/string-mcp))  
> **Current Version**: `v12.0` (Stable URL: `https://version-12-0.string-db.org`)  
> **Terms & Licensing**: CC-BY-4.0 for academic/research use ([https://string-db.org/cgi/access](https://string-db.org/cgi/access))  
> **Related Repository Guides**: [HIPPIE Database (`docs/hippie.md`)](hippie.md) | [PPI Hierarchy (`docs/ppi_evaluation_hierarchy.md`)](ppi_evaluation_hierarchy.md)

---

## 1. Overview & System Architecture

**STRING** (Search Tool for the Retrieval of Interacting Genes/Proteins) is one of the world's most widely used repositories of direct (physical) and indirect (functional) protein-protein associations. Developed by a consortium including the University of Zurich (von Mering group), EMBL, CPR/University of Copenhagen, and SIB Swiss Institute of Bioinformatics, STRING maps hundreds of millions of interactions across thousands of organisms.

### Core Distinctions
* **Functional vs. Physical Interactions**: STRING captures both *physical binding* (complexes, binary direct interactions) and *functional linkages* (proteins participating in the same metabolic pathway, shared operon, or regulatory cascade).
* **Cross-Species Evidence Transfer**: Benchmarked evidence from model organisms is orthology-transferred across the entire phylogenetic tree using SIMAP/Smith-Waterman alignments.

### Programmatic Access Channels

```
┌────────────────────────────────────────────────────────────────────────┐
│                        STRING Database Engine                          │
└──────┬──────────────────────┬──────────────────────┬───────────────────┘
       │                      │                      │
       ▼                      ▼                      ▼
┌──────────────┐       ┌──────────────┐       ┌──────────────┐
│  Synchronous │       │ Asynchronous │       │  Model       │
│  REST API    │       │ Values/Ranks │       │  Context     │
│  (/api/...)  │       │ GSEA API     │       │  Protocol    │
└──────┬───────┘       └──────┬───────┘       └──────┬───────┘
       │                      │                      │
  [TSV / JSON]         [API Key Queue]       [Agent Tools]
  - Map IDs            - Full logFC lists    - Anthropic Claude
  - Networks           - KS / AFC Tests      - Custom GPTs
  - Homology           - Bidirectional       - AI Pipelines
  - ORA Enrichment       GSEA Enrichment
```

1. **Synchronous REST API**: Fast, on-demand querying for identifier resolution, neighborhood network extraction, pairwise homology, and over-representation analysis (ORA).
2. **Asynchronous Values/Ranks API**: Background job queuing for full genome-wide differential expression / ranked dataset analysis (GSEA-equivalent using Kolmogorov-Smirnov and Aggregate Fold Change statistics).
3. **Model Context Protocol (MCP) Server**: Remote and self-hostable server (`https://mcp.string-db.org/`) providing structured, LLM-optimized tools for autonomous AI agents.
4. **SPARQL / RDF Endpoint**: Semantic web queries against RDF representations of STRING data.
5. **Bulk Flat Files**: Complete proteome-wide data dumps available from the [STRING Download Page](https://string-db.org/cgi/download) for large-scale offline matrix computations.

---

## 2. API Request Conventions & Communication Rules

### 2.1 Base URLs & Version Pinning
* **Development / Dynamic Endpoint**:  
  `https://string-db.org/api/[format]/[method]`
* **Production Version-Pinned Endpoint (Mandatory for Reproducibility)**:  
  `https://version-12-0.string-db.org/api/[format]/[method]`

> [!IMPORTANT]
> Always use the version-pinned URL (`https://version-12-0.string-db.org/api/...`) in analytical pipelines and published code. The unversioned root (`https://string-db.org`) will silently change when a new major STRING release is deployed.

### 2.2 Output Formats
The `[format]` path token controls the serialization format:

| Format Token | Output Type | Typical Application |
|---|---|---|
| `tsv` | Tab-separated text with headers | Standard data pipelines, pandas, R dataframes |
| `tsv-no-header`| Tab-separated text without header | Stream processing, CLI awk/grep filtering |
| `json` | Structured JSON array of objects | Web apps, REST services, automated parsing |
| `xml` | Standard XML document | Legacy enterprise pipelines |
| `psi-mi` | Proteomics Standards Initiative XML | Standardized proteomics interchange |
| `psi-mi-tab` | PSI-MITAB tabular format | Molecular interaction community standards |
| `image` | Standard PNG image (alpha transparency) | Web visualization, embedded graphics |
| `highres_image`| High-resolution PNG image | Publication figures, poster renders |
| `svg` | Scalable Vector Graphics | Interactive UI, vector graphics editing |

### 2.3 HTTP Methods & Delimiters
* **GET Requests**: Acceptable for quick interactive lookups or single proteins. Delimit multiple identifiers with URL-encoded carriage returns: `%0d` or `%0a` (e.g. `identifiers=TP53%0dCDK2`).
* **POST Requests**: **Strongly recommended** for any batch analysis or pipelines with $\ge 5$ proteins. Long query strings in GET requests easily exceed web server URL length limits (HTTP 414).
  * In POST bodies (`application/x-www-form-urlencoded`), separate identifiers using `\r`, `\n`, or `%0d`.

### 2.4 Mandatory & Core Parameters

| Parameter | Required? | Description & Constraints |
|---|---|---|
| `identifiers` | **Yes** (most endpoints) | One or more gene symbols, UniProt accessions, or STRING internal IDs separated by `\r` / `%0d`. |
| `species` | **Mandatory for >10 proteins** | NCBI Taxonomy ID (integer, e.g. `9606` for Human, `10090` for Mouse, `7227` for Fly, `4932` for Yeast). Always specify this parameter to prevent expensive multi-genome disambiguation. |
| `caller_identity` | **Required by etiquette** | A reverse-DNS or project identifier for your tool/script (e.g., `my_microprotein_pipeline.org`). Used by STRING operators to monitor usage without blocking IPs. |
| `required_score` | Optional | Interaction confidence threshold on a **0 to 1000** integer scale (default: `400`). |

### 2.5 Server Politeness & Rate Limits
1. **Pacing**: Introduce a delay of at least **1.0 second** between successive API requests (`sleep(1)`).
2. **Concurrency**: Avoid parallel multi-threaded hammering of the endpoints. STRING web servers will throttle or temporarily ban abusive IP addresses.
3. **404 Errors**: If none of the input identifiers can be resolved to valid proteins in the specified organism, STRING returns an HTTP **404 Not Found**. Always check HTTP status codes.

---

## 3. Interaction Scoring Scheme & Channels

STRING calculates a **combined score** ($S \in [0, 1000]$) by integrating distinct biological evidence channels. Scores represent the probability that an interaction is biologically meaningful relative to a random pair of proteins.

### Confidence Tiers
* **Low Confidence**: $150 \le \text{score} < 400$
* **Medium Confidence**: $400 \le \text{score} < 700$ (Standard baseline)
* **High Confidence**: $700 \le \text{score} < 900$
* **Highest Confidence**: $900 \le \text{score} \le 1000$

### Evidence Channel Breakdown

```
Combined Score [score]
  ├── Experimental Evidence    [escore] : In vitro / in vivo biochemical assays, affinity-MS, Y2H
  ├── Curated Databases        [dscore] : Annotated pathways (KEGG, Reactome, BioCyc, etc.)
  ├── Automated Textmining     [tscore] : Co-occurrence in PubMed abstracts and PMC full-text
  ├── Co-expression            [ascore] : Correlated gene expression across RNA-seq / microarrays
  ├── Conserved Neighborhood   [nscore] : Operon-like genomic proximity in prokaryotes
  ├── Gene Fusion              [fscore] : Orthologs fused into a single polypeptide chain
  └── Phylogenetic Profiles    [pscore] : Co-presence / co-absence across genomes
```

### Probabilistic Combination Formula
Assuming independent evidence channels, individual subscores are combined using a naive Bayes formulation corrected for the random background interaction probability ($p_0$):

$$S = 1 - (1 - p_0) \prod_{k \in \text{channels}} \left( 1 - \frac{S_k - p_0}{1 - p_0} \right)$$

---

## 4. Complete REST API Endpoints Reference

### 4.1 Identifier Mapping (`get_string_ids`)
Maps ambiguous user inputs (gene symbols, synonyms, UniProt accessions, Ensembl gene/protein IDs) to STRING's internal canonical protein identifiers (`stringId`).

* **Endpoint**: `/api/[format]/get_string_ids`
* **Method**: `GET` / `POST`
* **Key Parameters**:
  * `identifiers`: Input terms (separated by `\r` or `%0d`).
  * `species`: NCBI Taxon ID (e.g. `9606`).
  * `echo_query`: `1` to include the user's original input term in the output table, `0` otherwise (default: `0`).
* **Returned Fields**:
  * `queryItem`: User-supplied input name.
  * `queryIndex`: 0-indexed position in user's query list.
  * `stringId`: Canonical STRING identifier (e.g. `9606.ENSP00000269305`).
  * `ncbiTaxonId`: NCBI taxon ID.
  * `taxonName`: Species common/scientific name.
  * `preferredName`: Canonical gene symbol (e.g. `TP53`).
  * `annotation`: Short protein functional summary.

```bash
# Example: Map human TP53, CDK2, and UniProt P04637
curl -d "identifiers=TP53%0dCDK2%0dP04637&species=9606&echo_query=1&caller_identity=microprotein_doc" \
     "https://version-12-0.string-db.org/api/tsv/get_string_ids"
```

---

### 4.2 Interaction Network Retrieval (`network`)
Retrieves the subnetwork among query proteins, or expands the neighborhood around them.

* **Endpoint**: `/api/[format]/network`
* **Method**: `GET` / `POST`
* **Key Parameters**:
  * `identifiers`: Protein IDs (preferably canonical `stringId`).
  * `species`: NCBI Taxon ID.
  * `required_score`: Minimum combined score threshold (0-1000).
  * `network_type`: `functional` (default; associations + physical) or `physical` (physical binding only).
  * `add_nodes`: Number of additional high-confidence interaction partners to pull into the subnetwork (default: `0` for $\ge 2$ input proteins; automatically defaults to `10` if only `1` protein is submitted).
  * `network_term_id`: Optional functional term ID (e.g. GO or KEGG ID); when supplied, STRING builds the network of proteins bearing that term instead of an explicit protein list.
* **Returned Fields**:
  * `stringId_A`, `stringId_B`: Canonical STRING IDs for both interactors.
  * `preferredName_A`, `preferredName_B`: Gene symbols for interactors.
  * `ncbiTaxonId`: Species taxonomy ID.
  * `score`: Combined confidence score ($0.0 - 1.0$ in float format).
  * Channel scores: `escore`, `dscore`, `tscore`, `ascore`, `nscore`, `fscore`, `pscore`.

```bash
# Example: Physical interactions between human CDK1, CDK2, and Cyclins with high confidence (>=0.7)
curl -d "identifiers=CDK1%0dCDK2%0dCCNA2%0dCCNB1&species=9606&required_score=700&network_type=physical&caller_identity=microprotein_doc" \
     "https://version-12-0.string-db.org/api/tsv/network"
```

---

### 4.3 Interaction Partners of a Protein Set (`interaction_partners`)
Unlike `network` (which looks for edges *between* query proteins), `interaction_partners` finds the top interaction partners across the **entire proteome** for each submitted protein.

* **Endpoint**: `/api/[format]/interaction_partners`
* **Method**: `GET` / `POST`
* **Key Parameters**:
  * `identifiers`: Input protein identifiers.
  * `species`: NCBI Taxon ID.
  * `limit`: Maximum number of top interaction partners returned per query protein (sorted descending by confidence score).
  * `required_score`: Minimum score cutoff (0-1000).
  * `network_type`: `functional` or `physical`.
* **Returned Fields**:
  * Same columns as `network` (`stringId_A`, `stringId_B`, `preferredName_A`, `preferredName_B`, `score`, channel scores).

```bash
# Example: Find top 5 interactors for human BRCA1
curl -d "identifiers=BRCA1&species=9606&limit=5&caller_identity=microprotein_doc" \
     "https://version-12-0.string-db.org/api/tsv/interaction_partners"
```

---

### 4.4 Homology & Sequence Similarity (`homology` & `homology_best`)
STRING precomputes pairwise Smith-Waterman bit scores via SIMAP.

#### Pairwise Within-Species Homology (`homology`)
* **Endpoint**: `/api/[format]/homology`
* **Output**: `ncbiTaxonId_A`, `stringId_A`, `ncbiTaxonId_B`, `stringId_B`, `bitscore`.
* Returns half of the symmetric similarity matrix (plus self-hits) for scores $\ge 50$.

#### Cross-Species Best Hit Homology (`homology_best`)
* **Endpoint**: `/api/[format]/homology_best`
* **Key Parameters**:
  * `identifiers`: Query proteins in query species.
  * `species`: Query species taxon ID.
  * `species_b`: Comma- or `%0d`-delimited list of target organism taxon IDs (e.g. `10090%0d7227` for Mouse and Drosophila).
* **Output**: Highest-scoring ortholog/homolog in each target organism.

```bash
# Example: Best mouse (10090) homolog for human TP53
curl -d "identifiers=TP53&species=9606&species_b=10090&caller_identity=microprotein_doc" \
     "https://version-12-0.string-db.org/api/tsv/homology_best"
```

---

### 4.5 Functional Enrichment Analysis (`enrichment`)
Performs statistical Over-Representation Analysis (ORA) against functional annotation libraries:
* Gene Ontology: Biological Process (`Process`), Molecular Function (`Function`), Cellular Component (`Component`)
* Pathways: `KEGG`, Reactome (`RCTM`), `WikiPathways`
* Domains: `Pfam`, `InterPro`, `SMART`
* Phenotypes & Diseases: `HPO` (Human Phenotype Ontology), `DISEASES`, `COMPARTMENTS`, `TISSUES`
* Literature & Clusters: `PMID` (PubMed citations), `NetworkNeighborAL` (STRING network clusters)

* **Endpoint**: `/api/[format]/enrichment`
* **Method**: `GET` / `POST`
* **Key Parameters**:
  * `identifiers`: Protein set to test.
  * `species`: NCBI Taxon ID.
  * `background_string_identifiers`: Optional custom background proteome (e.g. only proteins detected in your cell type/mass-spec screen). Must be `%0d`-delimited canonical STRING IDs.
* **Returned Fields**:
  * `category`: Source database (e.g. `Process`, `KEGG`, `Pfam`).
  * `term`: Term identifier (e.g. `GO:0006281`).
  * `description`: Term text.
  * `number_of_genes`: Count of input proteins annotated with term.
  * `number_of_genes_in_background`: Total count of genes in proteome with term.
  * `p_value`: Raw hypergeometric p-value.
  * `fdr`: Benjamini-Hochberg False Discovery Rate.
  * `inputGenes`, `preferredNames`: List of matching input proteins.

```bash
# Example: Functional enrichment for human DNA repair factors
curl -d "identifiers=BRCA1%0dRAD51%0dATM%0dCHEK2%0dTP53&species=9606&caller_identity=microprotein_doc" \
     "https://version-12-0.string-db.org/api/tsv/enrichment"
```

---

### 4.6 Functional Annotations Retrieval (`functional_annotation`)
Retrieves all annotations mapped to the query proteins (regardless of statistical enrichment).

* **Endpoint**: `/api/[format]/functional_annotation`
* **Method**: `GET` / `POST`
* **Key Parameters**:
  * `identifiers`, `species`
  * `allow_pubmed`: `1` to include PubMed co-citation annotations (default: `0` due to large size).
  * `only_pubmed`: `1` to exclusively return PubMed literature links.
* **Returned Fields**: `category`, `term`, `description`, `number_of_genes`, `ratio_in_set`, `preferredNames`.

---

### 4.7 Reverse Term Lookup (`functional_terms`)
Finds all proteins associated with a keyword, disease, or ontology ID in a given species.

* **Endpoint**: `/api/[format]/functional_terms`
* **Method**: `GET` / `POST`
* **Key Parameters**:
  * `term_text`: Free-text phrase or ID (e.g. `"Melanoma"`, `"GO:0007049"`, `"MAP kinase"`).
  * `species`: Single NCBI Taxon ID (default: `9606`).
* **Returned Fields**: `category`, `term`, `description`, `proteinCount`, `preferredNames` (top 100), `stringIds` (top 100).

```bash
# Example: Find human proteins associated with "Ribosome biogenesis"
curl -d "term_text=Ribosome+biogenesis&species=9606&caller_identity=microprotein_doc" \
     "https://version-12-0.string-db.org/api/tsv/functional_terms"
```

---

### 4.8 Automated Gene Set Thematic Summary (`geneset_description`)
Performs automatic functional summarization of an unlabelled gene cluster (e.g. single-cell cluster or co-expression module) and returns three concise thematic descriptions.

* **Endpoint**: `/api/[format]/geneset_description`
* **Returned Fields**: `primary_description`, `secondary_description`, `tertiary_description`.

---

### 4.9 Protein-Protein Interaction Enrichment (`ppi_enrichment`)
Tests whether a group of proteins has significantly **more physical or functional edges** than expected by chance for a random subset of proteins with identical degree distributions.

* **Endpoint**: `/api/[format]/ppi_enrichment`
* **Returned Fields**:
  * `number_of_nodes`: Protein node count.
  * `number_of_edges`: Observed network edge count.
  * `average_node_degree`: Mean degree per node.
  * `local_clustering_coefficient`: Mean clustering coefficient.
  * `expected_number_of_edges`: Null-model expected edge count.
  * `p_value`: Statistical significance of PPI enrichment.

---

### 4.10 Network & Figure Visualization Endpoints

STRING can generate ready-to-render image files directly via HTTP:

1. **Network Diagram**:  
   `/api/[image|highres_image|svg]/network`
   * Visual options:
     * `network_flavor`: `evidence` (colored by channel), `confidence` (line thickness = score), or `actions` (arrows for activation/inhibition).
     * `hide_node_labels`: `0` or `1`.
     * `hide_disconnected_nodes`: `0` or `1`.
     * `flat_node_design`: `1` to disable 3D bubbles.
     * `add_color_nodes`: Add $N$ colored interaction partners.
     * `add_white_nodes`: Add $N$ white peripheral interactors.

2. **Enrichment Dot/Bubble Plot**:  
   `/api/[image|highres_image|svg]/enrichmentfigure`
   * Visual options:
     * `category`: E.g. `Process`, `KEGG`, `RCTM`, `Function`.
     * `x_axis`: Sorting and X-axis metric (`signal`, `strength`, `FDR`, or `gene_count`).
     * `group_by_similarity`: Clustering similarity cutoff ($0.1 - 1.0$).
     * `color_palette`: `mint_blue`, `lime_emerald`, `green_blue`, `yellow_pink`, etc.

3. **Interactive JavaScript Embedding**:
   Add STRING's standalone interactive network viewer directly into any HTML dashboard:
   ```html
   <script type="text/javascript" src="https://string-db.org/javascript/combined_embedded_network_v2.0.6.js"></script>
   <div id="stringEmbedded"></div>
   <script>
     getSTRING('https://string-db.org', {
       'species': '9606',
       'identifiers': ['TP53', 'MDM2', 'ATM', 'CHEK2'],
       'network_flavor': 'confidence',
       'caller_identity': 'microprotein_app'
     });
   </script>
   ```

---

## 5. Asynchronous Values/Ranks Enrichment API (GSEA-Style)

When you have a complete genome-wide experimental dataset (e.g. all 20,000 human genes with $\log_2 \text{FC}$, $t$-statistics, or negative $\log_{10} p$-values from RNA-seq / proteomics), traditional cut-off-based ORA misses subtle pathway-wide shifts. STRING provides an asynchronous, queue-based GSEA equivalent.

### Method Characteristics
* **Zero Thresholding**: Evaluates the full ranked distribution.
* **Dual Statistical Engines**: Combines **Kolmogorov-Smirnov (KS)** tests (fast, robust for large pathways) with **Aggregate Fold Change (AFC)** tests (high sensitivity for subtle shifts).
* **Bidirectional Detection**: Identifies pathways enriched at the top ($\uparrow$), bottom ($\downarrow$), or both extremes simultaneously (pathways featuring bidirectional perturbation).

### Execution Workflow

```
[ Step 1: Obtain Anonymous Key ]
          │
          ▼
   get_api_key  ───► Returns api_key (Active within 30 min, allows 1,000 jobs)
          │
          ▼
[ Step 2: Submit Dataset ]
          │
          ▼
   valuesranks_enrichment_submit  ───► Returns job_id
          │
          ▼
[ Step 3: Poll Status & Download ]
          │
          ▼
   valuesranks_enrichment_status (loop while status == 'queued' or 'running')
          │
          ├──► status: "success"
          └──► download_url (TSV table) & graph_url (Plot)
```

#### Step 1: Obtain API Key
```bash
curl "https://version-12-0.string-db.org/api/json/get_api_key"
# Response: [{"api_key": "YOUR_KEY", "note": "This key will be activated within 30 minutes."}]
```

#### Step 2: Prepare Input Data & Submit
Input file must be headerless, tab-separated `[Gene_Identifier]\t[Numeric_Value]`:
```text
TP53	3.42
CDK2	-1.85
BRCA1	2.11
MDM2	-2.90
```

Submit via POST:
```bash
curl -X POST "https://version-12-0.string-db.org/api/json/valuesranks_enrichment_submit" \
     -d "api_key=YOUR_KEY" \
     -d "species=9606" \
     -d "ge_fdr=0.05" \
     -d "caller_identity=microprotein_doc" \
     --data-urlencode "identifiers@my_experiment.tsv"
# Response: [{"job_id": "b0TyfmJPSkIs", "status": "queued"}]
```

#### Step 3: Check Status & Retrieve
```bash
curl "https://version-12-0.string-db.org/api/json/valuesranks_enrichment_status?api_key=YOUR_KEY&job_id=b0TyfmJPSkIs"
```
When `status == "success"`, the response contains:
* `download_url`: TSV table of all enriched pathways with FDR, direction, and gene counts.
* `page_url`: Direct interactive web view on STRING.
* `graph_url`: Auto-generated enrichment bubble plot.

---

## 6. AI Agent Integration: Model Context Protocol (MCP)

In addition to traditional REST endpoints, STRING provides a dedicated **Model Context Protocol (MCP)** server tailored specifically for autonomous AI agents (such as Google Antigravity, Anthropic Claude, and custom GPTs).

* **Live Public Endpoint**: `https://mcp.string-db.org/`
* **Open-Source Repository**: [https://github.com/meringlab/string-mcp](https://github.com/meringlab/string-mcp) (MIT License)

### Why MCP vs. Raw REST for AI Agents?
1. **Curated Token Footprint**: Raw STRING tables can return tens of thousands of rows, exhausting agent context windows. MCP tools paginate, filter, and truncate outputs to preserve context.
2. **Self-Describing Metadata**: Tool definitions contain extensive prompt engineering, schema validations, and biological hints so LLMs avoid hallucinations (e.g. reminding models to resolve taxon IDs first).
3. **Chained Multi-Step Workflows**: Built-in tools allow an agent to take a user prompt like *"Find why p53 is phosphorylated during DNA damage"*, map identifiers, extract kinase interaction partners, run GO enrichment, and return a synthesized summary with image links.

---

## 7. Production Python Client Implementation

Below is a robust, production-ready Python client demonstrating identifier mapping, subnetwork retrieval, functional enrichment, and rate-limiting.

```python
#!/usr/bin/env python3
"""STRING Database Python API Client (STRING v12.0).

Handles identifier disambiguation, network interactions, and enrichment
with polite rate-limiting and error handling.
"""

from __future__ import annotations

import io
import time
from typing import Any, Dict, List, Literal, Optional
import pandas as pd
import requests

STRING_VERSION = "12-0"
BASE_URL = f"https://version-{STRING_VERSION}.string-db.org/api"


class StringDbClient:
    """Client for querying the STRING database REST API."""

    def __init__(self, caller_identity: str = "microprotein_project", min_delay: float = 1.0):
        self.caller_identity = caller_identity
        self.min_delay = min_delay
        self._last_request_time = 0.0

    def _wait_for_rate_limit(self) -> None:
        elapsed = time.time() - self._last_request_time
        if elapsed < self.min_delay:
            time.sleep(self.min_delay - elapsed)

    def _post(self, endpoint: str, format_type: str, params: Dict[str, Any]) -> requests.Response:
        self._wait_for_rate_limit()
        url = f"{BASE_URL}/{format_type}/{endpoint}"
        params["caller_identity"] = self.caller_identity

        response = requests.post(url, data=params, timeout=60)
        self._last_request_time = time.time()

        if response.status_code == 404:
            raise ValueError(f"STRING could not resolve any provided identifiers for species {params.get('species')}.")
        response.raise_for_status()
        return response

    def map_identifiers(self, identifiers: List[str], species: int = 9606) -> pd.DataFrame:
        """Map common gene symbols / IDs to canonical STRING identifiers."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "echo_query": 1,
        }
        res = self._post("get_string_ids", "tsv", params)
        return pd.read_csv(io.StringIO(res.text), sep="\t")

    def get_network(
        self,
        identifiers: List[str],
        species: int = 9606,
        required_score: int = 400,
        network_type: Literal["functional", "physical"] = "functional",
        add_nodes: int = 0,
    ) -> pd.DataFrame:
        """Retrieve interaction edges between the provided proteins."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "required_score": required_score,
            "network_type": network_type,
            "add_nodes": add_nodes,
        }
        res = self._post("network", "tsv", params)
        return pd.read_csv(io.StringIO(res.text), sep="\t")

    def get_interaction_partners(
        self,
        identifiers: List[str],
        species: int = 9606,
        limit: int = 10,
        required_score: int = 700,
        network_type: Literal["functional", "physical"] = "physical",
    ) -> pd.DataFrame:
        """Retrieve top interaction partners across the proteome."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "limit": limit,
            "required_score": required_score,
            "network_type": network_type,
        }
        res = self._post("interaction_partners", "tsv", params)
        return pd.read_csv(io.StringIO(res.text), sep="\t")

    def get_enrichment(
        self,
        identifiers: List[str],
        species: int = 9606,
        background_ids: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """Perform Over-Representation Analysis (ORA) across GO, KEGG, Pfam, etc."""
        params: Dict[str, Any] = {
            "identifiers": "\r".join(identifiers),
            "species": species,
        }
        if background_ids:
            params["background_string_identifiers"] = "%0d".join(background_ids)

        res = self._post("enrichment", "tsv", params)
        return pd.read_csv(io.StringIO(res.text), sep="\t")

    def download_network_image(
        self,
        identifiers: List[str],
        output_path: str,
        species: int = 9606,
        high_res: bool = True,
        network_flavor: Literal["evidence", "confidence", "actions"] = "confidence",
    ) -> None:
        """Download rendered PNG network graphic."""
        fmt = "highres_image" if high_res else "image"
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "network_flavor": network_flavor,
        }
        res = self._post("network", fmt, params)
        with open(output_path, "wb") as f:
            f.write(res.content)


if __name__ == "__main__":
    client = StringDbClient()

    # 1. Map gene symbols
    genes = ["TP53", "MDM2", "CDKN1A", "ATM"]
    print("Mapping genes...")
    mapped_df = client.map_identifiers(genes, species=9606)
    print(mapped_df[["queryItem", "stringId", "preferredName", "annotation"]].head())

    # 2. Extract physical interaction subnetwork (high confidence >= 0.7)
    string_ids = mapped_df["stringId"].tolist()
    print("\nRetrieving physical interactions...")
    net_df = client.get_network(string_ids, species=9606, required_score=700, network_type="physical")
    print(net_df[["preferredName_A", "preferredName_B", "score", "escore", "dscore"]])

    # 3. Perform functional enrichment
    print("\nRunning pathway enrichment...")
    enrich_df = client.get_enrichment(string_ids, species=9606)
    sig_enrich = enrich_df[enrich_df["fdr"] < 0.01][["category", "term", "description", "fdr", "preferredNames"]]
    print(sig_enrich.head())
```

---

## 8. Best Practices, Pitfalls & Comparison Matrix

### 8.1 Critical Implementation Checklist
1. **Always Disambiguate First**: Never feed raw, user-typed gene symbols directly into downstream network or enrichment endpoints in a tight loop. Always pass them through `get_string_ids` first. Canonical IDs (`stringId`) bypass backend database lookup caches, drastically cutting query response times from $\sim 800\text{ms}$ to $<50\text{ms}$.
2. **Never Omit Species**: Queries with $>10$ proteins without a `species` argument are rejected with an error. Even for smaller queries, omitting species triggers expensive disambiguation across thousands of taxa.
3. **Delimiter Precision**: Do not separate identifiers with commas in GET requests. STRING expects `%0d` (carriage return `\r`) or `%0a` (newline `\n`). In Python `requests.post(data=...)`, pass `"\r".join(ids)`.
4. **Single-Protein Auto-Expansion Trap**: Note that submitting a single protein to `/network` or `/enrichment` automatically sets `add_nodes=10` or `add_white_nodes=10` to display its local neighborhood. If you want enrichment or interactions for *only* your query protein, check the returned `inputGenes` column.
5. **Bulk Data Boundaries**: If you need to construct a full proteome adjacency matrix or perform graph embedding (e.g. DeepWalk, Node2Vec, or MAPPIE interface models), **do not scrape the API**. Download `protein.links.v12.0.txt.gz` from [https://string-db.org/cgi/download](https://string-db.org/cgi/download).

### 8.2 Comparison: STRING vs. HIPPIE

| Feature / Dimension | STRING Database | HIPPIE Database (`docs/hippie.md`) |
|---|---|---|
| **Scope of Interactions** | Functional linkages + Direct physical bindings | Strictly direct physical interactions |
| **Taxonomic Coverage** | Thousands of organisms (Bacteria, Archaea, Eukaryota) | Strictly *Homo sapiens* (human only) |
| **Confidence Metric** | Combined probabilistic score ($0 - 1000$) | Experimentally weighted score ($0.00 - 1.00$) |
| **Evidence Channels** | 7 channels (experiments, databases, textmining, co-expression, synteny, fusions, phylogeny) | 3 channels (study count, experimental technique scores, orthology) |
| **Functional Enrichment** | Built-in ORA & GSEA (GO, KEGG, Reactome, Pfam, Monarch) | External GO / MeSH annotations |
| **Tissue Specificity** | Semi-quantitative via TISSUES channel | Quantitative GTEx RNA-seq RPKM $\ge 1$ filters |
| **AI Agent Interface** | Native MCP server (`https://mcp.string-db.org`) | REST API + local flat files |
