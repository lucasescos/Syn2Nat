# BioGRID: Biological General Repository for Interaction Datasets

> **Official Web Portal**: [https://thebiogrid.org](https://thebiogrid.org)  
> **REST API Documentation**: [https://wiki.thebiogrid.org/doku.php/biogridrest](https://wiki.thebiogrid.org/doku.php/biogridrest)  
> **API Key Registration**: [https://webservice.thebiogrid.org/](https://webservice.thebiogrid.org/)  
> **Formal WADL Specification**: [https://webservice.thebiogrid.org/application.wadl](https://webservice.thebiogrid.org/application.wadl)  
> **Downloads Portal**: [https://downloads.thebiogrid.org/BioGRID](https://downloads.thebiogrid.org/BioGRID)  
> **Python Client**: [`scripts/biogrid_client.py`](../scripts/biogrid_client.py)  
> **Related Projects**: [BioGRID ORCS](https://orcs.thebiogrid.org) (CRISPR Screens), [HIPPIE Reference (`docs/hippie.md`)](hippie.md)

---

## 1. Overview & System Architecture

**BioGRID** (Biological General Repository for Interaction Datasets) is a biomedical database that provides comprehensive curation of physical and genetic interactions, chemical interactions, post-translational modifications (PTMs), and phenotypic consequences from primary biomedical literature across major model organisms and humans.

Maintained by the **Tyers Lab** (Université de Montréal), the **Mike Tyers / Stark teams**, and collaborating international institutions, BioGRID is updated on the **4th of every month**.

```
Primary Literature & Curated Screens (PubMed, High-Throughput Assays, PDB)
                                 │
                                 ▼
                     [ BioGRID Curation Pipeline ]
         - Experimental System Annotation (PSI-MI mapping)
         - Multi-Taxon Harmonization (Human, Yeast, Mouse, Fly, Worm, etc.)
         - Throughput Stratification (Low vs. High Throughput)
         - Chemical / Drug / CRISPR Phenotype Linkage
                                 │
                                 ▼
                   [ BioGRID Central Database ]
             Interaction Records: Physical, Genetic, Chemical
                                 │
         ┌───────────────────────┴───────────────────────┐
         ▼                                               ▼
[ REST Web Service API ]                        [ Bulk Flat Files ]
  - Endpoint: webservice.thebiogrid.org           - Downloads: downloads.thebiogrid.org
  - Interactive queries via HTTP GET/POST         - TAB2, TAB3, PSI-MI 2.5, BioPAX
  - Real-time subnetwork retrieval                - Whole-interactome offline graphs
  - Formats: JSON, TAB2, Count, Extended
```

### Key Highlights
* **Breadth & Coverage**: Over 2.8 million curated interactions covering >85 model organisms (notably *Homo sapiens*, *Saccharomyces cerevisiae*, *Mus musculus*, *Drosophila melanogaster*, *Caenorhabditis elegans*, *Arabidopsis thaliana*).
* **Multi-Modal Interactions**: Curates direct physical interactions (e.g. two-hybrid, reconstituted complex, co-crystal), indirect associations (affinity capture-MS), genetic interactions (synthetic lethality, phenotypic enhancement/suppression), and chemical–protein associations.
* **Monthly Releases**: Updated on the 4th of every month with new interactions curated directly from PubMed literature.
* **Open REST Web Service**: Provides parameterized HTTP programmatic endpoints for direct integration into bioinformatics pipelines and computational frameworks.

---

## 2. Authentication & Access Keys

To use the BioGRID Web Service, all API requests require a valid, 32-character alphanumeric **Access Key** (`accesskey=[ACCESSKEY]`).

### How to Obtain an Access Key
1. Visit the web portal: [https://webservice.thebiogrid.org/](https://webservice.thebiogrid.org/)
2. Enter your:
   - **First Name**
   - **Last Name**
   - **Email Address**
   - **Project Name** (e.g. `Microprotein Interactome Profiling`)
3. Click **"Generate Access Key"**.
4. Your 32-character key will be displayed immediately.
5. Store the key in your local project environment:
   ```env
   # .env
   BIOGRID_ACCESS_KEY=abcdef0123456789abcdef0123456789
   ```

> [!NOTE]
> Access keys are free and instant for academic, non-profit, and research use. The email is collected only for administrative notices regarding service disruptions or updates.

---

## 3. REST API Architecture & Endpoints

The API is served over HTTPS from the base URL:
```
https://webservice.thebiogrid.org
```

Both `GET` and `POST` requests are supported. For pipeline scripts and reproducible queries, standard HTTP `GET` with encoded query parameters is recommended.

### Formal WADL Description
BioGRID provides a formal **Web Application Description Language (WADL)** specification describing all endpoints and parameters at:
```
https://webservice.thebiogrid.org/application.wadl
```

### Summary of Endpoints

| Endpoint | HTTP Method | Description |
| :--- | :--- | :--- |
| `/interactions/` | `GET`, `POST` | Primary search endpoint: fetches interactions matching gene lists, taxa, evidence types, PubMed IDs, and throughput. |
| `/interaction/{interactionId}` | `GET` | Fetches a single interaction by its unique BioGRID Interaction ID (integer). |
| `/organisms/` | `GET` | Lists all supported NCBI taxonomy IDs and corresponding scientific/common names. |
| `/identifiers/` | `GET` | Lists all supported external identifier types (e.g., UniProt, Ensembl, RefSeq, SGD). |
| `/evidence/` | `GET` | Lists all supported experimental evidence codes (e.g., Two-hybrid, Affinity Capture-MS). |
| `/version/` | `GET` | Returns the active BioGRID database and REST API version string. |

---

## 4. Query Parameters Reference

All query parameters are passed via query string (e.g., `?accesskey=...&geneList=TP53&...`).

### Complete Parameter Specification

| Parameter | Type | Default | Valid Values | Description |
| :--- | :--- | :--- | :--- | :--- |
| `accesskey` (or `accessKey`) | `string` | *None* | 32-char alphanumeric string | **Required**. Your personal API access key. |
| `geneList` | `string` | *Empty* | Pipe-delimited list (`\|`) | Target genes/proteins to search for (e.g. `TP53\|MDM2\|CDKN1A`). |
| `searchNames` | `boolean` | `false` | `true`, `false` | When `true`, matches entries in `geneList` against `OFFICIAL_SYMBOL`. |
| `searchIds` | `boolean` | `false` | `true`, `false` | When `true`, matches against `ENTREZ_GENE`, `ORDERED LOCUS`, and `SYSTEMATIC_NAME`. |
| `searchSynonyms` | `boolean` | `false` | `true`, `false` | When `true`, matches against known gene/protein aliases and synonyms. |
| `searchBiogridIds` | `boolean` | `false` | `true`, `false` | When `true`, matches against internal BioGRID IDs. |
| `additionalIdentifierTypes` | `string` | *Empty* | Pipe-delimited (e.g. `UNIPROT\|ENSEMBL`) | Matches against external database identifiers (e.g. UniProt accessions `P04637`, Ensembl IDs, RefSeq). |
| `excludeGenes` | `boolean` | `false` | `true`, `false` | When `true`, inverts the filter: excludes interactions involving genes in `geneList`. |
| `includeInteractors` | `boolean` | `true` | `true`, `false` | `true`: returns interactions where **at least one** partner is in `geneList` (1st-order neighborhood).<br>`false`: returns interactions strictly **between** genes in `geneList` (internal subnetwork / Layer 0). |
| `includeInteractorInteractions` | `boolean` | `false` | `true`, `false` | When `true` (and `includeInteractors=true`), also returns interactions between the interactors themselves (complete induced 1-hop subgraph). |
| `taxId` | `string` | `"All"` | NCBI Tax ID (e.g. `9606`, `10090`, `559292`) | Restricts interactor search to specific organism(s). Pipe-separate multiple tax IDs. |
| `interSpeciesExcluded` | `boolean` | `false` | `true`, `false` | When `true`, discards cross-species interactions (e.g. Human interactor interacting with a Mouse protein). |
| `selfInteractionsExcluded` | `boolean` | `false` | `true`, `false` | When `true`, discards homomers / self-interactions ($A \leftrightarrow A$). |
| `evidenceList` | `string` | *Empty* | Pipe-delimited list of evidence names | List of experimental systems (e.g. `Two-hybrid\|Affinity Capture-MS`). |
| `includeEvidence` | `boolean` | `false` | `true`, `false` | `false`: `evidenceList` acts as a **blacklist** (excludes listed systems).<br>`true`: `evidenceList` acts as a **whitelist** (returns only listed systems). |
| `pubmedList` | `string` | *Empty* | Pipe-delimited PubMed IDs | Filters by publication PMIDs (e.g. `18316726\|17662948`). |
| `excludePubmeds` | `boolean` | `false` | `true`, `false` | `false`: include only listed PMIDs.<br>`true`: exclude listed PMIDs. |
| `throughputTag` | `string` | `"any"` | `"any"`, `"low"`, `"high"` | Filters by experimental scale (`"low"`: low-throughput; `"high"`: high-throughput screens). |
| `htpThreshold` | `integer` | `2147483647` | Integer $\ge 0$ | Discards interactions from any publication reporting more than `htpThreshold` interactions (curtails mass screen noise). |
| `start` | `integer` | `0` | $0 \le n \le 2147483647$ | Result offset index (0-indexed) for pagination. |
| `max` | `integer` | `10000` | $1 \le n \le 10000$ | Maximum interactions per request (capped at 10,000 by the server). |
| `format` | `string` | `"tab2"` | `"tab1"`, `"tab2"`, `"extendedTab2"`, `"json"`, `"jsonExtended"`, `"count"` | Response payload format. |
| `includeHeader` | `boolean` | `false` | `true`, `false` | In tab-delimited formats, includes the `#` header comment line at line 1. |
| `translate` | `boolean` | `false` | `true`, `false` | When `true`, prepends a diagnostic summary showing how the server interpreted the query parameters. |

---

## 5. Network Topologies & Query Logic

Understanding how BioGRID interprets `includeInteractors` and `includeInteractorInteractions` is crucial for constructing network analyses:

```
Scenario A: includeInteractors = false
Query: [Gene A, Gene B, Gene C]
Result: Only edges strictly BETWEEN the query genes.
(Equivalent to HIPPIE Layer 0)
        [Gene A] ──────── [Gene B]
            \               /
             \             /
              ── [Gene C] ─

──────────────────────────────────────────────────────────────────

Scenario B: includeInteractors = true, includeInteractorInteractions = false
Query: [Gene A]
Result: Star graph of Gene A and all its immediate interaction partners.
(Equivalent to HIPPIE Layer 1)
           [Partner 1]
               │
   [Partner 4]─[Gene A]─[Partner 2]
               │
           [Partner 3]

──────────────────────────────────────────────────────────────────

Scenario C: includeInteractors = true, includeInteractorInteractions = true
Query: [Gene A]
Result: Full induced subnetwork (Star graph + edges connecting partners).
           [Partner 1] ──────── [Partner 2]
               │     \         /     │
               │      [Gene A]       │
               │     /         \     │
           [Partner 4] ──────── [Partner 3]
```

---

## 6. Output Formats & Schemas

BioGRID supports multiple response formats via the `format` parameter:

### 1. `tab2` (BioGRID Tab 2.0 Format) — Default
A standardized, 24-column tab-delimited table. Missing values are represented as hyphens (`-`).

| Col # | Field Name | Description | Example |
| :--- | :--- | :--- | :--- |
| 1 | `BIOGRID_INTERACTION_ID` | Unique interaction record identifier | `616539` |
| 2 | `ENTREZ_GENE_A` | NCBI Entrez Gene ID for Interactor A | `7157` |
| 3 | `ENTREZ_GENE_B` | NCBI Entrez Gene ID for Interactor B | `4193` |
| 4 | `BIOGRID_ID_A` | BioGRID identifier for Interactor A | `112999` |
| 5 | `BIOGRID_ID_B` | BioGRID identifier for Interactor B | `110482` |
| 6 | `SYSTEMATIC_NAME_A` | Systematic or ORF name | `TP53` |
| 7 | `SYSTEMATIC_NAME_B` | Systematic or ORF name | `MDM2` |
| 8 | `OFFICIAL_SYMBOL_A` | Official gene symbol | `TP53` |
| 9 | `OFFICIAL_SYMBOL_B` | Official gene symbol | `MDM2` |
| 10 | `SYNONYMS_A` | Pipe-separated aliases for Interactor A | `P53\|BCC7\|LFS1` |
| 11 | `SYNONYMS_B` | Pipe-separated aliases for Interactor B | `HDMX\|hdm2` |
| 12 | `EXPERIMENTAL_SYSTEM` | Evidence assay classification | `Two-hybrid` |
| 13 | `EXPERIMENTAL_SYSTEM_TYPE` | `physical` or `genetic` | `physical` |
| 14 | `PUBMED_AUTHOR` | First author of discovery publication | `Momand J (1992)` |
| 15 | `PUBMED_ID` | NCBI PubMed ID | `1533783` |
| 16 | `ORGANISM_A_ID` | NCBI Taxonomy ID for Interactor A | `9606` |
| 17 | `ORGANISM_B_ID` | NCBI Taxonomy ID for Interactor B | `9606` |
| 18 | `THROUGHPUT` | `Low Throughput`, `High Throughput`, or both | `Low Throughput` |
| 19 | `QUANTITATIVE_SCORE` | Reported affinity, p-value, or SGA score | `-` or `0.85` |
| 20 | `POST_TRANSLATIONAL_MODIFICATION` | PTM substrate modification (for biochemical assays) | `Ubiquitination` |
| 21 | `PHENOTYPES` | Observed genetic or cellular phenotypes | `lethal` |
| 22 | `QUALIFICATIONS` | Specific assay conditions or plain text notes | `-` |
| 23 | `TAGS` | Curatorial classification tags | `-` |
| 24 | `SOURCE_DATABASE` | Database origin of the record | `BIOGRID` |

### 2. `extendedTab2`
Includes all 24 columns from `tab2` plus 3 additional columns:
- Column 25: `SOURCE_DATABASE_IDENTIFIERS`
- Column 26: `NUMBER_OF_INTERACTIONS_PER_PUBLICATION`
- Column 27: `ADDITIONAL_IDENTIFIERS` (e.g. UniProt accessions, RefSeq IDs)

### 3. `json`
Returns an object whose keys are BioGRID Interaction IDs (as strings), mapping to interaction dictionaries:
```json
{
  "439969": {
    "BIOGRID_INTERACTION_ID": 439969,
    "BIOGRID_ID_A": 33287,
    "BIOGRID_ID_B": 31905,
    "ENTREZ_GENE_A": "852931",
    "ENTREZ_GENE_B": "851396",
    "OFFICIAL_SYMBOL_A": "KSS1",
    "OFFICIAL_SYMBOL_B": "STE7",
    "SYSTEMATIC_NAME_A": "YGR040W",
    "SYSTEMATIC_NAME_B": "YDL159W",
    "SYNONYMS_A": "FUS4",
    "SYNONYMS_B": "STE7",
    "EXPERIMENTAL_SYSTEM": "Affinity Capture-MS",
    "EXPERIMENTAL_SYSTEM_TYPE": "physical",
    "PUBMED_AUTHOR": "Breitkreutz A (2010)",
    "PUBMED_ID": "20489953",
    "ORGANISM_A": 559292,
    "ORGANISM_B": 559292,
    "THROUGHPUT": "High Throughput",
    "QUANTITATIVE_SCORE": "-",
    "MODIFICATION": "-",
    "PHENOTYPES": "-",
    "QUALIFICATIONS": "-",
    "TAGS": "-",
    "SOURCE_DATABASE": "BIOGRID"
  }
}
```

### 4. `count`
Returns a plain-text integer representing the total number of interactions matching the specified parameters.
* **Best Practice**: Use `format=count` before initiating bulk queries to verify result volume and plan pagination loops.

---

## 7. Experimental Evidence Codes

BioGRID curates evidence across physical and genetic categories:

### Physical Evidence Types
* `Affinity Capture-MS`: Mass spectrometry identification after epitope or affinity purification.
* `Affinity Capture-RNA`: Protein-RNA complex association.
* `Affinity Capture-Western`: Western blot detection after co-immunoprecipitation or pull-down.
* `Biochemical Activity`: In vitro enzymatic assays (kinase, phosphatase, ligase, etc.).
* `Co-crystal Structure`: High-resolution direct atomic interface from X-ray / Cryo-EM / NMR (PDB).
* `Co-fractionation`: Chromatography or gradient sedimentation profiling.
* `Co-localization`: Fluorescent or confocal microscopy co-occurrence.
* `Co-purification`: Native complex purification.
* `Far Western`: Overlay assay using labeled protein probes.
* `FRET`: Fluorescence Resonance Energy Transfer in living cells.
* `PCA`: Protein Fragment Complementation Assay.
* `Protein-peptide`: Interaction between a full-length protein and synthetic peptide.
* `Protein-RNA`: Direct protein-RNA binding.
* `Proximity Label-MS`: BioID / APEX proximity labeling followed by mass spec.
* `Reconstituted Complex`: In vitro binding of purified, individually expressed recombinant proteins.
* `Surface Plasmon Resonance`: Real-time biophysical binding kinetics (Biacore).
* `Two-hybrid`: Yeast or mammalian two-hybrid (Y2H / M2H) reporter assay.

### Genetic Evidence Types
* `Synthetic Lethality`: Double mutant exhibits lethality absent in single mutants.
* `Synthetic Growth Defect`: Double mutant shows impaired fitness.
* `Synthetic Rescue`: Second mutation or overexpression rescues mutant phenotype.
* `Phenotypic Enhancement`: Secondary mutation worsens specific mutant defect.
* `Phenotypic Suppression`: Secondary mutation mitigates mutant phenotype.
* `Dosage Growth Defect`: Overexpression of one gene causes fitness defect in specific genetic background.
* `Dosage Lethality`: Overexpression causes lethality.
* `Dosage Rescue`: Overexpression restores wild-type phenotype.
* `Positive Genetic` / `Negative Genetic`: Quantitative genetic interaction scores (e.g. SGA/E-MAP).

To dynamically fetch the live list of evidence codes:
```
GET https://webservice.thebiogrid.org/evidence/?accesskey=[KEY]&format=json
```

---

## 8. Supported Organisms & Taxon IDs

Common taxonomy IDs supported in BioGRID queries (`taxId`):

| Organism | Common Name | NCBI TaxID |
| :--- | :--- | :--- |
| *Homo sapiens* | Human | `9606` |
| *Mus musculus* | House mouse | `10090` |
| *Rattus norvegicus* | Norway rat | `10116` |
| *Drosophila melanogaster* | Fruit fly | `7227` |
| *Caenorhabditis elegans* | Nematode | `6239` |
| *Saccharomyces cerevisiae* | Baker's yeast | `559292` (or `4932`) |
| *Schizosaccharomyces pombe* | Fission yeast | `4896` |
| *Danio rerio* | Zebrafish | `7955` |
| *Arabidopsis thaliana* | Thale cress | `3702` |
| *Severe acute respiratory syndrome coronavirus 2* | SARS-CoV-2 | `2697049` |

To retrieve the complete, dynamic list of all supported organisms:
```
GET https://webservice.thebiogrid.org/organisms/?accesskey=[KEY]&format=json
```

---

## 9. Programmatic Usage Guide

### Method 1: Using the Project Python Client (`scripts/biogrid_client.py`)

A full-featured, zero-dependency Python client is included in [`scripts/biogrid_client.py`](../scripts/biogrid_client.py).

#### Python API Examples

```python
from scripts.biogrid_client import BioGridClient

# 1. Initialize client (automatically reads BIOGRID_ACCESS_KEY from .env)
client = BioGridClient()

# 2. Get interaction count for a gene (pre-flight check)
count = client.count(genes="TP53", tax_id=9606)
print(f"Total interactions for TP53: {count}")

# 3. Retrieve 1st-order physical interactors (low-throughput / high-confidence)
df_tp53 = client.query(
    genes=["TP53"],
    tax_id=9606,
    throughput="low",
    exp_type="physical",
    as_dataframe=True
)
print(df_tp53[["OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B", "EXPERIMENTAL_SYSTEM", "PUBMED_AUTHOR"]].head())

# 4. Check if two specific proteins physically interact
result = client.query_pair("TP53", "MDM2", exp_type="physical")
if result:
    print(f"Found {len(result)} reported physical assays between TP53 and MDM2")
    for hit in result:
        print(f"  - {hit['EXPERIMENTAL_SYSTEM']} ({hit['PUBMED_AUTHOR']}, PMID: {hit['PUBMED_ID']})")

# 5. Extract an internal subnetwork (edges strictly between query set)
candidates = ["TP53", "MDM2", "CDKN1A", "BAX", "ATM", "PML"]
subnetwork_df = client.query_subnetwork(candidates, tax_id=9606, as_dataframe=True)
print(f"Subnetwork contains {len(subnetwork_df)} internal edges")

# 6. Query by UniProt Accession
uniprot_df = client.query_by_identifier(
    identifiers=["P04637"],  # TP53 UniProt accession
    id_type="UNIPROT",
    tax_id=9606,
    as_dataframe=True
)
```

#### CLI Terminal Usage

```powershell
# Check total interaction count for human MDM2
python scripts/biogrid_client.py --count --genes MDM2 --taxid 9606

# Query low-throughput physical interactors for TP53
python scripts/biogrid_client.py --genes TP53 --taxid 9606 --throughput low --exp-type physical

# Check pairwise interaction between two proteins
python scripts/biogrid_client.py --pair TP53 MDM2 --taxid 9606

# Extract internal subnetwork (Layer 0) between candidate proteins
python scripts/biogrid_client.py --genes TP53,MDM2,CDKN1A,BAX,ATM --subnetwork --output ppi_subnetwork.tsv

# Export single interaction details by BioGRID Interaction ID
python scripts/biogrid_client.py --interaction-id 439969 --format json
```

---

### Method 2: Direct REST Requests with Python `requests`

```python
import os
import requests
import pandas as pd

API_KEY = os.getenv("BIOGRID_ACCESS_KEY", "YOUR_32_CHAR_ACCESS_KEY")
URL = "https://webservice.thebiogrid.org/interactions/"

params = {
    "accesskey": API_KEY,
    "format": "json",
    "geneList": "TP53|MDM2",
    "searchNames": "true",
    "includeInteractors": "false",  # Strictly interactions BETWEEN TP53 and MDM2
    "taxId": 9606,
    "interSpeciesExcluded": "true",
    "selfInteractionsExcluded": "true"
}

response = requests.get(URL, params=params)
response.raise_for_status()
data = response.json()

# Convert to DataFrame
df = pd.DataFrame.from_dict(data, orient="index")
print(f"Retrieved {len(df)} interactions:")
print(df[["OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B", "EXPERIMENTAL_SYSTEM", "PUBMED_ID"]])
```

---

### Method 3: Handling Large Queries (> 10,000 interactions) via Pagination

BioGRID limits single responses to **10,000 records** (`max=10000`). For whole-interactome screens or massive hubs, paginate using `start`:

```python
import os
import requests
import pandas as pd

API_KEY = os.getenv("BIOGRID_ACCESS_KEY")
URL = "https://webservice.thebiogrid.org/interactions/"

# Step 1: Pre-flight count
count_res = requests.get(URL, params={
    "accesskey": API_KEY,
    "geneList": "TP53",
    "searchNames": "true",
    "taxId": 9606,
    "format": "count"
})
total_interactions = int(count_res.text.strip())
print(f"Total available interactions: {total_interactions}")

# Step 2: Paginate in chunks of 10,000
all_records = {}
page_size = 10000
start = 0

while start < total_interactions:
    params = {
        "accesskey": API_KEY,
        "geneList": "TP53",
        "searchNames": "true",
        "taxId": 9606,
        "format": "json",
        "start": start,
        "max": page_size
    }
    r = requests.get(URL, params=params)
    chunk = r.json()
    if not chunk:
        break
    all_records.update(chunk)
    start += page_size
    print(f"Fetched {len(all_records)} / {total_interactions}...")

print(f"Successfully retrieved all {len(all_records)} interactions.")
```

---

### Method 4: Command-Line (cURL + jq)

```bash
# Set your API key
export BIOGRID_KEY="your_access_key_here"

# 1. Fetch JSON interactions for TP53
curl -s "https://webservice.thebiogrid.org/interactions/?accesskey=${BIOGRID_KEY}&geneList=TP53&searchNames=true&taxId=9606&max=10&format=json" \
  | jq 'to_entries[] | .value | {id: .BIOGRID_INTERACTION_ID, interactorA: .OFFICIAL_SYMBOL_A, interactorB: .OFFICIAL_SYMBOL_B, system: .EXPERIMENTAL_SYSTEM}'

# 2. Get interaction count directly
curl -s "https://webservice.thebiogrid.org/interactions/?accesskey=${BIOGRID_KEY}&geneList=TP53&searchNames=true&taxId=9606&format=count"

# 3. Check version
curl -s "https://webservice.thebiogrid.org/version/?accesskey=${BIOGRID_KEY}"
```

---

### Method 5: R Integration (`httr2` & `jsonlite`)

```r
library(httr2)
library(jsonlite)
library(dplyr)

api_key <- Sys.getenv("BIOGRID_ACCESS_KEY")

req <- request("https://webservice.thebiogrid.org/interactions/") %>%
  req_url_query(
    accesskey = api_key,
    geneList = "TP53|MDM2",
    searchNames = "true",
    taxId = "9606",
    format = "json"
  )

resp <- req_perform(req)
json_data <- resp_body_json(resp)

# Bind list of interactions into a tibble
ppi_df <- bind_rows(json_data)
print(select(ppi_df, OFFICIAL_SYMBOL_A, OFFICIAL_SYMBOL_B, EXPERIMENTAL_SYSTEM, PUBMED_ID))
```

---

## 10. Bulk Downloads vs. REST API

For tasks requiring the complete interactome (e.g., training graph neural networks, computing global topological centrality, or running whole-genome latent space embeddings), downloading flat files is much faster and eliminates network overhead:

| Feature | REST Web Service | Bulk Flat Files |
| :--- | :--- | :--- |
| **Best For** | Targeted gene lists, candidate validation, live pipeline queries, web servers | Interactome-wide training, graph deep learning, offline clusters |
| **Data Freshness** | Live release (updated monthly) | Synced monthly with release tarballs |
| **Latency** | Network latency per request | Zero latency once downloaded locally |
| **Record Limit** | 10,000 records per HTTP call | Entire multi-million edge dataset |
| **Standard File** | Returns on-the-fly TSV or JSON | `BIOGRID-ALL-latest.tab3.zip` (~150 MB compressed) |
| **Download URL** | N/A | `https://downloads.thebiogrid.org/BioGRID/Release-Archive/` |

---

## 11. Integration with Microprotein Research & PPI Hierarchy

In the context of the **Human Microprotein Characterization and Interactome Structural Discovery** pipeline:

1. **Orthogonal Validation**:
   - High-throughput mass-spec or two-hybrid interactions in BioGRID provide empirical confirmation for predicted microprotein–target complexes.
   - Cross-referencing against BioGRID physical assays distinguishes direct stoichiometric partners from non-specific co-complexes.

2. **Integration with HIPPIE & STRING**:
   - BioGRID is one of the primary data sources aggregated by **HIPPIE** ([`docs/hippie.md`](hippie.md)).
   - While HIPPIE pre-computes confidence scores ($0.0 \le S \le 1.0$) based on study count and technique weights, direct BioGRID API access enables filtering down to specific individual assay conditions (e.g., `Co-crystal Structure`, `Reconstituted Complex`), exact PubMed citations, or post-translational modification statuses (`MODIFICATION`).

3. **Induced Subnetwork Mapping**:
   - Setting `includeInteractors=false` retrieves the direct interaction web connecting candidate physiological partners of microprotein targets, facilitating downstream pathway enrichment and complex assembly modeling.
