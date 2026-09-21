# IntAct Molecular Interaction Database: API Reference & Programmatic Access Guide

> **Official Portal**: [https://www.ebi.ac.uk/intact/](https://www.ebi.ac.uk/intact/)  
> **API Base URL**: [https://www.ebi.ac.uk/intact/ws/](https://www.ebi.ac.uk/intact/ws/)  
> **Documentation & Technical Corner**: [https://www.ebi.ac.uk/intact/resources](https://www.ebi.ac.uk/intact/resources)  
> **GitHub Organization**: [https://github.com/intact-portal](https://github.com/intact-portal)  
> **Python Client**: [`scripts/intact_client.py`](../scripts/intact_client.py)  
> **Related PPI Resources**: [HIPPIE (`docs/hippie.md`)](hippie.md) | [MAPPIE (`docs/mappie.md`)](mappie.md)

---

## 1. Overview & System Architecture

**IntAct** is an open-source, open-access database system and analysis tool for molecular interaction data, developed and maintained by the Molecular Interactions Team at the **European Bioinformatics Institute (EMBL-EBI)**. IntAct is a founding member of the **International Molecular Exchange (IMEx) Consortium**.

### Key Facts & Metrics
* **Scale**: Over **1,780,000+ binary interactions** involving **145,000+ interactors** across thousands of organisms.
* **Curation Standard**: Full manual curation adhering to **IMEx standards** directly from peer-reviewed literature, as well as curated data exchanges with MINT, DIP, MatrixDB, and UniProtKB.
* **Data Conformity**: Completely standard-compliant with the **Proteomics Standards Initiative Molecular Interactions (PSI-MI)** community formats (PSI-MI XML 2.5/3.0, PSI-MITAB 2.5/2.6/2.7/2.8, and MI-JSON).
* **Confidence Scoring**: Each interaction is scored with **MIscore** ($0.00$ to $1.00$), an empirical confidence metric incorporating publication count, experimental detection methods, and interaction types.

### Architectural Stack
IntAct is architected as a set of decoupled backend microservices:

```
                      [ User / Client Applications ]
                                    │
           ┌────────────────────────┼────────────────────────┐
           ▼                        ▼                        ▼
      cURL / HTTP              Python Client            Cytoscape / Web
           │                        │                        │
           └────────────────────────┼────────────────────────┘
                                    │
                                    ▼
                 [ IntAct REST Gateway: /intact/ws/ ]
                                    │
       ┌───────────────┬────────────┴───┬────────────────┐
       ▼               ▼                ▼                ▼
 /ws/interactor  /ws/interaction   /ws/network       /ws/graph
  (Resolving &    (Binary Search &   (Cytoscape.js   (Neo4j Graph,
   Suggestions)    Faceting / Solr)   JSON format)    Details & Exports)
       │               │                                 │
       ▼               ▼                                 ▼
 [ Apache Solr Search Index ]               [ Neo4j Graph Database ]
 (Fast MIQL text & facet search)            (Deep node/edge traversal)
```

---

## 2. IntAct Microservice REST Endpoints

The modern IntAct platform exposes **four specialized REST microservices** under `https://www.ebi.ac.uk/intact/ws/`. Each service publishes an interactive **Swagger UI** and **OpenAPI 3.0** specification.

| Microservice | Base Path | Swagger UI | OpenAPI 3.0 Spec | Core Purpose |
| :--- | :--- | :--- | :--- | :--- |
| **Interactor** | `/ws/interactor` | [`/ws/interactor/swagger-ui.html`](https://www.ebi.ac.uk/intact/ws/interactor/swagger-ui.html) | [`/v3/api-docs`](https://www.ebi.ac.uk/intact/ws/interactor/v3/api-docs) | Interactor lookup, autocomplete, batch identifier resolution |
| **Interaction** | `/ws/interaction` | [`/ws/interaction/swagger-ui.html`](https://www.ebi.ac.uk/intact/ws/interaction/swagger-ui.html) | [`/v3/api-docs`](https://www.ebi.ac.uk/intact/ws/interaction/v3/api-docs) | Query binary interactions, faceted filtering, MIQL queries |
| **Network** | `/ws/network` | [`/ws/network/swagger-ui.html`](https://www.ebi.ac.uk/intact/ws/network/swagger-ui.html) | [`/v3/api-docs`](https://www.ebi.ac.uk/intact/ws/network/v3/api-docs) | Generate network graph JSON compatible with Cytoscape.js |
| **Graph** | `/ws/graph` | [`/ws/graph/swagger-ui.html`](https://www.ebi.ac.uk/intact/ws/graph/swagger-ui.html) | [`/v3/api-docs`](https://www.ebi.ac.uk/intact/ws/graph/v3/api-docs) | Deep entity details (participants, features, mutations) & multi-format exports |

---

### 2.1. Interactor Microservice (`/ws/interactor`)

Used to search, validate, and resolve molecular interactors (proteins, small molecules, nucleic acids, peptides).

#### Endpoints
* **`GET /findInteractor/{query}`**: Searches interactors matching a keyword, gene symbol, or accession.
  * `query`: Search string (e.g., `BRCA1`, `P38398`, `tumor suppressor`).
  * `page`: Zero-based page index (default: `0`).
  * `pageSize`: Number of items per page (default: `10`).
* **`POST /list/resolve`**: Batch resolves a list of terms/identifiers to canonical IntAct interactors.
  * Query parameters / payload:
    * `query`: Comma or space-separated list of identifiers.
    * `fuzzySearch`: Boolean (`true` / `false`).
    * `page`: `0`
    * `pageSize`: `10`
  * Response: A dictionary mapping `{ query_term: { content: [...] } }`.
* **`GET /countTotal`**: Returns the total number of interactors in IntAct (e.g., ~145,800).

#### Response Schema Example (`GET /findInteractor/BRCA1`)
```json
{
  "totalElements": 61,
  "totalPages": 7,
  "content": [
    {
      "interactorAc": "EBI-349905",
      "interactorName": "BRCA1",
      "interactorPreferredIdentifier": "P38398",
      "interactorDescription": "Breast cancer type 1 susceptibility protein",
      "interactorSpecies": "Homo sapiens",
      "interactorTaxId": 9606,
      "interactorType": "protein",
      "interactorTypeMiIdentifier": "MI:0326"
    }
  ]
}
```

---

### 2.2. Interaction Microservice (`/ws/interaction`)

The primary search engine for molecular interactions. Built on Apache Solr, it supports simple keyword queries, complex multi-field **MIQL** queries, and extensive faceting.

#### Endpoints
* **`GET /findInteractions/{query}`**: Simple GET search by token or accession.
  * Parameters: `query`, `page`, `pageSize`.
* **`POST /list`** or **`POST /list/body`**: Advanced filtered interaction search.
  * Accepts a JSON configuration object or query parameters:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `query` | `string` | `*` | Search query or MIQL expression |
| `advancedSearch` | `boolean` | `false` | Must be `true` when submitting MIQL syntax |
| `batchSearch` | `boolean` | `false` | Set `true` if querying a list of interactors |
| `minMIScore` | `float` | `0.0` | Lower bound of confidence score $[0.0, 1.0]$ |
| `maxMIScore` | `float` | `1.0` | Upper bound of confidence score $[0.0, 1.0]$ |
| `expansionFilter` | `boolean` | `false` | If `true`, excludes spoke-expanded co-complexes (returns **only** true binary experiments) |
| `negativeFilter` | `string` | `POSITIVE_ONLY` | Options: `POSITIVE_ONLY`, `NEGATIVE_ONLY`, `POSITIVE_AND_NEGATIVE` |
| `mutationFilter` | `boolean` | `false` | If `true`, returns only interactions affected by mutations |
| `intraSpeciesFilter`| `boolean` | `false` | If `true`, restricts to interactions within the same species |
| `interactorSpeciesFilter` | `list[string]` | `[]` | Filter by species (e.g., `["Homo sapiens"]` or `["9606"]`) |
| `interactionTypesFilter` | `list[string]` | `[]` | Filter by interaction type (e.g., `["direct interaction"]`) |
| `interactionDetectionMethodsFilter`| `list[string]` | `[]`| Filter by detection method (e.g., `["two hybrid"]`, `["pull down"]`) |
| `page` | `int` | `0` | Page index (0-based) |
| `pageSize` | `int` | `10` | Records per page (up to 100) |

* **`POST /findInteractionFacets`**: Computes facet statistics for a given query (taxa distributions, detection methods, MIscore breakdown, mutation counts).
* **`GET /countTotal`**: Returns total interactions count (over 1,780,000+).
* **`POST /countInteractionResult/{interactorAc}`**: Returns total interactions involving a specific accession.

#### Interaction Record Structure
Each interaction object inside `content` delivers rich metadata, including pre-formatted PSI-MI interchange representations:

```json
{
  "ac": "EBI-21978305",
  "binaryInteractionId": 5123594,
  "idA": "P38398 (uniprotkb)",
  "idB": "EBI-21978329 (intact)",
  "moleculeA": "BRCA1",
  "moleculeB": "ydisqvfpf",
  "typeA": "protein",
  "typeB": "peptide",
  "taxIdA": 9606,
  "speciesA": "Homo sapiens",
  "detectionMethod": "pull down",
  "detectionMethodMIIdentifier": "MI:0096",
  "type": "physical association",
  "typeMIIdentifier": "MI:0915",
  "intactMiscore": 0.4,
  "negative": false,
  "affectedByMutation": true,
  "expansionMethod": null,
  "publicationPubmedIdentifier": "15133502",
  "tab27Format": "uniprotkb:P38398\tintact:EBI-21978329\t...",
  "jsonFormat": "{\"object\":\"interactor\",\"id\":\"intact_EBI-21978329\",...}"
}
```

---

### 2.3. Network Microservice (`/ws/network`)

Generates interaction subgraphs optimized for graphical network visualizers (e.g., Cytoscape.js, ComplexViewer).

#### Endpoint
* **`POST /getInteractions`** (or **`POST /getInteractions/body`**)
  * Query parameters identical to `POST /interaction/list`.
  * Response contains two keys:
    * `data`: Array of Cytoscape.js nodes and edges.
    * `legend`: Color palettes, shapes, and taxonomy visual mappings.

#### Cytoscape Node / Edge Example
```json
{
  "data": [
    {
      "data": {
        "id": "EBI-349905",
        "label": "BRCA1 (P38398)",
        "interactor_name": "BRCA1",
        "preferred_id": "P38398",
        "species": "Homo sapiens",
        "tax_id": 9606,
        "interactor_type": "protein",
        "shape": "ellipse",
        "color": "#335e94",
        "mutation": true
      }
    },
    {
      "data": {
        "id": "edge_5123594",
        "source": "EBI-349905",
        "target": "EBI-21978329",
        "interaction_type": "physical association",
        "mi_score": 0.4
      }
    }
  ]
}
```

---

### 2.4. Graph Microservice (`/ws/graph`)

Interfaces with the IntAct Neo4j graph database to retrieve fine-grained biological annotations (molecular features, binding regions, residue mutations, kinetic parameters) and export data in standard formats.

#### Endpoints
* **`GET /interaction/details/{ac}`**: Complete interaction curation record (authors, publication, detection method, host organism, confidences, parameters).
* **`GET /participants/details/{ac}`**: Participant entities in interaction `ac` (biological and experimental roles, stoichiometry, expressed in host).
* **`GET /features/details/{ac}`**: Detailed sequence features on participants (binding domains, tags, mutations, PTMs, residue ranges).
* **`GET /export/interaction/{ac}?format={format}`**: Export single interaction.
* **`POST /export/interaction/list?query={query}&format={format}`**: Batch export interactions matching a query.

#### Supported Export Formats
| Format Identifier | MIME Type | Description |
| :--- | :--- | :--- |
| `miTab25` | `text/plain` | PSI-MITAB 2.5 (15 standard tab-delimited columns) |
| `miTab26` | `text/plain` | PSI-MITAB 2.6 (adds biological/experimental roles, participant identification) |
| `miTab27` | `text/plain` | PSI-MITAB 2.7 (42 columns, adds features, stoichiometry, participant types) |
| `featureTab` | `text/plain` | PSI-MI FeatureTab (detailed residue ranges and mutations) |
| `miJSON` | `application/json`| Compact interaction JSON used by ComplexViewer |
| `miXML25` | `application/xml` | Standard PSI-MI XML v2.5 |
| `miXML30` | `application/xml` | Standard PSI-MI XML v3.0 |

---

## 3. Molecular Interaction Query Language (MIQL)

IntAct implements **MIQL 2.8**, an extension of Apache Lucene query syntax designed for molecular interaction data. When calling `/ws/interaction/list` or `/findInteractions`, set `advancedSearch=true` to execute MIQL queries.

### 3.1. Syntax Rules
* **Operators**: `AND`, `OR`, `NOT` (must be capitalized).
* **Grouping**: Use parentheses `( ... )`.
* **Phrases**: Wrap strings containing spaces or hyphens in double quotes: `type:"direct interaction"`.
* **Ranges**: Use brackets `[lower TO upper]`: `intact-miscore:[0.60 TO 1.00]`.
* **Wildcards**: `*` (multiple characters) and `?` (single character) allowed at middle or end of words (not at beginning).
* **Null Check**: Use `-` as value to test for null or empty values.

### 3.2. MIQL Field Reference

| MIQL Field | Search Target | Example |
| :--- | :--- | :--- |
| `id` | Any identifier of Interactor A or B | `id:P38398` |
| `idA`, `idB` | Identifier of Interactor A or B specifically | `idA:P38398 AND idB:P51587` |
| `alias`, `identifier` | UniProt IDs, gene symbols, synonyms | `alias:(BRCA1 OR BRCA2)` |
| `geneName` | Gene name of Interactor A or B | `geneName:TP53` |
| `geneNameA`, `geneNameB` | Specific gene name for A or B | `geneNameA:MDM2 AND geneNameB:TP53` |
| `species` | Species name or NCBI Taxonomy ID | `species:human` or `species:9606` |
| `taxidA`, `taxidB` | Taxonomy ID of Interactor A or B | `taxidA:9606 AND taxidB:9606` |
| `taxidHost` | Host experimental expression organism | `taxidHost:9606` |
| `type` | PSI-MI interaction type | `type:"direct interaction"` |
| `detmethod` | Interaction detection method | `detmethod:"two hybrid"` |
| `intact-miscore` | MIscore confidence range | `intact-miscore:[0.70 TO 1.00]` |
| `complex` | Expansion algorithm used | `complex:spoke` or `-complex:*` (non-expanded) |
| `negative` | Negative interaction flag | `negative:false` |
| `mutation` | Interaction affected by mutation | `mutation:true` |
| `pubid` | PubMed ID or publication accession | `pubid:15133502` |
| `pubauth`, `pubyear` | First author or publication year range | `pubyear:[2020 TO 2026]` |
| `ftype` | Feature type (e.g., binding region, tag) | `ftype:"binding region"` |

### 3.3. Practical MIQL Query Examples
* **High-confidence direct human interactions for TP53**:
  ```lucene
  species:9606 AND geneName:TP53 AND type:"direct interaction" AND intact-miscore:[0.60 TO 1.00] AND negative:false
  ```
* **Binary physical interactions verified without complex expansion**:
  ```lucene
  id:P04637 AND -complex:* AND negative:false
  ```
* **Interactions disrupted or affected by mutations**:
  ```lucene
  geneName:BRCA1 AND mutation:true
  ```

---

## 4. Confidence Scoring: The MIscore Model

IntAct interactions are annotated with **MIscore** ($S_{\text{MI}}$), a heuristic scoring methodology developed by the HUPO-PSI consortium.

### Mathematical Formulation
MIscore normalizes multiple lines of experimental evidence into a single confidence value $S_{\text{MI}} \in [0.00, 1.00]$:

$$S_{\text{MI}} = \frac{K_p \cdot S_p(n) + K_m \cdot S_m(\text{cv}) + K_t \cdot S_t(\text{cv})}{K_p + K_m + K_t}$$

Where:
* $K_p, K_m, K_t \in [0, 1]$: Weight factors for publications, detection methods, and interaction types.
* $S_p, S_m, S_t \in [0, 1]$: Normalized sub-scores.

#### 1. Publication Score ($S_p$)
Logarithmic saturation based on independent publications ($n$):
$$S_p = \log_{(b + 1)}(n + 1)$$
*(Default saturation parameter $b = 7$, reaching $1.0$ at 7+ independent studies).*

#### 2. Detection Method Score ($S_m$) & Type Score ($S_t$)
Calculated from the cumulative diversity and quality weights of the PSI-MI Controlled Vocabulary terms:
$$S = \log_{(b + 1)}(a + 1) \quad \text{where} \quad a = \sum (scv_i \times n_i)$$

**Method Quality Weights ($scv_i$)**:
* Biophysical methods (`MI:0013`, X-ray, NMR, ITC, SPR): **$1.00$**
* Biochemical methods (`MI:0401`, pull-down, co-IP): **$1.00$**
* Protein complementation assays (`MI:0090`, Y2H, PCA): **$0.66$**
* Imaging techniques (`MI:0428`, FRET, microscopy): **$0.33$**
* Genetic interference (`MI:0254`): **$0.10$**

**Interaction Type Weights ($scv_i$)**:
* Direct interaction (`MI:0407`): **$1.00$**
* Physical association (`MI:0915`): **$0.66$**
* Association (`MI:0914`): **$0.33$**
* Colocalization (`MI:0403`): **$0.33$**
* Genetic interaction (`MI:0208`): **$0.10$**

#### UniProtKB / GOA Export Threshold
For export to UniProtKB and Gene Ontology annotation (GOA), IntAct applies an even more stringent integer scoring metric requiring:
1. Cumulative score $\ge 9$
2. At least two pieces of experimental evidence
3. At least one non-expanded binary physical evidence.

---

## 5. Complex Expansion: Spoke vs. Matrix

High-throughput affinity purification methods (e.g., AP-MS, TAP) isolate multi-protein co-complexes involving $N > 2$ molecules. To represent these in binary network structures, an expansion algorithm must be applied.

```
       [ Original Co-Complex: 1 Bait + 3 Preys ]
                       (Bait)
                      /   |   \
                   (P1)  (P2)  (P3)

Spoke Expansion (N-1 = 3 edges)       Matrix Expansion (N*(N-1)/2 = 6 edges)
          (Bait)                                  (Bait)
         /   |   \                               /  |   \
      (P1)  (P2)  (P3)                        (P1)──┼──(P2)
                                                \   |   /
(Links only Bait -> each Prey;                   \ (P3)/
 avoids fictitious Prey-Prey edges)     (Links all-to-all; generates
                                        excessive false positives)
```

* **IntAct Standard**: IntAct **strictly uses Spoke expansion** for $N$-ary complexes.
* **Filtering Warning**: A spoke-expanded binary indicates co-occurrence in a complex, but does **not** prove direct physical contact. To filter for direct physical contacts, set:
  * In API: `expansionFilter=true`
  * In MIQL: `-complex:*`

---

## 6. Programmatic Usage in Python, R & cURL

### 6.1. Python (Complete REST Workflow)

```python
import requests

BASE_URL = "https://www.ebi.ac.uk/intact/ws"

def get_high_confidence_interactions(gene_symbol: str, min_score: float = 0.6):
    """Fetch high-confidence binary interactions from IntAct."""
    endpoint = f"{BASE_URL}/interaction/list"
    params = {
        "query": f"geneName:{gene_symbol}",
        "advancedSearch": True,
        "minMIScore": min_score,
        "expansionFilter": True,  # Exclude spoke expansions
        "negativeFilter": "POSITIVE_ONLY",
        "page": 0,
        "pageSize": 50
    }
    
    response = requests.post(endpoint, params=params)
    response.raise_for_status()
    data = response.json()
    
    interactions = []
    for item in data.get("content", []):
        interactions.append({
            "interaction_ac": item.get("ac"),
            "interactor_a": item.get("moleculeA"),
            "interactor_b": item.get("moleculeB"),
            "type": item.get("type"),
            "detection_method": item.get("detectionMethod"),
            "miscore": item.get("intactMiscore"),
            "pubmed_id": item.get("publicationPubmedIdentifier")
        })
    return interactions

if __name__ == "__main__":
    results = get_high_confidence_interactions("TP53", min_score=0.6)
    print(f"Retrieved {len(results)} high-confidence direct interactions for TP53:")
    for r in results[:5]:
        print(f"  {r['interactor_a']} <-> {r['interactor_b']} | Score: {r['miscore']} | Method: {r['detection_method']}")
```

### 6.2. Exporting PSI-MITAB 2.7 via Python

```python
import requests

def download_mitab(gene_name: str, output_file: str):
    """Download interaction records in standard PSI-MITAB 2.7 format."""
    url = "https://www.ebi.ac.uk/intact/ws/graph/export/interaction/list"
    params = {
        "query": f"geneName:{gene_name} AND species:human",
        "format": "miTab27"
    }
    response = requests.post(url, params=params, stream=True)
    response.raise_for_status()
    with open(output_file, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"Saved MITAB 2.7 to {output_file}")
```

### 6.3. R Example (`httr` & `jsonlite`)

```R
library(httr)
library(jsonlite)

query_intact <- function(gene = "BRCA1", min_score = 0.5) {
  url <- "https://www.ebi.ac.uk/intact/ws/interaction/list"
  res <- POST(url, query = list(
    query = paste0("geneName:", gene, " AND species:human"),
    advancedSearch = "true",
    minMIScore = min_score,
    expansionFilter = "true",
    pageSize = 25
  ))
  stop_for_status(res)
  data <- fromJSON(content(res, "text", encoding = "UTF-8"))
  return(data$content[, c("moleculeA", "moleculeB", "intactMiscore", "type", "detectionMethod")])
}

interactions <- query_intact("BRCA1", min_score = 0.5)
head(interactions)
```

### 6.4. cURL & Bash Examples

```bash
# 1. Resolve gene identifier
curl -s "https://www.ebi.ac.uk/intact/ws/interactor/findInteractor/BRCA1?page=0&pageSize=1" | jq .

# 2. Query high-confidence interactions (MIscore >= 0.70)
curl -s -X POST "https://www.ebi.ac.uk/intact/ws/interaction/list?query=geneName:TP53&advancedSearch=true&minMIScore=0.7&pageSize=5" | jq '.content[] | {moleculeA, moleculeB, intactMiscore, detectionMethod}'

# 3. Export interaction EBI-21978305 in PSI-MITAB 2.7 format
curl -s "https://www.ebi.ac.uk/intact/ws/graph/export/interaction/EBI-21978305?format=miTab27"

# 4. Fetch Cytoscape network elements
curl -s -X POST "https://www.ebi.ac.uk/intact/ws/network/getInteractions?query=P38398&pageSize=10" | jq .data[0]
```

---

## 7. Status of PSICQUIC vs. Modern REST Microservices

Users familiar with older proteomics literature may look for **PSICQUIC** (Proteomics Standard Initiative Common QUery InterfaCe). 

| Feature | Legacy PSICQUIC | Modern IntAct REST Microservices (`/ws/`) |
| :--- | :--- | :--- |
| **Endpoint** | `ebi.ac.uk/Tools/webservices/psicquic/...` | `https://www.ebi.ac.uk/intact/ws/` |
| **Current Status** | Deprecated / Service outages common | **Actively maintained & high performance** |
| **Backend** | Legacy index | Apache Solr + Neo4j Graph DB |
| **Facets & Stats** | No | Full faceted search & aggregation |
| **Network Output**| MITAB / XML only | Native Cytoscape.js JSON, ComplexViewer JSON |
| **Molecular Features**| Partial | Full residue ranges, mutations, PTMs |
| **Recommendation**| Use only for federated queries across other databases | **Recommended for all IntAct programmatic queries** |

---

## 8. Summary Checklist for Integration

When integrating IntAct data into automated pipelines or protein representation frameworks:
1. **Always use `/ws/interaction/list` or `/ws/graph/`** instead of legacy PSICQUIC endpoints.
2. **Apply `expansionFilter=true`** (or `-complex:*`) if your model requires true physical binary contacts rather than co-complex memberships.
3. **Filter by `minMIScore`** ($s \ge 0.60$ for medium-to-high confidence, $s \ge 0.72$ for high confidence) to remove single-detection or unverified interactions.
4. **Use `negativeFilter=POSITIVE_ONLY`** unless actively building a negative interaction benchmark set.
5. **Set `pageSize` appropriately** (typically 50–100) and paginate using `page=0, 1, 2...` until `last=true`.
