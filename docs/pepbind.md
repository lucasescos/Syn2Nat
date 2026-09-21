# PepBind: Peptide Binding Protein Database

> **Original Portal**: [http://pepbind.bicpu.edu.in/](http://pepbind.bicpu.edu.in/)  
> **Prediction Server**: [http://pepbind.bicpu.edu.in/PepBind_prediction_beta.php](http://pepbind.bicpu.edu.in/PepBind_prediction_beta.php)  
> **Reference Publication**: Das AA, Sharma OP, Kumar MS, Krishna R, Mathur PP. *PepBind: A comprehensive database and computational tool for analysis of protein–peptide interactions*. **Genomics, Proteomics & Bioinformatics** (2013) 11(4):241–246. [DOI: 10.1016/j.gpb.2013.03.002](https://doi.org/10.1016/j.gpb.2013.03.002) | [PMCID: PMC4357787](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4357787/)  
> **Development Center**: Centre for Bioinformatics, Pondicherry University, Puducherry, India  
> **Related Implementations**: [`scripts/pepbind_pici.py`](../scripts/pepbind_pici.py) (Local PICI geometric engine)

---

## 1. Overview & Database Architecture

**PepBind** (Peptide Binding Protein Database) is a specialized structural bioinformatics database developed to systematically catalog, analyze, and predict **protein–peptide interactions**. It captures complexes where a globular protein or domain binds a flexible, short linear peptide (typically $\le 35$ amino acid residues).

### Core Features
* **Curated Dataset**: Contains **3,100+** protein–peptide complexes derived from the Protein Data Bank (PDB).
* **Classification Scheme**:
  * **Peptide Length**: Stratified into dipeptide, tripeptide, oligopeptide, and polypeptide categories ($\le 35$ residues).
  * **Cellular Activity / Pathway**: Categorized across **19 cellular functions**, including regulatory pathways (>40%), immune system (~20%), hydrolases/proteases (~30%), transferases, and membrane transport.
  * **Structure Determination Method**: Grouped by X-ray diffraction, NMR spectroscopy, and Electron Microscopy (EM).
* **Interface Computation (PICI Tool)**: Calculates residue-level contact interfaces (H-bonds, hydrophobic contacts, salt bridges/ionic interactions, and disulfide bonds).
* **Predictive Capability**: Provides a domain binding prediction server to suggest protein domains likely to interact with a user-supplied peptide sequence.

```
                         [ Protein Data Bank (PDB) ]
                                      │
                         [ Screening & Filtering ]
                        (Peptide length <= 35 aa)
                                      │
                                      ▼
                        [ PepBind Core Database ]
                 (3,100+ Curated Protein-Peptide Complexes)
                                      │
        ┌─────────────────────────────┼─────────────────────────────┐
        ▼                             ▼                             ▼
 [ Classification ]             [ PICI Tool ]            [ Search & Discovery ]
 - Length (di/tri/oligo)        - H-bonds (<=3.5 Å)      - Simple & Keyword Search
 - 19 Cellular Functions        - Hydrophobic (<=5.0 Å)  - Advanced Multi-Filter
 - Experimental Methods         - Ionic (<=6.0 Å)        - BLAST & FATCAT Search
        │                             │                             │
        └─────────────────────────────┼─────────────────────────────┘
                                      ▼
                          [ Data Access & Outputs ]
                    - PDB Coordinates (v3.30)
                    - PDBML XML Data Files (v4.0)
                    - FASTA Sequences
                    - PICI Interaction Tables & PDF Reports
```

---

## 2. API Status & Programmatic Access Reality

### The Native API Reality
PepBind was architected in 2012–2013 on an Apache 2.2 / PHP 5 / MySQL 5.1 stack. 
* **No Formal REST/OpenAPI Specification**: PepBind **does not** provide an official JSON REST API, GraphQL endpoint, or OpenAPI/Swagger specification.
* **Server Accessibility & Uptime**: The university host (`pepbind.bicpu.edu.in` / `bicpu.edu.in`) is an institutional academic server that experiences intermittent downtime, connection resets (`[WinError 10054]`), and network timeouts.

### Disambiguation with Homonymous Tools
When searching for "PepBind", bioinformaticians frequently encounter three distinct tools:
1. **PepBind (Pondicherry University, 2013)**: The original structural database and PICI interaction tool covered in this document.
2. **PepBind (Yang Lab / Shandong University, 2021)**: A machine-learning sequence predictor for peptide-binding residues combining SVMpep, S-SITE, and TM-SITE (`yanglab.qd.sdu.edu.cn/PepBind`).
3. **PepBDB (Huang Lab / HUST, 2018)**: Peptide-Protein Binding Database (~14,000+ complexes), which also lacks a REST API but provides downloadable flat files and machine learning datasets (`PepBDB-ML`).

---

## 3. Legacy Web Scraping & Scripted HTTP Access

For workflows querying the original server, programmatic access is conducted by mimicking browser form submissions (`GET` and `POST`) against the underlying PHP scripts:

### Key Script Endpoints

| Endpoint | Method | Parameters | Description |
| :--- | :--- | :--- | :--- |
| `http://pepbind.bicpu.edu.in/search.php` | `GET` / `POST` | `pdb_id`, `pname` | Simple search by 4-letter PDB code or protein name |
| `http://pepbind.bicpu.edu.in/keyword_search.php` | `POST` | `keyword` | Full-text search across all database tables |
| `http://pepbind.bicpu.edu.in/advanced_search.php` | `POST` | `pep_length`, `activity`, `method`, `author` | Multi-filter search joined by logical `AND` |
| `http://pepbind.bicpu.edu.in/blast.php` | `POST` | `sequence`, `database`, `evalue`, `matrix` | BLAST sequence similarity against PepBind/PDB |
| `http://pepbind.bicpu.edu.in/PepBind_prediction_beta.php` | `POST` | `pep_seq`, `threshold` | Peptide domain-binding prediction server |
| `http://pepbind.bicpu.edu.in/details.php` | `GET` | `id=<PDB_ID>` | Summary, sequences, Ramachandran plot, external links |
| `http://pepbind.bicpu.edu.in/download.php` | `GET` | `type=all` or `id=<PDB_ID>` | Coordinate files, PDBML XML, and PICI tables |

### Python Scripting Example (Legacy Server Client)

```python
import requests
from bs4 import BeautifulSoup

BASE_URL = "http://pepbind.bicpu.edu.in"

def query_pepbind_entry(pdb_id: str):
    """
    Query an individual entry from PepBind via HTTP GET.
    """
    url = f"{BASE_URL}/details.php"
    params = {"id": pdb_id.upper()}
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    
    try:
        response = requests.get(url, params=params, headers=headers, timeout=10)
        response.raise_for_status()
        
        soup = BeautifulSoup(response.text, "html.parser")
        # Extract title and interaction details from tables
        title = soup.find("title").text if soup.find("title") else "N/A"
        return {"status": "success", "title": title, "raw_html": response.text}
    except requests.exceptions.RequestException as e:
        return {"status": "error", "message": str(e)}
```

---

## 4. The PICI (Protein Inter-Chain Interaction) Engine

A principal contribution of PepBind is the **PICI** algorithm, which calculates interfacial contacts between the peptide chain and the receptor protein chain using 3D Euclidean distances:

$$D(A, B) = \sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2}$$

### Distance Cutoffs & Chemical Rules

1. **Hydrogen Bonds ($D \le 3.5\,\text{Å}$)**:
   * Computed between any donor/acceptor nitrogen (N) and oxygen (O) pairs between the peptide chain and the protein chain.
2. **Hydrophobic Contacts ($D \le 5.0\,\text{Å}$)**:
   * Computed between carbon atoms of non-polar hydrophobic residues:
   * **Residues**: Ala (A), Val (V), Leu (L), Ile (I), Met (M), Phe (F), Trp (W), Pro (P), Tyr (Y).
3. **Ionic / Electrostatic Interactions ($D \le 6.0\,\text{Å}$)**:
   * Computed between charged sidechain functional atoms (or centroid distances):
   * **Cationic residues**: Arg (R), Lys (K), His (H)
   * **Anionic residues**: Asp (D), Glu (E)
4. **Disulfide Bridges ($D \le 2.2\,\text{Å}$)**:
   * Computed between cysteine sulfur ($S_\gamma$) atoms.

---

## 5. Modern Programmatic Alternatives via REST APIs

Because PepBind is built strictly from Protein Data Bank structures and geometric interface rules, any modern pipeline can reproduce and expand PepBind's functionality using active, production-grade biological APIs.

### 5.1 RCSB PDB Search & Data API (RESTful)

The **RCSB PDB API** allows automated querying of all protein–peptide complexes using exact peptide length filters ($\le 35$ aa) and structured JSON responses.

#### Endpoint 1: Search for Protein-Peptide Complexes

* **URL**: `POST https://search.rcsb.org/rcsbsearch/v2/query`
* **Content-Type**: `application/json`

```json
{
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "entity_poly.rcsb_sample_sequence_length",
          "operator": "less_or_equal",
          "value": 35
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entry_info.polymer_entity_count_protein",
          "operator": "greater_or_equal",
          "value": 2
        }
      }
    ]
  },
  "return_type": "entry",
  "request_options": {
    "paginate": {
      "start": 0,
      "rows": 25
    }
  }
}
```

#### Endpoint 2: Retrieve Full Entry & Polymer Metadata

* **URL**: `GET https://data.rcsb.org/rest/v1/core/entry/{entry_id}`
* **URL**: `GET https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}`

---

### 5.2 PDBe REST API & PDBe PISA Interfaces

The **Protein Data Bank in Europe (PDBe)** provides specialized REST API endpoints for interface contacts, chemical properties, and quaternary structure assemblies.

#### Interface Analysis Endpoint
* **URL**: `GET https://www.ebi.ac.uk/pdbe/api/pisa/interactions/{pdb_id}`
* **Capabilities**: Returns calculated hydrogen bonds, salt bridges, disulfide bonds, and buried surface areas ($\Delta \text{G}$, $\text{BSA}$) between specific chains.

#### Entry Summary Endpoint
* **URL**: `GET https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}`

---

### 5.3 Local Python Execution of the PepBind PICI Engine

To perform PepBind-style interface calculations without relying on external web servers, use the local implementation [`scripts/pepbind_pici.py`](../scripts/pepbind_pici.py).

```python
import sys
from Bio.PDB import PDBParser, NeighborSearch

def calculate_pepbind_interfaces(pdb_file: str, protein_chain: str, peptide_chain: str):
    """
    Computes H-bonds, hydrophobic contacts, and ionic interactions 
    strictly according to PepBind's PICI specifications.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    model = next(iter(structure))  # first NMR model or X-ray crystal
    
    chain_prot = model[protein_chain]
    chain_pep = model[peptide_chain]
    
    pep_atoms = [a for a in chain_pep.get_atoms() if a.element != "H"]
    prot_atoms = [a for a in chain_prot.get_atoms() if a.element != "H"]
    
    ns = NeighborSearch(prot_atoms)
    
    h_bonds = []
    hydrophobic = []
    ionic = []
    
    hydrophobic_res = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}
    pos_res = {"ARG", "LYS", "HIS"}
    neg_res = {"ASP", "GLU"}
    
    for atom_p in pep_atoms:
        res_p = atom_p.get_parent()
        # Find all protein atoms within 6.0 Å
        neighbors = ns.search(atom_p.coord, 6.0)
        for atom_r in neighbors:
            res_r = atom_r.get_parent()
            dist = atom_p - atom_r
            
            # 1. Hydrogen bonds (N/O atoms <= 3.5 Å)
            if dist <= 3.5 and atom_p.element in ("N", "O") and atom_r.element in ("N", "O"):
                h_bonds.append((res_p.get_resname(), res_p.id[1], res_r.get_resname(), res_r.id[1], dist))
                
            # 2. Hydrophobic interactions (<= 5.0 Å)
            if dist <= 5.0 and res_p.get_resname() in hydrophobic_res and res_r.get_resname() in hydrophobic_res:
                if atom_p.element == "C" and atom_r.element == "C":
                    hydrophobic.append((res_p.get_resname(), res_p.id[1], res_r.get_resname(), res_r.id[1], dist))
                    
            # 3. Ionic interactions (<= 6.0 Å between opposite charges)
            if dist <= 6.0:
                is_salt_bridge = (
                    (res_p.get_resname() in pos_res and res_r.get_resname() in neg_res) or
                    (res_p.get_resname() in neg_res and res_r.get_resname() in pos_res)
                )
                if is_salt_bridge:
                    ionic.append((res_p.get_resname(), res_p.id[1], res_r.get_resname(), res_r.id[1], dist))
                    
    return {
        "hydrogen_bonds": len(h_bonds),
        "hydrophobic_contacts": len(hydrophobic),
        "ionic_interactions": len(ionic)
    }
```

---

## 6. Comparison of Peptide-Binding Database APIs

| Feature | PepBind (Pondicherry Univ.) | PepBDB (HUST) | RCSB PDB REST / GraphQL | PDBe REST API |
| :--- | :--- | :--- | :--- | :--- |
| **Status** | Historical / Intermittent | Active (Web/Files) | Active (99.9% Uptime) | Active (99.9% Uptime) |
| **Direct REST API** | ❌ No | ❌ No | ✅ Full REST & GraphQL | ✅ Full REST API |
| **Entries** | 3,100 complexes | 14,000+ complexes | 220,000+ entries | 220,000+ entries |
| **Peptide Length Filter** | Built-in ($\le 35$ aa) | Curated complexes | Configurable query | Configurable query |
| **Interface Contacts** | PICI (H-bond, Hydrophobic, Ionic) | Contact tables | Via validation sub-endpoints | PDBe PISA interface service |
| **Programmatic Access** | HTTP Scraping / Flat files | Bulk tarballs / GitHub | `requests` / `urllib` / Python SDK | `requests` / `urllib` |
| **Format** | HTML / PDBML / Text | `.tar.gz` / PDB files | Standard JSON / XML | Standard JSON |

---

## 7. Recommended Workflow for Microprotein / Peptide Research

When analyzing peptide-protein interactions or microprotein binding partners in this project:
1. **Query Complex Candidates**: Use the **RCSB Search API** to identify relevant structural homologues with peptide-length criteria ($\le 35$ aa or custom length).
2. **Retrieve High-Resolution Coordinates**: Fetch `.cif` or `.pdb` directly from RCSB or AlphaFold Protein Structure Database.
3. **Calculate Interfaces**: Execute [`scripts/pepbind_pici.py`](../scripts/pepbind_pici.py) or PDBe PISA API to extract full contact networks (H-bonds, salt bridges, hydrophobic packing).
4. **Score & Embed Interactions**: Feed interfacial contacts into **MAPPIE** (`docs/mappie.md`) or **ESM-3 / ESMFold2** (`docs/biohub.md`) for interaction scoring and binding motif generation.
