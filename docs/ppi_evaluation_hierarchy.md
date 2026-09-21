# Computational Protein–Protein Interaction Evaluation Framework

**Target Domain**: High-Throughput Microprotein & Protein–Protein Complex Structure Analysis  
**Predictive Engines**: Biohub ESMFold2 (`fold_all_atom`, `fold`), AlphaFold-Multimer (v2/v3)  
**Implementation Tool**: [`scripts/extract_interface_contacts.py`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/scripts/extract_interface_contacts.py)  
**Associated Repositories**:
- PDB Coordinates: [`structures/pdbs/`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/structures/pdbs/)
- PAE Matrices: [`structures/pae/`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/structures/pae/)
- Master Structural Atlas: [`reports/master_microprotein_ppi_structural_atlas.parquet`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/reports/master_microprotein_ppi_structural_atlas.parquet)

---

## 1. The 4-Pillar Binding Evaluation Hierarchy

Evaluating whether two proteins truly form a functional, physical complex—as opposed to a non-specific surface collision or in silico artifact—requires moving through four tiers of biophysical evidence, from macroscopic topological confidence down to atom-specific bond chemistry.

```
                      ESMFold2 / AlphaFold PPI Evaluation Hierarchy
 ┌──────────────────────────────────────────────────────────────────────────────┐
 │ PILLAR 1: Interface Topological Confidence (ipTM & Model Ranking Score)       │
 │   • Evaluates inter-chain packing accuracy independent of monomer size       │
 │   • Standard metric: 0.8·ipTM + 0.2·pTM                                     │
 └──────────────────────────────────────┬───────────────────────────────────────┘
                                        │
 ┌──────────────────────────────────────▼───────────────────────────────────────┐
 │ PILLAR 2: Inter-Chain Predicted Aligned Error (PAE Positional Certainty)     │
 │   • Off-diagonal block inspection (PAE_AB and PAE_BA)                        │
 │   • Distinguishes rigid binding (<5 Å) from floating ambiguity (>15 Å)       │
 └──────────────────────────────────────┬───────────────────────────────────────┘
                                        │
 ┌──────────────────────────────────────▼───────────────────────────────────────┐
 │ PILLAR 3: Structural Plasticity & Induced Folding Footprint                  │
 │   • Surface involvement (% of microprotein engaged, typically >50%)          │
 │   • Interface pLDDT shift (unstructured monomer -> ordered bound strand/helix)│
 └──────────────────────────────────────┬───────────────────────────────────────┘
                                        │
 ┌──────────────────────────────────────▼───────────────────────────────────────┐
 │ PILLAR 4: Contact Chemistry & Physical Bond Topology                         │
 │   • Direct H-Bonds (≤3.5 Å), Salt Bridges (≤4.0 Å), Hydrophobic Packing (≤4.5 Å)│
 │   • Solvent exclusion and complementary shape matching                        │
 └──────────────────────────────────────────────────────────────────────────────┘
```

---

### Pillar 1: Interface Confidence ($\text{ipTM}$) & Ranking Score

In heterodimers with severe size asymmetry (e.g. a 30 aa microprotein paired with a 500 aa canonical enzyme), the global TM-score ($\text{pTM}$) is dominated by the large partner. An algorithm can achieve $\text{pTM} > 0.85$ even if the microprotein is floating freely in empty space.

**Interface Predicted TM-Score ($\text{ipTM}$)** solves this by calculating the TM-score exclusively across residue pairs belonging to different chains:

$$\text{ipTM} = \frac{1}{|\text{Interface}|} \sum_{i \in \text{Chain A}, j \in \text{Chain B}} \frac{1}{1 + \left(\frac{d_{ij} - \hat{d}_{ij}}{d_0}\right)^2}$$

#### Standard Interpretation Thresholds:
| Score Tier | Threshold | Biophysical & Biochemical Meaning | Empirical Human Proteome Correlation |
| :--- | :---: | :--- | :--- |
| **High Confidence** | **$\text{ipTM} \ge 0.75$** | **Stoichiometric, Obligate Complex**: High shape complementarity, rigid docking, low dissociation constant ($K_d \sim \text{nM}$). | Stable multi-protein machines (e.g. proteasome rings, PCNA clamps, serpin complexes). |
| **Plausible Interface** | **$0.60 \le \text{ipTM} < 0.75$** | **Specific Macromolecular Interface**: Biologically authentic regulatory binding, substrate engagement, or allosteric tuning ($K_d \sim \mu\text{M}$). | Kinase/phosphatase adaptors, histone acetyltransferase recruitment, transcription factors. |
| **Transient / Weak** | **$0.45 \le \text{ipTM} < 0.60$** | **Transient or Induced Fit Interface**: Weak or dynamic interaction; requires cofactors, membranes, or post-translational modifications. | Signaling scaffolds, flexible regulatory tails, transient chaperones. |
| **Non-Interacting** | **$\text{ipTM} < 0.45$** | **Unbound / Spatial Collision**: Positional uncertainty too high to support physical complexation in isolation. | Independent cellular proteins or spatial non-binders. |

#### Composite Ranking Score:
Following DeepMind (Evans et al., *bioRxiv* 2022) and ESMFold benchmarks, complexes are ranked by:
$$\mathbf{\text{Model Score} = 0.8 \cdot \text{ipTM} + 0.2 \cdot \text{pTM}}$$

---

### Pillar 2: Inter-Chain Predicted Aligned Error ($\text{PAE}_{\text{inter}}$)

The Predicted Aligned Error (PAE) matrix is an $L \times L$ array where entry $\text{PAE}(i, j)$ estimates the expected positional error (in \AA) of residue $i$'s $C_\alpha$ when the structure is aligned on residue $j$.

For a complex of Chain A (length $L_A$) and Chain B (length $L_B$):
- **Diagonal blocks** ($L_A \times L_A$ and $L_B \times L_B$): Monomer folding confidence.
- **Off-diagonal blocks** ($L_A \times L_B$ and $L_B \times L_A$): **Relative chain orientation certainty**.

```
       Chain A (uProt)      Chain B (Partner)
      ┌──────────────────┬──────────────────┐
    A │  Monomer A PAE   │  Inter-Chain PAE │ <--- Focus here: Low error (<5 Å)
      │  (uProt folding) │   (A aligned B)  │      indicates locked orientation
      ├──────────────────┼──────────────────┤
    B │  Inter-Chain PAE │  Monomer B PAE   │
      │   (B aligned A)  │ (Partner folding)│
      └──────────────────┴──────────────────┘
```

#### Interpretation Criteria:
- **$\text{PAE} < 5.0\text{ \AA}$**: Sub-nanometer positional locking. The relative orientation and translation of the two chains is predicted with very high precision.
- **$\text{PAE} < 1.0\text{ \AA}$**: Atomic-resolution certainty. Sidechain rotamers and backbone hydrogen bonds can be inferred with high reliability.
- **$\text{PAE} > 15.0\text{ \AA}$**: High positional ambiguity. The model predicts both proteins may be folded, but cannot place them relative to each other.

---

### Pillar 3: Structural Plasticity & Induced Folding Footprint

Small non-canonical microproteins (15–100 aa) frequently exist as intrinsically disordered peptides in the cytoplasm or nucleus. When they encounter their physiological partner, they often undergo **coupled folding and binding** (induced fit).

#### Two Quantitative Metrics:
1. **Microprotein Surface Involvement (% Footprint)**:
   $$\text{Footprint \%} = \frac{\text{Number of Microprotein Residues with Atomic Distance} \le 4.5\text{ \AA}}{\text{Total Microprotein Length}} \times 100$$
   - Genuine microprotein adaptors typically bury **$\ge 50\text{--}70\%$** of their primary sequence in the binding groove of their partner.
2. **Interface pLDDT ($\text{pLDDT}_{\text{interface}}$)**:
   - Calculate the mean pLDDT restricted specifically to interface contact residues.
   - If a microprotein has bulk monomer $\text{pLDDT} \sim 35\text{--}45\%$ (disordered) but shifts to $\text{pLDDT}_{\text{interface}} > 70\text{--}80\%$ at the contact site, it indicates an **induced structural fold** (e.g. forming an amphipathic $\alpha$-helix or $\beta$-sheet upon binding).

---

### Pillar 4: Contact Chemistry & Physical Bond Topology

A genuine interface is governed by physical chemistry, not just geometric proximity. Atomic interactions must be classified according to strict biophysical distance constraints:

| Interaction Type | Participating Chemical Groups | Maximum Distance Cutoff | Structural Significance |
| :--- | :--- | :---: | :--- |
| **Hydrogen Bond** | Backbone or sidechain donor/acceptor (N–O, O–N, O–O) | **$\le 3.5\text{ \AA}$** | Provides directional binding specificity and secondary structure alignment (e.g. $\beta$-strand addition). |
| **Salt Bridge** | Acidic oxygens (Asp OD, Glu OE) to Basic nitrogens (Arg NH, Lys NZ, His ND/NE) | **$\le 4.0\text{ \AA}$** | High-energy electrostatic anchor; dictates binding orientation and pH sensitivity. |
| **Hydrophobic Packing** | Aliphatic/aromatic sidechains (Ala, Val, Leu, Ile, Met, Phe, Trp, Pro) | **$\le 4.5\text{ \AA}$** | Primary thermodynamic driving force ($\Delta G_{\text{solv}}$) via water exclusion. |
| **Aromatic Stacking** | Aromatic rings (Phe, Tyr, Trp, His) | **$\le 4.8\text{ \AA}$** | $\pi$-$\pi$ and cation-$\pi$ interactions stabilizing core packing. |
| **Van der Waals** | Any non-hydrogen heavy atom pair | **$\le 4.0\text{ \AA}$** | Steric surface complementarity and van der Waals attraction. |

---

## 2. Residue-Level Contact Extraction Protocol

To extract exact residue-by-residue interactions from any cofolded complex, the [`scripts/extract_interface_contacts.py`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/scripts/extract_interface_contacts.py) pipeline performs the following algorithm:

1. **3D Coordinate Parsing**:
   - Parses all Cartesian coordinates $(x, y, z)$ and B-factors (pLDDT) from the PDB file, partitioned into Chain A (Microprotein) and Chain B (Partner).
2. **All-Atom Distance Matrix Computation**:
   - For every residue pair $(i \in \text{Chain A}, j \in \text{Chain B})$, computes pairwise Euclidean distances across all heavy atoms:
     $$d(i, j) = \min_{a \in \text{atoms}(i), b \in \text{atoms}(j)} \sqrt{(x_a - x_b)^2 + (y_a - y_b)^2 + (z_a - z_b)^2}$$
3. **Chemistry Classification**:
   - Filters pairs with $d(i, j) \le 4.5\text{ \AA}$. Identifies atom names and residue types to assign the interaction class (Salt Bridge, Hydrogen Bond, Hydrophobic Packing, Van der Waals).
4. **Confidence Weighting via PAE Matrix**:
   - Looks up the corresponding entries $\text{PAE}(i, j)$ and $\text{PAE}(j, i)$ in the complex's PAE JSON matrix. Assigns the mean inter-chain PAE to the contact pair.
5. **Output Generation**:
   - Exports a structured TSV file listing every interacting residue pair with atom names, distance, bond type, and PAE.
   - Generates a standalone PyMOL (`.pml`) script for visual inspection.

---

## 3. Validated Case Studies from the Human Microprotein Atlas

### Case 1: High-Confidence Stoichiometric Complex (`SNX13` + `SERPINE2`)
- **Topological Scores**: $\text{ipTM} = 0.8525$, $\text{pTM} = 0.8226$, Partner $\text{pLDDT} = 82.0\%$
- **PDB Structure**: [`structures/pdbs/c7riboseqorf18_partner_01_SERPINE2.pdb`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/structures/pdbs/c7riboseqorf18_partner_01_SERPINE2.pdb)
- **Interface Chemistry**: 208 atomic contacts, **27 Hydrogen Bonds**, **5 Salt Bridges**, 33 Hydrophobic packings.
- **Atomic-Resolution Core**: Residues 40–52 of the 52 aa SNX13 microprotein dock into the SERPINE2 surface with **$\text{PAE} < 1.0\text{ \AA}$**:
  - `Thr45(OG1)` $\rightarrow$ `Gln341(NE2)`: $2.69\text{ \AA}$ (H-bond, **$\text{PAE} = 0.65\text{ \AA}$**)
  - `Val47(N)` $\rightarrow$ `Leu183(O)`: $2.77\text{ \AA}$ (H-bond, **$\text{PAE} = 0.65\text{ \AA}$**)
  - `Phe41(N)` $\rightarrow$ `Phe189(O)`: $2.77\text{ \AA}$ (H-bond, **$\text{PAE} = 0.99\text{ \AA}$**)
  - `Arg38(NE)` $\rightarrow$ `Asp350(OD2)`: $2.46\text{ \AA}$ (**Salt Bridge**)

### Case 2: Induced $\beta$-Strand Addition in NuA4 Complex (`MORF4L2` + `KAT5`)
- **Topological Scores**: $\text{ipTM} = 0.7119$, $\text{pTM} = 0.6905$, Partner $\text{pLDDT} = 81.2\%$
- **PDB Structure**: [`structures/pdbs/cXriboseqorf59_partner_02_KAT5.pdb`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/structures/pdbs/cXriboseqorf59_partner_02_KAT5.pdb)
- **Interface Chemistry**: 57 atomic contacts, **12 Hydrogen Bonds**, 9 Hydrophobic packings.
- **Mechanism**: The 32 aa Tier 1 peptidein forms an anti-parallel $\beta$-strand addition against the catalytic MYST acetyltransferase domain:
  - `Gln22(O)` $\rightarrow$ `Leu474(N)`: $2.90\text{ \AA}$ (H-bond, **$\text{PAE} = 3.05\text{ \AA}$**)
  - `Val23(N)` $\rightarrow$ `Ile437(O)`: $2.94\text{ \AA}$ (H-bond, **$\text{PAE} = 3.32\text{ \AA}$**)
  - `Thr24(N)` $\rightarrow$ `Leu474(O)`: $2.86\text{ \AA}$ (H-bond, **$\text{PAE} = 3.54\text{ \AA}$**)
  - `Phe21(N)` $\rightarrow$ `Ile439(O)`: $2.87\text{ \AA}$ (H-bond, **$\text{PAE} = 3.60\text{ \AA}$**)

### Case 3: Phosphatase Catalytic Cleft Docking (`CDC123` + `PPP1CA`)
- **Topological Scores**: $\text{ipTM} = 0.6500$, $\text{pTM} = 0.8617$, Partner $\text{pLDDT} = 90.1\%$
- **PDB Structure**: [`structures/pdbs/c10norep31_partner_05_PPP1CA.pdb`](file:///c:/Users/Lucas.Escosteguy/documentos/microproteinproject/structures/pdbs/c10norep31_partner_05_PPP1CA.pdb)
- **Interface Chemistry**: 62 atomic contacts, **9 Hydrogen Bonds**, **1 Salt Bridge** (`Arg9` $\leftrightarrow$ `Glu54`), 9 Hydrophobic packings.
- **Mechanism**: The 29 aa dual-coding intORF peptide docks directly into the PP1 phosphatase active groove:
  - `Ser22(O)` $\rightarrow$ `Cys291(N)`: $2.75\text{ \AA}$ (H-bond, **$\text{PAE} = 3.23\text{ \AA}$**)
  - `Ser22(N)` $\rightarrow$ `Leu289(O)`: $3.11\text{ \AA}$ (H-bond, **$\text{PAE} = 3.47\text{ \AA}$**)
  - `Val21(N)` $\rightarrow$ `Asp242(OD2)`: $2.87\text{ \AA}$ (H-bond, **$\text{PAE} = 3.67\text{ \AA}$**)

---

## 4. Operational CLI Usage

To execute this evaluation protocol on any coordinate file:

```bash
python scripts/extract_interface_contacts.py <path_to_pdb> --out-tsv <path_to_tsv> --out-pml <path_to_pml>
```

### CLI Arguments:
- `pdb`: Path to the `.pdb` coordinate file (e.g. `structures/pdbs/cXriboseqorf59_partner_02_KAT5.pdb`).
- `--pae`: *(Optional)* Path to the PAE JSON matrix. Defaults automatically to `structures/pae/<pdb_stem>_pae.json`.
- `--cutoff`: Distance cutoff in \AA for contacts (default: `4.5`).
- `--out-tsv`: Path to save the tabular residue contact table.
- `--out-pml`: Path to save the PyMOL 3D visualization script.
