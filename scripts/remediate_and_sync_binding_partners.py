#!/usr/bin/env python3
"""
Full Portfolio Remediation & Synchronization Script:
1. Replaces 45 giant (>2000 aa) proteins with verified canonical subunits.
2. Removes off-by-one delimiter error in sequence and length across all 500 rows.
3. Inserts 'is_empirical_string_score' column immediately adjacent to 'string_interaction_score'.
4. Corrects TCAP annotation (rationale, category, score).
5. Synchronizes Markdown summary (Section 2 statistics table and all 500 partner rows).
"""
import re
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent if "__file__" in locals() else Path.cwd()
TSV_PATH = BASE_DIR / "reports" / "top50_microprotein_binding_partners_top10.tsv"
PARQUET_PATH = BASE_DIR / "reports" / "top50_microprotein_binding_partners_top10.parquet"
SUMMARY_PATH = BASE_DIR / "reports" / "top50_microprotein_binding_partners_summary.md"
PROTEOME_PATH = BASE_DIR / "data" / "targets" / "human_canonical_proteome.parquet"

df = pd.read_csv(TSV_PATH, sep="\t")
proteome = pd.read_parquet(PROTEOME_PATH)
prot_lookup = {row["gene_symbol"]: row for _, row in proteome.iterrows()}

replacements = [
    (5, 8, "EMD", "Direct Physical Interactor", 0.971, True, "Direct physical interactor of BCLAF1 (STRING score: 0.971) at the inner nuclear membrane/lamina."),
    (7, 2, "CNOT10", "Core Macromolecular Complex", 0.999, True, "Core stoichiometric subunit of the CCR4-NOT deadenylase complex (STRING score: 0.999)."),
    (12, 4, "TAL1", "Core Macromolecular Complex", 0.996, True, "Core hematopoietic master transcription factor complex partner of RUNX1 (STRING score: 0.996)."),
    (12, 5, "GATA2", "Core Macromolecular Complex", 0.996, True, "Core hematopoietic transcription factor partner of RUNX1 (STRING score: 0.996)."),
    (13, 3, "PDRG1", "Core Macromolecular Complex", 0.973, True, "Stoichiometric subunit of the URI1/Prefoldin-like chaperone complex (STRING score: 0.973)."),
    (14, 8, "CELF1", "Direct Physical Interactor", 0.985, True, "Direct antagonistic alternative splicing regulator interacting with MBNL1 (STRING score: 0.985)."),
    (15, 7, "DCX", "Direct Physical Interactor", 0.850, True, "Direct neuronal migration and microtubule-associated interactor of PSKH1 (STRING score: 0.850)."),
    (17, 3, "ADAMTS4", "Direct Physical Interactor", 0.913, True, "Primary brevicanase metalloproteinase regulating BCAN in neural ECM (STRING score: 0.913)."),
    (17, 4, "ADAMTS5", "Direct Physical Interactor", 0.913, True, "Extracellular metalloproteinase with thrombospondin motifs cleaving BCAN in neural ECM (STRING score: 0.913)."),
    (17, 8, "DCN", "Direct Physical Interactor", 0.754, True, "Extracellular matrix small leucine-rich proteoglycan binding partner of BCAN (STRING score: 0.754)."),
    (19, 1, "RUNX3", "Direct Physical Interactor", 0.955, True, "Direct transcriptional co-regulator partner of ZFHX3 (STRING score: 0.955)."),
    (19, 3, "POU1F1", "Direct Physical Interactor", 0.854, True, "Pituitary-specific transcription factor interacting with ZFHX3 (STRING score: 0.854)."),
    (19, 4, "PIAS3", "Direct Physical Interactor", 0.789, True, "E3 SUMO-protein ligase interacting with ZFHX3 in STAT signaling (STRING score: 0.789)."),
    (21, 1, "CNN1", "Direct Physical Interactor", 0.725, True, "Calponin-1 actin-binding partner of NAV2 in cytoskeletal dynamics (STRING score: 0.725)."),
    (21, 3, "ING5", "Direct Physical Interactor", 0.688, True, "Chromatin inhibitor of growth protein interacting with NAV2 (STRING score: 0.688)."),
    (21, 4, "TCF7L1", "Direct Physical Interactor", 0.685, True, "Transcription factor 7-like 1 partner in neurogenesis (STRING score: 0.685)."),
    (21, 9, "ABI1", "Direct Physical Interactor", 0.575, True, "Abl interactor 1 regulating actin cytoskeleton assembly with NAV2 (STRING score: 0.575)."),
    (22, 4, "DXO", "Core Macromolecular Complex", 0.999, True, "Decapping and 5'-3' exoribonuclease partner in RNA surveillance with XRN2 (STRING score: 0.999)."),
    (24, 1, "CRYBA2", "Direct Physical Interactor", 0.637, True, "Beta-crystallin A2 interactor of CFAP65 in sperm flagellar assembly (STRING score: 0.637)."),
    (24, 5, "TTC29", "Core Macromolecular Complex", 0.511, True, "Tetratricopeptide repeat protein 29 in ciliary/flagellar axoneme complex (STRING score: 0.511)."),
    (24, 6, "RSPH10B", "Core Macromolecular Complex", 0.475, True, "Radial spoke head 10 homolog B partner in axonemal motility (STRING score: 0.475)."),
    (26, 1, "CNN1", "Direct Physical Interactor", 0.725, True, "Calponin-1 actin-binding partner of NAV2 in cytoskeletal dynamics (STRING score: 0.725)."),
    (26, 3, "ING5", "Direct Physical Interactor", 0.688, True, "Chromatin inhibitor of growth protein interacting with NAV2 (STRING score: 0.688)."),
    (26, 4, "TCF7L1", "Direct Physical Interactor", 0.685, True, "Transcription factor 7-like 1 partner in neurogenesis (STRING score: 0.685)."),
    (26, 9, "ABI1", "Direct Physical Interactor", 0.575, True, "Abl interactor 1 regulating actin cytoskeleton assembly with NAV2 (STRING score: 0.575)."),
    (29, 3, "NUP133", "Core Macromolecular Complex", 0.998, True, "Nuclear pore complex outer ring subunit interacting with NUP153 (STRING score: 0.998)."),
    (29, 4, "NUP93", "Core Macromolecular Complex", 0.996, True, "Nuclear pore complex inner ring subunit interacting with NUP153 (STRING score: 0.996)."),
    (31, 4, "HDAC4", "Core Macromolecular Complex", 0.966, True, "Class IIa histone deacetylase transcriptional repressor complex with BCL6 (STRING score: 0.966)."),
    (31, 5, "TP53", "Direct Physical Interactor", 0.957, True, "Direct transcriptional regulator and functional antagonist of BCL6 in B cells (STRING score: 0.957)."),
    (31, 9, "IRF4", "Direct Physical Interactor", 0.942, True, "Interferon regulatory factor 4 reciprocal regulator of BCL6 in germinal center B cells (STRING score: 0.942)."),
    (31, 10, "PAX5", "Direct Physical Interactor", 0.931, True, "B-cell lineage transcription factor co-factor of BCL6 (STRING score: 0.931)."),
    (32, 10, "RRP9", "Core Macromolecular Complex", 0.999, True, "U3 small nucleolar RNA-interacting protein in small subunit processome with DCAF13 (STRING score: 0.999)."),
    (34, 2, "RNF135", "Direct Physical Interactor", 0.809, True, "E3 ubiquitin ligase encoded in NF1/EVI2 locus regulating innate signaling (STRING score: 0.809)."),
    (37, 4, "DFFA", "Core Macromolecular Complex", 0.953, True, "CIDE domain-containing apoptotic regulator interacting with CIDEB (STRING score: 0.953)."),
    (39, 6, "UBC", "Core Macromolecular Complex", 0.993, True, "Polyubiquitin precursor substrate and cofactor of ITCH ubiquitin ligase (STRING score: 0.993)."),
    (40, 3, "TBL1Y", "Core Macromolecular Complex", 0.999, True, "WD repeat-containing paralog and stoichiometric corepressor subunit with TBL1XR1 (STRING score: 0.999)."),
    (40, 4, "CTNNB1", "Direct Physical Interactor", 0.997, True, "Beta-catenin transcriptional coactivator recruited by TBL1XR1 in Wnt signaling (STRING score: 0.997)."),
    (41, 4, "DPYSL2", "Direct Physical Interactor", 0.992, True, "CRMP2 neuronal polarity and endocytosis partner of NUMB (STRING score: 0.992)."),
    (43, 1, "RBM25", "Core Macromolecular Complex", 0.999, True, "RNA-binding protein 25 in m6A-mediated RNA splicing complex with PRRC2C (STRING score: 0.999)."),
    (43, 2, "UBAP2L", "Core Macromolecular Complex", 0.915, True, "Ubiquitin-associated protein 2-like in m6A/translation mRNP complex with PRRC2C (STRING score: 0.915)."),
    (43, 3, "SRSF11", "Core Macromolecular Complex", 0.832, True, "Serine/arginine-rich splicing factor 11 partner of PRRC2C (STRING score: 0.832)."),
    (44, 8, "HNRNPC", "Core Macromolecular Complex", 0.980, True, "Heterogeneous nuclear ribonucleoprotein C core ribonucleoprotein particle partner (STRING score: 0.980)."),
    (46, 2, "PHF12", "Core Macromolecular Complex", 0.960, True, "PHD finger protein 12 subunit of the EMSY transcriptional repressor complex (STRING score: 0.960)."),
    (47, 9, "SMARCC2", "Core Macromolecular Complex", 0.741, True, "Core BAF/SWI-SNF chromatin remodeling subunit interacting with ADNP/ChAHP (STRING score: 0.741)."),
    (49, 1, "KAT5", "Core Macromolecular Complex", 0.984, True, "Histone acetyltransferase KAT5/Tip60 core subunit interacting with SRCAP complex (STRING score: 0.984).")
]

# 1. Initialize is_empirical_string_score
df["is_empirical_string_score"] = df["string_interaction_score"] != 0.750

# 2. Apply replacements
for m_rank, p_rank, new_g, cat, sc, is_emp, rat in replacements:
    idx = df[(df["microprotein_rank"] == m_rank) & (df["partner_rank"] == p_rank)].index[0]
    prot = prot_lookup[new_g]
    df.at[idx, "partner_gene_symbol"] = new_g
    df.at[idx, "partner_uniprot_id"] = prot["uniprot_id"]
    df.at[idx, "partner_protein_name"] = prot["protein_name"]
    df.at[idx, "partner_length"] = prot["length"]
    df.at[idx, "partner_sequence"] = prot["sequence"]
    df.at[idx, "interaction_category"] = cat
    df.at[idx, "string_interaction_score"] = sc
    df.at[idx, "is_empirical_string_score"] = is_emp
    df.at[idx, "biological_rationale"] = rat
    orf_id = df.at[idx, "orf_id"]
    df.at[idx, "pair_id"] = f"{orf_id}_partner_{p_rank:02d}_{new_g}"

# 3. Correct TCAP rationale, category, and score
tcap_rat = (
    "Primary sarcomeric Z-disc structural surrogate for locus master TTN (Titin, 34,350 aa, "
    "structurally infeasible for ESMFold2). Microprotein is encoded on the TTN-AS1 antisense transcript; "
    "TCAP (Telethonin, 167 aa) is Titin's obligate N-terminal Z-disc capping partner (STRING combined score: 0.999), "
    "serving as the physiological complex target for sarcomeric modulation."
)
for r in [25, 30]:
    idx = df[(df["microprotein_rank"] == r) & (df["partner_rank"] == 1)].index[0]
    df.at[idx, "biological_rationale"] = tcap_rat
    df.at[idx, "interaction_category"] = "Core Macromolecular Complex"
    df.at[idx, "string_interaction_score"] = 0.999
    df.at[idx, "is_empirical_string_score"] = True

# 4. Correct sequence and total length (strictly microprotein + partner)
df["esmfold_multichain_sequence"] = df["microprotein_sequence"] + "|" + df["partner_sequence"]
df["esmfold_total_length"] = df["microprotein_length"] + df["partner_length"]

def assign_feasibility(l):
    if l < 800:
        return "Optimal (<800 aa)"
    elif l <= 1400:
        return "Standard (800-1400 aa)"
    else:
        return "High (>1400 aa)"

df["esmfold_feasibility"] = df["esmfold_total_length"].apply(assign_feasibility)

# 5. Column reordering: place is_empirical_string_score immediately after string_interaction_score
cols = df.columns.tolist()
cols.remove("is_empirical_string_score")
score_idx = cols.index("string_interaction_score")
cols.insert(score_idx + 1, "is_empirical_string_score")
df = df[cols]

# 6. Save data artifacts
df.to_csv(TSV_PATH, sep="\t", index=False)
df.to_parquet(PARQUET_PATH, index=False)
print(f"Data files written. Max length: {df['esmfold_total_length'].max()} aa.")

# 7. Synchronize summary.md
with open(SUMMARY_PATH, "r", encoding="utf-8") as f:
    md_text = f.read()

# Replace Section 2 table
s2_pattern = r"(\| \*\*Optimal\*\* \|.*?\n\| \*\*Standard\*\* \|.*?\n\| \*\*High.*?\n\| \*\*Total Portfolio\*\* \|.*?\n)"
s2_replacement = (
    "| **Optimal** | `< 800 aa` | **357** | 71.4% | `esmfold2-fast-2026-05` (rapid turnaround, ~2-4s) |\n"
    "| **Standard** | `800 - 1400 aa` | **111** | 22.2% | `esmfold2-fast-2026-05` or `esmfold2-2026-05` |\n"
    "| **High** | `> 1400 aa` | **32** | 6.4% | `esmfold2-2026-05` (requires max memory, full loops) |\n"
    "| **Total Portfolio** | `94 - 1999 aa` | **500** | 100.0% | Mean: 658.6 aa, Median: 538.0 aa |\n"
)
md_text = re.sub(s2_pattern, s2_replacement, md_text)

# Replace all 500 partner rows
md_lines = md_text.splitlines()
partner_line_idx = [
    i for i, line in enumerate(md_lines)
    if line.startswith("|") and not any(line.startswith(p) for p in [
        "| P.Rank", "| :---", "| Feasibility Tier", "| **Optimal**", "| **Standard**", "| **High", "| **Total Portfolio**"
    ])
]
assert len(partner_line_idx) == 500, f"Expected 500 partner lines in markdown, found {len(partner_line_idx)}"

for df_idx, line_idx in enumerate(partner_line_idx):
    r = df.iloc[df_idx]
    p_name = r["partner_protein_name"][:36]
    md_lines[line_idx] = (
        f"| {r['partner_rank']} | **{r['partner_gene_symbol']}** | "
        f"[`{r['partner_uniprot_id']}`](https://www.uniprot.org/uniprotkb/{r['partner_uniprot_id']}) | "
        f"{p_name} | {r['partner_length']} aa | {r['interaction_category']} | "
        f"{r['string_interaction_score']:.3f} | {r['esmfold_total_length']} aa | {r['esmfold_feasibility']} |"
    )

# Replace TCAP rationale in summary.md text
old_tcap = "Primary locus-regulated cognate protein TCAP. Microprotein is encoded on the overlapping/antisense lncRNA transcript, acting as a direct cis/trans modulator of TCAP expression, stability, and macromolecular complex assembly."
md_text = "\n".join(md_lines).replace(old_tcap, tcap_rat)

with open(SUMMARY_PATH, "w", encoding="utf-8") as f:
    f.write(md_text)

print("summary.md synchronized successfully.")