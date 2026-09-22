import sys
from pathlib import Path
BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from esm.sdk.forge import SequenceStructureForgeInferenceClient, FoldingConfig
from esm.utils.structure.input_builder import StructurePredictionInput, ProteinInput
from scripts.test_esm3_api import get_key

def main():
    api_key = get_key()
    client = SequenceStructureForgeInferenceClient(
        model="esmfold2-fast-2026-05",
        url="https://biohub.ai",
        token=api_key
    )

    micro_seq = "MGSRGGKGRLWLPGEREQLRIPYEVPG" # c1riboseqorf157 (27 aa)
    
    # Load PSMA1 sequence
    import pandas as pd
    df_targets = pd.read_parquet("data/targets/human_canonical_proteome.parquet")
    psma1_seq = df_targets[df_targets["gene_symbol"] == "PSMA1"]["sequence"].iloc[0]

    print(f"Cofolding c1riboseqorf157 ({len(micro_seq)} aa) + PSMA1 ({len(psma1_seq)} aa)...")

    spi = StructurePredictionInput(
        sequences=[
            ProteinInput(sequence=micro_seq, id="A"),
            ProteinInput(sequence=psma1_seq, id="B")
        ]
    )
    config = FoldingConfig(num_loops=20, num_sampling_steps=100)

    res = client.fold_all_atom(spi, config=config)
    print(f"[SUCCESS] Cofolding completed!")
    print(f"  ipTM: {res.iptm:.4f}" if hasattr(res, 'iptm') and res.iptm is not None else "  ipTM: N/A")
    print(f"  pTM:  {res.ptm:.4f}" if hasattr(res, 'ptm') and res.ptm is not None else "  pTM: N/A")
    print(f"  Mean pLDDT: {res.plddt.mean().item():.2f}" if hasattr(res, 'plddt') and res.plddt is not None else "  pLDDT: N/A")

if __name__ == "__main__":
    main()
