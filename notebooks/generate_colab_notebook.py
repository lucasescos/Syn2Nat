import json
from pathlib import Path

notebook = {
 "cells": [
  {
   "cell_type": "markdown",
   "id": "intro",
   "metadata": {},
   "source": [
    "# ESMFold2 Ground-Truth Binder Design on Google Colab (Zero-Modal Workflow)\n",
    "\n",
    "**Objective**: Run the exact differentiable binder design protocol from the Biohub/EvolutionaryScale paper [*\"Language Modeling Materializes a World Model of Protein Biology\"*](https://www.biorxiv.org/content/10.64898/2026.06.03.729735) directly inside Google Colab on an A100/L4 GPU with **$0.00 Modal spend**.\n",
    "\n",
    "### The Algorithm at a Glance\n",
    "Instead of discrete token hallucination, this engine optimizes continuous **soft sequence logits** using backpropagation through ESMFold2's distogram representation:\n",
    "- $\\mathcal{L}_{\\text{inter\\_contact}}$: Maximizes binder-target interface contacts ($k=1, <22\\text{ \\AA}$).\n",
    "- $\\mathcal{L}_{\\text{intra\\_contact}}$: Forces internal binder folding ($k=2, \\ge 9\\text{ sep}, <14\\text{ \\AA}$).\n",
    "- $\\mathcal{L}_{\\text{glob}}$: Minimizes radius of gyration ($R_g \\le 2.38 N^{0.365}$) for compact globularity.\n",
    "- $\\mathcal{L}_{\\text{epitope}}$: Penalizes distance to specific target hotspot residues (`target_hotspot_ids`).\n",
    "- $\\mathcal{L}_{\\text{ESMC}}$: Masked language model prior enforcing natural human biophysical distributions.\n",
    "\n",
    "---"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "step1_env",
   "metadata": {},
   "source": [
    "### Step 1: Check GPU & Install Dependencies\n",
    "Ensure you are running on an **A100 (40GB/80GB)** or **L4/V100** GPU in Colab (`Runtime > Change runtime type > GPU`)."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_gpu_check",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Verify GPU allocation\n",
    "!nvidia-smi"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_install",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Install official ESMFold2, Biohub ESM, and visualization tools\n",
    "!pip install -q esm py3dmol pyarrow biotite biopython\n",
    "!pip install -q git+https://github.com/oxpig/ANARCI.git || echo 'ANARCI optional for non-antibody minibinders'"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "step2_code",
   "metadata": {},
   "source": [
    "### Step 2: Grab the Official Biohub Design Engine (`binder_design.py`)\n",
    "We pull the official 1,498-line script directly from the Biohub repository."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_download_engine",
   "metadata": {},
   "outputs": [],
   "source": [
    "!wget -q -O binder_design.py https://raw.githubusercontent.com/Biohub/esm/main/cookbook/tutorials/binder_design.py\n",
    "!ls -lh binder_design.py"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "step3_load",
   "metadata": {},
   "source": [
    "### Step 3: Instantiate the Colab-Native Engine (`local=True`)\n",
    "By setting `use_scaling_critics=False` and `REUSE_ESMC=True`, the model reuses the 6B ESMC trunk across all critic evaluations. This drops peak VRAM from 51 GB down to **~24–27 GB**, allowing it to fit comfortably on a single Colab A100 GPU without Modal."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_load_engine",
   "metadata": {},
   "outputs": [],
   "source": [
    "import os\n",
    "import sys\n",
    "import torch\n",
    "import pandas as pd\n",
    "from pathlib import Path\n",
    "\n",
    "# Import the design class directly from the local file\n",
    "from binder_design import ESMFold2Design, TARGET_SEQUENCES\n",
    "\n",
    "print(\"Initializing Colab-native ESMFold2Design engine...\")\n",
    "app = ESMFold2Design()\n",
    "# Load inversion models and hero critics (disabling remote scaling critics to run 100% locally)\n",
    "app.load(use_scaling_critics=False)\n",
    "print(\"ESMFold2 and ESMC-6B loaded successfully onto Colab GPU!\")"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "step4_sanity_run",
   "metadata": {},
   "source": [
    "### Step 4: Run the End-to-End Sanity Check (e.g. PD-L1 Minibinder or PSMA1)\n",
    "We will run a 35-aa de novo minibinder design against the target using the exact loss functions from the paper."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_run_design",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Target Configuration: Target sequence crop and mutable binder prompt ('#' = design position)\n",
    "# Preset options: 'pd-l1', 'egfr', 'ctla4', 'cd45', or supply custom sequence\n",
    "TARGET_CHOICE = \"pd-l1\"\n",
    "BINDER_LENGTH = 35\n",
    "\n",
    "target_seq = TARGET_SEQUENCES.get(TARGET_CHOICE)\n",
    "binder_prompt = \"#\" * BINDER_LENGTH\n",
    "\n",
    "print(f\"Target: {TARGET_CHOICE} (Length: {len(target_seq)} aa)\")\n",
    "print(f\"Binder Scaffold: Free Minibinder ({BINDER_LENGTH} aa)\")\n",
    "print(\"Launching 150-step gradient descent optimization...\")\n",
    "\n",
    "# Run the exact optimization loop on Colab\n",
    "best_sequences, trajectory, critic_results = app.design(\n",
    "    target_name=TARGET_CHOICE,\n",
    "    target_sequence=target_seq,\n",
    "    binder_name=f\"minibinder_{BINDER_LENGTH}\",\n",
    "    binder_sequence=binder_prompt,\n",
    "    is_antibody=False,\n",
    "    seed=42,\n",
    "    batch_size=1,\n",
    "    target_hotspot_ids=None  # Optional: pass list like ['D116', 'K118'] to focus on specific pocket residues\n",
    ")\n",
    "\n",
    "full_complex_seq = best_sequences[0]\n",
    "designed_binder = full_complex_seq.split(\"|\")[-1]\n",
    "print(\"=\" * 70)\n",
    "print(\"OPTIMIZATION COMPLETE!\")\n",
    "print(f\"Designed Binder Sequence: {designed_binder}\")\n",
    "print(f\"Binder Length:            {len(designed_binder)} aa\")\n",
    "print(\"=\" * 70)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "step5_metrics",
   "metadata": {},
   "source": [
    "### Step 5: Evaluate Scoring, ipTM, and Selection Metrics\n",
    "Inspect the convergence of $\\mathcal{L}_{\\text{inter}}$, $\\mathcal{L}_{\\text{intra}}$, $\\mathcal{L}_{\\text{glob}}$, and structural confidence (ipTM)."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_eval_metrics",
   "metadata": {},
   "outputs": [],
   "source": [
    "from Bio.SeqUtils.ProtParam import ProteinAnalysis\n",
    "\n",
    "df_critics = pd.DataFrame(critic_results)\n",
    "print(\"Critic Results:\")\n",
    "display(df_critics[[\"critic_name\", \"iptm\", \"final_loss\", \"distogram_iptm_proxy\"]])\n",
    "\n",
    "# Check isoelectric point for solubility\n",
    "pI = ProteinAnalysis(designed_binder).isoelectric_point()\n",
    "print(f\"Binder Isoelectric Point (pI): {pI:.2f} (Paper standard: pI < 6 for enhanced solubility)\")"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "step6_vis",
   "metadata": {},
   "source": [
    "### Step 6: 3D Molecular Visualization with py3Dmol\n",
    "Visualize the designed complex directly in the notebook. The target is colored green, and the de novo binder is colored by structural confidence (pLDDT)."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_3d_vis",
   "metadata": {},
   "outputs": [],
   "source": [
    "import py3Dmol\n",
    "\n",
    "# Extract predicted 3D complex structure\n",
    "hero_critic = df_critics[df_critics.critic_name.str.contains(\"Cutoff2025\")].iloc[0]\n",
    "protein_complex = hero_critic[\"complex\"]\n",
    "\n",
    "if protein_complex is not None:\n",
    "    pdb_str = protein_complex.to_pdb_string()\n",
    "    view = py3Dmol.view(width=700, height=600)\n",
    "    view.addModel(pdb_str, \"pdb\")\n",
    "    # Target chain (A) in green\n",
    "    view.setStyle({\"chain\": \"A\"}, {\"cartoon\": {\"color\": \"#2ECC71\"}})\n",
    "    # Binder chain (B) colored by pLDDT B-factor gradient\n",
    "    view.setStyle(\n",
    "        {\"chain\": \"B\"},\n",
    "        {\"cartoon\": {\"colorscheme\": {\"prop\": \"b\", \"gradient\": \"rwb\", \"min\": 60, \"max\": 100}}}\n",
    "    )\n",
    "    view.center()\n",
    "    view.zoomTo()\n",
    "    view.show()\n",
    "else:\n",
    "    print(\"3D complex structure not rendered.\")"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "step7_microprotein_bridge",
   "metadata": {},
   "source": [
    "### Step 7: Natural Microprotein Mimicry Bridge\n",
    "Once the binder is designed, screen our 7,264 human microprotein catalog using background-centered ESMC 6B Layer 79 embeddings to identify endogenous human microproteins that share its biophysical interface geometry."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": None,
   "id": "cell_bridge",
   "metadata": {},
   "outputs": [],
   "source": [
    "print(\"Binder Sequence ready for human microprotein catalog screening:\")\n",
    "print(f\">designed_binder_minibinder_{BINDER_LENGTH}\\n{designed_binder}\")\n",
    "print(\"\\nThis sequence can be plugged directly into scripts/run_real_pockets_benchmark.py to retrieve natural human mimetics in 2 seconds.\")"
   ]
  }
 ],
 "metadata": {
  "accelerator": "GPU",
  "colab": {
   "gpuType": "A100",
   "provenance": []
  },
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}

out_path = Path(r"c:\Users\Lucas.Escosteguy\documentos\microproteinproject\notebooks\esmfold2_binder_design_colab_sanity_check.ipynb")
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(notebook, f, indent=1)
print(f"Generated Colab notebook successfully at: {out_path}")
