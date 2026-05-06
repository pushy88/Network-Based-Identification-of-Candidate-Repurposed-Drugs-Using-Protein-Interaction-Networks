"""
RUN ALL — Full Alzheimer's Drug Repurposing Pipeline

Executes all 4 steps in sequence:
  Step 1 → Fetch AD genes from DisGeNET
  Step 2 → Build PPI network via STRING API
  Step 3 → Network analysis, identify hub genes
  Step 4 → Query DGIdb for drug-gene interactions

Run from the project root:
    python run_pipeline.py
"""

import sys
import os
import time

# Add scripts/ to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scripts"))

def print_header(step: int, title: str):
    print("\n")
    print("╔" + "═" * 62 + "╗")
    print(f"║  STEP {step}: {title:<54}║")
    print("╚" + "═" * 62 + "╝")


def run_pipeline():
    print("\n" + "█" * 64)
    print("█                                                              █")
    print("█   ALZHEIMER'S DISEASE — NETWORK-BASED DRUG REPURPOSING      █")
    print("█   Computational Biology Pipeline v1.0                        █")
    print("█                                                              █")
    print("█" * 64)

    start_time = time.time()

    # ── Step 1: Disease genes ─────────────────────────────────────────────────
    print_header(1, "Fetch Alzheimer's Genes from DisGeNET")
    from step1_get_disease_genes import fetch_disgenet_genes
    genes = fetch_disgenet_genes()

    if not genes:
        print("✗ No genes found. Exiting.")
        sys.exit(1)

    time.sleep(1)

    # ── Step 2: Build PPI network ─────────────────────────────────────────────
    print_header(2, "Build PPI Network via STRING API")
    from step2_build_network import build_string_network
    G, gene_names = build_string_network()

    time.sleep(1)

    # ── Step 3: Network analysis ──────────────────────────────────────────────
    print_header(3, "Network Analysis — Hub Gene Identification")
    from step3_network_analysis import analyze_network
    hub_df, G = analyze_network()

    time.sleep(1)

    # ── Step 4: Drug-gene interactions ────────────────────────────────────────
    print_header(4, "Query DGIdb for Drug-Gene Interactions")
    from step4_drug_query import query_drug_interactions
    result = query_drug_interactions()

    # ── Summary ───────────────────────────────────────────────────────────────
    elapsed = time.time() - start_time
    print("\n\n" + "═" * 64)
    print("  PIPELINE COMPLETE")
    print(f"  Total runtime: {elapsed:.1f}s")
    print("═" * 64)
    print("\n  Output files:")
    for f in [
        "data/alzheimer_genes.csv",
        "data/string_interactions.csv",
        "data/ppi_network.json",
        "data/ppi_network_edges.csv",
        "results/centrality_scores.csv",
        "results/hub_genes.csv",
        "results/drug_gene_interactions.csv",
        "results/drug_candidates.csv",
    ]:
        exists = "✓" if os.path.exists(f) else "✗"
        print(f"    {exists}  {f}")

    print("\n  Next steps:")
    print("    1. Launch dashboard:  streamlit run dashboard.py")
    print("    2. Import ppi_network.json into Cytoscape for visualization")
    print("       (File > Import > Network from File > select JSON)")
    print("\n")


if __name__ == "__main__":
    run_pipeline()
