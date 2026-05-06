"""
STEP 1 — Fetch Type 2 Diabetes Genes from DisGeNET
"""
import requests
import pandas as pd
import os

DISEASE_NAME  = "Type 2 Diabetes"
AD_UMLS_CUI   = "C0011860"
OUTPUT_FILE   = "data/alzheimer_genes.csv"
MIN_SCORE     = 0.3
TOP_N         = 50

FALLBACK_GENES = [
    "INS", "INSR", "IRS1", "IRS2", "PIK3CA", "AKT2", "GLUT4",
    "SLC2A4", "PPARG", "PPARA", "FOXO1", "FOXO3", "SIRT1",
    "AMPK", "PRKAA1", "PRKAA2", "MTOR", "RAPTOR", "GCK",
    "HNF1A", "HNF4A", "HNF1B", "PDX1", "NEUROD1", "KCNJ11",
    "ABCC8", "TCF7L2", "CDKAL1", "CDKN2A", "CDKN2B", "IGF2BP2",
    "HHEX", "SLC30A8", "JAZF1", "CDC123", "CAMK1D", "TSPAN8",
    "THADA", "ADAMTS9", "NOTCH2", "VEGFA", "IL6", "TNF",
    "ADIPOQ", "LEP", "LEPR", "RETN", "FFA1", "GLP1R",
    "GCGR", "DPP4", "SGLT2", "SLC5A2"
]

def fetch_disgenet_genes():
    print("=" * 60)
    print(f"STEP 1: Fetching {DISEASE_NAME} Genes")
    print("=" * 60)
    os.makedirs("data", exist_ok=True)

    url = f"https://www.disgenet.org/api/gda/disease/{AD_UMLS_CUI}"
    params = {"source": "ALL", "format": "json", "limit": 100, "min_score": MIN_SCORE}
    print(f"\n→ Querying DisGeNET API for CUI: {AD_UMLS_CUI}")
    try:
        resp = requests.get(url, params=params, timeout=20)
        if resp.status_code == 200:
            data = resp.json()
            if isinstance(data, list) and len(data) > 0:
                df = pd.DataFrame(data)
                col_map = {"gene_symbol": "gene_symbol", "score": "gda_score", "disease_name": "disease_name"}
                df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})
                if "gene_symbol" in df.columns and "gda_score" in df.columns:
                    df = df.sort_values("gda_score", ascending=False).head(TOP_N)
                    df["source_method"] = "DisGeNET_API"
                    df.to_csv(OUTPUT_FILE, index=False)
                    print(f"  ✓ Saved {len(df)} genes from DisGeNET API")
                    return df["gene_symbol"].tolist()
    except Exception as e:
        print(f"  ✗ API failed: {e}")

    print(f"\n→ Using curated fallback list ({len(FALLBACK_GENES)} genes)")
    df = pd.DataFrame({
        "gene_symbol": FALLBACK_GENES,
        "gda_score": [1.0] * len(FALLBACK_GENES),
        "disease_name": [DISEASE_NAME] * len(FALLBACK_GENES),
        "source_method": ["Curated_Literature"] * len(FALLBACK_GENES)
    })
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"  ✓ Saved {len(df)} genes to {OUTPUT_FILE}")
    return FALLBACK_GENES

if __name__ == "__main__":
    genes = fetch_disgenet_genes()
    print(f"\n  STEP 1 COMPLETE — {len(genes)} genes ready\n")
