"""
STEP 3 — Network Analysis: Identify Hub Genes
===============================================
Loads the PPI network and computes multiple centrality measures
to rank and identify the most important (hub) genes.

Centrality measures used:
  - Degree centrality        : direct connections
  - Betweenness centrality   : bridge genes in network
  - Closeness centrality     : proximity to all other nodes
  - PageRank                 : importance by neighbor quality (like Google)
  - Eigenvector centrality   : connected to other important nodes

Output:
  results/hub_genes.csv          — ranked hub gene table
  results/centrality_scores.csv  — full centrality matrix
"""

import pandas as pd
import networkx as nx
import json
import os

GRAPH_JSON_FILE   = "data/ppi_network.json"
EDGES_CSV         = "data/ppi_network_edges.csv"
OUTPUT_HUB_GENES  = "results/hub_genes.csv"
OUTPUT_CENTRALITY = "results/centrality_scores.csv"
TOP_HUB_N         = 15   # Number of top hub genes to export


def load_network() -> nx.Graph:
    """Load the PPI network from the edge list CSV."""
    print("\n→ Loading PPI network...")
    try:
        df = pd.read_csv(EDGES_CSV)
        G  = nx.Graph()
        for _, row in df.iterrows():
            G.add_edge(row["gene_a"], row["gene_b"],
                       weight=float(row.get("score", 0.5)))
        print(f"  ✓ Network loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        return G
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Edge file not found: {EDGES_CSV}\n"
            "Please run step2_build_network.py first."
        )


def compute_centrality(G: nx.Graph) -> pd.DataFrame:
    """Compute multiple centrality measures for all nodes."""
    print("\n→ Computing centrality measures...")

    nodes = list(G.nodes())

    # Degree centrality (fraction of nodes connected to)
    degree_c = nx.degree_centrality(G)
    print("  ✓ Degree centrality")

    # Betweenness centrality (how often a node is on shortest path)
    between_c = nx.betweenness_centrality(G, weight="weight", normalized=True)
    print("  ✓ Betweenness centrality")

    # Closeness centrality (average distance to all nodes)
    close_c = nx.closeness_centrality(G)
    print("  ✓ Closeness centrality")

    # PageRank (recursive importance) — pure python, no scipy
    try:
        pagerank = nx.pagerank(G, weight="weight", max_iter=500, dangling=None)
        print("  ✓ PageRank")
    except Exception:
        pagerank = {n: 1/len(nodes) for n in nodes}
        print("  ✓ PageRank (uniform fallback)")

    # Eigenvector centrality
    try:
        eigen_c = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)
        print("  ✓ Eigenvector centrality")
    except (nx.PowerIterationFailedConvergence, nx.NetworkXError):
        eigen_c = {n: 0.0 for n in nodes}
        print("  ✓ Eigenvector centrality (fallback to 0)")

    # Assemble dataframe
    df = pd.DataFrame({
        "gene":                  nodes,
        "degree":                [G.degree(n) for n in nodes],
        "degree_centrality":     [degree_c.get(n, 0) for n in nodes],
        "betweenness_centrality":[between_c.get(n, 0) for n in nodes],
        "closeness_centrality":  [close_c.get(n, 0) for n in nodes],
        "pagerank":              [pagerank.get(n, 0) for n in nodes],
        "eigenvector_centrality":[eigen_c.get(n, 0) for n in nodes],
    })

    # Composite hub score (weighted average of all measures)
    # Normalize each column to [0,1] then take weighted mean
    def minmax(col):
        mn, mx = col.min(), col.max()
        return (col - mn) / (mx - mn) if mx > mn else col * 0

    df["hub_score"] = (
        0.25 * minmax(df["degree_centrality"]) +
        0.25 * minmax(df["betweenness_centrality"]) +
        0.20 * minmax(df["closeness_centrality"]) +
        0.20 * minmax(df["pagerank"]) +
        0.10 * minmax(df["eigenvector_centrality"])
    )

    df = df.sort_values("hub_score", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))

    return df


def analyze_network():
    print("=" * 60)
    print("STEP 3: Network Analysis — Hub Gene Identification")
    print("=" * 60)

    os.makedirs("results", exist_ok=True)

    G  = load_network()
    df = compute_centrality(G)

    # Save full centrality table
    df.to_csv(OUTPUT_CENTRALITY, index=False)
    print(f"\n  ✓ Full centrality table saved: {OUTPUT_CENTRALITY}")

    # Save top hub genes
    hub_df = df.head(TOP_HUB_N)[["rank", "gene", "degree", "hub_score",
                                   "degree_centrality", "betweenness_centrality",
                                   "pagerank"]]
    hub_df.to_csv(OUTPUT_HUB_GENES, index=False)
    print(f"  ✓ Top {TOP_HUB_N} hub genes saved: {OUTPUT_HUB_GENES}")

    # Print summary table
    print(f"\n  {'─'*65}")
    print(f"  {'Rank':<6} {'Gene':<14} {'Degree':<8} {'Betweenness':<14} {'PageRank':<12} {'Hub Score'}")
    print(f"  {'─'*65}")
    for _, row in hub_df.iterrows():
        print(f"  {int(row['rank']):<6} {row['gene']:<14} {int(row['degree']):<8} "
              f"{row['betweenness_centrality']:.4f}        "
              f"{row['pagerank']:.4f}       "
              f"{row['hub_score']:.4f}")
    print(f"  {'─'*65}")

    # Network statistics
    print(f"\n  Network Statistics:")
    print(f"    Average degree    : {sum(dict(G.degree()).values())/G.number_of_nodes():.2f}")
    print(f"    Network density   : {nx.density(G):.4f}")
    try:
        gcc = G.subgraph(max(nx.connected_components(G), key=len))
        print(f"    Clustering coeff  : {nx.average_clustering(G):.4f}")
        print(f"    Avg path length   : {nx.average_shortest_path_length(gcc):.4f}")
    except Exception:
        print(f"    Clustering coeff  : {nx.average_clustering(G):.4f}")

    return hub_df, G


if __name__ == "__main__":
    hub_df, G = analyze_network()
    print(f"\n{'='*60}")
    print(f"  STEP 3 COMPLETE — Top hub genes identified")
    print(f"{'='*60}\n")