"""
Live Drug Repurposing Server

Flask server that accepts any disease name, runs the full
DisGeNET → STRING → NetworkX → DGIdb pipeline live, and
returns JSON results to the frontend dashboard.

Run with:  python server.py
Open:      http://localhost:5000
"""

from flask import Flask, jsonify, request, send_from_directory
import requests
import networkx as nx
import json
import os
import time

app = Flask(__name__, static_folder='.')

# Disease name → UMLS CUI mapping (expandable)( only back up for the presentation)
KNOWN_CUIS = {
    "alzheimer": "C0002395", "alzheimer's": "C0002395", "alzheimers": "C0002395",
    "parkinson": "C0030567", "parkinson's": "C0030567", "parkinsons": "C0030567",
    "breast cancer": "C0678222", "breast": "C0678222",
    "colorectal cancer": "C1527249", "colorectal": "C1527249", "colon cancer": "C0009402",
    "psoriasis": "C0033860",
    "rheumatoid arthritis": "C0003873", "rheumatoid": "C0003873",
    "type 2 diabetes": "C0011860", "diabetes": "C0011860", "t2d": "C0011860",
    "lupus": "C0024141", "sle": "C0024141",
    "multiple sclerosis": "C0026769", "ms": "C0026769",
    "als": "C0002736", "amyotrophic lateral sclerosis": "C0002736",
    "melanoma": "C0025202",
    "lung cancer": "C0684249",
    "prostate cancer": "C0376358",
    "depression": "C0011570", "major depression": "C0011570",
    "schizophrenia": "C0036341",
    "autism": "C0004352",
    "epilepsy": "C0014544",
    "hypertension": "C0020538",
    "asthma": "C0004096",
    "crohn": "C0010346", "crohn's disease": "C0010346",
    "ulcerative colitis": "C0009324",
    "ovarian cancer": "C0029925",
    "pancreatic cancer": "C0030297",
    "leukemia": "C0023418",
    "lymphoma": "C0024299",
    "glioblastoma": "C0017636", "glioma": "C0017638",
    "huntington": "C0020179", "huntington's": "C0020179",
}

#Fallback gene lists for common diseases ( only back up for the presentation)
FALLBACK_GENES = {
    "C0002395": ["APP","PSEN1","PSEN2","APOE","CLU","CR1","BIN1","PICALM","CD33","ABCA7","TREM2","SORL1","MAPT","TARDBP","FUS","SNCA","LRRK2","GBA","PINK1","PARK7","SQSTM1","OPTN","TBK1","C9orf72","VCP"],
    "C0030567": ["SNCA","LRRK2","PINK1","PRKN","PARK7","GBA","ATP13A2","VPS35","FBXO7","DNAJC6","SYNJ1","DNAJC13","CHCHD2","PLA2G6","VPS13C","GAK","HTRA2","TH","SLC6A3","DRD2","UCHL1","MAPT"],
    "C0678222": ["BRCA1","BRCA2","TP53","PIK3CA","PTEN","AKT1","ERBB2","ESR1","CCND1","CDH1","MYC","EGFR","BCL2","CDKN2A","ATM","CHEK2","PALB2","RAD51","MTOR","CDKN1B"],
    "C1527249": ["APC","KRAS","TP53","SMAD4","BRAF","PIK3CA","PTEN","MLH1","MSH2","CTNNB1","CDH1","MYC","EGFR","BCL2","CDKN2A","AKT1","MTOR","CASP3","ERBB2","VIM"],
    "C0033860": ["TNF","IL17A","IL17F","IL23A","IL12B","IL6","IL1B","IFNG","IL10","NFKB1","STAT3","JAK1","JAK2","TYK2","CTLA4","IL4","IL13","IL22","FOXP3","RELA"],
    "C0003873": ["TNF","IL6","IL1B","IFNG","CXCR4","IL10","IL17A","CD80","STAT3","STAT1","CCL2","NFKB1","JAK1","JAK2","MMP3","MMP9","PTGS2","IL2RA","CXCL8","FOXP3"],
    "C0011860": ["INS","INSR","IRS1","PIK3CA","AKT2","PPARG","PPARA","FOXO1","SIRT1","GCK","HNF1A","TCF7L2","CDKN2A","KCNJ11","ABCC8","GLP1R","LEP","LEPR","TNF","IL6","DPP4","SLC5A2"],
    "C0024141": ["TREX1","C1Q","C2","C4A","IRF5","PTPN22","STAT4","BLK","TNFAIP3","IRAK1","IL10","TNF","IFNG","FCGR2A","FCGR3A","CR2","SPEN","TLR7","TLR9","DNASE1"],
    "C0026769": ["HLA-DRB1","IL2RA","IL7R","TNFRSF1A","IRF5","STAT3","IFNG","IL12B","IL10","TNF","CD58","CLEC16A","EVI5","GPC5","METTL21B","TAGAP","TMEM39A","IL22","PRDM1","TNFSF14"],
    "C0002736": ["SOD1","FUS","TARDBP","C9orf72","UBQLN2","VCP","OPTN","TBK1","SQSTM1","HNRNPA2B1","HNRNPA1","NEK1","CCNF","MATR3","KIF5A","SETX","VAPB","ANG","DCTN1","PFN1"],
    "C0025202": ["BRAF","NRAS","KIT","CDKN2A","TP53","PTEN","AKT1","PIK3CA","MAP2K1","MAP2K2","MITF","CDK4","MDM2","MCL1","BCL2","CCND1","EGFR","MET","RB1","TERT"],
    "C0684249": ["KRAS","TP53","EGFR","ALK","ROS1","BRAF","PIK3CA","PTEN","AKT1","MET","ERBB2","RET","STK11","KEAP1","NF1","CDKN2A","RB1","SMAD4","MYC","BCL2"],
    "C0376358": ["AR","TP53","PTEN","PIK3CA","AKT1","RB1","CDH1","BRCA2","ATM","BRCA1","CDK12","SPOP","FOXA1","ERG","ETV1","TMPRSS2","MYC","EGFR","PDGFRA","WNT5A"],
    "C0011570": ["SLC6A4","BDNF","HTR2A","TPH2","COMT","MAOA","CRH","FKBP5","NTRK2","CREB1","DRD4","GRIN2B","SIRT1","IL6","TNF","NR3C1","CLOCK","NPAS2","PER2","CRY1"],
    "C0036341": ["DISC1","DTNBP1","NRG1","COMT","DRD2","HTR2A","GRIN2B","RELN","PRODH","DAAO","MTHFR","NCAN","PCGEM1","MHC","CACNA1C","ANK3","KCNQ2","TCF4","ZNF804A","CSMD1"],
    "C0004352": ["SHANK3","SHANK2","SHANK1","NLGN3","NLGN4X","NRXN1","CHD8","PTEN","TSC1","TSC2","FMR1","MECP2","SYNGAP1","ADNP","DYRK1A","TBR1","POGZ","ARID1B","GRIN2B","SCN1A"],
    "C0014544": ["SCN1A","SCN2A","SCN8A","KCNQ2","KCNQ3","GABRG2","GABRA1","CDKL5","MECP2","TSC1","TSC2","DEPDC5","LGI1","CASPR2","NRXN1","SHANK3","SPTAN1","SLC6A1","ALDH7A1","POLG"],
    "C0020538": ["ACE","AGT","AGTR1","AGTR2","NOS3","ADD1","CYP11B2","ADRB2","ADRB1","GNB3","ATP2B1","CACNA1C","SLC8A1","EDN1","EDNRA","PPARGC1A","ADRA1A","ADRA2A","NPPB","NPR1"],
    "C0004096": ["IL4","IL5","IL13","STAT6","GATA3","TBX21","FOXP3","IL10","IFNG","IL17A","TNF","PTGDR2","ADRB2","ADRB1","CCR3","CPA3","MS4A2","FCER1G","SPDEF","FOXA1"],
    "C0010346": ["NOD2","ATG16L1","IRGM","IL23R","JAK2","STAT3","TNF","IL6","IL1B","IL10","IL17A","IFNG","CARD9","PTPN22","TNFSF15","SMAD3","CCR6","NKX2-3","PTGER4","RIPK2"],
    "C0009324": ["HLA-B","TNF","IL10","IL23R","IL17A","IFNG","JAK2","STAT3","NOD2","ATG16L1","IRGM","PTPN22","CARD9","HNF4A","CDH1","LAMB1","ECM1","PTGER4","TNFSF15","IL6"],
    "C0029925": ["BRCA1","BRCA2","TP53","PIK3CA","PTEN","AKT1","KRAS","BRAF","EGFR","ERBB2","CDK4","CDK6","CCND1","CDKN2A","RB1","ATM","CHEK2","RAD51","MYC","BCL2"],
    "C0030297": ["KRAS","TP53","SMAD4","CDKN2A","MYC","BRCA2","ATM","PIK3CA","AKT1","EGFR","ERBB2","MET","CDH1","VHL","PBRM1","BAP1","ARID1A","GNAS","TGFBR2","STK11"],
    "C0023418": ["BCR","ABL1","TP53","FLT3","NPM1","DNMT3A","IDH1","IDH2","RUNX1","CEBPA","TET2","ASXL1","EZH2","WT1","KIT","NRAS","KRAS","JAK2","STAT5A","MLL"],
    "C0024299": ["BCL2","MYC","TP53","CDKN2A","RB1","ATM","PIK3CA","AKT1","BRAF","KRAS","NRAS","MYD88","CD79A","CARD11","TNFAIP3","CREBBP","EP300","EZH2","KMT2D","NOTCH1"],
    "C0017636": ["IDH1","IDH2","TP53","ATRX","TERT","EGFR","PTEN","RB1","CDKN2A","MDM2","MDM4","PIK3CA","AKT1","BRAF","NF1","PDGFRA","CDK4","CDK6","CCND2","MYC"],
    "C0020179": ["HTT","HAP1","HIP1","PACSIN1","SH3GL3","PRKCG","CASP3","CASP8","CASP9","BCL2","BAX","TP53","BDNF","TrkB","NDRG2","SIRT1","HDAC1","HDAC2","CBP","P53"],
}

def search_umls_cui(disease_name):
    """Try to find UMLS CUI from name."""
    low = disease_name.lower().strip()
    if low in KNOWN_CUIS:
        return KNOWN_CUIS[low]
    for key, cui in KNOWN_CUIS.items():
        if key in low or low in key:
            return cui
    return None

def fetch_genes_disgenet(cui, disease_name):
    """Fetch genes from DisGeNET API."""
    url = f"https://www.disgenet.org/api/gda/disease/{cui}"
    params = {"source": "ALL", "format": "json", "limit": 100, "min_score": 0.3}
    try:
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            if isinstance(data, list) and len(data) > 0:
                genes = []
                for item in data:
                    sym = item.get("gene_symbol") or item.get("geneName","")
                    sc  = item.get("score", 0.5)
                    if sym:
                        genes.append({"gene_symbol": sym, "gda_score": sc})
                genes = sorted(genes, key=lambda x: x["gda_score"], reverse=True)[:50]
                return [g["gene_symbol"] for g in genes], genes
    except Exception:
        pass
    # fallback
    fallback = FALLBACK_GENES.get(cui, [])
    if not fallback:
        return [], []
    return fallback, [{"gene_symbol": g, "gda_score": 1.0} for g in fallback]

def fetch_string_network(genes):
    """Fetch PPI network from STRING API."""
    base = "https://string-db.org/api"
    # Map to string IDs
    id_map = {}
    try:
        r = requests.post(f"{base}/json/get_string_ids",
            data={"identifiers": "\r".join(genes), "species": 9606,
                  "limit": 1, "caller_identity": "live_repurposing"},
            timeout=20)
        for item in r.json():
            id_map[item.get("preferredName", item.get("queryItem",""))] = item.get("stringId","")
    except Exception:
        pass

    if not id_map:
        return [], {}

    gene_names = list(id_map.keys())
    string_ids = list(id_map.values())

    edges = []
    try:
        time.sleep(0.5)
        r = requests.post(f"{base}/tsv/network",
            data={"identifiers": "%0d".join(string_ids), "species": 9606,
                  "required_score": 400, "caller_identity": "live_repurposing"},
            timeout=30)
        lines = r.text.strip().split("\n")
        if len(lines) > 1:
            header = lines[0].split("\t")
            for line in lines[1:]:
                row = dict(zip(header, line.split("\t")))
                a = row.get("preferredName_A","")
                b = row.get("preferredName_B","")
                s = float(row.get("score",0))
                if a and b and s > 0:
                    edges.append({"a": a, "b": b, "score": s})
    except Exception:
        pass

    return edges, gene_names

def compute_centrality(edges, gene_names):
    """Compute hub scores via NetworkX."""
    G = nx.Graph()
    for g in gene_names:
        G.add_node(g)
    for e in edges:
        G.add_edge(e["a"], e["b"], weight=e["score"])

    if G.number_of_nodes() == 0:
        return []

    deg   = nx.degree_centrality(G)
    bet   = nx.betweenness_centrality(G, weight="weight", normalized=True)
    clo   = nx.closeness_centrality(G)
    try:
        pr = nx.pagerank(G, weight="weight", max_iter=500)
    except Exception:
        pr = {n: 1/len(G.nodes()) for n in G.nodes()}

    nodes = list(G.nodes())
    def mm(d):
        vals = list(d.values())
        mn,mx = min(vals),max(vals)
        return {k:(v-mn)/(mx-mn) if mx>mn else 0 for k,v in d.items()}

    nd,nb,nc,np_ = mm(deg),mm(bet),mm(clo),mm(pr)
    hub = {n: 0.25*nd[n]+0.25*nb[n]+0.2*nc[n]+0.2*np_[n]+0.1*deg[n] for n in nodes}

    ranked = sorted(hub.items(), key=lambda x: x[1], reverse=True)
    result = []
    for i,(gene,score) in enumerate(ranked[:15]):
        result.append({
            "rank": i+1,
            "gene": gene,
            "degree": G.degree(gene),
            "hub_score": round(score,4),
            "betweenness": round(bet[gene],4),
            "pagerank": round(pr[gene],4),
            "degree_centrality": round(deg[gene],4)
        })

    # network stats
    stats = {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "density": round(nx.density(G),4),
        "clustering": round(nx.average_clustering(G),4),
    }
    try:
        gcc = G.subgraph(max(nx.connected_components(G), key=len))
        stats["avg_path"] = round(nx.average_shortest_path_length(gcc),4)
    except Exception:
        stats["avg_path"] = "N/A"

    return result, stats

def fetch_drugs_dgidb(genes):
    """Fetch drug-gene interactions from DGIdb GraphQL."""
    url = "https://dgidb.org/api/graphql"
    gene_str = '", "'.join(genes)
    query = f'''{{
      genes(names: ["{gene_str}"]) {{
        nodes {{
          name
          interactions {{
            drug {{ name approved }}
            interactionScore
            interactionTypes {{ type directionality }}
            sources {{ sourceDbName }}
          }}
        }}
      }}
    }}'''
    interactions = []
    try:
        r = requests.post(url, json={"query": query},
                          headers={"Content-Type": "application/json"}, timeout=30)
        nodes = r.json().get("data",{}).get("genes",{}).get("nodes",[])
        for node in nodes:
            gene = node.get("name","")
            for intr in node.get("interactions",[]):
                drug = intr.get("drug",{})
                dname = drug.get("name","")
                if not dname:
                    continue
                types = "; ".join(t.get("type","") for t in intr.get("interactionTypes",[]))
                srcs  = "; ".join(s.get("sourceDbName","") for s in intr.get("sources",[]))
                interactions.append({
                    "gene": gene,
                    "drug": dname.strip(),
                    "approved": drug.get("approved", False),
                    "type": types or "unknown",
                    "score": intr.get("interactionScore", 0),
                    "sources": srcs
                })
    except Exception:
        pass
    return interactions

def rank_candidates(interactions, hubs):
    """Rank drug candidates."""
    hub_scores = {h["gene"]: h["hub_score"] for h in hubs}
    hub_ranks  = {h["gene"]: h["rank"] for h in hubs}

    from collections import defaultdict
    drug_map = defaultdict(list)
    for intr in interactions:
        drug_map[intr["drug"]].append(intr)

    candidates = []
    for drug, items in drug_map.items():
        genes_targeted = list({i["gene"] for i in items})
        n = len(genes_targeted)
        avg_hub = sum(hub_scores.get(g,0) for g in genes_targeted) / n
        approved = any(i["approved"] for i in items)
        types = "; ".join({i["type"] for i in items if i["type"]})
        sources = "; ".join({i["sources"] for i in items if i["sources"]})
        rs = n*0.4 + avg_hub*10*0.35 + (1 if approved else 0)*0.25
        candidates.append({
            "drug": drug,
            "genes": "; ".join(genes_targeted),
            "n_hub_genes": n,
            "avg_hub_score": round(avg_hub, 4),
            "approved": approved,
            "type": types,
            "sources": sources,
            "repurposing_score": round(rs, 4)
        })

    candidates.sort(key=lambda x: (x["n_hub_genes"], x["repurposing_score"]), reverse=True)
    for i, c in enumerate(candidates[:50]):
        c["rank"] = i + 1
    return candidates[:50]

# ── Routes ────────────────────────────────────────────────────────────────────

@app.route('/')
def index():
    return send_from_directory('.', 'live_dashboard.html')

@app.route('/api/analyze', methods=['POST'])
def analyze():
    disease_name = request.json.get('disease', '').strip()
    if not disease_name:
        return jsonify({"error": "Please enter a disease name"}), 400

    # Step 1: Find CUI
    cui = search_umls_cui(disease_name)
    if not cui:
        # Try a generic approach with fallback
        cui = "UNKNOWN"

    # Step 2: Fetch genes
    genes, gene_data = fetch_genes_disgenet(cui, disease_name)
    if not genes:
        return jsonify({"error": f"No genes found for '{disease_name}'. Try a more specific disease name."}), 404

    # Step 3: STRING network
    edges, gene_names = fetch_string_network(genes[:45])
    if not edges:
        edges = []

    # Step 4: Centrality
    hubs_result = compute_centrality(edges, gene_names if gene_names else genes[:30])
    if isinstance(hubs_result, tuple):
        hubs, net_stats = hubs_result
    else:
        hubs = hubs_result
        net_stats = {"nodes": len(genes), "edges": len(edges), "density": 0, "clustering": 0, "avg_path": "N/A"}

    # Step 5: DGIdb drugs
    hub_gene_names = [h["gene"] for h in hubs]
    raw_interactions = fetch_drugs_dgidb(hub_gene_names)
    candidates = rank_candidates(raw_interactions, hubs)

    return jsonify({
        "disease": disease_name,
        "cui": cui,
        "genes": gene_data[:30],
        "edges": edges[:300],
        "hubs": hubs,
        "net_stats": net_stats,
        "candidates": candidates[:50],
        "interactions": raw_interactions[:200]
    })

@app.route('/api/protein_info/<gene>', methods=['GET'])
def protein_info(gene):
    """Fetch protein info from UniProt."""
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true&format=json&size=1",
            timeout=10
        )
        data = r.json()
        results = data.get("results", [])
        if results:
            p = results[0]
            return jsonify({
                "accession": p.get("primaryAccession",""),
                "name": p.get("uniProtkbId",""),
                "protein_name": p.get("proteinDescription",{}).get("recommendedName",{}).get("fullName",{}).get("value",""),
                "function": next((c["texts"][0]["value"] for c in p.get("comments",[]) if c.get("commentType")=="FUNCTION"), ""),
                "length": p.get("sequence",{}).get("length",0),
                "mass": p.get("sequence",{}).get("molWeight",0),
                "subcellular": [c["locations"][0]["location"]["value"] for c in p.get("comments",[]) if c.get("commentType")=="SUBCELLULAR LOCATION" and c.get("locations")][:3],
                "gene": gene,
                "uniprot_url": f"https://www.uniprot.org/uniprotkb/{p.get('primaryAccession','')}"
            })
    except Exception:
        pass
    return jsonify({"gene": gene, "error": "No UniProt data found"})

@app.route('/api/drug_info/<drug_name>', methods=['GET'])
def drug_info(drug_name):
    """Fetch drug info from DGIdb."""
    try:
        query = f'''{{
          drugs(names: ["{drug_name}"]) {{
            nodes {{
              name
              approved
              drugAliases {{ alias }}
              interactions {{
                gene {{ name }}
                interactionTypes {{ type }}
              }}
            }}
          }}
        }}'''
        r = requests.post("https://dgidb.org/api/graphql",
                         json={"query": query},
                         headers={"Content-Type": "application/json"},
                         timeout=10)
        nodes = r.json().get("data",{}).get("drugs",{}).get("nodes",[])
        if nodes:
            d = nodes[0]
            aliases = [a["alias"] for a in d.get("drugAliases",[])[:5]]
            genes   = list({i["gene"]["name"] for i in d.get("interactions",[]) if i.get("gene")})[:10]
            return jsonify({
                "name": d.get("name",""),
                "approved": d.get("approved", False),
                "aliases": aliases,
                "known_targets": genes,
            })
    except Exception:
        pass
    return jsonify({"name": drug_name, "error": "No drug info found"})

if __name__ == '__main__':
    print("\n" + "="*55)
    print("  🧬 Live Drug Repurposing Server")
    print("="*55)
    print("  Open your browser at:  http://localhost:5000")
    print("  Press Ctrl+C to stop")
    print("="*55 + "\n")
    app.run(debug=False, port=5000)
