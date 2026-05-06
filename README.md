# Network-Based Identification of Candidate Repurposed Drugs Using Protein Interaction Networks

A computational bioinformatics pipeline for identifying candidate repurposed drugs using protein-protein interaction (PPI) networks, graph theory, and drug-gene interaction analysis.

---

# Project Overview

This project presents a fully automated network-based drug repurposing pipeline capable of analyzing human diseases using biological interaction networks.

The pipeline integrates multiple biological databases and computational tools to:

1. Retrieve disease-associated genes
2. Build protein-protein interaction networks
3. Identify hub genes using graph theory
4. Query drug-gene interaction databases
5. Rank candidate repurposed drugs

The workflow combines:
- DisGeNET
- STRING Database
- NetworkX
- DGIdb
- Flask
- Plotly

The project successfully reproduced known therapeutic landscapes across multiple diseases including:
- Alzheimer’s disease
- Parkinson’s disease
- Breast cancer
- Colorectal cancer
- Type 2 diabetes
- Psoriasis
- Rheumatoid arthritis

---

# Scientific Motivation

Traditional drug discovery:
- Requires 10–15 years
- Costs billions of dollars
- Has high failure rates

Drug repurposing offers a faster and more cost-effective alternative because approved drugs already possess known safety profiles.

Network pharmacology enables identification of biologically important targets through analysis of disease interaction networks instead of isolated genes.

---

# Pipeline Architecture

The computational workflow contains four major stages.

---

## Step 1 — Disease Gene Retrieval

Disease-associated genes are retrieved from DisGeNET using disease-specific UMLS CUIs.

Main script:

```python
step1_get_disease_genes.py

Features:

Dynamic API querying
Confidence filtering
Curated fallback gene sets
Automatic CSV export

Output files:

alzheimer_genes.csv
Step 2 — Protein-Protein Interaction Network Construction

The STRING API is queried to retrieve protein-protein interaction data.

Main script:

step2_build_network.py

Features:

STRING identifier mapping
Confidence-based interaction filtering
NetworkX graph generation
Cytoscape-compatible JSON export

Output files:

string_interactions.csv
ppi_network.json
ppi_network_edges.csv
Step 3 — Network Analysis and Hub Gene Identification

The PPI network is analyzed using graph-theoretic centrality measures.

Main script:

step3_network_analysis.py

Centrality measures:

Degree centrality
Betweenness centrality
Closeness centrality
PageRank
Eigenvector centrality

A composite hub score is calculated to identify biologically important genes.

Output files:

hub_genes.csv
centrality_scores.csv
Step 4 — Drug-Gene Interaction Query

Top hub genes are queried against DGIdb to identify candidate drugs.

Main script:

step4_drug_query.py

Features:

DGIdb GraphQL integration
Drug interaction ranking
Candidate prioritization
Approved drug filtering

Output files:

drug_gene_interactions.csv
drug_candidates.csv
Full Pipeline Execution

Run the complete workflow:

python run_pipeline.py

The pipeline automatically:

Retrieves genes
Builds the interaction network
Performs network analysis
Identifies hub genes
Queries drug interactions
Generates result files
Interactive Dashboard

The project includes an interactive web dashboard built using:

Flask
Plotly.js
HTML/CSS/JavaScript

Main files:

server.py
live_dashboard.html

Dashboard capabilities:

Real-time disease querying
Interactive network visualization
Drug candidate exploration
Hub gene analysis
Dynamic metrics display

Run the dashboard:

python server.py

Then open:

http://localhost:5000
Technologies Used
Programming Language
Python
Bioinformatics Databases
DisGeNET
STRING v12
DGIdb v5
UniProt
PubChem
Python Libraries
pandas
numpy
networkx
requests
plotly
pyvis
dgipy
streamlit

Install dependencies:

pip install -r requirements.txt
Repository Structure
File	Description
run_pipeline.py	Executes the full computational workflow
server.py	Flask backend server
live_dashboard.html	Interactive frontend dashboard
step1_get_disease_genes.py	Disease gene retrieval
step2_build_network.py	STRING network construction
step3_network_analysis.py	Centrality analysis
step4_drug_query.py	Drug interaction querying
requirements.txt	Python dependencies
alzheimer_genes.csv	Retrieved disease genes
string_interactions.csv	STRING interactions
ppi_network.json	Cytoscape network export
ppi_network_edges.csv	Network edge list
centrality_scores.csv	Centrality analysis results
hub_genes.csv	Top-ranked hub genes
drug_gene_interactions.csv	Drug-gene interaction results
drug_candidates.csv	Ranked candidate drugs
Example Workflow

Example disease:

Alzheimer’s Disease

Pipeline process:

Retrieve Alzheimer-associated genes from DisGeNET
Construct the PPI network using STRING
Identify hub genes such as:
APP
APOE
MAPT
TREM2
Query DGIdb for interacting drugs
Rank candidate repurposed drugs
Results

The pipeline successfully reproduced known therapeutic landscapes across multiple diseases.

Examples:

Alzheimer’s disease → Lecanemab
Parkinson’s disease → Prasinezumab
Colorectal cancer → Cetuximab

These results validate the biological relevance of the network-based methodology.

Biological Significance

Hub genes identified through network centrality often correspond to:

Key signaling regulators
Disease progression drivers
Therapeutic intervention points

The project demonstrates how graph theory and systems biology can support:

Precision medicine
Drug repurposing
Computational pharmacology
Future Improvements

Potential future developments:

Machine learning-based drug ranking
Multi-omics integration
Gene expression analysis
Clinical trial integration
AI-assisted target prioritization
Docker deployment
Cloud computing support
Author

Peter El Hajjar

Capstone Bioinformatics Project
Lebanese American University

References
DisGeNET
STRING Database
DGIdb
NetworkX Documentation
UniProt
PubChem
Repository Purpose

This repository demonstrates:

Computational biology workflows
Bioinformatics pipeline engineering
Network pharmacology analysis
Drug repurposing methodology
Biomedical data integration
