# Generating Sub Figure 4: Hub Identification of Networks

To generate **Sub Fig 4** (Hub Identification of the Networks), correlation files were created from OTU data. 

### Step 1: Preparing OTU Tables
- **Merge OTU tables** of **Bacteria, Fungi, and Oomycetes**.
- Normalize samples by calculating **OTU relative abundances**, stored in `1_OTU_correlation/*A.txt`.
- The tables were transformed using **log10(x + 1)**.
- Correlation files were computed using Spearman correlation within **CoNet**, a Cytoscape app.
- The resulting correlation networks are stored in: `1_OTU_correlation/MonthlyNetworks_Raw.cys`.

### Step 2: Filtering Correlation Edges
- **Edge tables were filtered** to retain only edges with **p-values < 0.001**.
- Filtering script: `2_0.001_edges/filterPval.R`.
- Filtered edge tables are saved in: `2_0.001_edges/*.csv`.
- Filtered networks are stored in: `2_0.001_edges/month_Networks.cys`.
- These filtered networks were analyzed in Cytoscape to obtain **network properties such as centralities**.
- Output node tables were saved in: `2_0.001_edges/nodetables/*.csv`.

### Hub Identification Process:
1. The filtered node tables were analyzed using the script: `scripts/HubIdentificationSpearmanCorr.R`.
2. **Hubs were identified per month** by selecting nodes with the **highest betweenness centrality and closeness centrality (top 5%)**.
3. The script generates the **Sup_Fig4a_HubIdentification.pdf** plot.

### Core Microbes as Hubs:
- A second plot is generated in the same script to illustrate in which months **core microbes act as hubs** using this input file 2_0.001_edges/nodetables/Hubs_Spearman.txt
- The resulting plot is stored as: **Sup_Fig4b_HubIdentification.pdf**.








