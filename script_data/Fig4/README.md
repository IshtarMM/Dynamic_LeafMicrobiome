# Creating Figure 4 Network and properties

## Step 1: Prepare OTU Tables
1. Merge OTU tables of Bacteria, Fungi, and Oomycetes.
2. Divide the table into five separate tables based on sample months: **November, December, January, February, and March**.
   - Sometimes, month numbers are used in file names (e.g., 11 for November, 12 for December, etc.).
3. For each OTU table, retain only OTUs that appear in at least **five samples** with **more than 10 reads**.
   - *Note: This step is not included in the scripts.*
4. Output files are stored in: `Fig4/Data/1_OTUTables_Network/`
   - Files: `novA.txt`, `decA.txt`, `janA.txt`, `febA.txt`, `marA.txt`

## Step 2: Generate Correlation Files
1. Generate correlation files using the **SparCC algorithm**:
   - Script: `Step1_OTUTables_Network/Script_RunSparcc.sh`
   - *Recommendation:* A faster version of SparCC, called **FastSpar**, is available: [FastSpar](https://github.com/scwatts/FastSpar)
2. After running, files are stored in: `Step2_Raw_Correlation_tables/`
   - File formats:
     - `*_corsparcc.txt` (correlation matrix)
     - `*_pvals.two_sided.txt` (p-values of correlation)
   - Files are named based on month numbers (e.g., 11 for November, 12 for December, etc.).
3. Filter correlations with **p-values < 0.001** using the script:
   - `Step2_Raw_Correlation_tables/SparccPvalfiltering.R`

## Step 3: Clean Correlation Tables
- Outputs are stored in: `Step3_Clean_Correlation_Tables/`
- Example file header:
  ```
  Node1   Node2   Pvalue  Cor
  BacV5_Otu_000003        BacV5_Otu_000002        0       0.557889817476
  ```
  - **Node1**: Name of the first correlated OTU
  - **Node2**: Name of the second correlated OTU
  - **Pvalue**: P-value of correlation
  - **Cor**: Correlation value (positive or negative)

### Generating Figure 4 Panel C: Degree of Nodes
1. Open each correlation table in **Cytoscape**.
2. Navigate to **Tools > Analyze Network (undirected)** to calculate network properties.
3. Export the **node table** and store it in: `Step3_Clean_Correlation_Tables/Fig4C_DegreeOfNodes/*.csv`
4. Use `Fig4_C_DegreeNetwork.R` to plot node degrees as a boxplot (`Fig4C.png`).

## Step 4: DyNet Analysis
- **DyNet** ([DOI: 10.1093/bioinformatics/btw187](https://doi.org/10.1093/bioinformatics/btw187)) integrates time-series molecular interaction data and visualizes network changes over time.
- Compare correlation files of one month with the next (e.g., `11_12`, `12_1`, `1_2`, `2_3`).
- **How-To Guide:** `Step4_DyNetwork/HowtoRunDyNet.docx`

### DyNet Outputs:
- **Cytoscape files**: `Step4_DyNetwork/*.cys`
- **Edge tables**: `Step4_DyNetwork/node_edge_tables_FromDyNet/*edge.csv`
- **Node tables**: `Step4_DyNetwork/node_edge_tables_FromDyNet/*node.csv`

### Generating Figures 4A, 4B, 4D, and 4E
#### **Fig 4B: Number of Nodes and Edges per Month**
1. Open node and edge tables from `Step4_DyNetwork/node_edge_tables_FromDyNet/`.
2. Count the number of nodes and edges.
3. Use `Fig4_B_D.R` to generate the plot.

#### **Fig 4D: Inherited Nodes and Edges Between Months**
- **Formula for inheritance calculation:**
  ```
  n1 = Number of nodes in Network1
  n2 = Number of nodes in Network2
  overlap = Number of nodes common in both networks
  inherited nodes = overlap / n2
  ```
- Similar calculation applies for edges.
- Script: `Step4_DyNetwork/node_edge_tables_FromDyNet/Fig4_B_D.R`

#### **Fig 4A: Connected Network Visualization**
1. Merge edges from all networks.
2. Modify node names by prefixing them with the corresponding month for better visualization.
3. Include edges showing connections between nodes present in consecutive networks.
4. Generate final edge table using `Fig4A_Node_Edgetable.R`.
5. Store the output in: `Fig4A_Network/edge.csv`.
6. Open `edge.csv` in **Cytoscape** and visualize as `Co-abundaceNetwork_Fig4a.cys`.
7. Import `nodeattribute.csv` for enhanced visualization (e.g., highlighting core microbes).
To create monthly clusters of nodes in Cytoscape, I first selected the network by month (e.g., searching for “nov” to select nodes from November).
For each month’s subset of nodes, I applied the Compound Spring Embedded layout to arrange them locally.
I repeated this process for each month, resulting in distinct, visually organized clusters of nodes grouped by month.
8. Export final figure: `Co-abundaceNetwork_Fig4a.png`.

---

This document outlines the steps required to generate **Figure 4**, including data preparation, correlation analysis, network visualization, and final plotting.

