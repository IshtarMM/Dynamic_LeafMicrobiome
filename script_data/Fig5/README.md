# Generating Figure 5: Hub Identification of Networks

To generate **Figure 5** (Hub Identification of the Networks), correlation networks were created and analyzed using **Cytoscape** to extract network properties such as **degree, centrality,** and more, as described in the **Figure 4 README file**.

### Data Files Used:
The same node files from the **Data** folder were used:
- `1.csv` (January)
- `2.csv` (February)
- `3.csv` (March)
- `11.csv` (November)
- `12.csv` (December)

### Hub Identification Process:
1. These files were opened in the script: `scripts/HubIdentification.R`.
2. **Hubs were identified per month** by selecting nodes with the **highest betweenness centrality and closeness centrality (top 5%)**.
3. The script generates the **Fig5a_HubIdentification.pdf** plot.

### Core Microbes as Hubs:
- A second plot is generated in the same script to show in which months **core microbes act as hubs**.
- The resulting plot is stored as: **Fig5b_HubIdentification.pdf**.

