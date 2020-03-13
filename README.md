# Integrating-multi-omics-data-in-yeast
# Multi-Omics Multiplex Analysis (MOMA)
## Whats Whys and Hows
Multi-Omics Multiplex Analysis (MOMA) is a analysis software we built for full-service multi-omics data analysis, including multi-source data integration, 3D pathway visualization,(subnetwork inference), and genome scale metabolic models for fluxomics simulation, network visualization, hypothese generation, follow-up exeriment guidance and network inference.

### Data integration
The most prominent problem faced by multi-omics data analysis is how to reconcile and integrate data derived from different sources (different data sets), which entails drastically different data distribution and properties. While independently tackle and analyze each data set is an option, though not ideal and missed the real meaning of multi-omics analysis, the quesiton still remains after stats generated for each data set, how to combine the stats? Furthermore, this way connections between different data sets are inevitablly lost. Therefore, a method to integrate data sets with different underlying distributions should be applied before stats analysis is needed.  
**Projection to Latent Structures models (PLS) with Discriminant Analysis (PLSDA)** is used to calculate the loadings (the contribution in driving the variation of the data) of each gene and metabolite, which is then used to determine whether this gene/metabolite bear significant change between samples.  
**Canonical Correlation Analysis** is used to calculate the weights for the connection between nodes from different data sources (such as gene vs metabolite and gene vs protein). These weights are later used in the network visualization module (edges) and subnetwork inference module.  

### Network visualization
Based on [KEGG pathway](https://www.genome.jp/kegg/pathway.html) and [The CellMap](https://thecellmap.org/), we can get the underlying network for gene x gene interaction, metabolite x metabolite interaction and metabolite x gene interaction. The three layers of netwrok will then be ploted in 3D in [Gephi](https://gephi.org/), with each node size denotes the PLSDA calculated loadings (contribution to separation of data between samples). Different from previous plotting, we add correlation between each node (calculated from CCA) as the weight of the edge, which would render a more straightforward view about how nodes within each layer and nodes across different layers are correlated. This sets the fundation for further network based analysis, such as infer the subnetwork unit that is significantly perturbed based on the node and edge loadings. In this way, a new way/unit to check and analyze the network can be generated, which is completely independent from the predefined pathway or GO term.     

### Genome scale metabolic models
Genme scale metabolic models (GEM) provide a powerful platform for cell metabolism simulation. With the extensive knowledge about yeast, the consesus GEMs have been well maintained and updated, such as the [Yeast8](https://www.nature.com/articles/s41467-019-11581-3). The 3000s well annotated and manually curated reactions can be efficiently solved with the constraint-based reconstruction and analysis (COBRA) toolbox based on Matlab. Yeast8 and COBRA offers a platform to simulate yeast under various envirnmental and genetical perturbations. The output fluxomics provides an extra layer of knowledge, which is essential in guiding reaction selection and hypothese formation for the follow-up experiment design.

## Workflow
Data from different omics are firstly organized to have each node (gene or protein or metabolite) associated with pathway information from KEGG database. PLS-DA and CCA analysis are then perfromed to analyze the contribution of each node toward the variation of the data and also the correlation between omics data sources. 3D plots is then performed with each omics data plotted on a separate layer with inner layer and inter layer connection drawn from KEGG and TheCellMap. Node color indicates its data source, node size indicates its power in driving sample variation, and edge width indicates correlation strength between two nodes. Correlation matrix produced by CCA could also be used to infer new edge if needed. On the other side, constraints can be set based on experimental settings and sample preparation for COBRA to simulate steady state Yeast8 models and produce fluxomics data. Reactions bear drastic flux changes can then work as the guidance for follow-up experiment design.  
![alt text](https://github.com/NCBI-Codeathons/Integrating-multi-omics-data-in-yeast/blob/master/plots/workflow.png)


## Installation
Please contact the authors for installation guidance for now. We are not trying to be mysterious and we are hard-working people. But we do want MOMA's debut experience be smooth and happy.  

## Example input and output for test run
### Example Input Data:  
[RNA-seq](https://github.com/NCBI-Codeathons/Integrating-multi-omics-data-in-yeast/blob/master/RNAseq_wtHL/seqdata_vstnorm.csv)  
[Metabolimics](https://github.com/NCBI-Codeathons/Integrating-multi-omics-data-in-yeast/blob/master/Metab_wtHL/wtHL_cpd_path_stats_data_mancur.csv)

### Example Output Plots:  
#### Multiomics Network  
![alt text](https://github.com/NCBI-Codeathons/Integrating-multi-omics-data-in-yeast/blob/master/plots/Network_sucrose.png)  
#### Fluxomics simulation  
![alt text](https://github.com/NCBI-Codeathons/Integrating-multi-omics-data-in-yeast/blob/master/plots/fluxplot_3.png)

## Additionals
Subnetwork inference based on the network model 
New edge inference based on the fluxomics simulation

## Team member
[Shuang Li](https://github.com/Shuang-Plum)  
[Yue Hu](https://github.com/jechia)  
[Samuel Ramirez](https://github.com/samuramirez)  

