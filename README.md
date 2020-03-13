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
Genme scale metabolic models (GEM) provide a powerful platform for cell metabolism simulation. With the advancement of Yeast 

## Workflow



## Installation

## Example input and output for test run

## Additionals

## Team member
[Shuang Li](https://github.com/Shuang-Plum)  
[Yue Hu](https://github.com/jechia)  
[Samuel Ramirez](https://github.com/samuramirez)  

