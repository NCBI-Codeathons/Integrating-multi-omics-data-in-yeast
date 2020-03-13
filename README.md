# Integrating-multi-omics-data-in-yeast
# XXXXXX
## Introduction
XXX is a shiny app we built for full-service multi-omics data analysis, including multi-source data integration, 3D pathway visualization, and genome scale metabolic models for pathway enrichment analysis, network visualization, hypothese generation, follow-up exeriment guidance and network inference.

### Data integration
The most prominent problem faced by multi-omics data analysis is how to reconcile and integrate data derived from different sources (different data sets), which entails drastically different data distribution and properties. While independently tackle and analyze each data set is an option, though not ideal and missed the real meaning of multi-omics analysis, the quesiton still remains after stats generated for each data set, how to combine the stats? Furthermore, this way connections between different data sets are inevitablly lost. Therefore, a method to integrate data sets with different underlying distributions should be applied before stats analysis is needed. Canonical Correlation Analysis has been utilized to find common structure existed across different data sets, which would fit perfectly for the purpose of data set combining.The advantage of CCA over other method is that paires of canonical variables can be easily collected to calculate weights for the connection between nodes from different data sources (such as gene vs metabolite and gene vs protein). Projection to Latent Structures models (PLS) with Discriminant Analysis (PLSDA)

### Network visualization
Based on [KEGG pathway](https://www.genome.jp/kegg/pathway.html) and [The CellMap](https://thecellmap.org/), we can get the underlying network for gene x gene interaction, metabolite x metabolite interaction and metabolite x gene interaction. The three layers of netwrok will then be ploted in 3D with each node color denoted the .Different from previous plotting, we add correlation between each node (calculated from CCA) as the weight of the edge, which would render a more straightforward view about how nodes within each layer and nodes across different layers are correlated. This will also sent the base for further network based analysis, such as infer the subnetwork that is significantly perturbed based on the node and edge loadings. This would give us a new way/unit to check and analyze the network, without the need to centralize/project data onto a predefined pathway or GO term.     

### genome scale metabolic models


## Workflow



## Installation

## Example input and output for test run

## Additionals

## Team member
[Shuang Li](https://github.com/Shuang-Plum)

[Yue Hu](https://github.com/jechia)

[Samuel Ramirez](https://github.com/samuramirez)

[Elton Rexhapaj](https://github.com/erexhepa)

[Matthew Berginski](https://github.com/mbergins)


