# Integrating-multi-omics-data-in-yeast
# XXXXXX
## Introduction
XXX is a shiny app we built for full-service multi-omics data analysis, including multi-source data integration, 3D pathway visualization, and genome scale metabolic models for pathway enrichment analysis, network visualization, hypothese generation, follow-up exeriment guidance and network inference.

### Data integration
The most prominent problem faced by multi-omics data analysis is how to reconcile and integrate data derived from different sources (different data sets), which entails drastically different data distribution and properties. While independently tackle and analyze each data set is an option, though not ideal and missed the real meaning of multi-omics analysis, the quesiton still remains after stats generated for each data set, how to combine the stats? Furthermore, this way connections between different data sets are inevitablly lost. Therefore, a method to integrate data sets with different underlying distributions should be applied before stats analysis is needed. Canonical Correlation Analysis has been utilized to find common structure existed across different data sets, which would fit perfectly for the purpose of data set combining.(regularized CCA for two data sets and generalized CCA for more than two data sets.) The advantage of CCA over other method (such as PLS) is that paires of canonical variables can be easily collected to calculate weights for the connection between nodes? from different data sources (such as gene vs metabolite and gene vs protein). 

### Network visualization


### genome scale metabolic models


## Workflow



## Installation

## Example input and output for test run

## Additionals

## Team member
[Shuang Li](https://github.com/Shuang-Plum)

[Samuel Ramirez](https://github.com/samuramirez)

[Matthew Berginski](https://github.com/mbergins)

[Yue Hu](https://github.com/jechia)
