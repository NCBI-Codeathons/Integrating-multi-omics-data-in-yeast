# yeast-GEM: The consensus genome-scale metabolic model of _Saccharomyces cerevisiae_

[![DOI](https://zenodo.org/badge/52777598.svg)](https://zenodo.org/badge/latestdoi/52777598) [![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem.svg)](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem) [![Join the chat at https://gitter.im/SysBioChalmers/yeast-GEM](https://badges.gitter.im/SysBioChalmers/yeast-GEM.svg)](https://gitter.im/SysBioChalmers/yeast-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

* Brief Model Description:

This repository contains the current consensus genome-scale metabolic model of _Saccharomyces cerevisiae_. It is the continuation of the legacy project [yeastnet](https://sourceforge.net/projects/yeast/). For the latest release please [click here](https://github.com/SysBioChalmers/yeast-GEM/releases).

* Model KeyWords:

**GEM Category:** species; **Utilisation:** experimental data reconstruction, multi-omics integrative analysis, _in silico_ strain design, model template; **Field:** metabolic-network reconstruction; **Type of Model:** reconstruction, curated; **Model Source:** YeastMetabolicNetwork; **Omic Source:** genomics, metabolomics; **Taxonomy:** _Saccharomyces cerevisiae_; **Metabolic System:** general metabolism; **Bioreactor**; **Strain:** S288C; **Condition:** aerobic, glucose-limited, defined media;

* Last update: 2019-05-21

* Main Model Descriptors:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
|:-------:|:--------------:|:---------:|:----------:|:-----:|
|_Saccharomyces cerevisiae_|[Yeast 7.6](https://sourceforge.net/projects/yeast/)|3991|2691|1149|

This repository is administered by Benjamín J. Sánchez ([@BenjaSanchez](https://github.com/benjasanchez)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Installation

### Required Software - User:

* Matlab user:
  * A functional Matlab installation (MATLAB 7.3 or higher).
  * The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
* Python user:
  * Python 2.7, 3.4, 3.5 or 3.6
  * [cobrapy](https://github.com/opencobra/cobrapy)

### Required Software - Contributor:

* Both of the previous Matlab requirements.
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).
* A [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Dependencies - Recommended Software:
* For Matlab, the [libSBML MATLAB API](https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface/) (version 5.17.0 is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.

### Installation Instructions
* For users: Clone it from [`master`](https://github.com/SysBioChalmers/yeast-GEM) in the Github repo, or just download [the latest release](https://github.com/SysBioChalmers/yeast-GEM/releases).
* For contributors: Fork it to your Github account, and create a new branch from [`devel`](https://github.com/SysBioChalmers/yeast-GEM/tree/devel).

## Usage

Make sure to load/save the model with the corresponding wrapper functions!
* In Matlab:
  * Loading: `complementaryScripts/loadYeastModel.m`
  * Saving: `complementaryScripts/saveYeastModel.m`
* In Python:
  * Loading: `complementaryScripts/loadYeastModel.py`
  * Saving: currently unavailable

## Model Files

The model is available in `.xml`, `.txt`, `.yml`, `.mat` and `.xlsx` (the last 2 extensions only in `master`). Additionally, the following 2 files are available:
* `dependencies.txt`: Tracks versions of toolboxes & SBML used for saving the model.
* `boundaryMets.txt`: Contains a list of all boundary metabolites in model, listing the id and name.

### Complementary Scripts

* `missingFields`: Folder with functions for adding missing fields to the model.
* `modelCuration`: Folder with curation functions.
* `otherChanges`: Folder with other types of changes.
* `increaseVersion.m`: Updates the version of the model in `version.txt` and as metaid in the `.xml` file. Saves the model as `.mat` and as `.xlsx`
* `loadYeastModel.m`: Loads the yeast model from the `.xml` file for Matlab.
* `loadYeastModel.py`: Loads the yeast model from the `.xml` file for Python.
* `saveYeastModel.m`: Saves yeast model as a `.xml`, `.yml` and `.txt` file, and updates `boundaryMets.txt` and `dependencies.txt`.

### Complementary Data

* `databases`: Yeast data from different databases (KEGG, SGD, swissprot, etc).
* `modelCuration`: Data files used for performing curations to the model. Mostly lists of new rxns, mets or genes added (or fixed) in the model.
* `physiology`: Data on yeast growth under different conditions, biomass composition, gene essentiality experiments, etc.

## Citation

* All yeast-GEM releases are archived in [Zenodo](https://zenodo.org/badge/latestdoi/52777598), for you to cite the specific version of yeast-GEM that you used in your study. Find the corresponding DOI [here](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%221494182%22&sort=-publication_date&all_versions=True).

* Additionally, if you would like to cite the yeast-GEM project, you may also refer to the [yeast 7 paper](https://doi.org/10.1089/ind.2013.0013), and point to the new link in the text, e.g.: *"The yeast consensus genome-scale model [Aung et al. 2013], which is currently being hosted at https://github.com/SysBioChalmers/yeast-GEM, [...]"*.

* Finally, if you would like to cite the idea of hosting a genome-scale model in Github, you may also refer to the [RAVEN 2 paper](https://doi.org/10.1371/journal.pcbi.1006541), which mentions this idea and exemplifies it on [Streptomyces_coelicolor-GEM](https://github.com/SysBioChalmers/Streptomyces_coelicolor-GEM).

## Contributing

Contributions are always welcome! Please read the [contributions guideline](https://github.com/SysBioChalmers/yeast-GEM/blob/master/.github/CONTRIBUTING.md) to get started.

## Contributors

* [Mihail Anton](https://www.chalmers.se/en/staff/Pages/mihail-anton.aspx) ([@mihai-sysbio](https://github.com/mihai-sysbio)), Chalmers University of Technology, Sweden
* [Moritz Beber](https://www.dtu.dk/english/service/phonebook/person?id=121829&tab=2&qt=dtupublicationquery) ([@Midnighter](https://github.com/Midnighter)), Danish Technical University, Denmark
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden
* [Dimitra Lappa](https://www.chalmers.se/en/Staff/Pages/lappa.aspx) ([@demilappa](https://github.com/demilappa)), Chalmers University of Technology, Sweden
* [Feiran Li](https://www.chalmers.se/en/staff/Pages/feiranl.aspx) ([@feiranl](https://github.com/feiranl)), Chalmers University of Technology, Sweden
* [Christian Lieven](https://www.dtu.dk/english/service/phonebook/person?id=103199&tab=2&qt=dtupublicationquery) ([@ChristianLieven](https://github.com/ChristianLieven)), Danish Technical University, Denmark
* [Hongzhong Lu](https://www.chalmers.se/en/Staff/Pages/luho.aspx) ([@hongzhonglu](https://github.com/hongzhonglu)), Chalmers University of Technology, Sweden
* [Simonas Marcišauskas](https://www.chalmers.se/en/Staff/Pages/simmarc.aspx) ([@simas232](https://github.com/simas232)), Chalmers University of Technology, Sweden
* [Thomas Pfau](https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/people/thomas_pfau) ([@tpfau](https://github.com/tpfau)), University of Luxembourg, Luxembourg
* [Benjamín J. Sánchez](https://www.chalmers.se/en/staff/Pages/bensan.aspx) ([@BenjaSanchez](https://github.com/benjasanchez)), Chalmers University of Technology, Sweden
