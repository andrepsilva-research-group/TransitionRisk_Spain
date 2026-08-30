# 🏦 Nature-related transition risk: Potential exposure of economic sectors to future protected area expansion in Spain

[![DOI](https://img.shields.io/badge/DOI-coming_soon-blue?logo=doi&logoColor=white)](https://doi.org/10.0000/placeholder)  [![Project Page](https://img.shields.io/badge/Project_Website-MISTRAFinBio-green?logo=leaflet&logoColor=white)](https://finbio.org/)

A research pipeline to explore nature-related trasition risks, by coupling systematic conservation prioritisation with asset-level economic data and using Spain's protected-area network as a proof-of-concept case. 

> 🚧 **Under active development** 🚧<br>


## Repository Overview

This repository contains all scripts and supporting materials used in the modeling and analysis pipeline.<br>

| Folder | Description |
|:-------------|:---------------------------------------------------------|
| **data** | Raw data, including trait dataset for 13 mammal species, sourced from multiple databases. |
| **src** | All R scripts used in data analysis, modeling, and visualization. |

<br>

 ## Scripts Description
 <br>

`selectSpecies.R`  
Selects and filters the mammal and bird species included classified within a threathned category and available from intersecting the IUCN Spatial Distribution files and Iberian peninsula spatial boundaries.

`taxaOccurrence_mammals.R`  
Retrieves and processes species occurrence records from GBIF to generate the occurrence dataset used for species distribution modelling.

`taxaOccurrence_birds.R`  
Retrieves and processes species occurrence records from GBIF to generate the occurrence dataset used for species distribution modelling.

`SDM.R`  
Builds species distribution models through [biomod2](https://biomodhub.github.io/biomod2/), using species occurrences from GBIF, together with climate and land-use as predictors, to estimate current and future species suitability ranges. A custom function is created to iterate through the whole process of building an SDM and do it across multiple species and scenarios.

`SDMfigures.R`  
Generates presence maps with used records vs total available ones, variable importance tables, evaluation metrics plots together with current and projected distribution maps for all species, commonly used for quality control and interpretation of results before moving forward.

`averageSuitbaility.R`  
Average SDMs raster outputs by conservation status and scenarios to achieve a layer with the average suitability for the chosen category (conservation status or future scenario) for available species.

 
