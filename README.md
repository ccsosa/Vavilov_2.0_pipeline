# Vavilov_2.0_pipeline

This repository contains the R code used to perform species distribution modelling (SDM) for edible plant species and to identify candidate species for Vision for Adapted Crops and Soils (VACS).

For more details, see the Excel file `VAVILOV2_PIPELINE_DESCRIPTION.xlsx` included in this repository.

## Main input data

- A file with species names extracted from literature and databases
- Bioclimatic and soil variables
- A file with plant part and food group classifications (use literature or the Food Plants International database)
- [Plants of the World Online (POWO)](https://powo.science.kew.org/)

> [!IMPORTANT]
> Please be careful with your curation of plant part and food group. We recommend using the Cook, 1995 food group standard.
> This pipeline was designed to be used for a target region — please adapt the code accordingly.
> Please modify `basedir` according to your needs.

> [!CAUTION]
> This pipeline uses parallel processing with at least eight cores for several steps. Please only use this approach on a server.
> Due to the large amount of disk space required, it is recommended to have at least 2 TB of hard disk space available.

## R code folder structure

- **000_Preprocessing**: Performs preprocessing, including taxonomic name resolution and acquisition of GBIF occurrences
- **001_observed_richness_buffer**: Performs a 10 km buffer analysis around occurrences
- **002_SDM**: Performs species distribution modelling
- **003_CONCAVE_HULL**: Computes concave hulls used to constrain each species' native range and derive its realized niche
- **004_realized_niche**: Crops SDM outputs to the concave hulls and produces summary richness maps
- **005_soil_outliers**: Performs quantile analyses to identify candidate species for Vision for Adapted Crops and Soils (VACS)
- **006_STATS**: Generates barplots summarizing species counts by plant part and food group

## Pipeline scheme

<img width="1059" height="681" alt="image" src="https://github.com/user-attachments/assets/b8366359-9434-458f-9bed-4e1789095433" />

### Relevant bibliography and resources

- Cook, Frances E. M. 1995. *Economic Botany Data Collection Standard*. Prepared for the International Working Group on Taxonomic Databases for Plant Sciences (TDWG) by the Royal Botanic Gardens, Kew. ISBN 0947643710. http://rs.tdwg.org/ebdc/doc/specification/
- https://foodplantsinternational.com/

## Authors

Main authors: Chrystian C. Sosa, Tobias Fremout, and Evert Thomas.
