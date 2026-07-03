# Vavilov_2.0_pipeline
This is a repo with R codes to perform edible plants species distribution models and identify possible VAC species
For more details please see the Excel file attached named VAVILOV2_PIPELINE_DESCRIPTION.xlsx 
# Structure:
- 000_Preprocessing: Performs preprocessing including taxa name resolution and GBIF occurrences obtention
- 001_observed_richness_buffer: Performs 10 km ratio around occurrences
- 002_SDM: Performs species distribution models
- 003_CONCAVE_HULL: Performs concave hull to be used as native range for species and obtain realized niche
- 004_realized_niche: Performs crops to obtain SDM cuts and summary maps into the concave hulls
- 005_soil_outliers: Performs quantile analyses to obtain Vision candidates species for Adapted Crops and Soils (VACS) visition
- 006_STATS: Performs barplots using a main archive with species name, plant parts and food groups

## Pipeline scheme
<img width="1059" height="681" alt="image" src="https://github.com/user-attachments/assets/b8366359-9434-458f-9bed-4e1789095433" />

### Relevant bibliography or resources
- Cook, Frances E. M. 1995. Economic Botany Data Collection Standard. Prepared for the International Working Group on Taxonomic Databases for Plant Sciences (TDWG) by the Royal Botanic Gardens, Kew. ISBN 0947643710. http://rs.tdwg.org/ebdc/doc/specification/
- https://foodplantsinternational.com/
