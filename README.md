# Asad-Tabeni_RiparianPlantMetacommunityStructure
R Script for the analyisis of the scientific paper: Bridging human drivers on riparian plant metacommunity structure along a dryland natural-rural-urban gradient

This repository contains the R scripts used in the analysis for the paper:
"Urban environmental filters shape plant functional composition and diversity in dryland riparian habitats"

##Contents

• statistical_analisis.R: Prepares the data, performs the calculations of the Beta-diversity (Bray-Curtis abundance based) index, and performs statistical analysis and generates plots for the main results presented in the paper.
• LDI_fragmentation_calc: Prepares the geospatial data already gathered and curated in a GIS environment and performs the calculations of the Landscape Division Index (LDI) of each sampling site.

##Requirements

Data availability for these analysis are presented in the supplementary matterial of the paper (STable1; STable2; STable3), and are specified in the script where to import them to the R environment. 

The scripts were written in RStudio 2024.12.1 and use standard ecological, geospatial and statistical packages. All required packages are loaded at the beginning of each script.

You can install any missing packages using the following example command in R:

install.packages("name_of_package")
