# Asad-Tabeni_RiparianPlantMetacommunityStructure
R Script for the analyisis of the scientific paper: Bridging human drivers on riparian plant metacommunity structure along a dryland natural-rural-urban gradient

This repository contains the R scripts used in the analysis for the paper:
"Bridging human drivers on riparian plant metacommunity structure along a dryland natural-rural-urban gradient"

## Contents and data structure

• statistical_analisis.R: Preprocesses datasets, computes β-diversity metrics (abundance-based Bray–Curtis dissimilarity), conducts statistical analyses, and generates figures corresponding to the main results reported in the manuscript.

• LDI_fragmentation_calc: Processes geospatial datasets previously curated in a GIS environment and computes the Landscape Division Index (LDI) for each sampling site.

## Data Requirements

Input datasets required for statistical_analysis.R are available in the Supplementary Materials of the paper:
STable1.xlsx
STable2.xlsx
STable3.xlsx

The script specifies the import paths and expected structure.

Important formatting notes: Files may need conversion from .xlsx to .csv depending on the user's R configuration.

Each spreadsheet contains: Sheet 1: Metadata; Sheet 2: Raw data

Ensure that the correct sheet is selected during import.

## Software environment
Scripts were devolped and tested using RStudio 2024.12.1 They rely on commonly used ecological, geospatial, and statistical libraries. All dependencies are explicitly loaded at the beginning of each script.

Install any missing package with the following command in R:

```r
install.packages("name_of_package")

