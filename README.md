# Characterisation of large transgene integrations in Chinese hamster ovary cells using a bioengineered mammalian transposase

This repository complements the publication by Marx et al (2024) (DOI link)

## Folders

(Not all of these folders appear in the gitHub repository.)

-   `data`: raw input files and tables
-   `analysis`: output tables generated by the scripts in this repository
-   `figures`: plots generated by the plotting scripts
-   `mab_sequence`: fasta file containing the sequence of the monoclonal antibody (mAb) trastuzumab

## Analysis description

This repository contains scripts for the following analysis:

1.  Quantification of the mAb subunits from the peak areas of UV-chromatograms

    -   Input: data/RelQuantIntact01.csv

    -   Script: [plot_subunits.R](plot_subunits.R) - Loads csv table with the subunits quantification and plots stacked barplot.

    -   Output: Saves stacked barplot figure_8.png in the folder `figures`
    ![Relative quantification of subunits](figures/stacked_bar/figure_8.png)

2.  Quantification and hexose bias correction of N-glycans attached to the Fc region of the mAb 2.1 Quantification of uncorrected N-glycans abundance 2.2 Quantification of hexose bias 2.3 Correction of hexose bias and final quantification of corrected N-glycans abudance

3.  Quantification of lysine variants and glycation on the intact mAb
