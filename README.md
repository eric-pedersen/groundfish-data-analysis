# groundfish-data-analysis
Code for analyzing groundfish data off newfoundland.

## To use this code

Ensure the working directory for R is set to the overall project folder

Run the scripts in the following order: 

1. compiling_data.R -> loads data, cleans it, and resaves data as files for other analyses
2. population_mean_and_synchrony.R -> generates figure 1
3. potential_drivers.R -> generates figure 2
4. ordinations_and_relative_composition.R -> generates figure 3
5. functional_diversity.R -> generates figure 4, figure S2, and secondarily estimates of FDis used in figure 7; these are saved in the data folder
6. clustering_and_plotting_maps.R -> Generates figure 5. 
7. distance_diversity_analyses.R -> generates figure 6 and figure S4
8. overall_community_change.R -> Generates figure 7.
9. conversion_factor_analysis.R -> Generates figures S5 and S6.

