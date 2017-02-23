# groundfish-data-analysis
Code for analyzing groundfish data off newfoundland.

Note: this code requires a cleaned community data file, available by request from Fisheries and Oceans Canada; contact Pierre Pepin (<pepinp@dfo-mpo.gc.ca>).


## To use this code

Place the following data files in the data folder:

1. DFO_SURVEYS_text_cleaner.csv
2. DFO_RV_Surveys_kg-per-tow_fish+shrimp+crab_Spring_Fall_2013.csv
3. Heike_Traits.csv
4. voronoi_shapes.Rdat

Ensure the working directory for R is set to the overall project folder

Run the scripts in the following order: 

1. compiling_data.R -> loads data, cleans it, and resaves data as files for other analyses
2. population_mean_and_synchrony.R -> generates figure 1
3. ordinations_and_relative_composition.R -> generates figure 2
4. functional_diversity.R -> generates figure 3, figure S2, and secondarily estimates of FDis used in figure 6; these are saved in the data folder
5. clustering_and_plotting_maps.R -> Generates figure 4. 
6. distance_diversity_analyses.R -> generates figure 5 and figure S4
7. overall_community_change.R -> Generates figure 6. 

