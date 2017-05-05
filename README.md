# groundfish-data-analysis
Code for analyzing groundfish data off newfoundland.

Note: this code requires a cleaned community data file, available by request from Fisheries and Oceans Canada; contact Pierre Pepin (<pepinp@dfo-mpo.gc.ca>).


## To use this code

Place the following data files in the data folder:

1. DFO_SURVEYS_text_cleaner.csv
2. DFO_RV_Surveys_kg-per-tow_fish+shrimp+crab_Spring_Fall_2013.csv
3. Heike_Traits.csv
4. voronoi_shapes.Rdat
5. fishing_effort_data.csv
6. composite_index.csv
7. In a subfolder labelled "NAFO effort data":
  * NAFO21B-80-89.txt
  * NAFO21B-90-99.txt
  * NAFO21B-2000-09.txt
  * NAFO21B-2010-14.csv

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

