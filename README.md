# avida-spatial-tools
Some useful Python scripts for working with spatial structure in Avida

#Dependencies:
* Python2.x
* Matplotlib
* Seaborn
* Numpy
* Scipy
* Pysal

#Workflow:

* Parse the data:
  * One environment file: parse_environment_file()
  * Multiple environment files: parse_environment_file_list()
  * One or more grid_task files: load_grid_data()
* Transform the data (optional):
  * Convert values to counts of tasks performed or resources available: make_count_grid()
  * Convert phenotype values to numbers representing deviations from then optimal phenotype: make_optimal_phenotype_grid()
  * Convert values to ranks indicating the complexity of a phenoytype or resource set relative to others in the environment: assign_ranks_by_cluster()
  * Convert values to lists representing the percentage of organisms in that cell doing each task: task_percentages()
* Aggregate the data (if you have multiple files):
  * agg_grid()
* Visualize the data:
  * Color grid
  * Color grid by hue mixing
  * Overlay circles representing phenotypes on colored grid representing environment
    * Represent phenotypes as circles with single color
    * Represent phenotypes as concentric circles representing tasks that individual can do
  * Make a movie showing overlaid phenotypes changing over time


