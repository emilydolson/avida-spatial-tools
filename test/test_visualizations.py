from visualizations import *

def test_visualize_environment():
    #env = ["distance0", "distance10", "distance21", "distance29"]
    env = ["/media/emily/hdd/resource-heterogeneity/environmentFiles/env50012.cfg", "/media/emily/hdd/resource-heterogeneity/environmentFiles/env50013.cfg"]#, "/media/emily/hdd/resource-heterogeneity/environmentFiles/env50014.cfg", "/media/emily/hdd/resource-heterogeneity/environmentFiles/env50015.cfg"]
    #grid = "/home/emily/repos/resource-heterogeneity/experiment/radius8_distance12_common/heterogeneity_replication_11857/grid_task.100000.dat"
    grid = visualize_environment(env, (59,59), "test_env.png")

def test_color_percentages():
    file_list = glob.glob("/media/emily/hdd/resource-heterogeneity/experiment/randomAnchors/*/grid_task.100000.dat")[3:4]
    color_percentages2(file_list)

def test_make_species_grid():
    grid = "/media/emily/hdd/resource-heterogeneity/experiment/inflow100_radius08_distance12_common/heterogeneity_replication_31212/grid_task.100000.dat"
    make_species_grid(grid)

def test_paired_environment_phenotype_grid():
    env = "/media/emily/hdd/resource-heterogeneity/environmentFiles/env11097.cfg"
    grid = "/media/emily/hdd/resource-heterogeneity/experiment/inflow100_radius24_commonresources/heterogeneity_replication_11097/grid_task.100000.dat"
    #grid = "/home/emily/repos/resource-heterogeneity/experiment/randomized_entropy/heterogeneity_replication_31204/grid_task.100000.dat"
    paired_environment_phenotype_grid_circles(grid, env)

def test_optimal_phenotypes():
    env = "/media/emily/hdd/resource-heterogeneity/environmentFiles/env50047.cfg"
    grid = "/media/emily/hdd/resource-heterogeneity/experiment/randomized_entropy/heterogeneity_replication_50047/grid_task.100000.dat"
    optimal_phenotypes(env, grid, (59,59))
