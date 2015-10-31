from visualizations import *

def test_visualize_environment():
    env = "test/example_environment.cfg"
    environment = parse_environment_file(env, (11,5))
    environment = convert_world_to_phenotype(environment)
    heat_map(environment, "visualize_environment")
    #grid = visualize_environment(env, (11,5), "test_env.png")

def test_color_percentages():
    file_list = glob.glob("test/grid_task.*.dat")
    color_percentages2(file_list)

def test_make_species_cluster_grid():
    grid = "test/grid_task.100000.dat"
    grid = load_grid_data(grid)
    grid, n = assign_ranks_by_cluster(grid, 12)
    grid = agg_grid(grid)
    heat_map(grid, "species_cluster_grid")

def test_make_species_hue_grid():
    grid = "test/grid_task.100000.dat"
    grid = load_grid_data(grid)
    grid = agg_grid(grid)
    heat_map(grid, "species_hue_grid")

def test_paired_environment_phenotype_grid():
    env = "test/example_environment.cfg"
    grid = "test/grid_task.100000.dat"

    world = parse_environment_file(env, (11,5))
    phenotypes = load_grid_data(grid)
    world, phenotypes, n = rank_environment_and_phenotypes(world, phenotypes, 20)
    phenotypes = agg_grid(phenotypes)
    
    paired_environment_phenotype_grid(world, phenotypes)
    plt.clf()

def test_paired_environment_phenotype_grid_circles():
    env = "test/example_environment.cfg"
    grid = "test/grid_task.100000.dat"
    world = parse_environment_file(env, (11,5))
    world = convert_world_to_phenotype(world)
    phenotypes = load_grid_data(grid)
    phenotypes = agg_grid(phenotypes)
    paired_environment_phenotype_grid_circles(world, phenotypes)
    plt.clf()


def test_optimal_phenotypes():
    env = "test/example_environment.cfg"
    grid = "test/grid_task.100000.dat"
    grid = load_grid_data(grid)
    env = parse_environment_file(env, (11,5))
    phenotypes = make_optimal_phenotype_grid(env, grid)
    phenotypes = agg_grid(phenotypes)
    heat_map(phenotypes, "optimal")

def test_paired_environment_phenotype_movie():
    phenotypes = load_grid_data(glob.glob("test/grid_task.*.dat"))
    env = parse_environment_file("test/example_environment.cfg", (11,5))

    env, phenotypes, n = rank_environment_and_phenotypes(env, phenotypes, 20)
    
    paired_environment_phenotype_movie(env, phenotypes)

if __name__ == "__main__":
    #test_paired_environment_phenotype_movie()
    #test_color_percentages()
    #test_visualize_environment()
    #test_paired_environment_phenotype_grid_circles()
    test_paired_environment_phenotype_grid()
