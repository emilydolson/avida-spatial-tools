from visualizations import *

def test_visualize_environment():
    env = "test/example_environment.cfg"

    grid = visualize_environment(env, (11,5), "test_env.png")

def test_color_percentages():
    file_list = glob.glob("test/grid_task.*.dat")
    color_percentages2(file_list)

def test_make_species_grid():
    grid = "test/grid_task.100000.dat"
    make_species_grid(grid)

def test_paired_environment_phenotype_grid():
    env = "test/example_environment.cfg"
    grid = "test/grid_task.100000.dat"

    paired_environment_phenotype_grid(grid, env)

def test_paired_environment_phenotype_grid_circles():
    env = "test/example_environment.cfg"
    grid = "test/grid_task.100000.dat"

    paired_environment_phenotype_grid_circles(grid, env)


def test_optimal_phenotypes():
    env = "test/example_environment.cfg"
    grid = "test/grid_task.100000.dat"
    optimal_phenotypes(env, grid)

def test_paired_environment_phenotype_movie():
    paired_environment_phenotype_movie(
        glob.glob("test/grid_task.*.dat")[:], 
        "test/example_environment.cfg", 15)

if __name__ == "__main__":
     test_paired_environment_phenotype_movie()
