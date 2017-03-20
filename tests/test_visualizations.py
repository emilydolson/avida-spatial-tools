from avidaspatial import *
import pytest
from matplotlib.testing.decorators import image_comparison


@pytest.mark.mpl_image_compare
def test_visualize_environment():
    fig = plt.figure()
    env = "tests/example_environment.cfg"
    environment = parse_environment_file(env, (11, 5))
    environment = convert_world_to_phenotype(environment)
    heat_map(environment, "visualize_environment", mask_zeros=True)
    # grid = visualize_environment(env, (11,5), "test_env.png")
    return fig


@pytest.mark.mpl_image_compare
def test_n_tasks_heatmap():
    fig = plt.figure()
    grid = "tests/grid_task.100000.dat"
    grid = load_grid_data(grid)
    grid = make_count_grid(grid)
    grid = agg_grid(grid, mean)
    heat_map(grid, "n_tasks")
    return fig


@pytest.mark.mpl_image_compare
def test_color_percentages():
    fig = plt.figure()
    file_list = glob.glob("tests/grid_task.*.dat")
    color_percentages(file_list, 2)
    return fig


@pytest.mark.mpl_image_compare
def test_color_percentages2():
    fig = plt.figure()
    file_list = glob.glob("tests/grid_task.*.dat")
    color_percentages2(file_list)
    return fig


@pytest.mark.mpl_image_compare
def test_make_species_cluster_grid():
    fig = plt.figure()
    grid = "tests/grid_task.100000.dat"
    grid = load_grid_data(grid)
    grid, n = assign_ranks_by_cluster(grid, 12)
    grid = agg_grid(grid, median)
    heat_map(grid, "species_cluster_grid", mask_zeros=True)
    return fig


@pytest.mark.mpl_image_compare
def test_make_species_hue_grid():
    fig = plt.figure()
    grid = "tests/grid_task.100000.dat"
    grid = load_grid_data(grid)
    grid = agg_grid(grid)
    heat_map(grid, "species_hue_grid")
    return fig


@pytest.mark.mpl_image_compare
def test_paired_environment_phenotype_grid():
    fig = plt.figure()
    env = "tests/example_environment.cfg"
    grid = "tests/grid_task.100000.dat"

    world = parse_environment_file(env, (11, 5))
    phenotypes = load_grid_data(grid)
    # world, phenotypes, n=rank_environment_and_phenotypes(world,phenotypes,20)
    world = convert_world_to_phenotype(world)
    phenotypes, n = assign_ranks_by_cluster(phenotypes, 12)
    phenotypes = agg_grid(phenotypes)
    length = get_pallete_length(phenotypes)
    palette = sns.hls_palette(length, s=1)
    world.task_palette = palette

    paired_environment_phenotype_grid(world, phenotypes, denom=n)
    return fig


@pytest.mark.mpl_image_compare
def test_paired_environment_phenotype_grid2():
    fig = plt.figure()
    env = "tests/example_environment2.cfg"
    grid = "tests/grid_task.100000.dat"

    world = parse_environment_file(env, (11, 5))
    phenotypes = load_grid_data(grid)
    world, phenotypes, n = \
        rank_environment_and_phenotypes(world, phenotypes, 20)
    phenotypes = agg_grid(phenotypes)
    palette = sns.hls_palette(n, s=1)
    world.task_palette = palette
    world.resource_palette = palette
    paired_environment_phenotype_grid(world, phenotypes, denom=n)
    return fig


@pytest.mark.mpl_image_compare
def test_paired_environment_phenotype_grid3():
    fig = plt.figure()
    env = "tests/conservation-9patches_10each_environment.cfg"
    grid = "tests/grid_task60by60.0.dat"

    world = parse_environment_file(env, (60, 60))
    phenotypes = load_grid_data(grid)

    world = convert_world_to_phenotype(world)
    phenotypes = agg_grid(phenotypes)

    paired_environment_phenotype_grid(world, phenotypes)
    return fig


# @image_comparison(baseline_images=['phenotype_niches_circlesexample_environment.png'])
@pytest.mark.mpl_image_compare
def test_paired_environment_phenotype_grid_circles():
    fig = plt.figure()
    env = "tests/example_environment.cfg"
    grid = "tests/grid_task.100000.dat"
    world = parse_environment_file(env, (11, 5))
    world = convert_world_to_phenotype(world)
    phenotypes = load_grid_data(grid)
    phenotypes = agg_grid(phenotypes)
    length = get_pallete_length(phenotypes)
    palette = sns.hls_palette(length, s=1)
    world.task_palette = palette
    paired_environment_phenotype_grid_circles(world, phenotypes)
    return fig


@pytest.mark.mpl_image_compare
def test_optimal_phenotypes():
    fig = plt.figure()
    env = "tests/example_environment.cfg"
    grid = "tests/grid_task.100000.dat"
    grid = load_grid_data(grid)
    env = parse_environment_file(env, (11, 5))
    phenotypes = make_optimal_phenotype_grid(env, grid)
    phenotypes = agg_grid(phenotypes)
    heat_map(phenotypes, "optimal")
    return fig


def test_paired_environment_phenotype_movie():
    fig = plt.figure()
    phenotypes = load_grid_data(glob.glob("tests/grid_task.*.dat"))
    env = parse_environment_file("tests/example_environment.cfg", (11, 5))
    env = convert_world_to_phenotype(env)
    # env, phenotypes, n = rank_environment_and_phenotypes(env, phenotypes, 20)
    # print env.grid
    phenotypes, n = assign_ranks_by_cluster(phenotypes, 12)
    paired_environment_phenotype_movie(env, phenotypes)


def test_make_movie():
    fig = plt.figure()
    print glob.glob("tests/positivegrid_task.*.dat")
    phenotypes = load_grid_data(glob.glob("tests/positivegrid_task.*.dat"),
                                "float")
    make_movie(phenotypes)

if __name__ == "__main__":
    # test_paired_environment_phenotype_movie()
    # test_color_percentages()
    # test_visualize_environment()
    # test_paired_environment_phenotype_grid_circles()
    # test_paired_environment_phenotype_grid()
    test_make_movie()
