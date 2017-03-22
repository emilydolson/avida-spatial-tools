from avidaspatial import *
import numpy as np


expected = [[1.4466166676282082, 0.9864267287308424, 0.9864267287308424,
            1.224394445405986, 1.224394445405986, 0.5032583347756456, 0.0,
            0.0, 0.0, 0.5032583347756456, 1.4466166676282082],
            [1.4466166676282082, 0.9864267287308424, 0.9864267287308424,
             1.3516441151533924, 1.3516441151533924, 0.7642045065086203,
             0.0, 0.0, 0.0, 0.5032583347756456, 1.4466166676282082],
            [1.224394445405986, 0.9864267287308424, 1.224394445405986,
             1.5304930567574826, 1.5304930567574826, 0.7642045065086203,
             0.0, 0.0, 0.0, 0.5032583347756456, 1.224394445405986],
            [0.5032583347756456, 0.5032583347756456, 0.9182958340544896,
             1.3516441151533924, 1.3516441151533924, 0.5032583347756456, 0.0,
             0.0, 0.0, 0.0, 0.5032583347756456],
            [0.9864267287308424, 0.9864267287308424, 0.7642045065086203,
             0.7642045065086203, 0.7642045065086203, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.9864267287308424]]


def test_diversity_map_spatial_weights():
    world = load_grid_data("tests/grid_task.10000.dat")
    world = agg_grid(world, mode)

    w = make_toroidal_weights(len(world), len(world[0]), rook=False)

    data = diversity_map_spatial_weights(world, w)
    for i, row in enumerate(data):
        for j, num in enumerate(row):
            assert(np.isclose(num, expected[i][j]))


def test_diversity_map():
    world = load_grid_data("tests/grid_task.10000.dat")
    world = agg_grid(world, mode)

    data = diversity_map(world)
    for i, row in enumerate(data):
        for j, num in enumerate(row):
            assert(np.isclose(num, expected[i][j]))


if __name__ == "__main__":
    test_diversity_map()
    test_diversity_map_spatial_weights()
