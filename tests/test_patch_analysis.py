from avidaspatial import *
import numpy as np


def test_area():
    assert(area([[0, 0], [0, 1], [0, 2], [0, 3]]) == 4)


def test_perimeter():
    assert(perimeter([[0, 0], [0, 1], [0, 2], [0, 3], [1, 1]]) == 12)


def test_get_edge_locations():
    patch = [[0, 0], [0, 1], [0, 2], [0, 3],
             [1, 0], [1, 1], [1, 2], [1, 3],
             [2, 0], [2, 1], [2, 2], [2, 3],
             [3, 1],
             [4, 0], [4, 1], [4, 2], [4, 3],
             [5, 0], [5, 1], [5, 2], [5, 3],
             [6, 0], [6, 1], [6, 2], [6, 3],
             [7, 0], [7, 1], [7, 2], [7, 3]]
    edge = get_edge_locations(patch)
    edge.sort()
    assert(edge == [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 3),
                    (2, 0), (2, 2), (2, 3),
                    (3, 1),
                    (4, 0), (4, 2), (4, 3),
                    (5, 0), (5, 3),
                    (6, 0), (6, 3),
                    (7, 0), (7, 1), (7, 2), (7, 3)]
           )


def test_is_edge():
    patch = [(0, 0), (0, 1), (0, 2), (0, 3),
             (1, 0), (1, 1), (1, 2), (1, 3),
             (2, 0), (2, 1), (2, 2), (2, 3),
             (3, 1),
             (4, 0), (4, 1), (4, 2), (4, 3),
             (5, 0), (5, 1), (5, 2), (5, 3),
             (6, 0), (6, 1), (6, 2), (6, 3),
             (7, 0), (7, 1), (7, 2), (7, 3)]
    assert(isedge((0, 0), patch))
    assert(not isedge((1, 1), patch))


def test_get_moore_neighbors():
    moore = get_moore_neighbors([0, 0], 60)
    moore.sort()
    expected = [[0, 59], [59, 0], [1, 0], [0, 1],
                [1, 59], [59, 59], [1, 1], [59, 1]]
    expected.sort()
    assert(moore == expected)


def test_get_rook_neighbors():
    rook = get_rook_neighbors([0, 0], 60)
    rook.sort()
    expected = [[0, 59], [59, 0], [1, 0], [0, 1]]
    expected.sort()
    assert(rook == expected)


def test_centroid():
    patch = [[0, 0], [0, 1], [0, 2], [1, 1]]
    c = centroid(patch)
    assert(np.isclose(c[0], .25))
    assert(np.isclose(c[1], 1))


def test_radius_of_gyration():
    patch = [[0, 1],
             [1, 0], [1, 2],
             [2, 1]]
    assert(np.isclose(radius_of_gyration(patch), 1))


def test_shape_index():
    patch = [[0, 0], [0, 1], [0, 2],
             [1, 0], [1, 1], [1, 2],
             [2, 0], [2, 1], [2, 2]]
    assert(shape_index(patch) == 1)
    patch = [[0, 0], [0, 1], [0, 2],
             [1, 0], [1, 1], [1, 2],
             [2, 0], [2, 1], [2, 2],
             [3, 1]]
    assert(np.isclose(shape_index(patch), 1.10679))


def test_fractal_dimension():
    patch = [[0, 0], [0, 1], [0, 2],
             [1, 0], [1, 1], [1, 2],
             [2, 0], [2, 1], [2, 2]]
    assert(fractal_dimension(patch) == 1)
    patch = [[0, 0], [0, 1], [0, 2],
             [1, 0], [1, 1], [1, 2],
             [2, 0], [2, 1], [2, 2],
             [3, 1]]
    assert(np.isclose(fractal_dimension(patch), 1.088136))


def test_related_circumscribing_circle():
    patch = [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5],
             [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [1, 5],
             [2, 0], [2, 1], [2, 2], [2, 3], [2, 4], [2, 5],
             [3, 0], [3, 1], [3, 2], [3, 3], [3, 4], [3, 5],
             [4, 0], [4, 1], [4, 2], [4, 3], [4, 4], [4, 5],
             [5, 0], [5, 1], [5, 2], [5, 3], [5, 4], [5, 5]]
    assert(np.isclose(related_circumscribing_circle(patch), 0.083267))


def test_number_core_areas():
    patch = [[0, 0], [0, 1], [0, 2], [0, 3],
             [1, 0], [1, 1], [1, 2], [1, 3],
             [2, 0], [2, 1], [2, 2], [2, 3],
             [3, 1],
             [4, 0], [4, 1], [4, 2], [4, 3],
             [5, 0], [5, 1], [5, 2], [5, 3],
             [6, 0], [6, 1], [6, 2], [6, 3],
             [7, 0], [7, 1], [7, 2], [7, 3]]

    assert(number_core_areas(patch, 1, 60) == 2)
    assert(number_core_areas(patch, 0, 60) == 1)
    assert(number_core_areas(patch, 2, 60) == 0)


def test_get_core_areas():
    patch = [[0, 0], [0, 1], [0, 2], [0, 3],
             [1, 0], [1, 1], [1, 2], [1, 3],
             [2, 0], [2, 1], [2, 2], [2, 3],
             [3, 1],
             [4, 0], [4, 1], [4, 2], [4, 3],
             [5, 0], [5, 1], [5, 2], [5, 3],
             [6, 0], [6, 1], [6, 2], [6, 3],
             [7, 0], [7, 1], [7, 2], [7, 3]]

    cores = get_core_areas(patch, 1, 60)
    cores.sort()
    for core in cores:
        core.sort()

    assert(cores == [[(1, 1), (1, 2)], [(5, 1), (5, 2), (6, 1), (6, 2)]])


def test_core_area_index():
    patch = [[0, 0], [0, 1], [0, 2], [0, 3],
             [1, 0], [1, 1], [1, 2], [1, 3],
             [2, 0], [2, 1], [2, 2], [2, 3]]
    expected = 1.0/6.0
    assert(np.isclose(core_area_index(patch, 1), expected))
