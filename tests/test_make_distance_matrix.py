from avidaspatial import *


def test_make_toroidal_weights():
    weights = make_toroidal_weights(5, 5)
    for el in weights.neighbors.values():
        assert(len(el) == 4)


def test_write_matrix():
    write_toroidal_contig_matrix_to_file(6, 6)
