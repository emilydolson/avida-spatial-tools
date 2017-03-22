from pysal.esda.getisord import G_Local
from collections import Counter
from .utils import *
from .patch_analysis import *
from .landscape_stats import *


def diversity_map_spatial_weights(world, weights):
    return make_diversity_map(world, weights=weights)


def diversity_map(world, neighbor_func=get_moore_neighbors_toroidal):
    return make_diversity_map(world, neighbor_func=neighbor_func)


def make_diversity_map(world, neighbor_func=None, weights=None):
    world_x = len(world[0])
    world_y = len(world)

    data = initialize_grid((world_x, world_y), -1)

    for y in range(world_y):
        for x in range(world_x):
            local_vals = Counter()

            if weights:
                neighbors = weights.neighbors[world_x*y + x]
                local_vals = Counter()
                local_vals[world[y][x]] += 1
                for i, n in enumerate(neighbors):
                    local_vals[world[n // world_x][n % world_x]] \
                        += weights.weights[world_x*y + x][i]

            elif neighbor_func:
                neighbors = neighbor_func([x, y], (world_x, world_y))
                neighbors.append([x, y])

                for n in neighbors:
                    local_vals[world[n[1]][n[0]]] += 1

            else:
                raise RuntimeError("""Diversity map needs a neighbor_func or
                                   weights matrix.""")

            data[y][x] = entropy(dict(local_vals))

    return data
