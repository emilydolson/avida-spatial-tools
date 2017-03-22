from pysal.esda.getisord import G_Local
from collections import Counter
from .utils import *
from .patch_analysis import *
from .landscape_stats import *


def compute_diversity_gradient_snazzy_gaussian(world, sd=5, recursive=False):
    world_x = len(world[0])
    world_y = len(world)
    if sd is None:
        sd = calculate_optimal_grain(world)

    focal = (0, 0)

    # Intialize phenotype count dicts
    dicts = deepcopy(world)
    for i in range(world_x):
        for j in range(world_y):
            dicts[i][j] = {}

    for i in range(world_x):
        for j in range(world_y):

            # Edge cases
            next_x = floor((i+1)/float(x_regions))
            x_prop = 1
            y_prop = 1
            if i < world_x - 1 and next_x != x:
                x_prop = (i/float(x_regions)-x)/((i+1)/float(x_regions) -
                                                 i/float(x_regions))
            if j < world_y - 1 and floor((j+1)/float(y_regions)) != y:
                y_prop = (j/float(y_regions)-y)/((j+1)/float(y_regions) -
                                                 j/float(y_regions))

            dict_increment(regions[x][y], world[i][j][0], x_prop*y_prop)

            if x_prop < 1 and y_prop < 1:
                dict_increment(regions[x+1][y+1], world[i][j][0],
                               (1-x_prop)*(1-y_prop))
            if x_prop < 1:
                dict_increment(regions[x+1][y], world[i][j][0],
                               (1-x_prop)*y_prop)
            if y_prop < 1:
                dict_increment(regions[x][y+1], world[i][j][0],
                               x_prop*(1-y_prop))

    entropies = []
    for i in range(grain):
        entropies.append([])
        for j in range(grain):
            # print regions[i][j]
            entropies[i].append(entropy(regions[i][j]))

    if not recursive:
        plt.imshow(entropies, interpolation="none", cmap="jet")
        plt.show()


def get_prop_in_region(var, index, regions):
    return (index/float(regions) - var)/(
        (index+1)/float(regions) - index/float(regions))


def compute_diversity_gradient(world, grain=None, recursive=False):
    world_x = len(world[0])
    world_y = len(world)
    if grain is None:
        grain = calculate_optimal_grain(world)

    regions = []
    for i in range(grain):
        regions.append([])
        for j in range(grain):
            regions[i].append({})

    x_regions = world_x/float(grain)
    y_regions = world_y/float(grain)

    for i in range(world_x):
        for j in range(world_y):
            x = int(floor(i/float(x_regions)))
            y = int(floor(j/float(y_regions)))

            # Edge cases
            next_x = floor((i+1)/float(x_regions))
            x_prop = 1
            y_prop = 1
            if i < world_x - 1 and next_x != x:
                x_prop = get_prop_in_region(x, i, x_regions)
            if j < world_y - 1 and floor((j+1)/float(y_regions)) != y:
                y_prop = get_prop_in_region(y, j, y_regions)

            dict_increment(regions[x][y], world[i][j][0], x_prop*y_prop)

            if x_prop < 1 and y_prop < 1:
                dict_increment(regions[x+1][y+1], world[i][j][0],
                               (1-x_prop)*(1-y_prop))
            if x_prop < 1:
                dict_increment(regions[x+1][y], world[i][j][0],
                               (1-x_prop)*y_prop)
            if y_prop < 1:
                dict_increment(regions[x][y+1], world[i][j][0],
                               x_prop*(1-y_prop))

    entropies = []
    for i in range(grain):
        entropies.append([])
        for j in range(grain):
            # print regions[i][j]
            entropies[i].append(entropy(regions[i][j]))

    if not recursive:
        plt.imshow(entropies, interpolation="none")
        plt.show()

    # getis_ord(entropies)
    return entropies


def diversity_map(world, neighbor_func=get_moore_neighbors_toroidal):
    world_x = len(world[0])
    world_y = len(world)

    data = initialize_grid((world_x, world_y), -1)

    for y in range(world_y):
        for x in range(world_x):
            neighbors = neighbor_func([x, y], (world_x, world_y))
            neighbors.append([x, y])
            local_vals = Counter()
            for n in neighbors:
                local_vals[world[n[1]][n[0]]] += 1
            data[y][x] = entropy(dict(local_vals))

    return data


def calculate_optimal_grain(world, increment=None):
    start = int(len(world)/2.0)
    if increment is None:
        increment = 1

    z_vals = []
    for grain in range(start, 3, increment*-1):
        entropies = compute_diversity_gradient(world, grain, True)
        w, data = convert_to_pysal(entropies)
        z_vals.append(pysal.Moran(data, w).z_rand)
    plt.plot(range(start, 3, increment*-1), z_vals)
    plt.show()
    return 21


def getis_ord(data):
    print(np.array(data))
    w, data = convert_to_pysal(data)
    print(G_Local(data, w).p_sim)
    # plt.imshow()
    plt.show()
