import random
import collections
import math
import patch_analysis as pa
from copy import deepcopy
import scipy.spatial.distance
import numpy as np
import matplotlib.pyplot as plt

POP_SIZE = 100
GENERATIONS = 100
THRESHOLD = .005
MUTATION_RATE = .005
CROSSOVER_PROB = .05
WORLD_X = 60
WORLD_Y = 60
K = 5
TOURNAMENT_SIZE = 10

ALL_CELLS = set([(x, y) for x in range(WORLD_X) for y in range(WORLD_Y)])

ShapeStats = collections.namedtuple("ShapeStats",
                                ["area", "perimeter", "perimeter_area_ratio",
                                 "radius_of_gyration", "shape_index",
                                 "fractal_dimension",
                                 "related_circumscribing_circle",
                                 "contiguity_index", "core_area1","core_area3",
                                 "core_area6", "n_core_areas1","n_core_areas3",
                                 "n_core_areas6", "core_area_index1",
                                 "core_area_index3", "core_area_index6"])

class Patch:
    def __init__(self, patch):
        self.cells = patch
        self.calc_stats()

    def __getitem__(self, b):
        return self.stats[b]

    def calc_stats(self):

        #Find biggest patch
        cores = pa.traverse_core(self.cells, WORLD_X)
        if len(cores) == 0:
            patch = []
        else:
            patch = max(cores, key=len)
        self.patch = patch
        self.grid = cells_to_grid(self.cells)

        #Calc metrics
        area = pa.area(patch)
        perimeter = pa.perimeter(patch)
        if area == 0:
            perimeter_area_ratio = 0
        else:
            perimeter_area_ratio = float(perimeter)/float(area)
        radius_of_gyration = pa.radius_of_gyration(patch)
        shape_index = pa.shape_index(patch, float(perimeter))
        fractal_dimension = pa.fractal_dimension(patch, float(perimeter))
        related_circumscribing_circle =\
                                    pa.related_circumscribing_circle(patch)
        contiguity_index = pa.contiguity_index(patch)

        #Sequentially removoe layers of edge cells
        #ca1  = pa.get_core_areas(patch, 1)
        #ca3  = pa.get_core_areas([cell for core in ca1 for cell in core], 2)
        #ca6  = pa.get_core_areas([cell for core in ca3 for cell in core], 3)

        core_area1 = 0#sum([len(i) for i in ca1])
        core_area3 = 0#sum([len(i) for i in ca3])
        core_area6 = 0#sum([len(i) for i in ca6])
        n_core_areas1 = 0#len(ca1)
        n_core_areas3 = 0#len(ca3)
        n_core_areas6 = 0#len(ca6)

        if area == 0:
            core_area_index1 = 0
            core_area_index3 = 0
            core_area_index6 = 0
        else:
            core_area_index1 = core_area1/float(area)
            core_area_index3 = core_area3/float(area)
            core_area_index6 = core_area6/float(area)

        perimeter = 0

        self.stats = ShapeStats(area, perimeter, perimeter_area_ratio,\
                                 radius_of_gyration, shape_index, \
                                 fractal_dimension,
                                 related_circumscribing_circle,
                                 contiguity_index, core_area1, core_area3,\
                                 core_area6, n_core_areas1, n_core_areas3,\
                                 n_core_areas6, core_area_index1,\
                                 core_area_index3, core_area_index6)

def sparseness(patch, archive):

    #Check to make sure this is a patch we even care about
    #print patch.stats.area
    if patch.stats.area < 25 or patch.stats.area > 900:
        return 0

    #print patch
    if not archive:
        #print "not in archive"
        return float("inf")

    dists = scipy.spatial.distance.pdist([patch.stats[1:]] + [el.stats[1:] for el in archive], "seuclidean")
    dists = dists[:len(archive)]
    dists.sort()

    if dists[0] == 0:
        #If this is identical to something in the archive
        #I don't care how far away it is from everything else
        return 0
    else:
        result =  sum(dists[:K])/float(K)

    if not np.isnan(result):
        return result

    dists = scipy.spatial.distance.pdist([patch.stats[1:]] + [el.stats[1:] for el in archive], "euclidean")
    dists = dists[:len(archive)]
    dists.sort()

    if dists[0] == 0:
        return 0
    else:
        return sum(dists[:K])/float(K)

def cells_to_grid(cells):
    image = np.zeros((WORLD_X, WORLD_Y))
    for cell in cells:
        image[tuple(cell)] = 1
    return image

def write_results(patches):
    lines = [", ".join(Patch([]).stats._fields) + ", id"]

    for i, patch in enumerate(patches):
        lines.append(", ".join([str(stat) for stat in patch.stats]))
        lines[-1] += ", " + str(i)
        outfile = open("cells_" + str(i), "w")
        outfile.write(", ".join([str(cell) for cell in patch.cells]))
        outfile.close()
        image = cells_to_grid(patch.patch)
        plt.matshow(image)
        plt.grid(False)
        plt.savefig("patch_"+str(i)+".png")
        plt.close()

    with open("patch_stats.csv", "w") as outfile:
        outfile.write("\n".join(lines))


seed_patch_2cores = Patch([[20, 6], [20, 7], [20, 8], [20, 9], [20, 10], [20, 11], [20, 12], [20, 13], [20, 14], [20, 15], [20, 16], [20, 17], [20, 18], [20, 19], [20, 20], [20, 21], [20, 22], [21, 6], [21, 7], [21, 8], [21, 9], [21, 10], [21, 11], [21, 12], [21, 13], [21, 14], [21, 15], [21, 16], [21, 17], [21, 18], [21, 19], [21, 20], [21, 21], [21, 22], [22, 6], [22, 7], [22, 8], [22, 9], [22, 10], [22, 11], [22, 12], [22, 13], [22, 14], [22, 15], [22, 16], [22, 17], [22, 18], [22, 19], [22, 20], [22, 21], [22, 22], [23, 6], [23, 7], [23, 8], [23, 9], [23, 10], [23, 11], [23, 12], [23, 13], [23, 14], [23, 15], [23, 16], [23, 17], [23, 18], [23, 19], [23, 20], [23, 21], [23, 22], [24, 6], [24, 7], [24, 8], [24, 9], [24, 10], [24, 11], [24, 12], [24, 13], [24, 14], [24, 15], [24, 16], [24, 17], [24, 18], [24, 19], [24, 20], [24, 21], [24, 22], [25, 6], [25, 7], [25, 8], [25, 9], [25, 10], [25, 11], [25, 12], [25, 13], [25, 14], [25, 15], [25, 16], [25, 17], [25, 18], [25, 19], [25, 20], [25, 21], [25, 22], [26, 6], [26, 7], [26, 8], [26, 9], [26, 10], [26, 11], [26, 12], [26, 13], [26, 14], [26, 15], [26, 16], [26, 17], [26, 18], [26, 19], [26, 20], [26, 21], [26, 22], [27, 6], [27, 7], [27, 8], [27, 9], [27, 10], [27, 11], [27, 12], [27, 13], [27, 14], [27, 15], [27, 16], [27, 17], [27, 18], [27, 19], [27, 20], [27, 21], [27, 22], [28, 6], [28, 7], [28, 8], [28, 9], [28, 10], [28, 11], [28, 12], [28, 13], [28, 14], [28, 15], [28, 16], [28, 17], [28, 18], [28, 19], [28, 20], [28, 21], [28, 22], [29, 6], [29, 7], [29, 8], [29, 9], [29, 10], [29, 11], [29, 12], [29, 13], [29, 14], [29, 15], [29, 16], [29, 17], [29, 18], [29, 19], [29, 20], [29, 21], [29, 22], [30, 6], [30, 7], [30, 8], [30, 9], [30, 10], [30, 11], [30, 12], [30, 13], [30, 14], [30, 15], [30, 16], [30, 17], [30, 18], [30, 19], [30, 20], [30, 21], [30, 22], [31, 6], [31, 7], [31, 8], [31, 9], [31, 10], [31, 11], [31, 12], [31, 13], [31, 14], [31, 15], [31, 16], [31, 17], [31, 18], [31, 19], [31, 20], [31, 21], [31, 22], [32, 6], [32, 7], [32, 8], [32, 9], [32, 10], [32, 11], [32, 12], [32, 13], [32, 14], [32, 15], [32, 16], [32, 17], [32, 18], [32, 19], [32, 20], [32, 21], [32, 22], [33, 6], [33, 7], [33, 8], [33, 9], [33, 10], [33, 11], [33, 12], [33, 13], [33, 14], [33, 15], [33, 16], [33, 17], [33, 18], [33, 19], [33, 20], [33, 21], [33, 22], [34, 6], [34, 7], [34, 8], [34, 9], [34, 10], [34, 11], [34, 12], [34, 13], [34, 14], [34, 15], [34, 16], [34, 17], [34, 18], [34, 19], [34, 20], [34, 21], [34, 22], [35, 6], [35, 7], [35, 8], [35, 9], [35, 10], [35, 11], [35, 12], [35, 13], [35, 14], [35, 15], [35, 16], [35, 17], [35, 18], [35, 19], [35, 20], [35, 21], [35, 22], [36, 6], [36, 7], [36, 8], [36, 9], [36, 10], [36, 11], [36, 12], [36, 13], [36, 14], [36, 15], [36, 16], [36, 17], [36, 18], [36, 19], [36, 20], [36, 21], [36, 22], [37, 6], [37, 7], [37, 8], [37, 9], [37, 10], [37, 11], [37, 12], [37, 13], [37, 14], [37, 15], [37, 16], [37, 17], [37, 18], [37, 19], [37, 20], [37, 21], [37, 22], [29, 23], [29, 24], [29, 25], [29, 26], [29, 27], [29, 28], [29, 29], [29, 30], [29, 31], [20, 32], [20, 33], [20, 34], [20, 35], [20, 36], [20, 37], [20, 38], [20, 39], [20, 40], [20, 41], [20, 42], [20, 43], [20, 44], [20, 45], [20, 46], [21, 32], [21, 33], [21, 34], [21, 35], [21, 36], [21, 37], [21, 38], [21, 39], [21, 40], [21, 41], [21, 42], [21, 43], [21, 44], [21, 45], [21, 46], [22, 32], [22, 33], [22, 34], [22, 35], [22, 36], [22, 37], [22, 38], [22, 39], [22, 40], [22, 41], [22, 42], [22, 43], [22, 44], [22, 45], [22, 46], [23, 32], [23, 33], [23, 34], [23, 35], [23, 36], [23, 37], [23, 38], [23, 39], [23, 40], [23, 41], [23, 42], [23, 43], [23, 44], [23, 45], [23, 46], [24, 32], [24, 33], [24, 34], [24, 35], [24, 36], [24, 37], [24, 38], [24, 39], [24, 40], [24, 41], [24, 42], [24, 43], [24, 44], [24, 45], [24, 46], [25, 32], [25, 33], [25, 34], [25, 35], [25, 36], [25, 37], [25, 38], [25, 39], [25, 40], [25, 41], [25, 42], [25, 43], [25, 44], [25, 45], [25, 46], [26, 32], [26, 33], [26, 34], [26, 35], [26, 36], [26, 37], [26, 38], [26, 39], [26, 40], [26, 41], [26, 42], [26, 43], [26, 44], [26, 45], [26, 46], [27, 32], [27, 33], [27, 34], [27, 35], [27, 36], [27, 37], [27, 38], [27, 39], [27, 40], [27, 41], [27, 42], [27, 43], [27, 44], [27, 45], [27, 46], [28, 32], [28, 33], [28, 34], [28, 35], [28, 36], [28, 37], [28, 38], [28, 39], [28, 40], [28, 41], [28, 42], [28, 43], [28, 44], [28, 45], [28, 46], [29, 32], [29, 33], [29, 34], [29, 35], [29, 36], [29, 37], [29, 38], [29, 39], [29, 40], [29, 41], [29, 42], [29, 43], [29, 44], [29, 45], [29, 46], [30, 32], [30, 33], [30, 34], [30, 35], [30, 36], [30, 37], [30, 38], [30, 39], [30, 40], [30, 41], [30, 42], [30, 43], [30, 44], [30, 45], [30, 46], [31, 32], [31, 33], [31, 34], [31, 35], [31, 36], [31, 37], [31, 38], [31, 39], [31, 40], [31, 41], [31, 42], [31, 43], [31, 44], [31, 45], [31, 46], [32, 32], [32, 33], [32, 34], [32, 35], [32, 36], [32, 37], [32, 38], [32, 39], [32, 40], [32, 41], [32, 42], [32, 43], [32, 44], [32, 45], [32, 46], [33, 32], [33, 33], [33, 34], [33, 35], [33, 36], [33, 37], [33, 38], [33, 39], [33, 40], [33, 41], [33, 42], [33, 43], [33, 44], [33, 45], [33, 46], [34, 32], [34, 33], [34, 34], [34, 35], [34, 36], [34, 37], [34, 38], [34, 39], [34, 40], [34, 41], [34, 42], [34, 43], [34, 44], [34, 45], [34, 46], [35, 32], [35, 33], [35, 34], [35, 35], [35, 36], [35, 37], [35, 38], [35, 39], [35, 40], [35, 41], [35, 42], [35, 43], [35, 44], [35, 45], [35, 46], [36, 32], [36, 33], [36, 34], [36, 35], [36, 36], [36, 37], [36, 38], [36, 39], [36, 40], [36, 41], [36, 42], [36, 43], [36, 44], [36, 45], [36, 46], [37, 32], [37, 33], [37, 34], [37, 35], [37, 36], [37, 37], [37, 38], [37, 39], [37, 40], [37, 41], [37, 42], [37, 43], [37, 44], [37, 45], [37, 46], [38, 32], [38, 33], [38, 34], [38, 35], [38, 36], [38, 37], [38, 38], [38, 39], [38, 40], [38, 41], [38, 42], [38, 43], [38, 44], [38, 45], [38, 46], [39, 32], [39, 33], [39, 34], [39, 35], [39, 36], [39, 37], [39, 38], [39, 39], [39, 40], [39, 41], [39, 42], [39, 43], [39, 44], [39, 45], [39, 46], [40, 32], [40, 33], [40, 34], [40, 35], [40, 36], [40, 37], [40, 38], [40, 39], [40, 40], [40, 41], [40, 42], [40, 43], [40, 44], [40, 45], [40, 46], [41, 32], [41, 33], [41, 34], [41, 35], [41, 36], [41, 37], [41, 38], [41, 39], [41, 40], [41, 41], [41, 42], [41, 43], [41, 44], [41, 45], [41, 46], [42, 32], [42, 33], [42, 34], [42, 35], [42, 36], [42, 37], [42, 38], [42, 39], [42, 40], [42, 41], [42, 42], [42, 43], [42, 44], [42, 45], [42, 46], [43, 32], [43, 33], [43, 34], [43, 35], [43, 36], [43, 37], [43, 38], [43, 39], [43, 40], [43, 41], [43, 42], [43, 43], [43, 44], [43, 45], [43, 46], [43, 47], [20, 47], [20, 48], [21, 47], [21, 48], [22, 47], [22, 48], [23, 47], [23, 48], [24, 47], [24, 48], [25, 47], [25, 48], [26, 47], [26, 48], [27, 47], [27, 48], [28, 47], [28, 48], [29, 47], [29, 48], [30, 47], [30, 48], [31, 47], [31, 48], [32, 47], [32, 48], [33, 47], [33, 48], [34, 47], [34, 48], [35, 47], [35, 48], [36, 47], [36, 48], [37, 47], [37, 48], [38, 47], [38, 48], [39, 47], [39, 48], [40, 47], [40, 48], [41, 47], [41, 48], [42, 47], [42, 48], [43, 48]])

seed_patch_square = Patch([[16, 14], [16, 15], [16, 16], [16, 17], [16, 18], [16, 19], [16, 20], [16, 21], [16, 22], [16, 23], [16, 24], [16, 25], [16, 26], [16, 27], [16, 28], [16, 29], [16, 30], [16, 31], [16, 32], [17, 14], [17, 15], [17, 16], [17, 17], [17, 18], [17, 19], [17, 20], [17, 21], [17, 22], [17, 23], [17, 24], [17, 25], [17, 26], [17, 27], [17, 28], [17, 29], [17, 30], [17, 31], [17, 32], [18, 14], [18, 15], [18, 16], [18, 17], [18, 18], [18, 19], [18, 20], [18, 21], [18, 22], [18, 23], [18, 24], [18, 25], [18, 26], [18, 27], [18, 28], [18, 29], [18, 30], [18, 31], [18, 32], [19, 14], [19, 15], [19, 16], [19, 17], [19, 18], [19, 19], [19, 20], [19, 21], [19, 22], [19, 23], [19, 24], [19, 25], [19, 26], [19, 27], [19, 28], [19, 29], [19, 30], [19, 31], [19, 32], [20, 14], [20, 15], [20, 16], [20, 17], [20, 18], [20, 19], [20, 20], [20, 21], [20, 22], [20, 23], [20, 24], [20, 25], [20, 26], [20, 27], [20, 28], [20, 29], [20, 30], [20, 31], [20, 32], [21, 14], [21, 15], [21, 16], [21, 17], [21, 18], [21, 19], [21, 20], [21, 21], [21, 22], [21, 23], [21, 24], [21, 25], [21, 26], [21, 27], [21, 28], [21, 29], [21, 30], [21, 31], [21, 32], [22, 14], [22, 15], [22, 16], [22, 17], [22, 18], [22, 19], [22, 20], [22, 21], [22, 22], [22, 23], [22, 24], [22, 25], [22, 26], [22, 27], [22, 28], [22, 29], [22, 30], [22, 31], [22, 32], [23, 14], [23, 15], [23, 16], [23, 17], [23, 18], [23, 19], [23, 20], [23, 21], [23, 22], [23, 23], [23, 24], [23, 25], [23, 26], [23, 27], [23, 28], [23, 29], [23, 30], [23, 31], [23, 32], [24, 14], [24, 15], [24, 16], [24, 17], [24, 18], [24, 19], [24, 20], [24, 21], [24, 22], [24, 23], [24, 24], [24, 25], [24, 26], [24, 27], [24, 28], [24, 29], [24, 30], [24, 31], [24, 32], [25, 14], [25, 15], [25, 16], [25, 17], [25, 18], [25, 19], [25, 20], [25, 21], [25, 22], [25, 23], [25, 24], [25, 25], [25, 26], [25, 27], [25, 28], [25, 29], [25, 30], [25, 31], [25, 32], [26, 14], [26, 15], [26, 16], [26, 17], [26, 18], [26, 19], [26, 20], [26, 21], [26, 22], [26, 23], [26, 24], [26, 25], [26, 26], [26, 27], [26, 28], [26, 29], [26, 30], [26, 31], [26, 32], [27, 14], [27, 15], [27, 16], [27, 17], [27, 18], [27, 19], [27, 20], [27, 21], [27, 22], [27, 23], [27, 24], [27, 25], [27, 26], [27, 27], [27, 28], [27, 29], [27, 30], [27, 31], [27, 32], [28, 14], [28, 15], [28, 16], [28, 17], [28, 18], [28, 19], [28, 20], [28, 21], [28, 22], [28, 23], [28, 24], [28, 25], [28, 26], [28, 27], [28, 28], [28, 29], [28, 30], [28, 31], [28, 32], [29, 14], [29, 15], [29, 16], [29, 17], [29, 18], [29, 19], [29, 20], [29, 21], [29, 22], [29, 23], [29, 24], [29, 25], [29, 26], [29, 27], [29, 28], [29, 29], [29, 30], [29, 31], [29, 32], [30, 14], [30, 15], [30, 16], [30, 17], [30, 18], [30, 19], [30, 20], [30, 21], [30, 22], [30, 23], [30, 24], [30, 25], [30, 26], [30, 27], [30, 28], [30, 29], [30, 30], [30, 31], [30, 32], [31, 14], [31, 15], [31, 16], [31, 17], [31, 18], [31, 19], [31, 20], [31, 21], [31, 22], [31, 23], [31, 24], [31, 25], [31, 26], [31, 27], [31, 28], [31, 29], [31, 30], [31, 31], [31, 32], [32, 14], [32, 15], [32, 16], [32, 17], [32, 18], [32, 19], [32, 20], [32, 21], [32, 22], [32, 23], [32, 24], [32, 25], [32, 26], [32, 27], [32, 28], [32, 29], [32, 30], [32, 31], [32, 32], [33, 14], [33, 15], [33, 16], [33, 17], [33, 18], [33, 19], [33, 20], [33, 21], [33, 22], [33, 23], [33, 24], [33, 25], [33, 26], [33, 27], [33, 28], [33, 29], [33, 30], [33, 31], [33, 32], [34, 14], [34, 15], [34, 16], [34, 17], [34, 18], [34, 19], [34, 20], [34, 21], [34, 22], [34, 23], [34, 24], [34, 25], [34, 26], [34, 27], [34, 28], [34, 29], [34, 30], [34, 31], [34, 32]])


def main():

    archive = []

    population = []
    for _ in range(POP_SIZE):
        # patch = [(random.randint(0, WORLD_X-1), random.randint(0, WORLD_Y-1))\
        #          for _ in range(random.randint(100, WORLD_X*WORLD_Y/2))]
        #
        # patch = list(set(patch)) #remove duplicates
        # patch = [list(cell) for cell in patch]

        cell_set = set()
        for _ in range(5):
            rad = random.randint(6, 30)
            center = [random.randint(0, WORLD_X), random.randint(0, WORLD_Y)]
            for i in range(rad+1):
                for j in range(rad+1):
                    if math.sqrt(i**2 + j**2) < rad:
                        cell_set.add(((center[0]+i)%WORLD_X, (center[1]+j)%WORLD_Y))
                        cell_set.add(((center[0]-i)%WORLD_X, (center[1]-j)%WORLD_Y))
                        cell_set.add(((center[0]+i)%WORLD_X, (center[1]-j)%WORLD_Y))
                        cell_set.add(((center[0]-i)%WORLD_X, (center[1]+j)%WORLD_Y))

        cells = [list(c) for c in cell_set]

        population.append(Patch(cells))


    for _ in range(GENERATIONS):
        print(_)
        for patch in population:
            fitness = sparseness(patch, archive)

            patch.fitness = fitness

            if fitness > THRESHOLD:
                archive.append(patch)

        new_pop = []
        for _ in range(POP_SIZE):
            winner = max(random.sample(population, TOURNAMENT_SIZE),
                         key = lambda a: a.fitness)

            # winner = deepcopy(winner)
            # for i in range(WORLD_X):
            #     for j in range(WORLD_Y):
            #         if random.random() < MUTATION_RATE/10.0:
            #             winner.grid[i,j] = int(not winner.grid[i,j])

            changed = False

            # if random.random() < MUTATION_RATE:
            #     #Add cell
            #     winner = deepcopy(winner)
            #     options = list(ALL_CELLS.difference(set([tuple(c) for c in winner.cells])))
            #     winner.cells.append(list(random.choice(options)))
            #
            #     changed = True
            #
            # if random.random() < MUTATION_RATE:
            #     #Remove cell
            #     winner = deepcopy(winner)
            #     winner.cells.pop(random.choice(range(len(winner.cells))))
            #     changed = True

            if random.random() < MUTATION_RATE:
                #Remove cell
                winner = deepcopy(winner)
                translation = deepcopy(winner.cells)
                dist_x = random.randint(0, WORLD_X)
                dist_y = random.randint(0, WORLD_Y)
                translation = [((c[0]+dist_x)%WORLD_X, (c[1]+dist_y)%WORLD_Y) for c in translation]
                winner.cells = [list(cell) for cell in set([tuple(c) for c in winner.cells]).union(set(translation))]
                changed = True

            if random.random() < MUTATION_RATE:
                winner = deepcopy(winner)
                center = random.choice(winner.cells)
                rad = random.randint(5, 20)
                cell_set = set([tuple(c) for c in winner.cells])
                for i in range(rad+1):
                    for j in range(rad+1):
                        if math.sqrt(i**2 + j**2) < rad:
                            cell_set.add(((center[0]+i)%WORLD_X, (center[1]+j)%WORLD_Y))
                            cell_set.add(((center[0]-i)%WORLD_X, (center[1]-j)%WORLD_Y))
                            cell_set.add(((center[0]+i)%WORLD_X, (center[1]-j)%WORLD_Y))
                            cell_set.add(((center[0]-i)%WORLD_X, (center[1]+j)%WORLD_Y))
                winner.cells = [list(c) for c in cell_set]
                changed = True

            if random.random() < MUTATION_RATE/2.0:
                winner = deepcopy(winner)
                center = random.choice(winner.cells)
                rad = random.randint(5, 20)
                for i in range(rad+1):
                    for j in range(rad+1):
                        if math.sqrt(i**2 + j**2) < rad:
                            up_right = ((center[0]+i)%WORLD_X, (center[1]+j)%WORLD_Y)
                            low_right = ((center[0]+i)%WORLD_X, (center[1]-j)%WORLD_Y)
                            low_left = ((center[0]-i)%WORLD_X, (center[1]-j)%WORLD_Y)
                            up_left = ((center[0]-i)%WORLD_X, (center[1]+j)%WORLD_Y)
                            if up_right in winner.cells:
                                winner.cells.remove(up_right)
                            if low_right in winner.cells:
                                winner.cells.remove(low_right)
                            if up_left in winner.cells:
                                winner.cells.remove(up_left)
                            if low_left in winner.cells:
                                winner.cells.remove(low_left)
                changed = True

            if random.random() < MUTATION_RATE/2.0:
                winner = deepcopy(winner)
                center = random.choice(winner.cells)
                rad = random.randint(5, 20)
                for i in range(rad+1):
                    for j in range(rad+1):
                        if i**2 >= j**2:
                            up_right = ((center[0]+i)%WORLD_X, (center[1]+j)%WORLD_Y)
                            low_right = ((center[0]+i)%WORLD_X, (center[1]-j)%WORLD_Y)
                            low_left = ((center[0]-i)%WORLD_X, (center[1]-j)%WORLD_Y)
                            up_left = ((center[0]-i)%WORLD_X, (center[1]+j)%WORLD_Y)
                            if up_right in winner.cells:
                                winner.cells.remove(up_right)
                            if low_right in winner.cells:
                                winner.cells.remove(low_right)
                            if up_left in winner.cells:
                                winner.cells.remove(up_left)
                            if low_left in winner.cells:
                                winner.cells.remove(low_left)
                changed = True


            if random.random() < MUTATION_RATE/2.0:
                winner = deepcopy(winner)
                center = random.choice(winner.cells)
                rad = random.randint(5, 20)
                for i in range(rad+1):
                    for j in range(rad+1):
                        up_right = ((center[0]+i)%WORLD_X, (center[1]+j)%WORLD_Y)
                        low_right = ((center[0]+i)%WORLD_X, (center[1]-j)%WORLD_Y)
                        low_left = ((center[0]-i)%WORLD_X, (center[1]-j)%WORLD_Y)
                        up_left = ((center[0]-i)%WORLD_X, (center[1]+j)%WORLD_Y)
                        if up_right in winner.cells:
                            winner.cells.remove(up_right)
                        if low_right in winner.cells:
                            winner.cells.remove(low_right)
                        if up_left in winner.cells:
                            winner.cells.remove(up_left)
                        if low_left in winner.cells:
                            winner.cells.remove(low_left)
                changed = True

            if random.random() < MUTATION_RATE:
                winner = deepcopy(winner)
                center = random.choice(winner.cells)
                rad = random.randint(5, 20)
                cell_set = set([tuple(c) for c in winner.cells])
                for i in range(rad+1):
                    for j in range(rad+1):
                        cell_set.add(((center[0]+i)%WORLD_X, (center[1]+j)%WORLD_Y))
                        cell_set.add(((center[0]-i)%WORLD_X, (center[1]-j)%WORLD_Y))
                        cell_set.add(((center[0]+i)%WORLD_X, (center[1]-j)%WORLD_Y))
                        cell_set.add(((center[0]-i)%WORLD_X, (center[1]+j)%WORLD_Y))
                winner.cells = [list(c) for c in cell_set]
                changed = True

            if random.random() < MUTATION_RATE:
                winner = deepcopy(winner)
                center = random.choice(winner.cells)
                rad = random.randint(5, 20)
                cell_set = set([tuple(c) for c in winner.cells])
                for i in range(rad+1):
                    for j in range(rad+1):
                        if i**2 >= j**2:
                            cell_set.add(((center[0]+i)%WORLD_X, (center[1]+j)%WORLD_Y))
                            cell_set.add(((center[0]+i)%WORLD_X, (center[1]-j)%WORLD_Y))
                            cell_set.add(((center[0]-i)%WORLD_X, (center[1]+j)%WORLD_Y))
                winner.cells = [list(c) for c in cell_set]
                changed = True

            if random.random() < CROSSOVER_PROB:
                other_parent = max(random.sample(population, TOURNAMENT_SIZE),
                         key = lambda a: a.fitness)
                winner = deepcopy(winner)
                winner.cells += other_parent.patch
                changed = True

            if changed:
                winner.calc_stats()
            new_pop.append(winner)


        population = new_pop
    for el in archive:
        print(el.stats)

    write_results(archive)


main()
