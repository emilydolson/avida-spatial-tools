from math import sqrt, pi, floor, ceil
from .utils import *
import sys
import glob
from collections import deque
import numpy as np
from scipy.spatial import ConvexHull


def area(patch):
    return len(patch)


def perimeter(patch, world_size=60):
    """
    Count cell faces in patch that do not connect to part of patch.
    This preserves various square geometry features that would not
    be preserved by merely counting the number of cells that touch
    an edge.
    """
    edge = 0
    patch = set([tuple(i) for i in patch])
    for cell in patch:
        neighbors = [tuple(i) for i in get_rook_neighbors(cell, world_size)]
        neighbors = [n for n in neighbors if n not in patch]
        edge += len(neighbors)

    return edge


def get_edge_locations(patch, world_size=60):

    edge = []  # list of cells on edge
    patch = set([tuple(i) for i in patch])

    for cell in patch:
        if isedge(cell, patch, world_size):
            edge.append(cell)

    return edge


def isedge(cell, patch, world_size=60):
    neighbors = [tuple(i) for i in get_rook_neighbors(cell, world_size)]
    neighbors = [n for n in neighbors if n not in patch]
    return bool(neighbors)


def get_rook_neighbors(cell, world_size=60):
    neighbors = [[cell[0]+1, cell[1]], [cell[0]-1, cell[1]],
                 [cell[0], cell[1]+1], [cell[0], cell[1]-1]]

    for neighbor in neighbors:
        for j in range(2):
            if neighbor[j] < 0:
                neighbor[j] += world_size
            elif neighbor[j] >= world_size:
                neighbor[j] -= world_size

    return neighbors


def get_moore_neighbors(cell, world_size=60):
    neighbors = [[x % world_size, y % world_size]
                 for x in range(cell[0]-1, cell[0]+2)
                 for y in range(cell[1]-1, cell[1]+2)]
    neighbors.remove(cell)
    return neighbors


def weighted_perimeter(patch):
    edge = 0
    for cell in patch:
        neighbors = get_moore_neighbors(cell)
        neighbors = [n for n in neighbors if n in patch]
        edge += (8.0 - len(neighbors))/8.0
    return edge


def centroid(patch):
    if len(patch) == 0:
        return (0, 0)

    x_avg = float(sum([cell[0] for cell in patch]))
    y_avg = float(sum([cell[1] for cell in patch]))

    x_avg /= float(len(patch))
    y_avg /= float(len(patch))

    return (x_avg, y_avg)


def radius_of_gyration(patch):
    if len(patch) == 0:
        return 0
    cent = centroid(patch)
    dist_sum = float(sum([dist(cell, cent) for cell in patch]))
    return dist_sum/len(patch)


def perimeter_area_ratio(patch):
    return float(perimeter(patch))/float(area(patch))


def shape_index(patch, perim=None):
    if perim is None:
        perim = float(perimeter(patch))
    patch_area = float(area(patch))
    if patch_area == 0:
        return 0

    # print "perim", perim, "area", patch_area, "ratio", perim/patch_area

    return (.25*perim)/sqrt(patch_area)


def fractal_dimension(patch, perim=None):
    if perim is None:
        perim = float(perimeter(patch))
    perim *= .25
    patch_area = float(area(patch))

    if patch_area == 1:
        # this a square, but will produce a divide by 0 error so we manually
        # return 1
        return 1
    elif patch_area == 0:
        return -1
    elif perim == 0:
        return -1

    return 2*log(perim)/log(patch_area)


def related_circumscribing_circle(patch, formula=True):
    """
    Formula (bool): True indicates that the area of the circumscribing circle
    should be calculated as pi*r^2. This is a more perfect circle than any
    circle composed of squares can be. As a result, the circumscribing circle
    loses the property of representing the maximum number of cells that could
    possibly be in the patch - this can produce negative numbers, particularly
    for small patches. Setting this value to False will calculate a
    circumscribing circle by directly counting grid cells that would be
    included in the circle. This retains the property of never being less than
    patch area and so will never return a value less than 0. However, it
    produces some strange artifacts for small patches and less precisely
    approximates the values reported in the original paper introducing this
    metric (Baker and Cai, 1992). It will also be slightly slower.
    """
    patch_area = float(area(patch))
    max_dist = 0.0
    cell_pair = (None, None)

    try:
        hull = ConvexHull(patch)
        edge = list(np.array(patch)[hull.vertices])
    except:
        edge = patch

    for i, cell1 in enumerate(edge):
        for j, cell2 in enumerate(edge[i+1:]):
            squared_dist = squared_toroidal_dist(cell1, cell2)
            if squared_dist > max_dist:
                max_dist = squared_dist
                cell_pair = (cell1, cell2)

    radius = sqrt(max_dist)/2.0  # only take sqrt once

    if radius == 0:
        # This is a 1-cell patch - manually return 0
        return 0

    if formula:
        return 1-(patch_area/((radius**2)*pi))

    center = ((cell_pair[0][0]+cell_pair[1][0])/2.0,
              ((cell_pair[0][1]+cell_pair[1][1])/2.0))

    # Calculating area of circumscrbing circle
    # by brute force. Turns out that this is the
    # Gauss circle problem, which is solved by an
    # infinite sum, so brute force will be more
    # precise

    circle_area = 0.0
    x_floor = int(floor(center[0]-radius))
    x_ceil = int(ceil(center[0]+radius)+1)

    y_floor = int(floor(center[1]-radius))
    y_ceil = int(ceil(center[1]+radius)+1)

    for x in range(x_floor, x_ceil):
        for y in range(y_floor, y_ceil):
            if dist((x, y), center) <= radius:
                circle_area += 1

    return 1 - (patch_area/circle_area)


def contiguity_index(patch):
    patch_area = float(area(patch))
    if patch_area == 0:
        return 0
    patch = set([tuple(i) for i in patch])
    contiguity = 0.0

    for cell in patch:
        for neighbor in [(cell[0]-1, cell[1]), (cell[0]+1, cell[1]),
                         (cell[0], cell[1]-1), (cell[0], cell[1]+1)]:
            if neighbor in patch:
                contiguity += 2

        for neighbor in [(cell[0]-1, cell[1]-1), (cell[0]+1, cell[1]-1),
                         (cell[0]-1, cell[1]+1), (cell[0]+1, cell[1]+1)]:
            if neighbor in patch:
                contiguity += 1

    # 13 is the maximum value that can be added for any given cell
    return (contiguity/patch_area - 1) / (13-1)


def core_area(patch, distance, world_size=60):
    edge = get_edge_locations(patch, world_size)
    core_area = 0

    for cell in patch:
        dist_to_edge = min([toroidal_dist(cell, other, world_size, world_size)
                            for other in edge])

        if dist_to_edge >= distance:
            core_area += 1
    return core_area


def number_core_areas(patch, distance, world_size=60):
    core_patch = []
    edge = get_edge_locations(patch, world_size)

    for cell in patch:
        dist_to_edge = min([toroidal_dist(cell, other, world_size, world_size)
                            for other in edge])
        if dist_to_edge >= distance:
            core_patch.append(cell)

    return len(traverse_core(core_patch, world_size))


def get_core_areas(patch, distance, world_size=60):
    core_patch = []
    # edge = get_edge_locations(patch, world_size)
    squared_distance = distance*distance

    patch = set([tuple(i) for i in patch])

    edges = set()

    for cell in patch:
        if cell in edges:
            if squared_distance == 0:
                core_patch.append(cell)
            continue

        dist_to_edge = -1
        queue = deque()
        curr = cell
        seen = set([])
        while True:

            if curr in edges:
                dist_to_edge = toroidal_dist(cell, curr,
                                             world_size, world_size)
                break

            else:
                neighbors = [tuple(i) for i in
                             get_moore_neighbors(list(curr), world_size)]

                for n in neighbors:
                    if not (n in patch):
                        dist_to_edge = toroidal_dist(cell, curr,
                                                     world_size, world_size)
                        edges.add(curr)
                        break

                    elif n not in seen:
                        queue.append(n)
                        seen.add(n)

                if dist_to_edge != -1:
                    break

                if len(queue) == 0:
                    dist_to_edge = float("inf")
                    break
                curr = queue.popleft()

        print(dist_to_edge, squared_distance, cell)
        if dist_to_edge >= squared_distance:
            # squared distance is more efficient
            core_patch.append(cell)
        print(core_patch)
    return traverse_core(core_patch, world_size)


def traverse_core(core_area, world_size=60):
    """
    Treat cells in core_area like a graph and traverse it to
    see how many connected components there are.
    """

    if not core_area:
        return []
    core_area = [tuple(i) for i in core_area]
    curr = core_area[0]
    core_area = set(core_area[1:])
    to_explore = []
    cores = [[curr]]

    while core_area:
        print(core_area)
        neighbors = [tuple(i) for i in
                     get_moore_neighbors(list(curr), world_size)]

        for n in neighbors:
            if n in core_area:
                core_area.remove(n)
                to_explore.append(n)
                cores[-1].append(n)

        if to_explore:
            curr = to_explore.pop()
        else:
            curr = core_area.pop()
            cores.append([curr])

    return cores


def core_area_index(patch, distance):
    core = float(core_area(patch, distance))
    patch_area = float(area(patch))
    return core/patch_area


def contrast_weighted_edge_density(patch):
    # For future projects - have a look at other contrast metrics too
    pass


def euclidean_nearest_neighbor(patch, other_patches):
    pass


def proximity_index(patch, other_patches):
    pass


def similarity_index(patch, other_patches):
    pass

if __name__ == "__main__":
    main()
