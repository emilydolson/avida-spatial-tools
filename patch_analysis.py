from math import sqrt, pi, floor, ceil
from utils import *
from avidaSpatialTools import *
import sys, glob

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
    for cell in patch:
        neighbors = get_rook_neighbors(cell, world_size)
        neighbors = [n for n in neighbors if n not in patch]
        edge += len(neighbors)

    return edge
        
def get_rook_neighbors(cell, world_size=60):
    neighbors = [ [cell[0]+1, cell[1]], [cell[0]-1, cell[1]], 
                  [cell[0], cell[1]+1], [cell[0], cell[1]-1] ]
    
    for neighbor in neighbors:
        for j in range(2):
            if neighbor[j] < 0:
                neighbor[j] += world_size
            elif neighbor[j] >= world_size:
                neighbor[j] -= world_size
        
    return neighbors

def get_moore_neighbors(cell, world_size=60):
    neighbors = []
    for x in range(cell[0]-1, cell[0]+2):
        for y in range(cell[1]-1, cell[1]+2):
            neighbors.append([x,y])
            for j in range(2):
                if neighbors[-1][j] < 0:
                    neighbors[-1][j] += world_size
                elif neighbors[-1][j] >= world_size:
                    neighbors[-1][j] -= world_size
        
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
    x_avg = float(sum([cell[0] for cell in patch]))
    y_avg = float(sum([cell[1] for cell in patch]))

    x_avg /= float(len(patch))
    y_avg /= float(len(patch))

    return (x_avg, y_avg)
    

def radius_of_gyration(patch):
    cent = centroid(patch)
    dist_sum = float(sum([dist(cell, cent) for cell in patch]))

    return dist_sum/len(patch)
    
def perimeter_area_ratio(patch):
    return float(perimeter(patch))/float(area(patch))

def shape_index(patch):
    perim = float(perimeter(patch))
    patch_area = float(area(patch))

    print "perim", perim, "area", patch_area, "ratio", perim/patch_area

    return (.25*perim)/sqrt(patch_area)
    
def fractal_dimension(patch):
    perim = .25 * float(perimeter(patch))
    patch_area = float(area(patch))

    if patch_area == 1:
        #this a square, but will produce a divide by 0 error so we manually
        #return 1
        return 1
    
    return 2*log(perim)/log(patch_area)

def related_circumscribing_circle(patch, formula=True):
    """
    Formula (bool): True indicates that the area of the circumscribing circle
    should be calculated as pi*r^2. This is a more perfect circle than any
    circle composed of squares can be. As a result, the circumscribing circle
    loses the property of representing the maximum number of cells that could
    possibly be in the patch - this can produce negative numbers, particularly
    for small patches. Setting this value to False will calculate a 
    circumscribing circle by directly counting grid cells that would be included
    in the circle. This retains the property of never being less than patch area
    and so will never return a value less than 0. However, it produces some
    strange artifacts for small patches and less precisely approximates the
    values reported in the original paper introducing this metric (Baker and 
    Cai, 1992). It will also be slightly slower.
    """
    patch_area = float(area(patch))
    max_dist = 0.0
    cell_pair = (None,None)
    for cell1 in patch:
        for cell2 in patch:
            max_dist = max(dist(cell1, cell2), max_dist)
            cell_pair = (cell1, cell2)

    radius = max_dist/2.0

    if radius == 0:
        #This is a 1-cell patch - manually return 0
        return 0

    if formula:
        return 1-(patch_area/((radius**2)*pi))

    center = ((cell_pair[0][0]+cell_pair[1][0])/2.0, ((cell_pair[0][1]+cell_pair[1][1])/2.0))

    #Calculating area of circumscrbing circle
    #by brute force. Turns out that this is the
    #Gauss circle problem, which is solved by an
    #infinite sum, so brute force will be more
    #precise

    circle_area = 0.0

    for x in range(int(floor(center[0]-radius)), int(ceil(center[0]+radius))+1):
        for y in range(int(floor(center[1]-radius)), int(ceil(center[1]+radius))+1):
            if dist((x,y),center) <= radius:
                circle_area += 1

    return 1 - (patch_area/circle_area)

def contiguity_index(patch):
    patch_area = float(area(patch))
    contiguity = 0.0

    for cell in patch:
        for neighbor in [[cell[0]-1,cell[1]], [cell[0]+1,cell[1]], 
                         [cell[0], cell[1]-1], [cell[0], cell[1]+1]]:
            if neighbor in patch:
                contiguity += 2

        for neighbor in [[cell[0]-1,cell[1]-1], [cell[0]+1,cell[1]-1], 
                         [cell[0]-1, cell[1]+1], [cell[0]+1, cell[1]+1]]:
            if neighbor in patch:
                contiguity += 1

    #13 is the maximum value that can be added for any given cell
    return (contiguity/patch_area - 1) / (13-1)

def core_area(patch, distance):
    pass

def number_core_areas(patch, distance):
    pass

def core_area_index(patch, distance):
    core = float(core_area(patch, distance))
    patch_area = float(area(patch))
    return core/patch_area

def contrast_weighted_edge_density(patch):
    #For future projects - have a look at other contrast metrics too
    pass

def euclidean_nearest_neighbor(patch, other_patches):
    pass

def proximity_index(patch, other_patches):
    pass

def similarity_index(patch, other_patches):
    pass

if __name__ == "__main__":
    main()
