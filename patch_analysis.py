from math import sqrt, pi, floor, ceil
from avidaSpatialTools import *
import sys, glob

#mapping of environment files to conditions because I planned badly
awful_dict = {"600cellsquare_env.cfg": "600square",
              "600_cell_star_thick_center_env.cfg": "600thickcenterstar",
              "star_env.cfg": "star",
              "300cellsquare_env.cfg": "300square",
              "300cellline_env.cfg": "300line",
              "devolab_env.cfg": "devolab",
              "420cellline_env.cfg":"420line",
              "square_test_env.cfg": None,
              "180cellline_env.cfg":"180line",
              "2cell_disperse_env.cfg":"2disperse",
              "540cellline_env.cfg": "540line",
              "480cellline_env.cfg":"480line",
              "3_20_patch_env.cfg":"3by20patch",
              "180cell_circle_env.cfg":"180circle",
              "480cellsquare_env.cfg":"480square",
              "4_15_patch_env.cfg":"4by15patch",
              "1cell_disperse_env.cfg":"1disperse",
              "120cellsquare_env.cfg":"120square",
              "4cell_disperse_env.cfg":"4disperse",
              "180cellsquare_env.cfg":"180square",
              "2_patches_env.cfg":"2patch",
              "oval_env.cfg":"oval",
              "thick_x_env.cfg":"thickx",
              "1_60_patch_env.cfg":"1by60patch",
              "360cellline_env.cfg":"360line",
              "120cellline_env.cfg":"120line",
              "diamond_env.cfg":"diamond",
              "600_snake_env.cfg":"600snake",
              "cross_square_env.cfg":None,
              "x_env.cfg":"xpatch",
              "240cellsquare_env.cfg":"240square",
              "checker_env.cfg":"checker",
              "2_30_patch_env.cfg":"2by30patch",
              "540cellsquare_env.cfg":"540square",
              "600cellline_env.cfg":"600line",
              "360cellsquare_env.cfg":"360square",
              "60cell_circle_env.cfg":"60circle",
              "240cellline_env.cfg":"240line",
              "120cell_circle_env.cfg":"120circle",
              "3cell_disperse_env.cfg":"3disperse",
              "6_10_patch_env.cfg":"6by10patch",
              "420cellsquare_env.cfg":"420square",
              "5_12_patch_env.cfg":"5by12patch",
              "600_cell_filledcircle_env.cfg":"600filledcircle",
              "5cell_disperse_env.cfg":"5disperse",
              "600_cell_star_env.cfg":"600star",
              "empty__patch_env.cfg":None}


for i in range(10):
    filename = "random_10cells_"+str(i)+"_env.cfg"
    awful_dict[filename] = "10random"+str(i)

for i in range(10,20):
    filename = "random_100cells_"+str(i)+"_env.cfg"
    awful_dict[filename] = "100random"+str(i)

for i in range(20,30):
    filename = "random_500cells_"+str(i)+"_env.cfg"
    awful_dict[filename] = "500random"+str(i)

for i in range(30, 40):
    filename = "random_1000cells_"+str(i)+"_env.cfg"
    awful_dict[filename] = "1000random"+str(i)

for i in range(40,50):
    filename = "random_2000cells_"+str(i)+"_env.cfg"
    awful_dict[filename] = "2000random"+str(i)

for i in range(50,60):
    filename = "random_3000cells_"+str(i)+"_env.cfg"
    awful_dict[filename] = "3000random"+str(i)

awful_dict["random_100cells_i_env.cfg"] = None
awful_dict["random_100cells__env.cfg"] = None
awful_dict["all_patch_env.cfg"] = None



def main():
    environments_to_analyze = glob.glob(sys.argv[1])
    print sys.argv[1], environments_to_analyze
    world_size = 60
    outstr = ",".join(["environment", "area", "perimeter", "radius_of_gyration", "perimeter_area_ratio", "shape_index", "fractal_dimension", "related_circumscribing_circle", "contiguity_index"]) + "\n"

    for env in environments_to_analyze:
        envname = env.split("/")[-1]
        if awful_dict[envname] is None:
            continue
        world = parse_environment_file(env, world_size)
        
        niches = make_niche_dictionary(world, world_size, mode="cells")
        del niches[frozenset([])]
        if len(niches) == 0:
            continue
        
        visualize_environment(env, 60, awful_dict[envname])

        patch = niches.values()[0]
        outstr += str(awful_dict[envname]) + ","
        outstr += str(area(patch)) + ","
        outstr += str(perimeter(patch)) + ","
        outstr += str(radius_of_gyration(patch)) + ","
        #outstr += str(centroid(patch)) + ","
        outstr += str(perimeter_area_ratio(patch)) + ","
        outstr += str(shape_index(patch)) + ","
        outstr += str(fractal_dimension(patch)) + ","
        outstr += str(related_circumscribing_circle(patch)) + ","
        outstr += str(contiguity_index(patch)) + "\n"
        
    outfile = open("environment_stats.csv", "w")
    outfile.write(outstr)
    outfile.close()

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



def dist(p1, p2):
    """
    Returns the distance between the two given tuples.
    """
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

if __name__ == "__main__":
    main()
