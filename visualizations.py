__author__ = "Emily Dolson"
__copyright__ = "Copyright 2014, Emily Dolson"
__version__ = "0.9"
__maintainer__ = "Emily Dolson"
__email__ = "EmilyLDolson@gmail.com"
__status__ = "Development"

import random, glob, re, string
import scipy.stats as stats
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as hierarchicalcluster
from matplotlib import animation
from matplotlib.collections import PatchCollection
from pysal.esda.getisord import G_Local
from utils import *
from matplotlib import pyplot as plt
import matplotlib
from parse_files import *

#RES_SET = ["safe"]
RES_SET = ["equ", "xor", "nor", "andn", "or", "orn", "and", "nand", "not"]

hues = [.01, .1, .175, .375, .475, .575, .71, .8, .9]
#hues = [.01, .1, .175, .375, .475, .575, .71, .8, .9]
#hues = [0, .075, .175, .2, .425, .575, .01, .5]
#random.shuffle(hues)

def assign_ranks_by_cluster(grid, n, ranks=None):
    if ranks is None:
        ranks = generate_ranks(grid, n)
    return assign_ranks_to_grid(grid, ranks), len(ranks)

def generate_ranks(grid, n):
    phenotypes = deepcopy(grid)
    if type(phenotypes) is list and type(phenotypes[0]) is list:
        phenotypes = flatten_array(phenotypes)
    
    #Remove duplicates from types
    types = list(frozenset(phenotypes))
    if len(types) < n:
        ranks = rank_types(types)
    else:
        ranks = cluster_types(types, n)
    
    return ranks


def assign_ranks_to_grid(grid, ranks):
    assignments = deepcopy(grid)
    ranks["0b0"] = 0
    ranks["-0b1"] = -1
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            assignments[i][j] = ranks[grid[i][j]]

    return assignments

def prepend_zeros_to_lists(ls):
    longest = max([len(l) for l in ls])

    for i in range(len(ls)):
        while len(ls[i]) < longest:
            ls[i].insert(0, "0")


def cluster_types(types, max_clust=12):
    """
    Generates a dictionary mapping each binary number in types to an integer
    from 0 to max_clust. Hierarchical clustering is used to determine which
    which binary numbers should map to the same integer.
    """
    
    #Fill in leading zeros to make all numbers same length.
    ls = [list(t[2:]) for t in types]
    prepend_zeros_to_lists(ls)

    #Do actual clustering
    dist_matrix = pdist(ls, weighted_hamming)
    clusters = hierarchicalcluster.complete(dist_matrix)
    if len(types) < max_clust:
        max_clust = len(types)
    clusters = hierarchicalcluster.fcluster(clusters, max_clust, \
                                            criterion="maxclust")

    #Group members of each cluster together
    cluster_dict = dict((c, []) for c in set(clusters))
    for i in range(len(types)):
        cluster_dict[clusters[i]].append(types[i])

    #Figure out the relative rank of each cluster
    cluster_ranks = dict.fromkeys(cluster_dict.keys())
    for key in cluster_dict:
        cluster_ranks[key] = eval(string_avg(cluster_dict[key], binary=True))

    i = len(cluster_ranks)
    for key in sorted(cluster_ranks, key=cluster_ranks.get):
        cluster_ranks[key] = i
        i -= 1

    #Create a dictionary mapping binary numbers to indices
    ranks = {}
    for key in cluster_dict:
        for typ in cluster_dict[key]:
            ranks[typ] = cluster_ranks[key]
   
    return ranks

def rank_types(types):
    """
    Takes a list of binary numbers and returns a dictionary mapping each
    binary number to an integer indicating it's rank within the list.
    """
    include_null = '0b0' in types
    sorted_types = deepcopy(types)
    for i in range(len(sorted_types)):
        sorted_types[i] = int(sorted_types[i], 2)
    sorted_types.sort()

    ranks = {}
    for t in types:
        ranks[t] = sorted_types.index(eval(t)) + int(not include_null)

    return ranks




#~~~~~~~~~~~~~~~~~~~~VISUALIZATION FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def make_n_tasks_grid(file_list, agg=mean):
    """
    Load a grid_task style file (specified in filename) and visualize in a
    heat map such that cells with phenotypes that perform more tasks are
    cooler colors.
    """
    data = load_grid_data([file_list])

    for i in range(len(data[0])):
        for j in range(len(data)):
            for k in range(len(data[i][j])):
                data[i][j][k] = data[i][j][k].count("1")

    data = agg_grid(data, agg)

    #TODO: Is "nTasks" option actually necessary?
    return color_grid(data, "nTasks")



def optimal_phenotypes(env_file, grid_file, agg=mean):

    phenotypes = load_grid_data([grid_file])
    world_size = (len(phenotypes[0]), len(phenotypes))
    world = parse_environment_file(env_file, world_size)


    for i in range(world_size[1]):
        for j in range(world_size[0]):
            for k in range(len(phenotypes[i][j])):
                phenotype = phenotype_to_res_set(phenotypes[i][j][k])
                diff = len(world[i][j].symmetric_difference(phenotype))
                phenotypes[i][j][k] = diff

    phenotypes = agg_grid(phenotypes, agg)

    grid = color_grid(phenotypes)
    make_imshow_plot(grid, "test3")
    return grid

def paired_environment_phenotype_movie(species_files, env_file, k=15):
    """
    Makes an animation overlaying colored circles representing phenotypes over
    an imshow() plot indicating the resources present in each cell. Colors
    are determined by clustering all represented combinations of resources/tasks
    (in phenotypes or niches) using complete link UPGMA. Distance is calculated
    using a weighted Hamming distance that prioritizes left-most bits but
    excludes EQU (this is kind of focused on the resource heterogeneity project
    and should probably be changed).

    Inputs:
          species_files - a list of strings indicating the names of all of the
          files describing phenotype locations to be included in the animation.
          These should all be of the grid_tasks format and order should be 
          specified by the number between "grid_tasks" and the file extension 
          (as is created by Avida by default).

          env_file - The name of an Avida environment file (string).

          k (default: 15) - the number of clusters of phenotypes/niches to make.
               
    Outputs:
         Returns a matplotlib animation object. 
         Saves animation in the file: 
            [environment_file_identifier]_phenotype_overlay.mp4
    """

    #Put grid_tasks files in correct order
    species_files.sort(key=lambda f: int(re.sub("[^0-9]", "", f)))
    
    #load data
    data = load_grid_data(species_files)
    world_size = (len(data[0]), len(data)) #extract size from grid_tasks
    world = parse_environment_file(env_file, world_size)

    #Create list of all niches and all phenotypes, in phenotype format
    niches = flatten_array(world)
    
    niches = [res_set_to_phenotype(i, world.resources) for i in niches]
    phenotypes = [phen for col in data for row in col for phen in row]
                
    #TODO: FIX THIS!
    types = set(phenotypes+niches)

    types.discard("-0b1") #We'll handle this specially
    types.discard("0b0") #We'll handle this specially

    #Do all clustering ahead of time so colors remain consistent.
    types = generate_ranks(list(types), k)
    
    types["-0b1"] = -1 # The empty phenotype/niche should always be rank 0
    types["0b0"] = 0 # The empty phenotype/niche should always be rank 0

    #So we don't have to keep initializig new arrays or screw up original data
    phen_grid = deepcopy(data)

    #Create figure to do plotting
    fig = plt.figure(figsize=(20,20))

    #Create list of circles at every place in environment
    patches = []

    for i in range(len(phen_grid)):
        for j in range(len(phen_grid[i])):
            patches.append(plt.Circle((j,i), radius=.3, lw=2, ec="black", facecolor=None))
        
    #This will be called to color niches, which are always in background
    def init():
        plot_world(world, k, types)
        for p in patches:
            fig.gca().add_patch(p)

    #Change colors of circles as appropriate for new time step
    def animate(n):
        print n
        #Load in data from time step n
        for i in range(len(data)):
            for j in range(len(data[i])):
                phen_grid[i][j] = data[i][j][n]

        #Recolor circles
        plot_phens_blits(phen_grid, k, types, patches)

        return patches,

    #Do actual animation
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(species_files), blit=True, interval=750)
    
    anim.save(world.name + "_phenotype_overlay.mov")
    return anim

def plot_phens(phen_grid, k, types):
    
    assignments, n = assign_ranks_by_cluster(phen_grid, k, types)
    grid = color_grid(assignments, k+1, True)
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] != [0,0,0]:
                plt.gca().add_patch(plt.Circle((j,i), \
                            radius=.3, lw=1, ec="black", facecolor=grid[i][j]))

def plot_phens_circles(phen_grid):
    grid = phen_grid
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] != "0b0":
                first = True
                for k in range(2, len(grid[i][j])):
                    if int(grid[i][j][k]) == 1:
                        arr = np.zeros((1,1,3))
                        arr[0,0,0] = hues[k-2]
                        arr[0,0,1] = 1
                        arr[0,0,2] = 1
                        rgb = matplotlib.colors.hsv_to_rgb(arr)
                        rgb = rgb.flatten()
                        plt.gca().add_patch(plt.Circle((j,i), radius=(11-k)*.05, lw=.1 if first else 0, ec="black", facecolor=rgb))
                        first = False

def plot_phens_blits(phen_grid, k, types, patches):
    
    assignments, n = assign_ranks_by_cluster(phen_grid, k, types)
    grid = color_grid(assignments, k+1, True)

    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] == [0,0,0]:
                patches[i * len(grid[i]) + j].set_visible(False)
            else:
                patches[i*len(grid[i])+j].set_facecolor(grid[i][j])
                patches[i*len(grid[i])+j].set_visible(True)

    return patches

def plot_world(world, k, types, p=None):

    #Convert to binary numbers if needed
    try:
        int(world[0][0], 2)
    except:
        conversion_func = function_with_args(res_set_to_phenotype, world.resources)
        grid = agg_grid(deepcopy(world), conversion_func)

    assignments, n = assign_ranks_by_cluster(grid, k, types)
    world = color_grid(assignments, k+1, True)
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.tight_layout()
    plt.imshow(world, interpolation="none", hold=True)
    axes = plt.gca()
    axes.autoscale(False)

    if p != None:
        axes.add_collection(p)

def paired_environment_phenotype_grid(species_files, env_file, agg=mode):

    phen_grid = agg_grid(load_grid_data(species_files), agg)
    world_size = (len(phen_grid[0]), len(phen_grid))
    world = parse_environment_file(env_file, world_size)

    phenotypes = flatten_array(phen_grid)
    niches = deepcopy(world)
    niches = flatten_array(niches)
    niches = [res_set_to_phenotype(i, world.resources) for i in niches]

    types = set(phenotypes+niches)
    types.discard("0b0")
    k = 25
    types = generate_ranks(list(types), k)
    types["0b0"] = 0

    plot_world(world, len(types.keys()), types)
    plot_phens(phen_grid, len(types.keys()), types)

    plt.savefig("phenotype_niches_"+world.name, dpi=1000)

def paired_environment_phenotype_grid_circles(species_files, env_files, agg=mode, name=""):
    plt.gcf().set_size_inches(40,40)
    phen_grid = agg_grid(load_grid_data(species_files), agg)
    world_size = (len(phen_grid[0]), len(phen_grid))
    world = parse_environment_file(env_files, world_size)

    for i in range(world_size[1]):
        for j in range(world_size[0]):
            world[i][j] = res_set_to_phenotype(world[i][j], world.resources)

    #plot_world(world, len(types.keys()), types)
    world_grid = color_by_phenotype(world, 9.0, True)
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.imshow(world_grid, interpolation="none", hold=True)
    #plt.show()
    axes = plt.gca()
    axes.autoscale(False)
    plot_phens_circles(phen_grid)

    if name == "":
        plt.savefig("phenotype_niches_"+world.name, dpi=500)
    else:
        plt.savefig("phenotype_niches_"+name, dpi=500)
    return plt.gcf()

def make_species_grid(file_list, agg=mode, name="speciesgrid"):
    """
    Makes a heat map of the most common phenotypes (by phenotype) in each
    cell in specified in file (string).
    """
    data = agg_grid(load_grid_data(file_list), agg)
    data, k = assign_ranks_by_cluster(data, 27)
    grid = color_grid(data, k, True)
    #grid = color_by_phenotype(data, 9.0, True)
    return make_imshow_plot(grid, name)

def color_grid(data, denom=9.0, mask_zeros = False):
    """
    Loads specified data into a grid to create a heat map of phenotypic
    complexity and location.
    """
    grid = []
    for row in range(len(data)):
        grid.append([])
        for col in range(len(data[row])):
            arr = np.zeros((1,1,3))
            
            if float(data[row][col]) > 0:
                arr[0,0,0] = (float(data[row][col])/denom)
                arr[0,0,1] = 1
                arr[0,0,2] = 1

            elif float(data[row][col]) == 0:
                arr[0,0,0] = int(not mask_zeros)
                arr[0,0,1] = int(not mask_zeros)
                arr[0,0,2] = 1
            else: #-1
                arr[0,0,0] = int(not mask_zeros)
                arr[0,0,1] = int(not mask_zeros)
                arr[0,0,2] = int(not mask_zeros)

            rgb = matplotlib.colors.hsv_to_rgb(arr)
            grid[row].append([rgb[0,0,0], rgb[0,0,1], rgb[0,0,2]])

    return grid

def color_percentages(file_list, file_name="color_percent.png", \
                      intensification_factor=1.2):
    """
    Creates an image in which each cell in the avida grid is represented as
    a square of 9 sub-cells. Each of these 9 sub-cells represents a different
    task, and is colored such that cooler colors represent more complex tasks.
    The saturation of each sub-cell indicates the percentage of grids in the
    given data-set in which the organism in that cell could perform the
    corresponding task.

    Inputs: file_list - list of names of of avida task grid files to be used
            in making figure.

            intensification_factor (default 1.2): A number to multiply
            the percentage of organisms doing a task by in order to increase
            visibility. This can be useful in cases where a lot of the
            percentages are too low to be easily visualized.

    Returns: Grid indicating appropriate color values for images.
    """
    #Load data
    data = task_percentages(load_grid_data(file_list))

    #Initialize grid
    grid = [[]] * len(data)*3
    for i in range(len(grid)):
        grid[i] = [[]]*len(data[0])*3

    #Color grid
    for i in range(len(data[0])):
        for j in range(len(data)):
            for k in range(3): #create grid of sub-cells
                for l in range(3):
                    #build a color in matplotlib's preferred hsv format
                    arr = np.zeros((1, 1, 3))
                    arr[0, 0, 1] = float(data[i][j][k*3 + l]) \
                           *intensification_factor #saturation, based on data
                    arr[0, 0, 0] = (k*3 +l)/9.0 #hue based on task
                    arr[0, 0, 2] = 1 #value is always 1
                    rgb = matplotlib.colors.hsv_to_rgb(arr) #convert to rgb

                    grid[i*3+k][j*3+l] = ([rgb[0,0,0], rgb[0,0,1], rgb[0,0,2]])

    make_imshow_plot(grid, "colorpercentages")

    return grid

def color_percentages2(file_list):
    """
    Super experimental
    """
    data = task_percentages(load_grid_data(file_list))

    grid = [[]] * len(data)
    for i in range(len(grid)):
        grid[i] = [[]]*len(data[0])

    for i in range(len(data)):
        for j in range(len(data[0])):
            r = sum(data[i][j][:3])/3.0
            g = sum(data[i][j][3:6])/3.0
            b = sum(data[i][j][6:])/3.0
            grid[i][j] =(r, g, b)

    make_imshow_plot(grid, "colorpercentages2")
    return grid


def visualize_environment(filename, world_size=(60,60), outfile=""):
    worlds = parse_environment_file_list(filename, world_size)
    niches = []
    for world in worlds:
        temp_niches = [world[i][j] for i in range(len(world)) for j in range(len(world[i]))]
        niches += [res_set_to_phenotype(i, world.resources) for i in temp_niches]

    types = set(niches)
    types.discard("0b0")
    k = 30
    if len(types) > 1:
        types = generate_ranks(list(types), k)
    else:
        types = {list(types)[0]:1}
    types["0b0"] = 0
    
    for world in worlds:
        #grid, d = cluster(world[seed], len(types), types)
        grid = world.grid
        for i in range(len(grid)):
             grid[i] = [res_set_to_phenotype(grid[i][j], world.resources) for j in range(len(grid[i]))]
    
        grid = color_by_phenotype(grid,  9, True)

        print "making plot", "niches_"+ (world.name if outfile == "" else outfile)
        make_imshow_plot(grid, "niches_" + (world.name if outfile == "" else outfile))

    return grid


def make_imshow_plot(grid, name):
  plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
  plt.imshow(grid, interpolation="none", aspect=1)
  plt.tight_layout()
  plt.savefig(name, dpi=500, bbox_inches="tight")

def color_by_phenotype(data, denom=9.0, mask_zeros = False):
    """
    Loads specified data into a grid to create a heat map of phenotypic
    complexity and location.

    TODO: If there are only two colors, make them black and white
    """
    grid = []

    for row in range(len(data)):
        grid.append([])
        for col in range(len(data[row])):
            arr = np.zeros((1,1,3))
            if int(data[row][col], 2) > 0:

                #since this is a 1D array, we need the zeroth elements
                #of np.nonzero.
                locs = np.nonzero(data[row][col][2:])[0]
                color = sum([hues[i] for i in locs])/float(len(locs))
                
                arr[0,0,0] = color
                arr[0,0,1] = .9 + .1*((data[row][col].count("1"))/8.0)
                arr[0,0,2] = .9 + .1*((data[row][col].count("1")+2)**2/100.0) 

            else:
                arr[0,0,0] = int(not mask_zeros)
                arr[0,0,1] = int(not mask_zeros) * data[row][col].count("1")/8.0
                arr[0,0,2] = int(not mask_zeros)

            rgb = matplotlib.colors.hsv_to_rgb(arr)
            grid[row].append([rgb[0,0,0], rgb[0,0,1], rgb[0,0,2]])

    return grid


"""
#NOT WORKING YET
def compute_diversity_gradient_snazzy_gaussian(world, sd=5, recursive=False):
    world_x = len(world[0])
    world_y = len(world)
    if sd == None:
        sd = calculate_optimal_grain(world)
    
    focal = (0,0)
    
    #Intialize phenotype count dicts
    dicts = deepcopy(world)
    for i in range(world_x):
        for j in range(world_y):
            dicts[i][j] = {}

    for i in range(world_x):
        for j in range(world_y):

            #Edge cases
            next_x = floor((i+1)/float(x_regions))
            x_prop = 1
            y_prop = 1
            if i < world_x - 1 and next_x != x:
                x_prop = (i/float(x_regions)-x)/((i+1)/float(x_regions) - i/float(x_regions))
            if j < world_y - 1 and floor((j+1)/float(y_regions)) != y:
                y_prop = (j/float(y_regions)-y)/((j+1)/float(y_regions) - j/float(y_regions))
            
            dict_increment(regions[x][y], world[i][j][0], x_prop*y_prop)

            if x_prop < 1 and y_prop < 1:
                dict_increment(regions[x+1][y+1], world[i][j][0], (1-x_prop)*(1-y_prop))
            if x_prop < 1:
                dict_increment(regions[x+1][y], world[i][j][0], (1-x_prop)*y_prop)
            if y_prop < 1:
                dict_increment(regions[x][y+1], world[i][j][0], x_prop*(1-y_prop))
    
    entropies = []
    for i in range(grain):
        entropies.append([])
        for j in range(grain):
            #print regions[i][j]
            entropies[i].append(entropy(regions[i][j]))

    if not recursive:
        plt.imshow(entropies, interpolation="none", cmap="jet")
        plt.show()
"""

def compute_diversity_gradient(world, grain=None, recursive=False):
    world_x = len(world[0])
    world_y = len(world)
    if grain == None:
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

            #Edge cases
            next_x = floor((i+1)/float(x_regions))
            x_prop = 1
            y_prop = 1
            if i < world_x - 1 and next_x != x:
                x_prop = (i/float(x_regions)-x)/((i+1)/float(x_regions) - i/float(x_regions))
            if j < world_y - 1 and floor((j+1)/float(y_regions)) != y:
                y_prop = (j/float(y_regions)-y)/((j+1)/float(y_regions) - j/float(y_regions))
            
            dict_increment(regions[x][y], world[i][j][0], x_prop*y_prop)

            if x_prop < 1 and y_prop < 1:
                dict_increment(regions[x+1][y+1], world[i][j][0], (1-x_prop)*(1-y_prop))
            if x_prop < 1:
                dict_increment(regions[x+1][y], world[i][j][0], (1-x_prop)*y_prop)
            if y_prop < 1:
                dict_increment(regions[x][y+1], world[i][j][0], x_prop*(1-y_prop))
    
    entropies = []
    for i in range(grain):
        entropies.append([])
        for j in range(grain):
            #print regions[i][j]
            entropies[i].append(entropy(regions[i][j]))

    if not recursive:
        plt.imshow(entropies, interpolation="none", cmap="jet")
        plt.show()

    #getis_ord(entropies)
    return entropies

def calculate_optimal_grain(world, increment=None):
    start = int(len(world)/2.0)
    if increment == None:
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
    print np.array(data)
    w, data = convert_to_pysal(data)
    print G_Local(data, w).p_sim
    #plt.imshow()
    plt.show()

def task_percentages(data, n_tasks=9):
    """
    Calculates the percentage of organisms in each cell (across multiple files)
    that were doing a given task.
    """
    pdata = deepcopy(data)
    for i in range(len(data)):
        for j in range(len(data[0])):
            percentages = [0.0]*n_tasks
            for k in range(len(data[i][j])):
                b_ind = data[i][j][k].find("b")
                for l in range(b_ind+1, len(data[i][j][k])):
                    percentages[l-2] += int(data[i][j][k][l])
            for p in range(len(percentages)):
                percentages[p] /= len(data[i][j])
            pdata[i][j] = percentages

    return pdata


if __name__ == "__main__":
    main()
