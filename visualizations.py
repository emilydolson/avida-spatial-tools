__author__ = "Emily Dolson"
__copyright__ = "Copyright 2014, Emily Dolson"
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Emily Dolson"
__email__ = "EmilyLDolson@gmail.com"
__status__ = "Development"

from math import sqrt, log, floor, ceil
import copy
from copy import deepcopy
import random, glob, re, string
import numpy as np
import scipy.stats as stats
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as hierarchicalcluster
from matplotlib import animation
from matplotlib.collections import PatchCollection
import pysal
from pysal.esda.getisord import G_Local
from utils import *

try:
    import matplotlib
    from matplotlib import pyplot as plt
except:
    print "Matplotlib import failed"
    print "AvidaSpatialTools is running in HPCC Compatability Mode"

RES_SET = ["safe"]
#RES_SET = ["not", "nand", "and", "orn", "or", "andn", "nor", "xor", "equ"]

hues = [.01, .1, .175, .375, .475, .575, .71, .8, .9]
#hues = [.01, .1, .175, .375, .475, .575, .71, .8, .9]
#hues = [0, .075, .175, .2, .425, .575, .01, .5]
#random.shuffle(hues)

def main():
    #data = load_grid_data("/media/emily/hdd/resource-heterogeneity/experiment/inflow100_radius24_commonresources/heterogeneity_replication_11097/grid_task.100000.dat")
    #compute_diversity_gradient(data, 15)
    #test_color_percentages()
    #test_optimal_phenotypes()
    #test_visualize_environment()
    #test_make_species_grid()
    #test_paired_environment_phenotype_grid()
    #paired_environment_phenotype_movie(glob.glob("/home/emily/repos/resource-heterogeneity/experiment/randomized_entropy/heterogeneity_replication_50047/grid_task.*.dat"), "/home/emily/repos/resource-heterogeneity/environmentFiles/env50047.cfg")
    paired_environment_phenotype_movie(glob.glob("/home/emily/hpcc/conservation/round_2_results/4_patches_15_cells_100_killed-pop1_10901/data/grid_task.*.dat"), "/home/emily/hpcc/conservation/configs/conservation-4patches_15each-environment.cfg", 15, ["safe"])

#~~~~~~~~~~~~~~~~~~~~~~AGGREGATION FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def cluster(grid, n, ranks=None):
    grid = deepcopy(grid)
    phenotypes = [grid[i][j] for i in range(len(grid)) for j in range(len(grid[i]))]

    #check if conversion to phenotype is necessary
    try:
        if phenotypes[0][:2] == "0b":
            pass
        else:
            phenotypes = [res_set_to_phenotype(i, RES_SET) for i in phenotypes]
    except:
        phenotypes = [res_set_to_phenotype(i, RES_SET) for i in phenotypes]

    assignments = deepcopy(grid)

    if ranks == None:
        #Remove duplicates from types
        types = list(frozenset(phenotypes))
        if len(types) < n:
            ranks = rank_types(types)
        else:
            ranks = cluster_types(types, n)

    ranks["0b000000000"] = 0
    ranks["0b0"] = 0
    ranks["0b1"] = -1
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            assignments[i][j] = ranks[phenotypes[i*len(grid[i])+j]]

    return assignments, len(ranks.keys())

def cluster_types(types, max_clust=12):
    ls = [list(t[2:]) for t in types]
    dist_matrix = pdist(ls, weighted_hamming)
    clusters = hierarchicalcluster.complete(dist_matrix)
    #hierarchicalcluster.dendrogram(clusters)
    #plt.show()
    if len(types) < max_clust:
        max_clust = len(types)
    clusters = hierarchicalcluster.fcluster(clusters, max_clust, criterion="maxclust")
    #print clusters
    
    cluster_dict = dict((c, []) for c in set(clusters))
    for i in range(len(types)):
        cluster_dict[clusters[i]].append(types[i])

    cluster_ranks = dict.fromkeys(cluster_dict.keys())
    for key in cluster_dict:
        cluster_ranks[key] = eval(string_agg(cluster_dict[key]))

    i = len(cluster_ranks)
    for key in sorted(cluster_ranks, key=cluster_ranks.get):
        cluster_ranks[key] = i
        i -= 1

    ranks = {}
    for key in cluster_dict:
        for typ in cluster_dict[key]:
            ranks[typ] = cluster_ranks[key]
   
    return ranks

def rank_types(types):
    include_null = '0b000000000' in types
    sorted_types = deepcopy(types)
    for i in range(len(sorted_types)):
        sorted_types[i] = eval(sorted_types[i])
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
            data[i][j] = agg(data[i][j])

    #TODO: Is "nTasks" option actually necessary?
    return color_grid(data, "nTasks")

def test_optimal_phenotypes():
    env = "/home/emily/repos/resource-heterogeneity/environmentFiles/env50047.cfg"
    grid = "/home/emily/repos/resource-heterogeneity/experiment/randomized_entropy/heterogeneity_replication_50047/grid_task.100000.dat"
    optimal_phenotypes(env, grid, 59)

def optimal_phenotypes(env_file, grid_file, world_size=60, agg=mean):
    world = parse_environment_file(env_file, world_size)
    phenotypes = load_grid_data([grid_file])

    for i in range(len(world)):
        for j in range(len(world)):
            for k in range(len(phenotypes[i][j])):
                phenotype = phenotype_to_res_set(phenotypes[i][j][k])
                diff = len(world[i][j].symmetric_difference(phenotype))
                if "equ" in phenotype:
                   diff -= 1
                phenotypes[i][j][k] = diff

            phenotypes[i][j] = agg(phenotypes[i][j])

    grid = color_grid(phenotypes)
    make_imshow_plot(grid, "test3")
    return grid

def paired_environment_phenotype_movie(species_files, env_file, k=15, res_list=["equ", "xor", "nor", "andn", "or", "orn", "and", "nand", "not"]):
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
    world_size = int(len(data)) #extract size from grid_tasks
    world_dict = parse_environment_file_list(env_file, world_size)
    seed = world_dict.keys()[0] #extract env file ID
    world = world_dict[seed]
    seed = seed.strip("abcdefghijklmnopqrstuvwxyzF/_-")#TODO: Allow non-digit ID

    #Create list of all niches and all phenotypes, in phenotype format
    niches = [world[i][j] for i in range(len(world)) \
              for j in range(len(world[i]))]
    
    niches = [res_set_to_phenotype(i, res_list) for i in niches]
    phenotypes = [phen for col in data for row in col for phen in row]

    #TODO: FIX THIS!
    if res_list != ["safe"]:
        types = set(phenotypes+niches)
    else:
        types = set(phenotypes)

    types.discard("-0b000000001") #We'll handle this specially
    types.discard("-0b000000000") #We'll handle this specially

    #Do all clustering ahead of time so colors remain consistent.
    types = cluster_types(list(types), k)
    types["-0b000000001"] = -1 # The empty phenotype/niche should always be rank 0
    types["-0b1"] = -1 # The empty phenotype/niche should always be rank 0
    types["0b000000000"] = 0 # The empty phenotype/niche should always be rank 0
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
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(species_files), blit=True)
    anim.save(str(seed) + "_phenotype_overlay.mp4", writer="mencoder", fps=2)

    return anim

def plot_phens(phen_grid, k, types):
    grid = color_grid(cluster(phen_grid, k, types), False, k+1, True)
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] != [0,0,0]:
                plt.gca().add_patch(plt.Circle((j,i), radius=.3, lw=1, ec="black", facecolor=grid[i][j]))

def plot_phens_circles(phen_grid):
    grid = phen_grid
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] != "0b000000000":
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
    grid = color_grid(cluster(phen_grid, k, types)[0], False, k+1, True)
    
    for i in range(len(phen_grid)):
        for j in range(len(phen_grid[i])):
            if grid[i][j] == [0,0,0]:
                patches[i*len(phen_grid)+j].set_visible(False)
            else:
                patches[i*len(phen_grid)+j].set_facecolor(grid[i][j])
                patches[i*len(phen_grid)+j].set_visible(True)
    return patches

def plot_world(world, k, types, p=None):
    world = color_grid(cluster(world, k, types)[0], False, k+1, True)
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.tight_layout()
    plt.imshow(world, interpolation="none", hold=True)
    axes = plt.gca()
    axes.autoscale(False)

    if p != None:
        axes.add_collection(p)

def test_paired_environment_phenotype_grid():
    env = "/media/emily/hdd/resource-heterogeneity/environmentFiles/env11097.cfg"
    grid = "/media/emily/hdd/resource-heterogeneity/experiment/inflow100_radius24_commonresources/heterogeneity_replication_11097/grid_task.100000.dat"
    #grid = "/home/emily/repos/resource-heterogeneity/experiment/randomized_entropy/heterogeneity_replication_31204/grid_task.100000.dat"
    paired_environment_phenotype_grid_circles(grid, env)

def paired_environment_phenotype_grid(species_files, env_files, agg=mode):

    phen_grid = agg_grid(load_grid_data(species_files), agg)
    world_size = len(phen_grid)
    world_dict = parse_environment_file_list(env_files, world_size)
    seed = world_dict.keys()[0]
    world = world_dict[seed]
    seed = seed.strip("abcdefghijklmnopqrstuvwxyzF/_-")

    phenotypes = [phen_grid[i][j] for i in range(len(phen_grid)) for j in range(len(phen_grid[i]))]
    niches = [world[i][j] for i in range(len(world)) for j in range(len(world[i]))]
    niches = [res_set_to_phenotype(i) for i in niches]

    types = set(phenotypes+niches)
    types.discard("0b000000000")
    k = 25
    types = cluster_types(list(types), k)
    types["0b000000000"] = 0

    plot_world(world, len(types.keys()), types)
    plot_phens(phen_grid, len(types.keys()), types)

    plt.savefig("phenotype_niches_"+str(seed), dpi=1000)

def paired_environment_phenotype_grid_circles(species_files, env_files, agg=mode, name=""):
    plt.gcf().set_size_inches(40,40)
    phen_grid = agg_grid(load_grid_data(species_files), agg)
    world_size = len(phen_grid)
    world_dict = parse_environment_file_list(env_files, world_size)
    seed = world_dict.keys()[0]
    world = world_dict[seed]
    seed = seed.strip(string.ascii_letters+"/_-")
    #if seed.contains("/"):
    #    seed = seed.split("/")[-1]

    #phenotypes = [phen_grid[i][j] for i in range(len(phen_grid)) for j in range(len(phen_grid[i]))]
    #niches = [world[i][j] for i in range(len(world)) for j in range(len(world[i]))]
    #niches = [res_set_to_phenotype(i) for i in niches]

    for i in range(world_size):
        for j in range(world_size):
            world[i][j] = res_set_to_phenotype(world[i][j])

    #types = set(niches)
    #types.discard("0b000000000")
    #k = 25
    #types = cluster_types(list(types), k)
    #types["0b000000000"] = 0

    #plot_world(world, len(types.keys()), types)
    world = color_by_phenotype(world, 9.0, True)
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.imshow(world, interpolation="none", hold=True)
    #plt.show()
    axes = plt.gca()
    axes.autoscale(False)
    plot_phens_circles(phen_grid)

    if name == "":
        plt.savefig("phenotype_niches_"+str(seed), dpi=2000)
    else:
        plt.savefig("phenotype_niches_"+name, dpi=2000)
    return plt.gcf()

def test_make_species_grid():
    grid = "/home/emily/repos/resource-heterogeneity/experiment/radius8_distance12_common/heterogeneity_replication_11857/grid_task.100000.dat"
    make_species_grid(grid)

def make_species_grid(file_list, agg=mode, name="speciesgrid"):
    """
    Makes a heat map of the most common phenotypes (by phenotype) in each
    cell in specified in file (string).
    """
    data = agg_grid(load_grid_data(file_list), agg)
    data, k = cluster(data, 27)
    grid = color_grid(data, False, k, True)
    #grid = color_by_phenotype(data, 9.0, True)
    return make_imshow_plot(grid, name)

def color_grid(data, mode="", denom=9.0, mask_zeros = False):
    """
    Loads specified data into a grid to create a heat map of phenotypic
    complexity and location.
    """
    grid = []
    for row in range(len(data)):
        grid.append([])
        for col in range(len(data[row])):
            arr = np.zeros((1,1,3))
            if mode == "nTasks":
                data[row][col] = n_tasks(data[row][col])
            #print data[row][col]
            
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

    for i in range(len(data[0])):
        for j in range(len(data)):
            r = sum(data[i][j][:3])/3.0
            g = sum(data[i][j][3:6])/3.0
            b = sum(data[i][j][6:])/3.0
            grid[i][j] =(r, g, b)

    make_imshow_plot(grid, "colorpercentages2")
    return grid

def test_color_percentages():
    file_list = glob.glob("/home/emily/repos/resource-heterogeneity/experiment/randomAnchors/*/grid_task.100000.dat")[3:4]
    color_percentages2(file_list)

def visualize_environment(filename, world_size=60, outfile=""):
    world = parse_environment_file_list(filename, world_size)
    seeds = world.keys()
    niches = []
    for seed in seeds:
        temp_niches = [world[seed][i][j] for i in range(len(world[seed])) for j in range(len(world[seed][i]))]
        niches += [res_set_to_phenotype(i) for i in temp_niches]

    types = set(niches)
    types.discard("0b000000000")
    k = 30
    two_color = False
    if len(types) > 1:
        types = cluster_types(list(types), k)
    else:
        types = {list(types)[0]:1}
        two_color = True
    types["0b000000000"] = 0
    
    for seed in seeds:
        #grid, d = cluster(world[seed], len(types), types)
        grid = world[seed]
        for i in range(len(grid)):
             grid[i] = [res_set_to_phenotype(grid[i][j]) for j in range(len(grid[i]))]
    
        grid = color_by_phenotype(grid,  9, True, two_color)

        print "making plot", "niches_"+ (str(seed) if outfile == "" else outfile)
        make_imshow_plot(grid, "niches_" + (str(seed) if outfile == "" else outfile))

    return grid

def test_visualize_environment():
    #env = ["distance0", "distance10", "distance21", "distance29"]
    env = ["/media/emily/hdd/resource-heterogeneity/environmentFiles/env50012.cfg", "/media/emily/hdd/resource-heterogeneity/environmentFiles/env50013.cfg"]#, "/media/emily/hdd/resource-heterogeneity/environmentFiles/env50014.cfg", "/media/emily/hdd/resource-heterogeneity/environmentFiles/env50015.cfg"]
    #grid = "/home/emily/repos/resource-heterogeneity/experiment/radius8_distance12_common/heterogeneity_replication_11857/grid_task.100000.dat"
    grid = visualize_environment(env, 59, "test_env.png")

def make_imshow_plot(grid, name):
  plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
  plt.imshow(grid, interpolation="none", aspect=1)
  plt.tight_layout()
  plt.savefig(name, dpi=1000, bbox_inches="tight")

def color_by_phenotype(data, denom=9.0, mask_zeros = False, two_color=False):
    """
    Loads specified data into a grid to create a heat map of phenotypic
    complexity and location.
    """
    grid = []

    for row in range(len(data)):
        grid.append([])
        for col in range(len(data[row])):
            arr = np.zeros((1,1,3))
            if data[row][col] != "0b000000000":
                locs = []
                curr = 1
                while curr >= 1 and curr < len(data[row][col]):
                    locs.append(data[row][col].find("1", curr)-2)
                    if curr-3  == locs[-1]:
                        locs = locs[:-1]
                        curr = -1
                    else:
                        curr = locs[-1] + 3
                if locs[-1] <= -1:
                    locs = locs[:-1]
                    curr = -1
                
                color = sum([hues[i] for i in locs])/float(len(locs))
                
                arr[0,0,0] = color
                arr[0,0,1] = .9 + .1*((data[row][col].count("1"))/8.0)
                if two_color:
                    arr[0,0,2] = 0
                else:
                    arr[0,0,2] = .9 + .1*((data[row][col].count("1")+2)**2/100.0) 

            else:
                arr[0,0,0] = int(not mask_zeros)
                
                if two_color:
                    arr[0,0,1] = 0
                else:
                    arr[0,0,1] = int(not mask_zeros) * data[row][col].count("1")/8.0 
                if two_color:
                    arr[0,0,2] = 1
                else:
                    arr[0,0,2] = int(not mask_zeros)

            rgb = matplotlib.colors.hsv_to_rgb(arr)
            grid[row].append([rgb[0,0,0], rgb[0,0,1], rgb[0,0,2]])

    return grid

#~~~~~~~~~~~~~~~~~~~~~~LANSCAPE-LEVEL CALCULATIONS~~~~~~~~~~~~~~~~~~~~~~#

#NOT WORKING YET
def compute_diversity_gradient_snazzy_gaussian(world, sd=5, recursive=False):
    world_x = len(world[0])
    world_y = len(world)
    if sd == None:
        sd = calculate_optimal_grain(world)
    
    focal = (0,0)
    
    #Intialize phenotype count dicts
    dicts = copy.deepcopy(world)
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

def convert_to_pysal(data):
    w = pysal.lat2W(len(data[0]), len(data))
    data = np.array(data)
    data = np.reshape(data, (len(data)*len(data[0]), 1))
    return w, data

def calc_environment_entropy(res_dict, world_size = (60,60), exclude_desert=False):
    """
    Calculate the Shannon entropy of a given environment, treating each niche
    (where niches are defined by regions in which different sets of resources
    are rewarded) as a category. The environment is specified with the
    following inputs:

    res_dict - a dictionary in which keys are resources in the environment
    and values are list of tuples representing the cells they're in.

    world_size - a tuple indicating the dimensions of the world.
           Default = 60x60, because that's the default Avida world siz

    excludeDesert - an optional argument which defaults to False. If True is
          specific, niches in which no tasks are rewarded
          will not be considered in the calculation.
    """

    #Initialize list of list of sets to record which niches are where
    world = make_niche_grid(res_dict, world_size)

    niches = make_niche_dictionary(world, world_size)

    if exclude_desert and niches.has_key(frozenset([])):
        del niches[frozenset([])]

    #Calculate entropy
    return entropy(niches)

def make_niche_dictionary(world, world_size, mode="freq"):
    #loop through world, counting frequency of each niche
    niches = {}
    for i in range(world_size):
        for j in range(world_size):
            #use frozensets because they are hashable
            if niches.has_key(frozenset(world[i][j])):
                if mode == "freq":
                    niches[frozenset(world[i][j])] += 1
                elif mode == "cells":
                    niches[frozenset(world[i][j])].append([i,j])
            else:
                if mode == "freq":
                    niches[frozenset(world[i][j])] = 1
                elif mode == "cells":
                    niches[frozenset(world[i][j])] = [[i,j]]
                else:
                    print "Unrecognized mode for make_niche_dictionary"
                    return
                    
    return niches

def sqrt_shannon_entropy(filename):
    """
    Calculates Shannon entropy based on square root of phenotype count.
    This might account for relationship between population size and
    evolvability.
    """
    data = load_grid_data(filename, "raw")
    phenotypes = {}
    for r in data:
        for c in r:
            if phenotypes.has_key(c):
                phenotypes[c] += 1
            else:
                phenotypes[c] = 1

    for key in phenotypes.keys():
        phenotypes[key] = sqrt(phenotypes[key])

    return entropy(phenotypes)

#~~~~~~~~~~~~~~~~~~~~~~~METRIC CALCULATION FUNCTIONS~~~~~~~~~~~~~~~~~~~~#


def patch_richness(res_dict, world_size=60):
    world = make_niche_grid(res_dict, world_size)
    niches = {}
    for i in range(world_size):
        for j in range(world_size):

            #use frozensets because they are hashable
            if niches.has_key(frozenset(world[i][j])):
                niches[frozenset(world[i][j])] += 1
            else:
                niches[frozenset(world[i][j])] = 1

    return len(niches.keys())

def entropy(dictionary):
    """
    Helper function for entropy calculations.
    Takes a frequency dictionary and calculates entropy of the keys.
    """
    total = 0.0
    entropy = 0
    for key in dictionary.keys():
        total += dictionary[key]

    for key in dictionary.keys():
        entropy += dictionary[key]/total * log(1.0/(dictionary[key]/total), 2)
    return entropy

def find_edges(res_dict, world_size=60):
    world = make_niche_grid(res_dict, anchors, rad, world_size)
    edge_count = 0
    for i in range(world_size):
        for j in range(world_size):
            if i >= 1:
                if world[i][j] != world[i-1][j]:
                    edge_count += 1
                elif j >= 1 and world[i][j] != world[i-1][j-1]:
                    edge_count += 1
                elif j < args.worldSize - 1 and world[i][j] != world[i-1][j+1]:
                    edge_count += 1
            elif j >= 1:
                if world[i][j] != world[i][j-1]:
                    edge_count += 1
                elif i < world_size - 1 and world[i][j] !=  world[i+1][j-1]:
                    edge_count += 1
            elif i < world_size - 1:
                if world[i][j] != world[i+1][j]:
                    edge_count += 1
                elif j < world_size - 1 and world[i][j] != world[i+1][j+1]:
                    edge_count += 1
            elif j < world_size - 1:
                if world[i][j] != world[i][j+1]:
                    edge_count += 1
    return edge_count


def n_tasks(dec_num):
    """
    Takes a decimal number as input and returns the number of ones in the
    binary representation.
    This translates to the number of tasks being done by an organism with a
    phenotype represented as a decimal number.
    """
    bitstring = ""
    try:
        bitstring = dec_num[2:]
    except:
        bitstring = bin(int(dec_num))[2:] #cut off 0b
    #print bin(int(dec_num)), bitstring
    return bitstring.count("1")


def weighted_hamming(b1, b2):
    assert(len(b1) == len(b2))
    hamming = 0
    for i in range(len(b1)):
        if b1[i] != b2[i]:
            #differences at more significant (leftward) bits are more important
            if i > 0:
                hamming += 1 + 1.0/i
    return hamming

def task_percentages(data):
    pdata = deepcopy(data)
    for i in range(len(data[0])):
        for j in range(len(data)):
            percentages = [0.0]*9
            for k in range(len(data[i][j])):
                for l in range(2, len(data[i][j][k])):
                    percentages[l-2] += int(data[i][j][k][l])
            for p in range(len(percentages)):
                percentages[p] /= len(data[i][j])
            pdata[i][j] = percentages

    return pdata


if __name__ == "__main__":
    main()
