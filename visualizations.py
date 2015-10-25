__author__ = "Emily Dolson"
__copyright__ = "Copyright 2014, Emily Dolson"
__version__ = "0.9"
__maintainer__ = "Emily Dolson"
__email__ = "EmilyLDolson@gmail.com"
__status__ = "Development"

import random, glob, re, string
from matplotlib.collections import PatchCollection
from transform_data import *
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.animation
from parse_files import *

#RES_SET = ["safe"]
RES_SET = ["equ", "xor", "nor", "andn", "or", "orn", "and", "nand", "not"]

hues = [.1, .7]
hues = [.01, .1, .175, .375, .475, .575, .71, .8, .9]
#hues = [.01, .1, .175, .375, .475, .575, .71, .8, .9]
#hues = [0, .075, .175, .2, .425, .575, .01, .5]
#random.shuffle(hues)


def make_visualization(env_files, grid_files, grid_transform, grid_agg, vis_func, name=""):

    phenotypes = load_grid_data(grid_files)
    world_size = (len(phenotypes[0]), len(phenotypes))

    worlds = parse_environment_file_list(env_files, world_size)[0]

    phenotypes = grid_transform(worlds, phenotypes)
    phenotypes = agg_grid(phenotypes, grid_agg)
    
    name = name + "_" + worlds.name + "_" + vis_func.__name__ + ".png" 

    vis_func(phenotypes, name)

def heat_map(grid, name):
    grid = color_grid(grid)
    make_imshow_plot(grid, name)

def optimal_phenotypes(env_file, grid_file, agg=mean):

    make_visualization(env_file, grid_file, make_optimal_phenotype_grid, mean, heat_map, "optimal_phenotypes")

def paired_environment_phenotype_movie(environment, phenotypes, k=15):
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

    #So we don't have to keep initializig new arrays or screw up original data
    phen_grid = deepcopy(phenotypes)
    
    types = get_ranks_for_environment_and_phenotypes(environment, phenotypes)

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
    anim = matplotlib.animation.FuncAnimation(
        fig, animate, init_func=init, 
        frames=len(species_files), blit=True, interval=750)
    
    anim.save(world.name + "_phenotype_overlay.mov")
    return anim

def plot_phens(phen_grid, denom=9):
    
    #assignments, n = assign_ranks_by_cluster(phen_grid, k, types)
    grid = color_grid(phen_grid, k+1, True)
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
            curr_patch = patches[i * len(grid[i]) + j]
            if grid[i][j] == [0,0,0]:
               curr_patch.set_visible(False)
            else:
                curr_patch.set_facecolor(grid[i][j])
                curr_patch.set_visible(True)

    return patches

def plot_world(world, denom=9.0, p=None):

    world = color_grid(assignments, denom, True)
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.tight_layout()
    plt.imshow(world, interpolation="none", hold=True)
    axes = plt.gca()
    axes.autoscale(False)

    if p != None:
        axes.add_collection(p)

def paired_environment_phenotype_grid(environment, phenotypes):

    ranks = get_ranks_for_environment_and_phenotypes(environment, phenotypes)
    environment_assignments, n = assign_ranks_by_cluster(environment, k, ranks)
    phenotype_assignments, n = assign_ranks_by_cluster(phenotypes, k, ranks)

    plot_world(environment_assignments, n)
    plot_phens(phenotype_assignments, n)

    plt.savefig("phenotype_niches_" + environment.name, dpi=1000)

def paired_environment_phenotype_grid_circles(environment, phenotypes):
    plt.gcf().set_size_inches(40,40)

    plot_world(environment)
    plot_phens_circles(phenotypes)

    plt.savefig("phenotype_niches_"+world.name, dpi=500)
    
    return plt.gcf()

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
            arr[0,0,0] = int(not mask_zeros)
            arr[0,0,1] = int(not mask_zeros)
            arr[0,0,2] = int(not mask_zeros)

            if type(data[row][col]) is str:
                arr = color_array_by_hue_mix(data[row][col], arr, denom)
            else:
                arr = color_array_by_value(data[row][col], arr, denom)
            
            rgb = matplotlib.colors.hsv_to_rgb(arr)
            grid[row].append(list(rgb[0][0]))

    return grid

def color_array_by_value(value, arr, denom):
    
    if float(value) > 0:
        arr[0,0,0] = (float(value)/denom)
        arr[0,0,1] = 1
        
    if float(value) >= 0:
        arr[0,0,2] = 1

    return arr

def color_array_by_hue_mix(value, arr, denom):

    if int(value, 2) > 0:
        
        #since this is a 1D array, we need the zeroth elements
        #of np.nonzero.
        locs = np.nonzero(value[2:])[0]
        color = sum([hues[i] for i in locs])/float(len(locs))
        
        arr[0,0,0] = color
        arr[0,0,1] = .9 + .1*((value.count("1"))/denom)
        arr[0,0,2] = .9 + .1*((value.count("1")+2)**2/100.0) 

    else:
        arr[0,0,1] *= value.count("1")/denom
        arr[0,0,2] = int(int(value, 2) == 0)

    return arr

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

                    grid[i*3+k][j*3+l] = list(rgb[0][0])

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
            grid[i][j] = (r, g, b)

    make_imshow_plot(grid, "colorpercentages2")
    return grid


def make_imshow_plot(grid, name):
  plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
  plt.imshow(grid, interpolation="none", aspect=1)
  plt.tight_layout()
  plt.savefig(name, dpi=500, bbox_inches="tight")



if __name__ == "__main__":
    main()
