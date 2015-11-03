import random, glob, re, string
from transform_data import *
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.animation
from parse_files import *
import seaborn as sns


def heat_map(grid, name, **kwargs):
    denom, palette = get_kwargs(grid, kwargs)
    if "mask_zeros" in kwargs:
        mask_zeros = kwargs["mask_zeros"]
    else:
        mask_zeros = False

    grid = color_grid(grid, palette, denom, mask_zeros)
    make_imshow_plot(grid, name)

def paired_environment_phenotype_movie(environment, phenotypes, **kwargs):
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

    #Create figure to do plotting
    fig = plt.figure(figsize=(20,20))

    #Create list of circles at every place in environment
    patches = []

    for i in range(len(phenotypes)):
        for j in range(len(phenotypes[i])):
            patches.append(plt.Circle((j,i), radius=.3, lw=2, ec="black", facecolor=None, zorder=2))
        
    #This will be called to color niches, which are always in background
    def init():
        plot_world(environment, **kwargs)
        for p in patches:
            fig.gca().add_patch(p)

    #Change colors of circles as appropriate for new time step
    def animate(n):
        phen_grid = slice_3d_grid(phenotypes, n)
        #Recolor circles
        plot_phens_blits(phen_grid, patches, **kwargs)

        return patches,

    #Do actual animation
    anim = matplotlib.animation.FuncAnimation(
        fig, animate, init_func=init, 
        frames=len(phenotypes[0][0]), blit=True, interval=750)
    
    anim.save(environment.name + "_phenotype_overlay.mov")
    return anim

def get_kwargs(grid, kwargs, phenotypes=False):
    
    denom = None
    if "denom" in kwargs:
        denom = kwargs["denom"]

    if "palette" in kwargs:
        palette = kwargs["palette"]
        if denom is None:
            denom = len(palette)
    elif "environment" in kwargs or isinstance(grid, EnvironmentFile):
        if "environment" in kwargs:
            env = kwargs["environment"]
        else:
            env = grid

        if phenotypes:
            palette = env.task_palette
            if denom is None:
                denom = len(env.tasks)
        else:
            palette = env.resource_palette
            if denom is None:
                denom = len(env.resources)

    else:
        length = get_pallete_length(elements)
        palette = sns.hls_palette(length, s=1)
        denom = length
    return denom, palette

def plot_phens(phen_grid, **kwargs):

    denom, palette = get_kwargs(phen_grid, kwargs, True)
     
    grid = color_grid(phen_grid, palette, denom)
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if tuple(grid[i][j]) != -1:
                plt.gca().add_patch(plt.Circle((j,i), \
                radius=.3, lw=1, ec="black", facecolor=grid[i][j], zorder=2))

def plot_phens_circles(phen_grid, **kwargs):

    denom, palette = get_kwargs(phen_grid, kwargs, True)

    n_tasks = len(palette)
    grid = phen_grid
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if int(grid[i][j], 2) != 0:
                first = True
                b_ind = grid[i][j].find("b")
                phen = grid[i][j][b_ind+1:]
                for k in range(len(phen)):
                    if int(phen[k]) == 1:
                        plt.gca().add_patch(
                            plt.Circle(
                                (j,i), radius=(n_tasks - k)*.05, \
                                lw=.1 if first else 0, ec="black",\
                                facecolor=palette[k], zorder=2+k))
                        first = False

def plot_phens_blits(phen_grid, patches, **kwargs):

    denom, palette = get_kwargs(phen_grid, kwargs)    
    grid = color_grid(phen_grid, palette, denom)

    for i in range(len(grid)):
        for j in range(len(grid[i])):
            curr_patch = patches[i * len(grid[i]) + j]
            if grid[i][j] == -1:
               curr_patch.set_visible(False)
            else:
                curr_patch.set_facecolor(grid[i][j])
                curr_patch.set_visible(True)

    return patches

def plot_world(world, **kwargs):
    denom, palette = get_kwargs(world, kwargs, False)
    world = color_grid(world, palette, denom, True)
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.tight_layout()
    plt.imshow(world, interpolation="none", hold=True, zorder=1)
    axes = plt.gca()
    axes.autoscale(False)

def paired_environment_phenotype_grid(environment, phenotypes, **kwargs):

    plot_world(environment, **kwargs)
    plot_phens(phenotypes, **kwargs)

    plt.savefig("phenotype_niches_" + environment.name, dpi=1000)

def paired_environment_phenotype_grid_circles(environment, phenotypes,**kwargs):
    plot_world(environment, **kwargs)
    plot_phens_circles(phenotypes)

    plt.savefig("phenotype_niches_circles"+environment.name, dpi=1000)
    
    return plt.gcf()

def color_grid(data, palette, denom=9.0, mask_zeros=True):
    """
    Loads specified data into a grid to create a heat map of phenotypic
    complexity and location.
    """
    grid = []

    try: 
        #If this isn't numeric, don't bother with this block
        float(data[0][0])

        #This is continuous data - we need a colormap rather than palette
        palette = matplotlib.colors.LinearSegmentedColormap.from_list(
                "color_grid", palette)
        palette.set_bad(alpha=0)
    except:
        pass

    for row in range(len(data)):
        grid.append([])
        for col in range(len(data[row])):
            
            if type(data[row][col]) is str:
                rgb = color_array_by_hue_mix(data[row][col], palette)
            else:
                rgb = color_array_by_value(data[row][col], palette, denom, mask_zeros)
            
            grid[row].append(rgb)
    
    return grid

def color_array_by_value(value, palette, denom, mask_zeros):
    if value == -1:
        return -1
    if value == 0 and mask_zeros:
        if type(palette) is list:
            return (1, 1, 1)
        return (1, 1, 1, 1)
    if type(palette) is list:
        return palette[value]
    return palette(float(value)/float(denom))

def color_array_by_hue_mix(value, palette):
    
    if int(value, 2) > 0:
        
        #Convert bits to list and reverse order to avoid issues with
        #differing lengths
        int_list = [int(i) for i in list(value[2:])]
        int_list.reverse()

        #since this is a 1D array, we need the zeroth elements
        #of np.nonzero.
        locs = np.nonzero(int_list)[0]
        rgb_vals = [palette[i] for i in locs]
    
        rgb = [0, 0, 0]    
        for val in rgb_vals:
            for index in range(len(val)):
                rgb[index] += val[index]
            
        for i in range(len(rgb)):
            rgb[i] /= len(locs)

        return tuple(rgb)
        
    if int(value, 2) == 0:
        return (1, 1, 1)

    return -1

def color_percentages(file_list, n_tasks=9, file_name="color_percent.png", \
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
    for i in range(len(data)):
        for j in range(len(data[i])):
            for k in range(3): #create grid of sub-cells
                for l in range(3):
                    if len(data[i][j]) > k*3+l:
                        #build a color in matplotlib's preferred hsv format
                        arr = np.zeros((1, 1, 3))
                        arr[0, 0, 1] = float(data[i][j][k*3 + l]) \
                            *intensification_factor #saturation, based on data
                        arr[0, 0, 0] = (k*3 +l)/9.0 #hue based on task
                        arr[0, 0, 2] = 1 #value is always 1
                        rgb = matplotlib.colors.hsv_to_rgb(arr) #convert to rgb

                        grid[i*3+k][j*3+l] = list(rgb[0][0])
                    else:
                        grid[i*3+k][j*3+l] = (1, 1, 1, 1)

    return make_imshow_plot(grid, "colorpercentages")

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
  plt.imshow(grid, interpolation="nearest", aspect=1, zorder=1)
  plt.tight_layout()
  plt.savefig(name, dpi=500, bbox_inches="tight")



if __name__ == "__main__":
    main()
