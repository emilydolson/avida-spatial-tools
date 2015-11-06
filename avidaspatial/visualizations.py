import random, glob, re, string
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.animation
from utils import *
from transform_data import *
import seaborn as sns


def heat_map(grid, name, **kwargs):
    """
    Generic function for making a heat map based on the values in a grid.
    
    Arguments: grid - the grid of numbers or binary strings to be visualized.
               name - string indicating what the file storing the image
                      should be called.

    kwargs:
               palette - a seaborn palette (list of RGB values) indicating
                         how to color values. Will be converted to a continuous
                         colormap if necessary
               denom   - the maximum value of numbers in the grid (only used
                         if the grid actually contains numbers). This is used
                         to normalize values and use the full dynamic range of
                         the color pallete.
    """
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
    an imshow() plot indicating the resources present in each cell. By default,
    color is determined using the palettes in the EnvironmentFile object
    passed as the first parameter. The easiest way to change color palettes
    is to assign new palettes to environment.task_palette and 
    environment.resource_palette before calling this function. If either the
    environment or phenotypes grids contain integers greater than 1, you should
    pass a `denom` keyword argument indicating how to normalize them. Using
    differnet denoms for the environment and phenotypes is not currently
    supported (if you need to, you should probably just divide everything by
    the appropraite denoms before passing them to this funciton).

    Inputs:
          environment - an EnvironmentFile object indicatng the distribution
                         of resources and the appropriate palettes to use.
          phenotypes  - a 2d array of numbers or binary strings representing
                        the placement of phenotypes across the environment

    kwargs:
          denom - an integer indicating how to normalize numbers in the
                  environment and phenotype grids if neccesary.
    Outputs:
         Returns a matplotlib animation object. 
         Saves animation in the file: 
            [environment_file_identifier]_phenotype_overlay.mp4
    """
    denom, palette = get_kwargs(environment, kwargs)

    #Create figure to do plotting
    fig = plt.figure(figsize=(20,20))

    #Create list of circles at every place in environment
    patches = []

    for i in range(len(phenotypes)):
        for j in range(len(phenotypes[i])):
            patches.append(plt.Circle((j,i), radius=.3, \
                            lw=2, ec="black", facecolor=None, zorder=2))
        
    #This will be called to color niches, which are always in background
    def init():
        plot_world(environment, palette = environment.resource_palette, 
                   denom = denom)
        for p in patches:
            fig.gca().add_patch(p)

    #Change colors of circles as appropriate for new time step
    def animate(n):
        phen_grid = slice_3d_grid(phenotypes, n)
        #Recolor circles
        plot_phens_blits(phen_grid, patches, 
                         palette = environment.task_palette, denom = denom)
        return patches,

    #Do actual animation
    anim = matplotlib.animation.FuncAnimation(
        fig, animate, init_func=init, 
        frames=len(phenotypes[0][0]), blit=True, interval=750)
    
    anim.save(environment.name + "_phenotype_overlay.mov")
    return anim

def plot_phens(phen_grid, **kwargs):
    """
    Plots circles colored according to the values in phen_grid.

    -1 serves as a sentinel value, indicating that a circle should not be
    plotted in that location.
    """
    denom, palette = get_kwargs(phen_grid, kwargs, True)
     
    grid = color_grid(phen_grid, palette, denom)
    
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if tuple(grid[i][j]) != -1:
                plt.gca().add_patch(plt.Circle((j,i), \
                radius=.3, lw=1, ec="black", facecolor=grid[i][j], zorder=2))

def plot_phens_circles(phen_grid, **kwargs):
    """
    Plots phenotypes represented as concentric circles. Each circle represents
    one task that the phenotype can perform, with larger circles representing
    more complex tasks.

    Arguments: phen_grid - a 2D array of strings representing binary numbers
  
    kwargs:
               palette - a seaborn palette (list of RGB values) indicating
                         how to color values. Will be converted to a continuous
                         colormap if necessary
               denom   - the maximum value of numbers in the grid (only used
                         if the grid actually contains numbers). This is used
                         to normalize values and use the full dynamic range of
                         the color pallete.

    TODO: come up with way to represent organisms that don't do any tasks.

    """
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
    """
    A version of plot_phens designed to be used in animations. Takes a 2D array
    of phenotypes and a list of matplotlib patch objects that have already
    been added to the current axes and recolors the patches based on the array.
    """

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
    """
    Addes a heat-map representing the data in world (an EnvironmentFile object)
    to the current plot.

    kwargs:
               palette - a seaborn palette (list of RGB values) indicating
                         how to color values. Will be converted to a continuous
                         colormap if necessary
               denom   - the maximum value of numbers in the grid (only used
                         if the grid actually contains numbers). This is used
                         to normalize values and use the full dynamic range of
                         the color pallete.
    """
    denom, palette = get_kwargs(world, kwargs, False)
    world = color_grid(world, palette, denom, True)
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.tight_layout()
    plt.imshow(world, interpolation="none", hold=True, zorder=1)
    axes = plt.gca()
    axes.autoscale(False)

def paired_environment_phenotype_grid(environment, phenotypes, **kwargs):
    """
    Plots the given environment (EnvironmentFile object) and phenotypes
    (2d array of numbers or binary strings) onto the same image and saves 
    the image based on the name of the environment file. The environment file
    will be represented by coloring square cells, while the phenotypes are
    circles overlaid on top.

    By default, color is determined using the palettes in the EnvironmentFile 
    object passed as the first parameter. The easiest way to change color 
    palettes is to assign new palettes to environment.task_palette and 
    environment.resource_palette before calling this function. If either the
    environment or phenotypes grids contain integers greater than 1, you should
    pass a `denom` keyword argument indicating how to normalize them. Using
    differnet denoms for the environment and phenotypes is not currently
    supported (if you need to, you should probably just divide everything by
    the appropraite denoms before passing them to this funciton).

    Inputs:
          environment - an EnvironmentFile object indicatng the distribution
                         of resources and the appropriate palettes to use.
          phenotypes  - a 2d array of numbers or binary strings representing
                        the placement of phenotypes across the environment

    kwargs:
          denom - an integer indicating how to normalize numbers in the
                  environment and phenotype grids if neccesary.

    
    """

    denom, palette = get_kwargs(environment, kwargs)

    plot_world(environment, palette = environment.resource_palette, denom=denom)
    plot_phens(phenotypes, palette = environment.task_palette, denom=denom)

    plt.savefig("phenotype_niches_" + environment.name + ".png", dpi=2000)

def paired_environment_phenotype_grid_circles(environment, phenotypes,**kwargs):
    """
    Plots the given environment (EnvironmentFile object) and phenotypes
    (2d array of binary strings) onto the same image and saves 
    the image based on the name of the environment file. The environment file
    will be represented by coloring square cells, while the phenotypes are
    represented as concentric circles indicating the set of tasks the organism
    at that location can perform.

    By default, color is determined using the palettes in the EnvironmentFile 
    object passed as the first parameter. The easiest way to change color 
    palettes is to assign new palettes to environment.task_palette and 
    environment.resource_palette before calling this function. If either the
    environment or phenotypes grids contain integers greater than 1, you should
    pass a `denom` keyword argument indicating how to normalize them. Using
    differnet denoms for the environment and phenotypes is not currently
    supported (if you need to, you should probably just divide everything by
    the appropraite denoms before passing them to this funciton).

    Inputs:
          environment - an EnvironmentFile object indicatng the distribution
                         of resources and the appropriate palettes to use.
          phenotypes  - a 2d array of binary strings representing
                        the placement of phenotypes across the environment

    kwargs:
          denom - an integer indicating how to normalize numbers in the
                  environment and phenotype grids if neccesary.

    """
    denom, palette = get_kwargs(environment, kwargs)

    plot_world(environment, palette = environment.resource_palette, denom=denom)
    plot_phens_circles(phenotypes, palette = environment.task_palette)

    plt.savefig("phenotype_niches_circles"+environment.name, dpi=1000)
    
    return plt.gcf()

def color_grid(data, palette, denom=9.0, mask_zeros=True):
    """
    Convert the given data (2d array of numbers or binary strings) to a 2d
    array of RGB or RGBA values which can then be visualized as a heat map.

    Arguments:
    data - 2d array of numbers or binary strings
    palette - a seaborn palette (list of RGB values) indicating how to convert 
              data to colors. Will be converted to a continuous colormap if 
              necessary. This should generally be the length of the longest 
              binary string or the highest possible number
    denom - if the data is composed of numbers rather than binary strings,
            this number will indicate how to normalize the data to [0, 1] should
            it be neccessary.
    mask_zeros - Boolean indicating whether 0s should be colored white rather
                 than the color specified by the palette. -1s always yield
                 -1 so that missing data can be handled appropriately.

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
                rgb = color_array_by_value(data[row][col], palette, denom, \
                                           mask_zeros)
            
            grid[row].append(rgb)
    
    return grid

def color_array_by_value(value, palette, denom, mask_zeros):
    """
    Figure out the appropriate RGB or RGBA color for the given numerical
    value based on the palette, denom, and whether zeros should be masked.
    """
    if value == -1: # sentinel value
        return -1

    if value == 0 and mask_zeros: #This value is masked
        if type(palette) is list:
            return (1, 1, 1)
        return (1, 1, 1, 1)

    if type(palette) is list: #This is a palette
        return palette[value]
    
    #This is continuous data so the palette is actually a colormap
    return palette(float(value)/float(denom))

def color_array_by_hue_mix(value, palette):
    """
    Figure out the appropriate color for a binary string value by averaging
    the colors corresponding the indices of each one that it contains. Makes
    for visualizations that intuitively show patch overlap.
    """
    if int(value, 2) > 0:
        
        #Convert bits to list and reverse order to avoid issues with
        #differing lengths
        int_list = [int(i) for i in list(value[2:])]
        int_list.reverse()

        #since this is a 1D array, we need the zeroth elements
        #of np.nonzero.
        locs = np.nonzero(int_list)[0]
        rgb_vals = [palette[i] for i in locs]
    
        rgb = [0]*len(rgb_vals[0]) #We don't know if it's rgb or rgba    
        for val in rgb_vals:
            for index in range(len(val)):
                rgb[index] += val[index]
            
        for i in range(len(rgb)):
            rgb[i] /= len(locs)

        return tuple(rgb)
        
    if int(value, 2) == 0:
        return (1, 1, 1) if len(palette[0])==3 else (1,1,1,1)

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
    """
    Takes a grid of RGB or RGBA values and a filename to save the figure into.
    Generates a figure by coloring all grid cells appropriately.
    """
    plt.tick_params(labelbottom="off", labeltop="off", labelleft="off", \
            labelright="off", bottom="off", top="off", left="off", right="off")
    plt.imshow(grid, interpolation="nearest", aspect=1, zorder=1)
    plt.tight_layout()
    plt.savefig(name, dpi=1000, bbox_inches="tight")



if __name__ == "__main__":
    main()
