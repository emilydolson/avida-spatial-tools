#This file contains functions that are used throuhgout avida-spatial-tools

from math import sqrt, log, floor, ceil
from copy import deepcopy
import pysal
import numpy as np
from environment_file import *
import seaborn as sns

def get_kwargs(grid, kwargs, phenotypes=False):
    """
    Helper function to figure out what denom and palette to use, based on the
    kwargs and the grid being plotted. The optional (default: false) argument
    indicates whether the grid contains phenotypes, as opposed to resources.
    """
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
        length = get_pallete_length(grid)
        palette = sns.hls_palette(length, s=1)
        denom = length
        
    return denom, palette

def get_pallete_length(grid):
    """
    Takes a 2d grid and figures out how many different elements are in it, so
    that we know how big to make the palette. Also avoids the unfortunate
    red/green palette that results from too few elements.

    Returns int indicating the length the palette should have.
    """
    elements = list(set(flatten_array(grid)))
    length = len(elements)
    
    if type(elements[0]) is str:
        lengths = [len(el) for el in elements if not el.startswith("-")]
        if max(lengths) < 5: #Mixing red and green 
            length += 2 #is not pretty so let's avoid it
    return length

def agg_grid(grid, agg=None):
    """
    Many functions return a 2d list with a complex data type in each cell.
    For instance, grids representing environments have a set of resources, while
    reading in multiple data files at once will yield a list containing the
    values for that cell from each file. In order to visualize these data types
    it is helpful to summarize the more complex data types with a single number.
    For instance, you might want to take the length of a resource set to see how
    many resource types are present. Alternately, you might want to take the
    mode of a list to see the most common phenotype in a cell. 

    This function facilitates this analysis by calling the given aggregation
    function (agg) on each cell of the given grid and returning the result.

    agg - A function indicating how to summarize grid contents. Default: len.
    """
    grid = deepcopy(grid)

    if agg is None:
        if type(grid[0][0]) is list and type(grid[0][0][0]) is str:
            agg = string_avg
        else:
            agg = mode

    for i in range(len(grid)):
        for j in range(len(grid[i])):
            grid[i][j] = agg(grid[i][j])

    return grid

def slice_3d_grid(grid, n):
    """
    Takes a three dimensional array and an integer (n) and returns a 2d array
    containing the Nth value from the 3rd dimension at each location in the 
    grid.
    """
    phen_grid = initialize_grid((len(grid[0]), len(grid)), 0)

    for i in range(len(grid)):
        for j in range(len(grid[i])):
            phen_grid[i][j] = grid[i][j][n]

    return phen_grid

def flatten_array(grid):
    """
    Takes a multi-dimensional array and returns a 1 dimensional array with the
    same contents.
    """
    grid = [grid[i][j] for i in range(len(grid)) for j in range(len(grid[i]))]
    while type(grid[0]) is list:
        grid = flatten_array(grid)
    return grid


def prepend_zeros_to_lists(ls):
    """
    Takes a list of lists and appends 0s to the beggining of each sub_list
    until they are all the same length. Used for sign-extending binary numbers.
    """
    longest = max([len(l) for l in ls])

    for i in range(len(ls)):
        while len(ls[i]) < longest:
            ls[i].insert(0, "0")

def dict_increment(d, key, amount):
    if key in d:
        d[key] += amount
    else:
        d[key] = amount
     
def dist(p1, p2):
    """
    Returns the distance between the two given tuples.
    """
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def function_with_args(func, *args):
    """
    Returns a function that calls a function with the specified arguments.
    The returned function still takes one argument representing the first
    positional argument. 

    This is mostly a helper function for using agg_grid with functions
    requiring more information than the cell contents.
    """
    def inner(arg):
        return func(arg, *args)
    return inner

def convert_world_to_phenotype(world):
    """
    Converts sets indicating the resources present in a single cell to binary
    strings (bit order is based on the order of resources in world.resources).

    TODO: Figure out how to handle relationship between resources and tasks

    Inputs: world - an EnvironmentFile object with a grid of resource sets
    Returns: an EnvironmentFile object with a grid of binary strings
    """
    if set(world.resources) != set(world.tasks):
        print "Warning: world phenotypes don't correspond to phenotypes"
    if set(world.resources).issubset(set(world.tasks)):
        conversion_func = function_with_args(res_set_to_phenotype, world.tasks)
    else:
        conversion_func = function_with_args(res_set_to_phenotype, world.resources)
    grid = agg_grid(deepcopy(world), conversion_func)
    return grid


def phenotype_to_res_set(phenotype, resources):
    """
    Converts a binary string to a set containing the resources indicated by
    the bits in the string.

    Inputs: phenotype - a binary string
            resources - a list of string indicating which resources correspond
                        to which indices of the phenotype

    returns: A set of strings indicating resources
    """
    assert(phenotype[0:2] == "0b")
    phenotype = phenotype[2:]
    #Fill in leading zeroes
    while len(phenotype) < len(resources):
        phenotype = "0" + phenotype

    res_set = set()

    for i in range(len(phenotype)):
        if phenotype[i] == "1":
            res_set.add(resources[i])
    
    assert(phenotype.count("1") == len(res_set))
    return res_set

def res_set_to_phenotype(res_set, full_list):
    """
    Converts a set of strings indicating resources to a binary string where
    the positions of 1s indicate which resources are present.

    Inputs: res_set - a set of strings indicating which resources are present
            full_list - a list of strings indicating all resources which could
                        could be present, and the order in which they should
                        map to bits in the phenotype
    returns: A binary string
    """

    phenotype = len(full_list) * ["0"]

    for i in range(len(full_list)):
        if full_list[i] in res_set:
            phenotype[i] = "1"
   
    assert(phenotype.count("1") == len(res_set))
    
    #Remove uneceesary leading 0s
    while phenotype[0] == "0" and len(phenotype) > 1:
        phenotype = phenotype[1:]

    return "0b"+"".join(phenotype)


def weighted_hamming(b1, b2):
    """
    Hamming distance that emphasizes differences earlier in strings.
    """
    assert(len(b1) == len(b2))
    hamming = 0
    for i in range(len(b1)):
        if b1[i] != b2[i]:
            #differences at more significant (leftward) bits are more important
            if i > 0:
                hamming += 1 + 1.0/i
                #This weighting is completely arbitrary
    return hamming


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

def convert_to_pysal(data):
    """
    Pysal expects a distance matrix, and data formatted in a numpy array.
    This functions takes a data grid and returns those things.
    """
    w = pysal.lat2W(len(data[0]), len(data))
    data = np.array(data)
    data = np.reshape(data, (len(data)*len(data[0]), 1))
    return w, data


#~~~~~~~~~~~~~~~~~~~~~~AGGREGATION FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Provided for easy use with agg_grid

def mode(ls):
    """
    Takes a list as an argument and returns the mode of (most common item in)
    that list.
    """
    return max(set(ls), key=ls.count)

def mean(ls):
    """
    Takes a list and returns the mean.
    """
    return float(sum(ls))/len(ls)

def median(ls):
    """
    Takes a list and returns the median.
    """
    ls = sorted(ls)
    return ls[int(floor(len(ls)/2.0))]

def string_avg(strings, binary=True):
    """
    Takes a list of strings of equal length and returns a string containing
    the most common value from each index in the string.

    Optional argument: binary - a boolean indicating whether or not to treat
    strings as binary numbers (fill in leading zeros if lengths differ).
    """

    if binary: #Assume this is a binary number and fill leading zeros
        strings = deepcopy(strings)
        longest = len(max(strings, key=len))

        for i in range(len(strings)):
            while len(strings[i]) < longest:
                split_string = strings[i].split("b") 
                strings[i] = split_string[0] + "b0" + split_string[1]

    avg = ""
    for i in (range(len(strings[0]))):
        opts = []
        for s in strings:
            opts.append(s[i])
        avg += max(set(opts), key=opts.count)

    return avg

def get_world_dimensions(gridfile):
    """
    This function takes the name of a file in grid_task format and returns
    the dimensions of the world it represents.
    """
    infile = open(gridfile)
    lines = infile.readlines()
    infile.close()
    world_x = len(lines[0].split())
    world_y = len(lines)
    return (world_x, world_y)

def initialize_grid(world_size, inner):
    """
    Creates an empty grid (2d list) with the dimensions specified in
    world_size. Each element is initialized to the inner argument.
    """
    data = []

    for i in range(world_size[1]):
        data.append([])
        for j in range(world_size[0]):
            data[i].append(deepcopy(inner))
    return data

