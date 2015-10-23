#This file contains functions that are used throuhgout avida-spatial-tools

from math import sqrt, log, floor, ceil
from copy import deepcopy
import pysal
import numpy as np

def flatten_array(grid):
    return [grid[i][j] for i in range(len(grid)) for j in range(len(grid[i]))]

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
    def inner(arg):
        return func(arg, *args)
    return inner

def phenotype_to_res_set(phenotype, resources = ["equ", "xor", "nor", "andn", "or", "orn", "and", "nand", "not"]):

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

def res_set_to_phenotype(res_set, full_list = ["equ", "xor", "nor", "andn", "or", "orn", "and", "nand", "not"]):

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
    ls.sort()
    return ls[int(floor(len(ls)/2.0))]

def string_avg(strings, binary=False):
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

