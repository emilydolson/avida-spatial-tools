#This file contains functions that are used throuhgout avida-spatial-tools

from math import sqrt
from copy import deepcopy


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
    return "0b"+"".join(phenotype)

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

def string_avg(strings):
    """
    Takes a list of strings of equal length and returns a string containing
    the most common value from each index in the string.
    """
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
    infile = open(file_list[0])
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
