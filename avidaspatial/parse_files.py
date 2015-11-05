#This file contains functions for parsing Avida environment files and spatial
#data output files.

import re
from utils import *
from copy import deepcopy
from environment_file import EnvironmentFile

def load_grid_data(file_list, data_type="binary", sort=True):
    """
    Loads data from one or multiple grid_task files.

    Arguments:
        file_list - either a string or a list of strings indicating files to
                    load data from. Files are assumed to be in grid_task.dat 
                    format (space delimited values, one per cell).

        data_type - a string representing what type of data is in the file.
                    Either "binary", "int", "float", or "string".

        sort      - If you're making a movie, you want the files to be in
                    chronological order. By default, they will be sorted.
                    If for some reason you don't want them in chronological
                    order, set sort to False.

    Returns: A three-dimensional array. The first dimension is columns, the
    second is rows. At each row,column index in the array is another list
    which holds the values that each of the requested files has at that location
    in the grid. If you want this list collapsed to a single representative
    number, you should use agg_niche_grid.
    """

    #If there's only one file, we pretend it's a list
    if not type(file_list) is list:
        file_list = [file_list]
    elif sort:
        #put file_list in chronological order
        file_list.sort(key=lambda f: int(re.sub("[^0-9]", "", f)))

    world_size = get_world_dimensions(file_list[0])

    #Initialize empty data array
    data = initialize_grid(world_size, [])

    #Loop through file list, reading in data
    for f in file_list:
        infile = open(f)
        lines = infile.readlines()
        for i in range(world_size[1]):
            lines[i] = lines[i].split()
            for j in range(world_size[0]):
                if data_type == "binary":
                    val = bin(int(lines[i][j]))
                elif data_type == "float":
                    val = float(lines[i][j])
                elif data_type == "int":
                    val = int(lines[i][j])
                elif data_type == "string":
                    val = str(lines[i][j])
                else:
                    print "Unsupported data_type passed to load_grid"
                    return
                data[i][j].append(val)

        infile.close()

    return data


def make_niche_grid(res_dict, world_size=(60,60)):
    """
    Converts dictionary specifying where resources are to nested lists
    specifying what sets of resources are where.

    res_dict - a dictionary in which keys are resources in the environment
    and values are list of tuples representing the cells they're in.

    world_size - a tuple indicating the dimensions of the world.
           Default = 60x60, because that's the default Avida world size

    Returns a list of lists of sets indicating the set of resources
    available at each x,y location in the Avida grid.
    """

    #Initialize array to represent world
    world = initialize_grid(world_size, set())

    #Fill in data on niches present in each cell of the world
    for res in res_dict:
        for cell in res_dict[res]:
            world[cell[1]][cell[0]].add(res)
 
    return world

def parse_environment_file_list(names, world_size=(60,60)):
    """
    Extract information about spatial resources from all environment files in
    a list.

    Arguments: 
    names - a list of strings representing the paths to the environment files.
    world_size - a tuple representing the x and y coordinates of the world.
                 (default: 60x60)

    Returns a dictionary in which the keys are filenames and the values are
    list of lists of sets indicating the set of resources
    available at each x,y location in the Avida grid for that environment.
    """

    #Convert single file to list if necessary
    try:
        names[0] = names[0]
    except:
        names = [names]

    envs = []
    for name in names:
        envs.append(parse_environment_file(name, world_size))

    return envs

def reduce_resource_name_to_task(res_name):
    """
    Assuming that the convention of naming resources associated with tasks as
    res[TASK][number], reduces such resource names to just the name of the task.
    This ensures that multiple copies of the same resource are treated the
    same. Resource names of different formats will be left untouched.
    """
    #Reduce resource names to tasks being rewarded
    if res_name[:3].lower() != "res":
        return res_name
    res_name = res_name[3:].lower()
    while res_name[-1].isdigit():
        res_name = res_name[:-1]
    return res_name

def parse_environment_file(filename, world_size=(60,60)):
    """
    Extract information about spatial resources from an environment file.

    Arguments: 
    filename - a string representing the path to the environment file.
    world_size - a tuple representing the x and y coordinates of the world.
                 (default: 60x60)

    Returns a list of lists of sets indicating the set of resources
    available at each x,y location in the Avida grid.
    """

    infile = open(filename)
    lines = infile.readlines()
    infile.close()

    tasks = []

    #Find all spatial resources and record which cells they're in
    res_dict = {}
    for line in lines:
        if line.startswith("GRADIENT_RESOURCE"):
            name, cells = parse_gradient(line, world_size)
        elif line.startswith("CELL"):
            name, cells = parse_cell(line, world_size)
        elif line.startswith("REACTION"):
            task = parse_reaction(line)
            if task not in tasks:
                tasks.append(task)
        else:
            continue
        
        dict_increment(res_dict, name, cells)
            
    #Create a map of niches across the environment and return it
    grid = make_niche_grid(res_dict, world_size)    

    return EnvironmentFile(grid, res_dict.keys(), world_size, filename, tasks)

def parse_reaction(line):
    """
    Takes a string declaring a new reaction in an Avida environment file and
    returns the name of the associated task.
    """
    sline = line.split()
    return sline[2].strip()

def parse_gradient(line, world_size):
    """
    Takes a string representing a GRADIENT_RESOURCE (as specified in Avida
    environment files) and a tuple representing the x and y dimensions of the
    world, and returns the name of the gradient resource and a list of
    tuples representing the cells it's in.
    """
    #remove "GRADIENT_RESOURCE"
    line = line[18:]

    sline = [el.strip() for el in line.split(":")]
    name = sline[0]
    name = reduce_resource_name_to_task(name)
    radius = None
    x = None
    y = None

    #Extract data
    for item in sline:
        if item.startswith("height"):
            radius = int(item.split("=")[1])
        elif item.startswith("peakx"):
            x = int(item.split("=")[1])
        elif item.startswith("peaky"):
            y = int(item.split("=")[1])
            
    #Translate circle to cells)
    cells = []
    for i in range(world_size[0]):
        for j in range(world_size[1]):
            if (dist((i,j), (x,y))) <= radius-1:
                cells.append((i,j))

    return (name, cells)

def parse_cell(line, world_size):
    """
    Takes a string representing a CELL resource (as specified in Avida
    environment files) and a tuple representing the x and y dimensions of the
    world, and returns the name of the resource and a list of
    tuples representing the cells it's in.
    """
    #Remove "CELL "
    line = line[4:]

    #Extract information
    sline = [i.strip() for i in line.split(":")]
    name = sline[0]
    name = reduce_resource_name_to_task(name)
    cell_components = sline[1].split(",")
    cells = []

    #List all cells
    for i in range(len(cell_components)):
        if ".." in cell_components[i]:
            cell_range = [int(j) for j in cell_components[i].split("..")]
            for j in range(int(cell_range[0]), int(cell_range[1]+1)):
                cells.append(j)
        else:
            cells.append(int(cell_components[i]))

    if cells == [""]:
        return (name, [])

    xy_pairs = [(int(c)%world_size[0], int(c)//world_size[0]) for c in cells]

    return (name, xy_pairs)
