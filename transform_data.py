from utils import *
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as hierarchicalcluster

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


def cluster_types(types, max_clust=12):
    """
    Generates a dictionary mapping each binary number in types to an integer
    from 0 to max_clust. Hierarchical clustering is used to determine which
    which binary numbers should map to the same integer.
    """
    
    if len(types) < max_clust:
        max_clust = len(types)
 
    #Do actual clustering
    cluster_dict = do_clustering(types, max_clust)
    cluster_ranks = rank_clusters(cluster_dict, types)

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

def make_count_grid(data):
    """
    Load a grid_task style file (specified in filename) and visualize in a
    heat map such that cells with phenotypes that perform more tasks are
    cooler colors.
    """
    data = deepcopy(data)

    for i in range(len(data[0])):
        for j in range(len(data)):
            for k in range(len(data[i][j])):
                try:
                    data[i][j][k] = data[i][j][k].count("1")
                except:
                    data[i][j][k] = len(data[i][j][k])

    return data

def make_optimal_phenotype_grid(environment, phenotypes):
    world_size = environment.world_size
    phenotypes = deepcopy(phenotypes)

    for i in range(world_size[1]):
        for j in range(world_size[0]):
            for k in range(len(phenotypes[i][j])):
                phenotype = phenotype_to_res_set(phenotypes[i][j][k])
                diff = len(world[i][j].symmetric_difference(phenotype))
                phenotypes[i][j][k] = diff

    return phenotypes

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