from utils import *
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as hierarchicalcluster

def rank_environment_and_phenotypes(environment, phenotypes, k=15):
    environment = convert_world_to_phenotype(environment)
    ranks = get_ranks_for_environment_and_phenotypes(environment, phenotypes)
    environment, n = assign_ranks_by_cluster(environment, k, ranks)
    phenotypes, n = assign_ranks_by_cluster(phenotypes, k, ranks)
    return environment, phenotypes, n


def do_clustering(types, max_clust):
    #Fill in leading zeros to make all numbers same length.
    ls = [list(t[t.find("b")+1:]) for t in types]
    prepend_zeros_to_lists(ls)

    dist_matrix = pdist(ls, weighted_hamming)
    clusters = hierarchicalcluster.complete(dist_matrix)
    clusters = hierarchicalcluster.fcluster(clusters, max_clust, \
                                            criterion="maxclust")

    #Group members of each cluster together
    cluster_dict = dict((c, []) for c in set(clusters))
    for i in range(len(types)):
        cluster_dict[clusters[i]].append(types[i])

    return cluster_dict

def rank_clusters(cluster_dict, types):
    #Figure out the relative rank of each cluster
    cluster_ranks = dict.fromkeys(cluster_dict.keys())
    for key in cluster_dict:
        cluster_ranks[key] = eval(string_avg(cluster_dict[key], binary=True))

    i = len(cluster_ranks)
    for key in sorted(cluster_ranks, key=cluster_ranks.get):
        cluster_ranks[key] = i
        i -= 1

    return cluster_ranks


def get_ranks_for_environment_and_phenotypes(environment, phenotypes, k=15):
    """
    Environment is expected to already have been converted to binary numbers
    (generally because this is being called by rank_environment_and_phenotypes)
    """
    #Create list of all niches and all phenotypes, in phenotype format
    niches = flatten_array(environment)    
    phenotypes = flatten_array(phenotypes)
                
    types = set(phenotypes+niches)

    types.discard("-0b1") #We'll handle this specially
    types.discard("0b0") #We'll handle this specially

    #Do all clustering ahead of time so colors remain consistent.
    ranks = generate_ranks(list(types), k)
    
    ranks["-0b1"] = -1 # The empty phenotype/niche should always be rank 0
    ranks["0b0"] = -1 # The empty phenotype/niche should always be rank 0

    return ranks
    

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
            if type(grid[i][j]) is list:
                for k in range(len(grid[i][j])):
                    assignments[i][j][k] = ranks[grid[i][j][k]]
            else:
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
                if type(data[i][j][k]) is list:
                    for l in range(len(data[i][j][k])):
                        try:
                            data[i][j][k] = data[i][j][k][l].count("1")
                        except:
                            data[i][j][k] = len(data[i][j][k][l])      
                else:
                    try:
                        data[i][j][k] = data[i][j][k].count("1")
                    except:
                        data[i][j][k] = len(data[i][j][k])

    return data

def make_optimal_phenotype_grid(environment, phenotypes):
    world_size = environment.size
    phenotypes = deepcopy(phenotypes)

    for i in range(world_size[1]):
        for j in range(world_size[0]):
            for k in range(len(phenotypes[i][j])):
                phenotype = phenotype_to_res_set(phenotypes[i][j][k], environment.tasks)
                diff = len(environment[i][j].symmetric_difference(phenotype))
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
