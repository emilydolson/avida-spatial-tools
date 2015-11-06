from avidaspatial import *
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as hierarchicalcluster

def rank_environment_and_phenotypes(environment, phenotypes, k=15):
    """
    Clusters sets of resources/tasks using a weighted hamming distance such that
    you can have few enough values to give each group of similar things a
    different color. This function is designed for cases when you want to
    color both an environment and a set of phenotypes such that the colors
    corespond to each other.

    Takes an EnvironmentFile object, a 2d array of phenotypes, and, optionally,
    a number indicating the maximum number of clusters (default 15).

    Returns: 
     - An EnvironmentFile in which the grid has been replaced with integers
       indicating which cluster a cell is a member of. Integers are assigned
       such that cells containing more or more complex resources have higher
       numbers.

     - A 2D grid of numbers representing the clusters each phenotype was
       assigned to.

     - An integer representing the total number of clusters.
    """
    environment = convert_world_to_phenotype(environment)
    ranks = get_ranks_for_environment_and_phenotypes(environment, phenotypes)
    environment, n = assign_ranks_by_cluster(environment, k, ranks)
    phenotypes, n = assign_ranks_by_cluster(phenotypes, k, ranks)
    return environment, phenotypes, n


def do_clustering(types, max_clust):
    """
    Helper method for clustering that takes a list of all of the things being
    clustered (which are assumed to be binary numbers represented as strings),
    and an int representing the maximum number of clusters that are allowed.

    Returns: A dictionary mapping cluster ids to lists of numbers that are part
    of that cluster.
    """
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

def rank_clusters(cluster_dict):
    """
    Helper function for clustering that takes a dictionary mapping cluster
    ids to lists of the binary strings that are part of that cluster and returns
    a dictionary mapping cluster ids to integers representing their "rank".
    Ranks provide an ordering for the clusters such that each cluster has
    its own rank, and clusters are ordered from simplest to most complex.
    """
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
    Takes an EnvironmentFile and a 2d array represemtimg phenotypes at each
    location. Optionally also takes an integer indicating the maximum number
    of clusters allowed to be created (default 15).

    Environment is expected to already have been converted to binary numbers
    (generally because this is being called by rank_environment_and_phenotypes)

    Return a dictionary mapping binary strings representing groups of
    resources/tasks that are present/performed in a given cell to integers
    indicating the ranked order of the cluster they're part of.
    """
    #Create list of all niches and all phenotypes, in phenotype format
    niches = flatten_array(environment)    
    phenotypes = flatten_array(phenotypes)
                
    types = set(phenotypes+niches)

    types.discard("-0b1") #We'll handle this specially
    types.discard("0b0") #We'll handle this specially

    #Do all clustering ahead of time so colors remain consistent.
    ranks = generate_ranks(list(types), k)
    
    ranks["-0b1"] = -1 # The empty phenotype/niche should always be rank -1
    ranks["0b0"] = 0 # The empty phenotype/niche should always be rank 0

    return ranks
    

def assign_ranks_by_cluster(grid, n, ranks=None):
    """
    Takes a 2D array representing phenotypes or resource sets across the world,
    and integer rpresenting the maximum number of clusters allowed, and 
    optionally a dictionary indicating the rank of the cluster of each
    phenotype/resource set. If this dictionary is not provided, one will be 
    generated.

    Returns: - A 2d array of numbers indicating the ranks of the clusters
               of the resource set/phenotype in each cell
             - An integer representing the number of clusters created.
    """
    if ranks is None:
        ranks = generate_ranks(grid, n)
    return assign_ranks_to_grid(grid, ranks), len(ranks)

def generate_ranks(grid, n):
    """
    Takes a grid of phenotypes or resource sets representing as strings
    representing binary numbers, and an integer indicating the maximum number
    of clusters to generated.

    Clusters the data in grid into a maximum of n groups, ranks each group by
    the complexity and length of its "average" member, and returns a dictionary
    mapping binary numbers to integers representing the rank of the cluster
    they're part of.
    """
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
    """
    Takes a 2D array of binary numbers represented as strings and a dictionary
    mapping binary strings to integers representing the rank of the cluster
    they belong to, and returns a grid in which each binary number has been
    replaced with the rank of its cluster.
    """
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
    cluster_ranks = rank_clusters(cluster_dict)

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

    This is basically the better alternative to cluster_types, that works
    in that perfect world where we have few enough types to represent each
    as its own color.
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
    Takes a 2 or 3d grid of strings representing binary numbers.

    Returns a grid of the same dimensions in which each binary number has been
    replaced by an integer indicating the number of ones that were in that
    number.
    """
    data = deepcopy(data)

    for i in range(len(data)):        
        for j in range(len(data[i])):
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
    """
    Takes an EnvironmentFile object and a 2d array of phenotypes and returns
    a 2d array in which each location contains an index representing the
    distance between the phenotype in that location and the optimal phenotype
    for that location.

    This is acheived by using the task list in the EnvironmentFile to convert
    the phenotypes to sets of tasks, and comparing them to the sets of
    resources in the environment. So if the environment file that you created
    the EnvironmentFile object from for some reason doesn't contain all of the
    tasks, or doesn't contain them in the right order this won't work. If this
    is the environment file that you used for the run of Avida that generated
    this data, you should be fine.
    """
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
    Takes a 3D array of strings representing binary numbers and calculates 
    the percentage of organisms in each cell (across multiple files)
    that were doing a given task.

    Returns an m x n x n_tasks array indicating the percentages of organisms
    at each location (across the 3rd dimension) that were doing each task.
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
