def patch_richness(res_dict, world_size=60):
    world = make_niche_grid(res_dict, world_size)
    niches = {}
    for i in range(world_size):
        for j in range(world_size):

            #use frozensets because they are hashable
            if niches.has_key(frozenset(world[i][j])):
                niches[frozenset(world[i][j])] += 1
            else:
                niches[frozenset(world[i][j])] = 1

    return len(niches.keys())

def calc_environment_entropy(res_dict, world_size = (60,60), exclude_desert=False):
    """
    Calculate the Shannon entropy of a given environment, treating each niche
    (where niches are defined by regions in which different sets of resources
    are rewarded) as a category. The environment is specified with the
    following inputs:

    res_dict - a dictionary in which keys are resources in the environment
    and values are list of tuples representing the cells they're in.

    world_size - a tuple indicating the dimensions of the world.
           Default = 60x60, because that's the default Avida world siz

    excludeDesert - an optional argument which defaults to False. If True is
          specific, niches in which no tasks are rewarded
          will not be considered in the calculation.
    """

    #Initialize list of list of sets to record which niches are where
    world = make_niche_grid(res_dict, world_size)

    niches = make_niche_dictionary(world, world_size)

    if exclude_desert and niches.has_key(frozenset([])):
        del niches[frozenset([])]

    #Calculate entropy
    return entropy(niches)

def make_niche_dictionary(world, world_size, mode="freq"):
    #loop through world, counting frequency of each niche
    niches = {}
    for i in range(world_size):
        for j in range(world_size):
            #use frozensets because they are hashable
            if niches.has_key(frozenset(world[i][j])):
                if mode == "freq":
                    niches[frozenset(world[i][j])] += 1
                elif mode == "cells":
                    niches[frozenset(world[i][j])].append([i,j])
            else:
                if mode == "freq":
                    niches[frozenset(world[i][j])] = 1
                elif mode == "cells":
                    niches[frozenset(world[i][j])] = [[i,j]]
                else:
                    print "Unrecognized mode for make_niche_dictionary"
                    return
                    
    return niches


def entropy(dictionary):
    """
    Helper function for entropy calculations.
    Takes a frequency dictionary and calculates entropy of the keys.
    """
    total = 0.0
    entropy = 0
    for key in dictionary.keys():
        total += dictionary[key]

    for key in dictionary.keys():
        entropy += dictionary[key]/total * log(1.0/(dictionary[key]/total), 2)
    return entropy

def sqrt_shannon_entropy(filename):
    """
    Calculates Shannon entropy based on square root of phenotype count.
    This might account for relationship between population size and
    evolvability.
    """
    data = load_grid_data(filename, "raw")
    phenotypes = {}
    for r in data:
        for c in r:
            if phenotypes.has_key(c):
                phenotypes[c] += 1
            else:
                phenotypes[c] = 1

    for key in phenotypes.keys():
        phenotypes[key] = sqrt(phenotypes[key])

    return entropy(phenotypes)

