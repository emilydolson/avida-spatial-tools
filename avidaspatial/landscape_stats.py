from .utils import *
from .parse_files import *


def patch_richness(world, world_size=(60,60)):
    niches = {}
    for i in range(world_size[1]):
        for j in range(world_size[0]):

            # use frozensets because they are hashable
            if frozenset(world[i][j]) in niches:
                niches[frozenset(world[i][j])] += 1
            else:
                niches[frozenset(world[i][j])] = 1

    return len(niches.keys())


def calc_environment_entropy(world, world_size=(60, 60),
                             exclude_desert=False):
    """
    Calculate the Shannon entropy of a given environment, treating each niche
    (where niches are defined by regions in which different sets of resources
    are rewarded) as a category. The environment is specified with the
    following inputs:

    world - a list of lists of sets of resources (strings) indicating
            the set of resources in every cell in the world.

    world_size - a tuple indicating the dimensions of the world.
           Default = 60x60, because that's the default Avida world siz

    excludeDesert - an optional argument which defaults to False. If True is
          specific, niches in which no tasks are rewarded
          will not be considered in the calculation.
    """

    niches = make_niche_dictionary(world, world_size)

    if exclude_desert and frozenset([]) in niches:
        del niches[frozenset([])]

    # Calculate entropy
    return entropy(niches)


def make_niche_dictionary(world, world_size, mode="freq"):
    # loop through world, counting frequency of each niche
    niches = {}
    for i in range(world_size[1]):
        for j in range(world_size[0]):
            # use frozensets because they are hashable
            if frozenset(world[i][j]) in niches:
                if mode == "freq":
                    niches[frozenset(world[i][j])] += 1
                elif mode == "cells":
                    niches[frozenset(world[i][j])].append([i, j])
            else:
                if mode == "freq":
                    niches[frozenset(world[i][j])] = 1
                elif mode == "cells":
                    niches[frozenset(world[i][j])] = [[i, j]]
                else:
                    print("Unrecognized mode for make_niche_dictionary")
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
    data = load_grid_data(filename, "int")
    data = agg_grid(data, mode)
    phenotypes = {}
    for r in data:
        for c in r:
            if c in phenotypes:
                phenotypes[c] += 1
            else:
                phenotypes[c] = 1

    for key in phenotypes.keys():
        phenotypes[key] = sqrt(phenotypes[key])

    return entropy(phenotypes)
