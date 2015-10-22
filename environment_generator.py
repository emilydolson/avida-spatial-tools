#!/usr/bin/python

#This program allows you to automatically generate environment files with
#complex spatial resource layouts

import argparse
import random
from math import sqrt, log, floor, ceil
from utils import *

parser = argparse.ArgumentParser(description="Generate environment files for environmental heterogeneity experiment.", add_help=False)

modes = parser.add_argument_group("","")

#Ways of running
mode_group = parser.add_mutually_exclusive_group()
mode_group.add_argument("--environmentPicker", action="store_true", help = "Use the environmentPicker GUI to select cells to fill with resources")
mode_group.add_argument("--twoCircles", action="store_true", help="Triggers twocircles mode, which generates a set of environment files designed for the two circle overlaps experiments.")
mode_group.add_argument("--graphMode", action="store_true", help="Do simulation rather than generate environment file.")

#Global settings
parser.add_argument("--worldSize", default=60, type=int, help="Specifies size of one side (assumes square world)")
parser.add_argument("--randomSeed", default=-1, type=int, help="The random seed to use. By default, a random one will be selected.")
parser.add_argument("--outfile", default="environment.cfg", help="The name of the environment file to be created")
parser.add_argument("--inflow", default=100, type=int, help="Sets amount of inflow.")
parser.add_argument("--outflow", default=.01, type=float, help="Sets amount of outflow.")


#Stuff you can add to the environment
parser.add_argument("--resources", default=["resNOT", "resAND", "resOR", "resNOR", "resNAND", "resORN", "resANDN", "resXOR", "resEQU"], nargs="*", type=str, help = "The set of resources to be used. Defaults to logic-9.")
parser.add_argument("--randomPatch", nargs="*", help = "Place a randomly generated patch in the environment. By default, this will include all resources. If additional command-line resources are specified, then only this set will be included in the patch. These resources will not be included with well-mixed resources.")
parser.add_argument("--evenSquares", type=int, help = "Place evenly spaced squares of the specified size in environment")
parser.add_argument("--gradientResources", nargs="*", help = "Designates which resources are to be implemented as gradient resources (circles in space). If this flag is not used, no resources will be gradients. If this flag is used and no additional arguments are passed, all resources will be gradients. Additional arguments can be passed to this flag to specify a subset of the resources to represent as cells.")
parser.add_argument("--cellResources", nargs="*", help = "Designates which resources are to be implemented as cell resources. If this flag is not used, no resources will be cells. If this flag is used and no additional arguments are passed, all resources will be cells. Additional arguments can be passed to this flag to specify a subset of the resources to represent as cells.")

#~~~~~~~~~~Settings for random patches~~~~~~~~~~~~~~~~~~~
parser.add_argument("--patchSize", default=600, type=int, help = "Size of random patch to be generated.")
parser.add_argument("--cellInflow", default=100, type=int, help="Inflow of cell resources")
parser.add_argument("--cellOutflow", default=.01, type=int, help="Outflow of cell resources")
parser.add_argument("--initial", default=10, type=int, help="Initial resource levels of cell resources")

#~~~~~~~~~~Settings for gradient resources~~~~~~~~~~~~~~~

#Choosing anchors:
anchors_group = parser.add_mutually_exclusive_group()
anchors_group.add_argument("--randAnchors", action="store_true", help="Generates random anchor points")
anchors_group.add_argument("--distance", default=None, type=int, help="Distance between anchor points of central patches")
anchors_group.add_argument("--evenAnchors", action="store_true", help="Places gradient resoruces as evenly as possible through environment.")
parser.add_argument("--boundedAnchors", action="store_true", help = "Forces patches to be within world.")

#Chossing radius
radius_group  = parser.add_mutually_exclusive_group()
parser.add_argument("--patchRadius", default=10, type=int, help="Radius of gradient resources")
parser.add_argument("--randRadius", nargs="*", type=int, help = "Sets radius of each patch to a random integer between two specified numebers.")

#Number of patches
parser.add_argument("--randNPatches", action="store_true", help = "Use a random number of patches within 10 of the number specified.")
parser.add_argument("--patchesPerSide", default=4, type=int, help="number of patches in one row or column of world (assume square layout)")

#Other options
parser.add_argument("--notcommon", action="store_true", help = "Included for backwards compatability - makes plateau gradient resources not be common. This should generally not be specified, although there isn't much evidence that it matters.")

#~~~~~~~~~Reactions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parser.add_argument("--reaction", nargs="*", help="A list of resources to make corresponding reactions for. Default if flag is set is to add all.")
parser.add_argument("--infinite", action="store_true", help="Don't deplete resources when reaction are performed.")
parser.add_argument("--frac", default = .0025, type=float, help="The fraction of available resource to use for reactions.")
parser.add_argument("--resMax", default = 25, type=float, help="Maximum units of resource an organism can use per reaction.")
parser.add_argument("--rxnType", default = "pow", help = "How is task reward granted?")
parser.add_argument("--maxCount", default = 1, type=int, help="maximum number of times a task can be completed")


#~~~~~~~Graph Mode Options~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parser.add_argument("--xvar", default="distance", type=str, help = "X variable for graph mode")
parser.add_argument("--yvar", default="entropy", type=str, help = "Y variable for graph mode")
parser.add_argument("--stepSize", default=1, type=int, help = "Step size for graph mode")
parser.add_argument("--reps", default=30, type=int, help = "Replications per step for graph mode")
parser.add_argument("--upLim", default=30, type=int, help = "Upper limit of range for graph mode")


args = parser.parse_args()#[],namespace=args)

taskValDict = {"not":1.0,
               "nand":1.0,
               "and":2.0,
               "orn":2.0,
               "or":3.0,
               "andn":3.0,
               "nor":4.0,
               "xor":4.0,
               "equ":5.0} #rewards for each task

def main():

    if args.randomSeed != -1:
        random.seed(args.randomSeed)

    args.nPatches = args.patchesPerSide**2
    if args.randNPatches:
        args.nPatches += random.randrange(-60,60)

    if args.twoCircles:
        print "Generating", len(args.resources)*len(args.resources), "environment files."
        genTwoCircles()
        print "Success!"
        exit(0)

    if args.graphMode:
        #graphRichnessVsDistance(14, genRandResources())
        #graphRichnessVsRadius(genRandResources(), calcEvenAnchors())
        #entropyHistogramRandomPlacements(genRandResources())
        graphEntropyVsRadius(genRandResources(), calcEvenAnchors())
        exit(0)

    outfile = open(args.outfile, "w")

    #print args

    if args.environmentPicker:
        from cell_picker import cell_picker
        cells = cell_picker(args.worldSize, args.worldSize).run()
        #print args.worldSize
        #cells = ep.cell_picker(60, 60, [(args.worldSize,args.worldSize)]).run()
        cells = list(set(cells))
        cells.sort()
        for resource in args.resources:
            outfile.write(genCell(resource, cells))

        outfile.write("\n")

        for resource in args.resources:
            #print resource
            outfile.write(genReaction(resource))
            
        outfile.close()
        exit(0)

    #No special mode. Just build a standard environment file

    spatial_resources = []
    if args.gradientResources is not None:
        spatial_resources += args.gradientResources if args.gradientResources != [] else args.resources
    if args.cellResources is not None:
        spatial_resources += args.cellResources if args.cellResources != [] else args.resources
    if args.randomPatch is not None:
        spatial_resources += args.randomPatch if args.randomPatch != [] else args.resources

    #Add non-spatial resources:
    well_mixed_resources = args.resources[:]
    for res in spatial_resources:
        if res in well_mixed_resources:
            well_mixed_resources.remove(res)

    for res in well_mixed_resources:
        outfile.write(genRes(res, args.inflow, args.outflow))

    #Add random patch
    if args.randomPatch is not None:
        cells = random_patch(args.patchSize)
        cells.sort()
        for resource in args.resources:
            outfile.write(genCell(resource, cells))

        outfile.write("\n")

        for resource in args.resources if len(args.randomPatch)==0 else args.randomPatch:
            #print resource
            outfile.write(genReaction(resource))
            
    outfile.write("\n")
    
    #Add gradient resources if appropriate
    if args.gradientResources is not None:
        anchors = []
        if args.randAnchors: #random anchor placement
            anchors = calcRandomAnchors(args.boundedAnchors)
        elif args.distance is None: #Evenly spaced anchor points
            anchors = calcEvenAnchors()
        else: #Anchors placed a specific distance from each other
            anchors = calcTightAnchors(args.distance, args.patchesPerSide)

        #Generate random resource order
        randResources = genRandResources(args.gradientResources if args.gradientResources != [] else args.resources)
        random.shuffle(randResources)


        #Write confiugration to specifed output file
        for i in range(len(anchors)):
            if args.randRadius is not None:
                radius = random.randrange(args.randRadius[0], args.randRadius[1])
                outfile.write(genGradient(randResources[i], args.inflow, radius, anchors[i],\
                                      common=(not args.notcommon)))
            else:
                outfile.write(genGradient(randResources[i], args.inflow, args.patchRadius, anchors[i], common=(not args.notcommon)))
                rads.append(args.patchRadius)
        
    if args.evenSquares is not None:
        cells = place_even_squares(args.evenSquares)
        for res in args.cellResources:
            outfile.write(genCell(res, cells))

    outfile.write("\n")

    for res in args.resources:
        outfile.write(genReaction(res, not args.infinite))

    outfile.close()

def place_even_squares(side):
    cells = []
    for y in range(0, args.worldSize, side*2):
        for x in range(0, args.worldSize, side*2):
            for z in range(side):
                for z2 in range(side):
                    cells.append((y+z)*args.worldSize + x + z2)
    cells.sort()
    return cells

def random_patch(size):
    start_point = [random.randrange(0,args.worldSize), 
                   random.randrange(0,args.worldSize)]
    curr_patch = [start_point]
    perimeter = [start_point]
    while len(curr_patch) < size:
        to_expand = random.choice(perimeter)
        neighbors = get_moore_neighbors(to_expand)
        neighbors = [n for n in neighbors if n not in curr_patch]
        if len(neighbors) == 0:
            perimeter.remove(to_expand)
            continue
        elif len(neighbors) == 1:
            perimeter.remove(to_expand)
        next_cell = random.choice(neighbors)
        perimeter.append(next_cell)
        curr_patch.append(next_cell)

    index_cells = []
    for cell in curr_patch:
        index_cells.append(cell[0]%args.worldSize + cell[1]*args.worldSize)

    return index_cells

def get_moore_neighbors(cell):
    neighbors = []
    for x in range(cell[0]-1, cell[0]+2):
        for y in range(cell[1]-1, cell[1]+2):
            neighbors.append([x,y])
            for j in range(2):
                if neighbors[-1][j] < 0:
                    neighbors[-1][j] += args.worldSize
                elif neighbors[-1][j] >= args.worldSize:
                    neighbors[-1][j] -= args.worldSize
        
    neighbors.remove(cell)
    return neighbors

def calcEvenAnchors():
    """
    Calculates anchor points evenly spaced across the world, given 
    user-specified parameters.
    """
    anchors = []
    dist = (args.worldSize+1)/(args.patchesPerSide+1)
    for i in range(dist-1, args.worldSize, dist):
        for j in range(dist-1, args.worldSize, dist):
            anchors.append((i,j))
    return anchors

def genRandResources(resources):
    """
    Generates a list of the appropriate length containing a roughly equal 
    number of all resources in a random order
    """
    randResources = []
    nEach = args.nPatches / len(resources)
    extras = args.nPatches % len(resources)
    for i in range(nEach):
        for res in args.resources:
            randResources.append(res + str(i))

    additional = random.sample(resources, extras)
    for res in additional:
        randResources.append(res + str(nEach))

    random.shuffle(randResources)
    return randResources

def genTwoCircles():
    """
    Generates all 64 potential 2 circles comparison environments.
    """
    if args.distance == None:
        print "Error: genTwoCircles requires a distance"
        return
    
    centerPoint = (int(args.worldSize/2), int(args.worldSize/2))
    anchors = [(centerPoint[0]-int(floor(args.distance/2.0)), centerPoint[1]),
               (centerPoint[0]+int(ceil(args.distance/2.0)), centerPoint[1])]
    for res1 in arg.resources:
        for res2 in args.resources:
            outfile = open(res1.strip("res").lower()+"X"+
                res2.strip("res").lower()+str(args.distance)+".cfg", "w")
            outfile.write(genGradient(res1+"1", args.inflow, 
                                args.patchRadius, anchors[0]))
            outfile.write(genGradient(res2+"2", args.inflow, args.patchRadius,
                                      anchors[1]))
            outfile.write("\n")
            
            outfile.write(genReaction(res1+"1", not args.infinite))
            outfile.write(genReaction(res2+"2", not args.infinite))

            outfile.close()            

def calcRandomAnchors(inworld=True):
    """
    Generates a list of random anchor points such that all circles will fit 
    in the world, given the specified radius and worldsize.
    The number of anchors to generate is determined by squaring the specified 
    number of patches per side.
    """
    anchors = []
    rng = (args.patchRadius, args.worldSize - args.patchRadius)
    if not inworld:
        rng = (0, args.worldSize)
    for i in range(args.nPatches):
        anchors.append((random.randrange(rng[0], rng[1]), 
                        random.randrange(rng[0], rng[1])))

    return anchors

def calcTightAnchors(d, patches):
    """
    Recursively generates the number of anchor points specified in the 
    patches argument, such that all patches are d cells away
    from their nearest neighbors.
    """
    centerPoint = (int(args.worldSize/2), int(args.worldSize/2))
    anchors = []
    if patches == 0:
        return anchors
    elif patches == 1:
        anchors.append(centerPoint)
        return anchors
   
    elif patches == 2:
        xs = [centerPoint[0]+int(floor(d/2)), centerPoint[0]-int(ceil(d/2))]
        ys = [centerPoint[1]+int(floor(d/2)), centerPoint[1]-int(ceil(d/2))]
        for i in xs:
            for j in ys:
                anchors.append((i,j))
        return anchors
    
    elif patches == 3:
        xs = [centerPoint[0] + int(d/2), centerPoint[0], 
              centerPoint[0] - int(d/2)]
        ys = [centerPoint[1] + int(d/2), centerPoint[1], 
              centerPoint[1] - int(d/2)]
        for i in xs:
            for j in ys:
                if not(i==centerPoint[0] and j==centerpoint[0]):
                    anchors.append((i,j))
        return anchors
    
    elif patches%2 == 0:
        dsout = (patches-2)/2 + 1

        xs = [centerPoint[0] + int(floor(d/2.0))+d*i for i in range(dsout)]\
             + [centerPoint[0] - (int(ceil(d/2.0))+d*i) for i in range(dsout)]
        ys = [centerPoint[1] + int(floor(d/2.0))+d*i for i in range(dsout)]\
             + [centerPoint[1] - (int(ceil(d/2.0))+d*i) for i in range(dsout)]

        for i in xs:
            anchors.append((i, max(ys)))
            anchors.append((i, min(ys)))
        for i in ys:
            anchors.append((max(xs), i))
            anchors.append((min(xs), i))

        if d != 0:
            anchors = list(set(anchors))
        anchors.sort()
        return (anchors + calcTightAnchors(d, patches-2))[:patches*patches] 
        #to cut off the extras in the case where d=0

    else:
        #Note - an odd number of args.patchesPerSide requires that there be 
        #a patch at the centerpoint
        dsout = (patches-1)/2
        xs = [centerPoint[0] + d*i for i in range(dsout)] + [centerPoint[0] \
                                        - d*i for i in range(dsout)]
        ys = [centerPoint[1] + d*i for i in range(dsout)] + [centerPoint[1] \
                                        - d*i for i in range(dsout)]
        for i in xs:
            anchors.append((i, max(ys)))
            anchors.append((i, min(ys)))
        for i in ys:
            anchors.append((max(xs), i))
            anchors.append((min(xs), i))

        return anchors + calcTightAnchors(d, patches-2)


def genGradient(resource, inflow, radius, loc, common=True):
    """
    Returns a line of text to add to an environment file, initializing a 
    gradient resource with the specified 
    name (string), inflow(int), radius(int), and location (tuple of ints)
    """
    return "GRADIENT_RESOURCE " + str(resource) + ":height=" + str(radius) + \
        ":plateau=" + str(inflow) +":spread=" + str(radius-1) + ":common=" + \
        str(int(common)) + ":updatestep=1000000:peakx=" + str(loc[0]) + \
        ":peaky=" + str(loc[1]) + ":plateau_inflow=" + str(inflow) + \
                  ":initial=" + str(inflow) + "\n"

def genRes(resource, inflow, outflow):
    """
    Returns a line of text to add to an environment file, initializing a 
    standard resource with the specified name (string), inflow(int), and 
    outflow(int)
    """
    return "RESOURCE " + resource + ":inflow=" + str(inflow) + \
                ":outflow=" + str(outflow) + "\n"

def genCell(resource, cells):
    return "CELL " + resource + ":" + ",".join([str(i) for i in cells]) \
        + ":inflow=" + str(args.cellInflow) + ":outflow=" + str(args.cellOutflow) + ":initial=" + str(args.inflow) + "\n"

def genReaction(resource, depletable=0):
    """
    Returns a line of text to add to an environment file, initializing a 
    reaction that uses the resource specified in the first
    argument to perform the associated task (resource names are expected to 
    be of the form "resTASK#" where "TASK" corresponds
    to the task the resource is associated with and # is an integer uniquely 
    identifying that specific gradient resource. For 
    example, the first AND resource would be named resAND0). An optional 
    second argument (int) specifies whether or not the reaction
    should deplete the resource (by default it will not).
    """
    task = resource.lower()
    if task[:3] == "res":
        task = task[3:]
    while task[-1].isdigit():
        task = task[:-1]

    name = resource[3:]
    return "REACTION " + name + " " + task + " process:resource=" + \
            resource + ":value=" + str(taskValDict[task]) + ":type=" \
            + args.rxnType + ":frac=" + str(args.frac) + ":max=" + str(args.resMax) + \
            ":depletable=" + str(int(depletable)) + " requisite:max_count=" \
            + str(args.maxCount) + "\n"

if __name__ == "__main__":
    main()
