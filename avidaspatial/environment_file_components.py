def gen_gradient(args, resource, inflow, radius, loc, common=True):
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

def gen_res(args, resource, inflow, outflow):
    """
    Returns a line of text to add to an environment file, initializing a 
    standard resource with the specified name (string), inflow(int), and 
    outflow(int)
    """
    return "RESOURCE " + resource + ":inflow=" + str(inflow) + \
                ":outflow=" + str(outflow) + "\n"

def gen_cell(args, resource, cells):
    return "CELL " + resource + ":" + ",".join([str(i) for i in cells]) \
        + ":inflow=" + str(args.cellInflow) + ":outflow=" + str(args.cellOutflow) + ":initial=" + str(args.inflow) + "\n"

def gen_reaction(args, resource, depletable=0):
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
            resource + ":value=" + str(args.taskValDict[task]) + ":type=" \
            + args.rxnType + ":frac=" + str(args.frac) + ":max=" + str(args.resMax) + \
            ":depletable=" + str(int(depletable)) + " requisite:max_count=" \
            + str(args.maxCount) + "\n"

def calcEvenAnchors(args):
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

def calcRandomAnchors(args, inworld=True):
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

def get_surrounding_coord_even(coord, d, dsout):
    return  [coord + int(floor(d/2.0))+d*i for i in range(dsout)]\
             + [coord - (int(ceil(d/2.0))+d*i) for i in range(dsout)]

def get_surrounding_coord_odd(coord, d, dsout):
    return [coord + d*i for i in range(dsout)] +\
           [coord - d*i for i in range(dsout)]

def add_anchors(centerPoint, d, dsout, anchors, even=True):

    if even:
        xs = get_surrounding_coord_even(centerPoint[0], d, dsout)
        ys = get_surrounding_coord_even(centerPoint[0], d, dsout)

    else:
        xs = get_surrounding_coord_odd(centerPoint[0], d, dsout)
        ys = get_surrounding_coord_odd(centerPoint[0], d, dsout)

    pairwise_point_combination(xs, ys, anchors)
    

def calcTightAnchors(args, d, patches):
    """
    Recursively generates the number of anchor points specified in the 
    patches argument, such that all patches are d cells away
    from their nearest neighbors.
    """
    centerPoint = (int(args.worldSize/2), int(args.worldSize/2))
    anchors = []
    if patches == 0:
        pass

    elif patches == 1:
        anchors.append(centerPoint)
           
    elif patches%2 == 0:
        dsout = (patches-2)/2 + 1
        add_anchors(centerPoint, d, dsout, anchors, True)
        if d != 0:
            anchors = list(set(anchors))
        anchors.sort()
        if dsout != 1:
            return (anchors + calcTightAnchors(d, patches-2))[:patches*patches] 
            #to cut off the extras in the case where d=0

    else:
        #Note - an odd number of args.patchesPerSide requires that there be 
        #a patch at the centerpoint
        dsout = (patches-1)/2
        add_anchors(centerPoint, d, dsout, anchors, False)
        if dsout != 1:
            return anchors + calcTightAnchors(d, patches-2)

    return anchors

def pairwise_point_combinations(xs, ys, anchors):
    """
    Does an in-place addition of the four points that can be composed by
    combining coordinates from the two lists to the given list of anchors
    """
    for i in xs:
        anchors.append((i, max(ys)))
        anchors.append((i, min(ys)))
    for i in ys:
        anchors.append((max(xs), i))
        anchors.append((min(xs), i))

def random_patch(args, size):
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

def get_moore_neighbors(args, cell):
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

def genRandResources(args, resources):
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

def place_even_squares(args, side):
    cells = []
    for y in range(0, args.worldSize, side*2):
        for x in range(0, args.worldSize, side*2):
            for z in range(side):
                for z2 in range(side):
                    cells.append((y+z)*args.worldSize + x + z2)
    cells.sort()
    return cells
