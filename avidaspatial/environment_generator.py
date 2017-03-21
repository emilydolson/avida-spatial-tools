#!/usr/bin/python

# WORK IN PROGRESS - do not use!

# This program allows you to automatically generate environment files with
# complex spatial resource layouts

import argparse
import random
from math import sqrt, log, floor, ceil


def get_args(cli=[]):
    parser = argparse.ArgumentParser(description="""Generate environment files
                                     for environmental heterogeneity
                                     experiment.""",
                                     add_help=False)

    modes = parser.add_argument_group("", "")

    # Ways of running
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--environmentPicker", action="store_true",
                            help="""Use the environmentPicker GUI to select
                            cells to fill with resources""")

    # Global settings
    parser.add_argument("--worldSize", default=60, type=int,
                        help="Specifies size of 1 side (assumes square world)")
    parser.add_argument("--randomSeed", default=-1, type=int,
                        help="""The random seed to use. By default, a random
                        one will be selected.""")
    parser.add_argument("--outfile", default="environment.cfg",
                        help="The name of the environment file to be created")
    parser.add_argument("--inflow", default=100, type=int,
                        help="Sets amount of inflow.")
    parser.add_argument("--outflow", default=.01, type=float,
                        help="Sets amount of outflow.")

    # Stuff you can add to the environment
    parser.add_argument("--resources",
                        default=["resNOT", "resAND", "resOR", "resNOR",
                                 "resNAND", "resORN", "resANDN", "resXOR",
                                 "resEQU"],
                        nargs="*",
                        type=str,
                        help="""The set of resources to be used.
                        Defaults to logic-9.""")
    parser.add_argument("--randomPatch", nargs="*",
                        help="""Place a randomly generated patch in the
                        environment. By default, this will include all
                        resources. If additional command-line resources are
                        specified, then only this set will be included in the
                        patch. These resources won't be included with
                        well-mixed resources.""")
    parser.add_argument("--evenSquares", type=int,
                        help="""Place evenly spaced squares of the specified
                        size in environment""")
    parser.add_argument("--gradientResources", nargs="*",
                        help="""Designates which resources are to be
                        implemented as gradient resources (circles in space).
                        If this flag is not used, no resources will be
                        gradients. If this flag is used and no additional
                        arguments are passed, all resources will be gradients.
                        Additional arguments can be passed to this flag to
                        specify a subset of the resources to represent as
                        cells.""")
    parser.add_argument("--cellResources", nargs="*",
                        help="""Designates which resources are to be implemented
                        as cell resources. If this flag is not used, no
                        resources will be cells. If this flag is used and no
                        additional arguments are passed, all resources will be
                        cells. Additional arguments can be passed to this flag
                        to specify a subset of the resources to represent as
                        cells.""")

    # ~~~~~~~~~~Settings for random patches~~~~~~~~~~~~~~~~~~~
    parser.add_argument("--patchSize", default=600, type=int,
                        help="Size of random patch to be generated.")
    parser.add_argument("--cellInflow", default=100, type=int,
                        help="Inflow of cell resources")
    parser.add_argument("--cellOutflow", default=.01, type=int,
                        help="Outflow of cell resources")
    parser.add_argument("--initial", default=10, type=int,
                        help="Initial resource levels of cell resources")

    # ~~~~~~~~~~Settings for gradient resources~~~~~~~~~~~~~~~

    # Choosing anchors:
    anchors_group = parser.add_mutually_exclusive_group()
    anchors_group.add_argument("--randAnchors", action="store_true",
                               help="Generates random anchor points")
    anchors_group.add_argument("--distance", default=None, type=int,
                               help="""Distance between anchor points of
                               central patches""")
    anchors_group.add_argument("--evenAnchors", action="store_true",
                               help="""Places gradient resoruces as evenly as
                               possible through environment.""")
    parser.add_argument("--boundedAnchors", action="store_true",
                        help="Forces patches to be within world.")

    # Choosing radius
    radius_group = parser.add_mutually_exclusive_group()
    parser.add_argument("--patchRadius", default=10, type=int,
                        help="Radius of gradient resources")
    parser.add_argument("--randRadius", nargs="*", type=int,
                        help="""Sets radius of each patch to a random integer
                        between two specified numebers.""")

    # Number of patches
    parser.add_argument("--randNPatches", action="store_true",
                        help="""Use a random number of patches within 10 of
                        the number specified.""")
    parser.add_argument("--patchesPerSide", default=4, type=int,
                        help="""number of patches in one row or column of
                        world (assume square layout)""")

    parser.add_argument("--nPatches", default=4, type=int,
                        help="""Number of patches.""")

    # Other options
    parser.add_argument("--notcommon", action="store_true",
                        help="""Included for backwards compatability - makes
                        plateau gradient resources not be common. This should
                        generally not be specified, although there isn't much
                        evidence that it matters.""")

    # ~~~~~~~~~Reactions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    parser.add_argument("--reaction", nargs="*",
                        help="""A list of resources to make corresponding
                        reactions for. Default if flag is set is to add all."""
                        )
    parser.add_argument("--infinite", action="store_true",
                        help="""Don't deplete resources when reaction
                        are performed.""")
    parser.add_argument("--frac", default=.0025, type=float,
                        help="""The fraction of available resource to use
                        for reactions.""")
    parser.add_argument("--resMax", default=25, type=float,
                        help="""Maximum units of resource an organism can
                        use per reaction.""")
    parser.add_argument("--rxnType", default="pow",
                        help="How is task reward granted?")
    parser.add_argument("--maxCount", default=1, type=int,
                        help="maximum number of times a task can be completed")

    args = parser.parse_args(cli)
    taskValDict = {"not": 1.0,
                   "nand": 1.0,
                   "and": 2.0,
                   "orn": 2.0,
                   "or": 3.0,
                   "andn": 3.0,
                   "nor": 4.0,
                   "xor": 4.0,
                   "equ": 5.0}  # rewards for each task

    args.taskValDict = taskValDict
    return args


def main():  # pragma: no cover

    args = get_args(sys.argv)

    if args.randomSeed != -1:
        random.seed(args.randomSeed)

    args.nPatches = args.patchesPerSide**2
    if args.randNPatches:
        args.nPatches += random.randrange(-60, 60)

    outfile = open(args.outfile, "w")

    # print args

    if args.environmentPicker:
        from cell_picker import cell_picker
        cells = cell_picker(args.worldSize, args.worldSize).run()
        # print args.worldSize
        # cells = ep.cell_picker
        # (60, 60, [(args.worldSize,args.worldSize)]).run()
        cells = list(set(cells))
        cells.sort()
        for resource in args.resources:
            outfile.write(gen_cell(resource, cells))

        outfile.write("\n")

        for resource in args.resources:
            # print resource
            outfile.write(gen_reaction(resource))

        outfile.close()
        exit(0)

    spatial_resources = []
    if args.gradientResources is not None:
        spatial_resources += args.gradientResources if args.gradientResources \
            != [] else args.resources
    if args.cellResources is not None:
        spatial_resources += args.cellResources if args.cellResources \
            != [] else args.resources
    if args.randomPatch is not None:
        spatial_resources += args.randomPatch if args.randomPatch \
            != [] else args.resources

    # Add non-spatial resources:
    well_mixed_resources = args.resources[:]
    for res in spatial_resources:
        if res in well_mixed_resources:
            well_mixed_resources.remove(res)

    for res in well_mixed_resources:
        outfile.write(gen_res(res, args.inflow, args.outflow))

    # Add random patch
    if args.randomPatch is not None:
        cells = random_patch(args.patchSize)
        cells.sort()
        for resource in args.resources:
            outfile.write(gen_cell(resource, cells))

        outfile.write("\n")

        for resource in args.resources \
                if len(args.randomPatch) == 0 else args.randomPatch:
            # print resource
            outfile.write(gen_reaction(resource))

    outfile.write("\n")

    # Add gradient resources if appropriate
    if args.gradientResources is not None:
        anchors = []
        if args.randAnchors:  # random anchor placement
            anchors = calcRandomAnchors(args.boundedAnchors)
        elif args.distance is None:  # Evenly spaced anchor points
            anchors = calcEvenAnchors()
        else:  # Anchors placed a specific distance from each other
            anchors = calcTightAnchors(args.distance, args.patchesPerSide)

        # Generate random resource order
        randResources = genRandResources(args.gradientResources
                                         if args.gradientResources != []
                                         else args.resources)
        random.shuffle(randResources)

        # Write confiugration to specifed output file
        for i in range(len(anchors)):
            if args.randRadius is not None:
                radius = random.randrange(args.randRadius[0],
                                          args.randRadius[1])
                outfile.write(gen_gradient(randResources[i], args.inflow,
                                           radius, anchors[i],
                                           common=(not args.notcommon)))
            else:
                outfile.write(gen_gradient(randResources[i], args.inflow,
                                           args.patchRadius, anchors[i],
                                           common=(not args.notcommon)))
                rads.append(args.patchRadius)

    if args.evenSquares is not None:
        cells = place_even_squares(args.evenSquares)
        for res in args.cellResources:
            outfile.write(gen_cell(res, cells))

    outfile.write("\n")

    for res in args.resources:
        outfile.write(gen_reaction(res, not args.infinite))

    outfile.close()

if __name__ == "__main__":
    main()
