import pysal
import numpy as np


def make_toroidal_weights(nrows, ncols, rook=True):
    x = ncols
    y = nrows
    w = pysal.lat2W(y, x, rook=rook).neighbors

    for i in range(x*y):
        if i == 0:
            # Upper left
            w[i].append(x-1)
            w[i].append(x*(y-1))

            if not rook:
                w[i].append((x*y)-1)
                w[i].append((2*x)-1)
                w[i].append(x*(y-1)+1)

        elif i == x-1:
            # Upper right
            w[i].append(0)
            w[i].append((x*y)-1)

            if not rook:
                w[i].append(x*(y-1))
                w[i].append(x)
                w[i].append((x*y)-2)

        elif i == x*(y-1):
            # Lower left
            w[i].append((x*y)-1)
            w[i].append(0)

            if not rook:
                w[i].append(x-1)
                w[i].append(1)
                w[i].append(i-1)

        elif i == (x*y)-1:
            # Lower right
            w[i].append(((x*y)-1))
            w[i].append(x-1)

            if not rook:
                w[i].append(0)
                w[i].append(x)
                w[i].append(1)

        elif i % x == 0:
            # we are on the left edge
            w[i].append(i+x-1)

            if not rook:
                w[i].append(i+(2*x)-1)
                w[i].append(i-1)

        elif i % x == x-1:
            # we are on the right edge
            w[i].append(i-x+1)

            if not rook:
                w[i].append(i+1)
                w[i].append(i-(2*x)+1)

        elif i // x == 0:
            # we are on the top edge
            w[i].append((x*(y-1)) + i)

            if not rook:
                w[i].append(((x*(y-1)) + i % x)-1)
                w[i].append(((x*(y-1)) + i % x)+1)

        elif i // x == y-1:
            w[i].append(i % x)

            if not rook:
                w[i].append((i % x)-1)
                w[i].append((i % x)+1)

    weights = pysal.weights.weights.W(w)
    return weights


def write_toroidal_contig_matrix_to_file(x, y,
                                         filename="toroidal_dist_matrix.csv"):
    weights = make_toroidal_weights(x, y)
    outfile = open(filename, "w")
    for i in weights.full()[0]:
        outfile.write(",".join([str(j) for j in i])+"\n")

    outfile.close()
