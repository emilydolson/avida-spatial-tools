import pysal
import numpy as np

x = 5
y = 5

w = pysal.lat2W(x, y).neighbors

for i in range(x*y):

    if i == 0:
        # Upper left
        w[i].append(x-1)
        w[i].append(x*(y-1))

    elif i == x-1:
        # Upper right
        w[i].append(0)
        w[i].append((x*y)-1)

    elif i == x*(y-1):
        # Lower left
        w[i].append((x*y)-1)
        w[i].append(0)

    elif i == (x*y)-1:
        # Lower right
        w[i].append(((x*y)-1))
        w[i].append(x-1)

    elif i % x == 0:
        # we are on the left edge
        w[i].append(i+x-1)

    elif i % x == x-1:
        # we are on the right edge
        w[i].append(i-x+1)

    elif i // x == 0:
        # we are on the top edge
        w[i].append((x*(y-1)) + i % x)

    elif i // x == y-1:
        w[i].append(i % x)


weights = pysal.weights.weights.W(w)
outfile = open("toroidal_dist_matrix.csv", "w")
for i in weights.full()[0]:
    outfile.write(",".join([str(j) for j in i])+"\n")

outfile.close()
