import pickle
import numpy as np
from IPython import embed
import matplotlib.pyplot as plt


if __name__=="__main__":
    with open("metadata.pkl", "rb") as f:
        fpin = pickle.load(f)

    # exract all node ids and coordinates
    coords = {}
    coord_data = fpin['COOR']
    line_ind = 0
    n_lines = len(fpin['COOR'])
    while True:
        line1 = coord_data[line_ind]
        line1 = [float(x) for x in line1.split()]

        coord, gen_incr, x, y = line1

        # determine generation increment
        if gen_incr == 1:
            line_ind += 1
            line2 = coord_data[line_ind]
            line2 = [float(x) for x in line2.split()]

            coord2, gen_incr2, x2, y2 = line2

            xs = np.linspace(x, x2, int(coord2 - coord + 1))
            ys = np.linspace(y, y2, int(coord2 - coord + 1))

            coords_ = np.arange(coord, coord2 + 1)

            for coord_, x_, y_ in zip(coords_, xs, ys):
                coords[coord_.item()] = [x_.item(), y_.item()]
            line_ind += 1
        else:
            coords[coord.item()] = [x.item(), y.item()]
            line_ind += 1

        if line_ind == n_lines:
            break

    # extract all elements
    ele_dict = {}
    elem_data = fpin["ELEM"]
    n_lines = len(elem_data)
    line_ind = 0
    while True:
        line1 = elem_data[line_ind]
        line1 = [float(x) for x in line1.split()]
        ele, mat_set, n1, n2, n3, n4, gen_incr = line1

        if gen_incr == 1:
            line_ind += 1
            line2 = elem_data[line_ind]
            line2 = [float(x) for x in line2.split()]

            ele2, mat_set2, n1_2, n2_2, n3_2, n4_2, gen_incr_2 = line2

            for i in range(int(ele), int(ele2) + 1):
                ele_dict[i] = [n1, n2, n3, n4]

                n1 += gen_incr
                n2 += gen_incr
                n3 += gen_incr
                n4 += gen_incr

            line_ind += 1
        else:
            ele_dict[ele.item()] = [n1, n2, n3, n4]

            line_ind += 1

        if line_ind == n_lines:
            break

    line_list = [[(x1, x2), (x2, x3), (x3, x4), (x4, x1)] for x1, x2, x3, x4 in ele_dict.values()]
    line_list = [x for ele in line_list for x in ele]

    f, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(*np.array(list(coords.values())).T, marker='o', color='r')
    for line in line_list:
        n1 = coords[line[0]]
        n2 = coords[line[1]]
        ax.plot([n1[0], n2[0]], [n1[1], n2[1]], 'k')

    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()

    f.savefig('q4a_mesh.pdf', bbox_inches='tight')

    embed()

