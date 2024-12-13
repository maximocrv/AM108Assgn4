import pickle
import numpy as np
from IPython import embed
import matplotlib.pyplot as plt


if __name__ == "__main__":
    thetas = np.linspace(0, np.pi / 2, 9)

    fpin = {"COOR": []}
    node = 1
    for i, theta in enumerate(thetas):
        # node string
        node_str = f"{node}"
        ele1 = " " * (5 - len(node_str)) + node_str

        # generation increment
        ele2 = " " * 4 + "1"

        # coordinates
        x = np.sin(theta)
        y = np.cos(theta)

        x_str = f"{x:.4f}"
        y_str = f"{y:.4f}"
        ele3 = " " * (10 - len(x_str)) + x_str + " " * (10 - len(y_str)) + y_str

        line1 = ele1 + ele2 + ele3
        fpin["COOR"].append(line1)

        # node string
        node += 10
        node_str = f"{node}"
        ele1 = " " * (5 - len(node_str)) + node_str

        # generation increment
        ele2 = " " * 4 + "0"

        # coordinates
        x_str = f"{2 * x:.4f}"
        y_str = f"{2 * y:.4f}"
        ele3 = " " * (10 - len(x_str)) + x_str + " " * (10 - len(y_str)) + y_str

        line2 = ele1 + ele2 + ele3

        fpin["COOR"].append(line2)

        node += 1

    # defining elements
    fpin["ELEM"] = []
    # define nodes in element 1
    n1 = 1
    n2 = 12
    n3 = 13
    n4 = 2
    elem = 1
    for i in range(8):
        # element id
        elem_str = f"{elem}"
        ele1 = " " * (5 - len(elem_str)) + elem_str

        # material set
        ele2 = " " * 4 + "1"

        # nodes
        n1_str = f"{n1}"
        n2_str = f"{n2}"
        n3_str = f"{n3}"
        n4_str = f"{n4}"

        ele3 = (
            " " * (5 - len(n1_str))
            + n1_str
            + " " * (5 - len(n2_str))
            + n2_str
            + " " * (5 - len(n3_str))
            + n3_str
            + " " * (5 - len(n4_str))
            + n4_str
        )

        # generation increment
        ele4 = " " * 4 + "1"

        # construct line
        line1 = ele1 + ele2 + ele3 + ele4
        fpin["ELEM"].append(line1)

        elem += 9

        # element id
        elem_str = f"{elem}"
        ele1 = " " * (5 - len(elem_str)) + elem_str

        # material set
        ele2 = " " * 4 + "1"

        # nodes
        n1 += 9
        n2 += 9
        n3 += 9
        n4 += 9

        n1_str = f"{n1}"
        n2_str = f"{n2}"
        n3_str = f"{n3}"
        n4_str = f"{n4}"

        ele3 = (
            " " * (5 - len(n1_str))
            + n1_str
            + " " * (5 - len(n2_str))
            + n2_str
            + " " * (5 - len(n3_str))
            + n3_str
            + " " * (5 - len(n4_str))
            + n4_str
        )

        # generation increment
        ele4 = " " * 4 + "0"

        # construct line
        line2 = ele1 + ele2 + ele3 + ele4
        fpin["ELEM"].append(line2)

        elem += 1
        n1 += 2
        n2 += 2
        n3 += 2
        n4 += 2

    # defining consistent nodal loads
    fpin["FORC"] = []
    dtheta = np.diff(thetas)[-1]
    node = 1
    for i, theta in enumerate(thetas):
        # node string
        node_str = f"{node}"
        ele1 = " " * (5 - len(node_str)) + node_str

        # generation increment
        ele2 = " " * 4 + "0"

        # proportional table numbers
        ele3 = " " * 4 + "1" + " " * 4 + "1"

        line1 = ele1 + ele2 + ele3

        # consistent nodal loads
        if node == 1:
            f1 = dtheta * np.sin(thetas[i]) / 2
            f2 = dtheta * np.cos(thetas[i]) / 2

            # using straight line segments
            # f = 2 * np.sin(dtheta / 2) / 2
            # f1 = f * np.sin(thetas[i] + dtheta / 2)
            # f2 = f * np.cos(thetas[i] + dtheta / 2)
        elif node == 89:
            f1 = dtheta * np.sin(thetas[i]) / 2
            f2 = dtheta * np.cos(thetas[i]) / 2
            # using straight line segments
            # f = 2 * np.sin(dtheta / 2) / 2
            # f1 = f * np.sin(thetas[i] - dtheta / 2)
            # f2 = f * np.cos(thetas[i] - dtheta / 2)
        else:
            f1 = dtheta * np.sin(thetas[i])
            f2 = dtheta * np.cos(thetas[i])

            # using straight line segments
            # f = 2 * np.sin(dtheta / 2) / 2
            # f1 = f * np.sin(thetas[i] - dtheta / 2) + f * np.sin(thetas[i] + dtheta / 2)
            # f2 = f * np.cos(thetas[i] - dtheta / 2) + f * np.cos(thetas[i] + dtheta / 2)

        f1_str = f"{f1:.4f}"
        f2_str = f"{f2:.4f}"

        line2 = " " * (10 - len(f1_str)) + f1_str + " " * (10 - len(f2_str)) + f2_str

        fpin["FORC"].append(line1)
        fpin["FORC"].append(line2)

        # node string
        node += 11

    with open("metadata.pkl", "wb") as f:
        pickle.dump(fpin, f)

    embed()
