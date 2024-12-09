import numpy as np
from IPython import embed


if __name__=="__main__":
    with open("fpprts.dat", "r") as f:
        lines = f.readlines()

    sigma_dict = {}

    for line in lines:
        line = line.split()
        if line[0] == "STRESS":
            ele = int(line[2])
            sigma_dict[ele] = {}

            sigma_xx = float(line[5])
            sigma_yy = float(line[6])
            sigma_dict[ele]['sigma_xx'] = sigma_xx
            sigma_dict[ele]['sigma_yy'] = sigma_yy

            if ele in np.arange(1, 12):
                theta = np.pi / 2 / 8 / 2
            elif ele in np.arange(31, 41):
                theta = 3 * np.pi / 2 / 8 + np.pi / 2 / 8 / 2
            elif ele in np.arange(71, 81):
                theta = np.pi / 2 - np.pi / 2 / 8 / 2

            sigma_dict[ele]['sigma_rr'] = sigma_yy * np.cos(theta) + sigma_xx * np.sin(theta)
            sigma_dict[ele]['sigma_tt'] = sigma_yy * np.sin(theta) - sigma_xx * np.cos(theta)

    for ele, sigmas in sigma_dict.items():
        print(f"Element: {ele}")
        print(f"sigma_tt = {sigmas['sigma_tt']}")
        print(f"sigma_rr = {sigmas['sigma_rr']}")
        print("\n")
    embed()
