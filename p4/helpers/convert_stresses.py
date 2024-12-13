import numpy as np
from IPython import embed


if __name__ == "__main__":
    with open("../p4a/fpprts.dat", "r") as f:
        lines = f.readlines()

    sigma_dict = {}

    for line in lines:
        line = line.split()
        if line[0] == "STRESS":
            ele = int(line[2])
            sigma_dict[ele] = {}

            sigma_xx = float(line[5])
            sigma_yy = float(line[6])
            sigma_xy = float(line[7])
            sigma_dict[ele]["sigma_xx"] = sigma_xx
            sigma_dict[ele]["sigma_yy"] = sigma_yy

            if ele in np.arange(1, 11):
                theta = np.pi / 2 - np.pi / 2 / 8 / 2
            elif ele in np.arange(21, 31):
                theta = 5 * np.pi / 2 / 8 + np.pi / 2 / 8 / 2
            elif ele in np.arange(31, 41):
                theta = 4 * np.pi / 2 / 8 + np.pi / 2 / 8 / 2
            elif ele in np.arange(41, 51):
                theta = 3 * np.pi / 2 / 8 + np.pi / 2 / 8 / 2
            elif ele in np.arange(71, 81):
                theta = np.pi / 2 / 8 / 2

            A = np.array(
                [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]
            )
            sigma_B = np.array([[sigma_xx, sigma_xy], [sigma_xy, sigma_yy]])
            sigma_Bprime = A @ sigma_B @ A.T
            sigma_rr = sigma_Bprime[0, 0]
            sigma_tt = sigma_Bprime[1, 1]
            # theta *= -1

            sigma_dict[ele]["sigma_rr"] = sigma_rr
            sigma_dict[ele]["sigma_tt"] = sigma_tt

            # sigma_dict[ele]['sigma_rr'] = sigma_xx * np.cos(theta) + sigma_yy * np.sin(theta)
            # sigma_dict[ele]['sigma_tt'] = - sigma_xx * np.sin(theta) + sigma_yy * np.cos(theta)

    for ele, sigmas in sigma_dict.items():
        print(
            f"Element: {ele} -> sigma_tt = {sigmas['sigma_tt']:.4f} || sigma_rr = {sigmas['sigma_rr']:.4f}"
        )
        if ele in [10, 30, 40, 50, 80]:
            print("\n")
    embed()
