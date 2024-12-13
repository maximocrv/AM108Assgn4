import numpy as np
from IPython import embed


if __name__ == "__main__":
    with open("../p4b/fpprts.dat", "r") as f:
        lines = f.readlines()

    sigma_dict = {}

    for line in lines:
        line = line.split()
        if line[0] == "STRESS":
            ele = int(line[2])
            sigma_dict[ele] = {}

            sigma_tt = float(line[4])
            sigma_rr = float(line[5])
            sigma_dict[ele]["sigma_tt"] = sigma_tt
            sigma_dict[ele]["sigma_rr"] = sigma_rr

    for k, v in sigma_dict.items():
        print(f"& {v['sigma_tt']:.4f} & {v['sigma_rr']:.4f}")

    embed()
