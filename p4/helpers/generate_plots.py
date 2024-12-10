import numpy as np
import matplotlib.pyplot as plt
from IPython import embed


if __name__=="__main__":
    rs = np.linspace(1.05, 1.95, 10)

    # reading pl strain results
    with open("../p4a/fpprts.dat", "r") as f:
        lines = f.readlines()

    pl_str_sigma_tt_1 = []
    pl_str_sigma_rr_1 = []

    pl_str_sigma_tt_2 = []
    pl_str_sigma_rr_2 = []

    pl_str_sigma_tt_3 = []
    pl_str_sigma_rr_3 = []
    for line in lines:
        line = line.split()
        if line[0] == "STRESS":
            ele = int(line[2])

            sigma_xx = float(line[5])
            sigma_yy = float(line[6])
            sigma_xy = float(line[7])

            if ele in np.arange(1, 11):
                theta = np.pi / 2 - np.pi / 2 / 8 / 2
            elif ele in np.arange(41, 51):
                theta = 3 * np.pi / 2 / 8 + np.pi / 2 / 8 / 2
            elif ele in np.arange(71, 81):
                theta = np.pi / 2 / 8 / 2
            else:
                continue

            A = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
            sigma_B = np.array([[sigma_xx, sigma_xy], [sigma_xy, sigma_yy]])
            sigma_Bprime = A @ sigma_B @ A.T
            sigma_tt = sigma_Bprime[1,1]
            sigma_rr = sigma_Bprime[0,0]

            if ele in np.arange(1, 11):
                pl_str_sigma_tt_1.append(sigma_tt)
                pl_str_sigma_rr_1.append(sigma_rr)
            elif ele in np.arange(41, 51):
                pl_str_sigma_tt_2.append(sigma_tt)
                pl_str_sigma_rr_2.append(sigma_rr)
            elif ele in np.arange(71, 81):
                pl_str_sigma_tt_3.append(sigma_tt)
                pl_str_sigma_rr_3.append(sigma_rr)

    pl_str_sigma_tt_1 = np.array(pl_str_sigma_tt_1)
    pl_str_sigma_rr_1 = np.array(pl_str_sigma_rr_1)
    pl_str_sigma_tt_2 = np.array(pl_str_sigma_tt_2)
    pl_str_sigma_rr_2 = np.array(pl_str_sigma_rr_2)
    pl_str_sigma_tt_3 = np.array(pl_str_sigma_tt_3)
    pl_str_sigma_rr_3 = np.array(pl_str_sigma_rr_3)

    # reading axisymmetric results
    with open("../p4b/fpprts.dat", "r") as f:
        lines = f.readlines()

    axi_sigma_tt = []
    axi_sigma_rr = []
    for line in lines:
        line = line.split()
        if line[0] == "STRESS":
            ele = int(line[2])

            sigma_tt = float(line[4])
            sigma_rr = float(line[5])

            axi_sigma_tt.append(sigma_tt)
            axi_sigma_rr.append(sigma_rr)

    axi_sigma_tt = np.array(axi_sigma_tt)
    axi_sigma_rr = np.array(axi_sigma_rr)

    # compute analytic results
    f_sigma_tt = lambda r: 1 / 3 * (1 + 4 / r ** 2)
    f_sigma_rr = lambda r: 1 / 3 * (1 - 4 / r ** 2)

    rs_ = np.linspace(1.05, 1.95, 256)
    analytic_sigma_tt = []
    analytic_sigma_rr = []
    for r in rs_:
        analytic_sigma_tt.append(f_sigma_tt(r))
        analytic_sigma_rr.append(f_sigma_rr(r))

    analytic_sigma_tt = np.array(analytic_sigma_tt)
    analytic_sigma_rr = np.array(analytic_sigma_rr)

    # sigma_tt plots
    f, ax = plt.subplots(figsize=(16,9))
    ax.set_title(r"$\sigma_{\theta\theta}$ vs $r$", fontsize=24)
    ax.set_xlabel(r'$r$', fontsize=18)
    ax.set_ylabel(r'$\sigma_{\theta \theta}$', fontsize=18)
    ax.plot(rs_, analytic_sigma_tt, 'k--', label="Analytic", linewidth=2, alpha=0.9)
    ax.plot(rs, axi_sigma_tt, 'g-.', label="Axisymmetric", linewidth=2, alpha=0.7)
    ax.plot(rs, pl_str_sigma_tt_1, 'r-.', label=r"Plane-$\varepsilon$ ($\theta = \pi / 32$)", linewidth=2, alpha=0.7)
    ax.plot(rs, pl_str_sigma_tt_2, 'b:', label=r"Plane-$\varepsilon$ ($\theta = 7 \pi / 32$)", linewidth=2, alpha=0.7)
    # ax.plot(rs, pl_str_sigma_tt_3, 'g--', label=r"Plane-$\varepsilon$ ($\theta = 15 \pi / 32$)")
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)
    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.grid()
    f.savefig("sigma_tt.pdf")
    plt.close()

    # sigma_rr plots
    f, ax = plt.subplots(figsize=(16,9))
    ax.set_title(r"$\sigma_{rr}$ vs $r$", fontsize=24)
    ax.set_xlabel(r'$r$', fontsize=18)
    ax.set_ylabel(r'$\sigma_{rr}$', fontsize=18)
    ax.plot(rs_, analytic_sigma_rr, 'k--', label="Analytic", linewidth=2, alpha=0.9)
    ax.plot(rs, axi_sigma_rr, 'g-.', label="Axisymmetric", linewidth=2, alpha=0.7)
    ax.plot(rs, pl_str_sigma_rr_1, 'r-.', label=r"Plane-$\varepsilon$ ($\theta = \pi / 32$)", linewidth=2, alpha=0.7)
    ax.plot(rs, pl_str_sigma_rr_2, 'b:', label=r"Plane-$\varepsilon$ ($\theta = 7 \pi / 32$)", linewidth=2, alpha=0.7)
    # ax.plot(rs, pl_str_sigma_rr_3, 'g--', label=r"Plane-$\varepsilon$ ($\theta = 15 \pi / 32$)")
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)
    ax.legend(fontsize=14)
    plt.tight_layout()
    plt.grid()
    f.savefig("sigma_rr.pdf")
    plt.close()

    embed()
