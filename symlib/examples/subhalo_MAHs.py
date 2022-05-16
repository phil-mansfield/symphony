import numpy as np
import matplotlib.pyplot as plt
import symlib
import sys
import os
import os.path as path

try:
    import palette
    palette.configure(False)
except:
    pass

def get_args():
    sim_dir = sys.argv[1] # directory containing the simulation suites
    suite_name = sys.argv[2] # Name of the simulation suite of the halo
    host_name = sys.argv[3] # ID of the host halo
    out_base_dir = sys.argv[4] # base dir where the frames are stored

    base_dir = path.join(sim_dir,  suite_name, host_name)
    out_dir = path.join(out_base_dir, suite_name, host_name)
    os.makedirs(out_dir, exist_ok=True)

    return base_dir, suite_name, out_dir

def main():
    base_dir, suite_name, out_dir = get_args()

    fig, ax = plt.subplots()

    scales = symlib.scale_factors(base_dir)
    snaps = np.arange(len(scales))
    param = symlib.parameter_table[suite_name]
    halos, hists = symlib.read_subhalos(param, base_dir)

    colors = ["k", "tab:red", "tab:orange", "tab:green",
              "tab:blue", "tab:purple"]

    for i in range(len(colors)):
        halo = halos[i]
        ok = halo["ok"]

        if i == 0:
            plt.plot(scales[ok], halo["mvir"][ok]/param["h100"],
                     lw=3, c=colors[i])
            continue

        ax.plot(scales[ok], halo["mvir"][ok]/param["h100"],
                "--", lw=2, color=colors[i % len(colors)])
        ok = halo["ok"] & (snaps < hists[i]["merger_snap"])
        ax.plot(scales[ok], halo["mvir"][ok]/param["h100"],
                lw=3, color=colors[i % len(colors)])

        first_infall_snap = hists[i]["first_infall_snap"]
        ax.plot([scales[first_infall_snap]],
                [halo[first_infall_snap]["mvir"]/param["h100"]],
                "o", color=colors[i % len(colors)])

    plt.plot([scales[0], scales[-1]], 2*[300*param["mp"]/param["h100"]],
             "--", lw=2, c="k")

    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_xlim(scales[0], scales[-1])

    ax.set_xlabel(r"$a(z)$")
    ax.set_ylabel(r"$M_{\rm vir}$")

    plt.savefig(path.join(out_dir, "subhalo_MAHs.png"))

if __name__ == "__main__": main()
