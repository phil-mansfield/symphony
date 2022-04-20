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

def main():
    sim_dir = sys.argv[1] # directory containing the simulation suites
    suite_name = sys.argv[2] # Name of the simulation suite of the halo
    host_id = int(sys.argv[3]) # ID of the host halo
    out_base_dir = sys.argv[4] # base dir where the frames are stored

    base_dir = path.join(sim_dir,  suite_name, "Halo%03d" % host_id)
    out_dir = path.join(out_base_dir, suite_name, "Halo%03d" % host_id)
    os.makedirs(out_dir, exist_ok=True)

    fig, ax = plt.subplots(1, 2, figsize=(16, 8))

    scales = symlib.scale_factors()
    param = symlib.parameter_table[suite_name]
    halos, _ = symlib.read_subhalos(param, base_dir)

    # present day
    snap0 = len(halos[0]) - 1
    # z=1
    snap1 = np.searchsorted(scales, 0.5)
    snaps = [snap1, snap0]

    for i in range(len(halos)):
        for j, snap in enumerate(snaps):
            halo = halos[i,snap]
            host = halos[0,snap]
            scale = scales[snap]

            if not halo["ok"] or not host["ok"]: continue
            dx = (halo["x"][0] - host["x"][0])*scale/param["h100"]*1e3
            dy = (halo["x"][1] - host["x"][1])*scale/param["h100"]*1e3
            r = halo["rvir"]*scale/param["h100"]*1e3

            if i == 0:
                lw = 3
                color = "tab:red"
                rvir_max = r
            else:
                lw = 1
                color = "tab:blue"

            symlib.plot_circle(ax[j], dx, dy, r, color=color, lw=lw)

    ax[0].set_xlim(-1.5*rvir_max, 1.5*rvir_max)
    ax[0].set_ylim(-1.5*rvir_max, 1.5*rvir_max)
    ax[1].set_xlim(-1.5*rvir_max, 1.5*rvir_max)
    ax[1].set_ylim(-1.5*rvir_max, 1.5*rvir_max)

    ax[0].set_title(r"$z=1$")
    ax[1].set_title(r"$z=0$")

    ax[0].set_xlabel(r"$X\ ({\rm kpc})$")
    ax[0].set_ylabel(r"$Y\ ({\rm kpc})$")
    ax[1].set_xlabel(r"$X\ ({\rm kpc})$")
    ax[1].axes.yaxis.set_visible(False)

    plt.savefig(path.join(out_dir, "subhalo_positions.png"))
    
if __name__ == "__main__": main()
