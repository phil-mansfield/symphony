import numpy as np
import matplotlib.pyplot as plt
import symlib

""" tree_mass_func.py constructs an Mpeak mass function of the host and the
field halos around it.
"""

# Import a library that makes plots prettier, if it's downloaded.
try:
    import palette
    palette.configure(True)
except:
    pass

def main():
    # Get file locations.
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"
    halo_num = 0
    sim_dir = symlib.get_host_directory(base_dir, suite, halo_num)

    # Read in simulation data and convert units.
    param = symlib.simulation_parameters(suite)

    # Read in tree data
    b = symlib.read_branches(sim_dir)
    # Tree variables are always returned as a list, so if you only specify one,
    # unpack it as a length-1 tuple.
    mvir, = symlib.read_tree(sim_dir, ["mvir"])
    # Convert units
    mvir, mp = mvir/param["h100"], param["mp"]/param["h100"]

    # Flag halo branches which are probably not artifacts.
    ok = b["is_real"] & (~b["is_disappear"])

    mpeak = np.zeros(len(b))

    for i in range(len(mpeak)):
        if not ok[i]: continue
        start, end = b[i]["start"], b[i]["end"]
        mpeak[i] = np.max(mvir[start: end])
    
    # Find host subhalos
    mpeak_host = mpeak[b["is_main_sub"]]

    # Calculate the mass function of both groups of halos
    bins = np.logspace(np.log10(param["mp"]/param["h100"]), 13, 200)
    n_host, _ = np.histogram(mpeak_host, bins=bins)
    n_all, _ = np.histogram(mpeak, bins=bins)
    N_host = np.cumsum(n_host[::-1])[::-1]
    N_all = np.cumsum(n_all[::-1])[::-1]

    # Plot
    fig, ax = plt.subplots()
    left_bins = bins[:-1]
    plt.plot(left_bins, N_host, c="tab:red", label=r"${\rm Host\ subhalos}$")
    plt.plot(left_bins, N_all, c="tab:blue", label=r"${\rm All}$")

    # Some plotting code to make things look nice.
    ax.set_xscale("log")
    ax.set_yscale("log")
    ylo, yhi = ax.get_ylim()
    ax.set_ylim(ylo, yhi)
    plt.plot(2*[300*param["mp"]/param["h100"]],
             [ylo, yhi], "--", c="k", lw=2.5)
    ax.legend(loc="upper right")
    ax.set_xlabel(r"$M_{\rm peak}\ (M_\odot)$")
    ax.set_ylabel(r"$N(>M_{\rm peak})$")
    fig.savefig("plots/tree_mass_func.png")

if __name__ == "__main__": main()
