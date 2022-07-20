import numpy as np
import matplotlib.pyplot as plt
import symlib

""" tree_mergers.py counts the number of mergers per snapshot as a function of
 a(z) for the main halo.
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

    param = symlib.simulation_parameters(suite)

    # We'll need the normal subhalos so we can find the host.
    scale = symlib.scale_factors(sim_dir)
    h, hist = symlib.read_subhalos(sim_dir)

    # Read in tree data
    b = symlib.read_branches(sim_dir)
    # Tree variables are always returned as a list, so if you only specify one,
    # unpack it as a length-1 tuple.
    dfid, next_co_prog, snap = symlib.read_tree(
        sim_dir, ["dfid", "next_co_prog", "snap"]
    )

    host_branch = b[hist["branch_idx"][0]]
    host_start, host_end = host_branch["start"], host_branch["end"]

    # Flag halo branches which are probably not artifacts.
    ok = b["is_real"] & (~b["is_disappear"])
    
    # Step 1: counting the number of mergers. Requires a lookup table, which
    # we construct from the branch information and the depth-first IDs ("dfid")
    table = symlib.merger_lookup_table(b, dfid)
    n_mergers = np.zeros(host_end - host_start, dtype=int)
    n_artifacts = np.zeros(host_end - host_start, dtype=int)
    for i in range(host_start, host_end):
        branch_idx = symlib.find_all_merger_branches(
            b, table, next_co_prog, i)
        n_mergers[i - host_start] = np.sum(ok[branch_idx])
        n_artifacts[i - host_start] = np.sum(~ok[branch_idx])

    # Step 2: getting the scale factor of each snapshot.
    host_snap = snap[host_start: host_end]
    host_scale = scale[host_snap]

    fig, ax = plt.subplots()
    ax.plot(host_scale, n_mergers, "tab:blue",
            label=r"$N_{\rm merger}$")
    ax.plot(host_scale, n_artifacts, "tab:red",
            label=r"$N_{\rm artifact}$")

    # Additional plotting code to make things prettier
    ax.set_xlabel(r"$a(z)$")
    ax.set_ylabel(r"$N$")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    ax.legend(loc="upper left")

    fig.savefig("plots/tree_mergers.png")

if __name__ == "__main__": main()
