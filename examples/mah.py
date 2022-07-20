import numpy as np
import matplotlib.pyplot as plt
import symlib

""" positions.py plots a host halo and its subhalos as circles.
"""

# Import a library that makes plots prettier, if it's downloaded.
try:
    import palette
    palette.configure(True)
except:
    pass

def main():
    # Get file locations. You could also point to the directory directly.
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"
    halo_num = 0
    sim_dir = symlib.get_host_directory(base_dir, suite, halo_num)

    # Read in simulation data.
    h, hist = symlib.read_subhalos(sim_dir)
    scale = symlib.scale_factors(sim_dir)

    # Snapshots, for the purpose of making cuts
    snaps = np.arange(len(h[0]))

    fig, ax = plt.subplots()
    colors = ["k", "tab:red", "tab:orange", "tab:green",
              "tab:blue", "tab:purple"]
    for i in range(6):
        ok = h[i,:]["ok"] # Snapshots where the halo exists
        if i == 0:
            # Plot the host halo
            plt.plot(scale[ok], h[i,ok]["mvir"], c=colors[i])
        else:
            # Plot the full history of the subhalo as a dahsed line
            plt.plot(scale[ok], h[i,ok]["mvir"], "--", c=colors[i], lw=2)
            # Plot its history inside the host halo as a solid line
            is_sub = (snaps >= hist["merger_snap"][i]) & ok
            plt.plot(scale[is_sub], h[i,is_sub]["mvir"], c=colors[i])
        
    # Some plotting code to make things look nice.
    param = symlib.simulation_parameters(suite)
    ax.set_xlim(0, 1)
    ax.set_ylim(30*param["mp"]/param["h100"], 2*h[0,-1]["mvir"])
    ax.set_yscale("log")
    ax.set_xlabel(r"$a(z)$")
    ax.set_ylabel(r"$M_{\rm vir}\ (M_\odot)$")
    fig.savefig("plots/mah.png")

if __name__ == "__main__": main()
