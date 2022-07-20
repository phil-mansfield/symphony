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
    # Read in simulation data.
    sim_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/SymphonyMilkyWay/Halo023"
    h, hist = symlib.read_subhalos(sim_dir)

    fig, ax = plt.subplots()

    host = h[0,-1] # First halo, last snapshot.
    symlib.plot_circle(ax, host["x"][0], host["x"][1],
                       host["rvir"], c="tab:red")

    for i in range(1, len(h)):
        sub = h[i,-1] # i-th halo, last snapshot.
        if sub["ok"]:
            symlib.plot_circle(ax, sub["x"][0], sub["x"][1],
                               sub["rvir"], c="tab:blue", lw=1.5)

        
    # Some plotting code to make things look nice.
    ax.set_xlim(-1.75*host["rvir"], 1.75*host["rvir"])
    ax.set_ylim(-1.75*host["rvir"], 1.75*host["rvir"])
    ax.set_xlabel(r"$X\ ({\rm kpc})$")
    ax.set_ylabel(r"$Y\ ({\rm kpc})$")
    fig.savefig("plots/positions.png")

if __name__ == "__main__": main()
