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
    # Get file locations.
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    suite = "SymphonyMilkyWay"

    fig, ax = plt.subplots()
    
    # Set up histogram bins
    param = symlib.simulation_parameters(suite)
    m_min = param["mp"]/param["h100"]*300
    m_max = 1e12 
    n_bins = 200
    bins = 10**np.linspace(np.log10(m_min), np.log10(m_max), n_bins+1)

    # Set up cumulative histograms for mass functions
    N_infall = np.zeros(n_bins)
    N_splashback = np.zeros(n_bins)
    N_vir = np.zeros(n_bins)

    n_hosts = symlib.n_hosts(suite)
    for halo_num in range(n_hosts):
        sim_dir = symlib.get_host_directory(base_dir, suite, halo_num)

        # Read in simulation data and convert units.
        scale = symlib.scale_factors(sim_dir)
        h, hist = symlib.read_subhalos(param, sim_dir)
        h = symlib.set_units_halos(h, scale, param)
        hist = symlib.set_units_histories(hist, scale, param)

        # All suriving subhalos which are currently within Rvir of the host.
        r = np.sqrt(np.sum(h["x"][:,-1]**2, axis=1)) # z=0 distance to the host
        host_rvir = h["rvir"][0,-1]
        ok = h["ok"][:,-1] & (r < host_rvir)
        n_vir, _ = np.histogram(hist["mpeak"][ok][1:], bins=bins)
        
        # Add to the cumulative histograms
        N_vir += np.cumsum(n_vir[::-1])[::-1]/n_hosts

    # plot
    left_bins = bins[:-1]
    plt.plot(left_bins, N_vir, c="tab:blue")

    # Some additional plotting code to make things look nice.
    ax.set_xlim(1e8, 1e12)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{\rm sub,peak}$")
    ax.set_ylabel(r"$N(>M_{\rm sub,peak})$")
    fig.savefig("plots/mass_func.png")

if __name__ == "__main__": main()
