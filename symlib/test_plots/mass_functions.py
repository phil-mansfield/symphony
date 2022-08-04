import numpy as np
import symlib
import os
import os.path as path
import matplotlib.pyplot as plt
try:
    import palette
    from palette import pc
    palette.configure(False)
except:
    pc = lambda x: "tab:" + x

def main():
    # Set up data locations
    base = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns/"
    out_base = "example_plots"

    # plotting globals
    n_suite = len(symlib.SUITE_NAMES)
    fig, ax = plt.subplots()
    bins = 10**np.linspace(5, 14, 200)

    colors = [pc("p"), pc("b"), pc("g"), pc("o"), pc("r")]
    for i_suite in range(n_suite):
        # Figure out directory names and initialize variables
        suite_name = symlib.SUITE_NAMES[i_suite]
        sub_mvir, sub_mpeak, sub_minfall, n_host = [], [], [], 0

        # loop over all hosts
        for i_host in range(symlib.n_halos(suite_name)):
            halo_dir = symlib.get_halo_dir(base, suite_name, i_host)

            # Read subhalo information
            param = symlib.parameter_table[suite_name]
            halos, hists = symlib.read_subhalos(param, halo_dir)

            # Remove disrupted subhaloes and the host.
            ok = ((halos["mvir"][:,-1] >= 0) &
                  (halos["mvir"][:,-1] < halos["mvir"][0,-1]))
            sub_mvir.append(halos["mvir"][ok,-1])
            sub_mpeak.append(hists["mpeak"][ok])

            # The infall mass function
            ok_infall = halos["mvir"][:,-1] < halos["mvir"][0,-1]
            minfall = halos["mvir"][np.arange(len(halos)),hists["merger_snap"]]
            sub_minfall.append(minfall[ok_infall])

            n_host += 1

        # Combine into one big array.
        sub_mvir, sub_mpeak = np.hstack(sub_mvir), np.hstack(sub_mpeak)
        sub_minfall = np.hstack(sub_minfall)

        # Histograming and stacking
        n_mvir, edges = np.histogram(sub_mvir, bins=bins)
        n_mpeak, _ = np.histogram(sub_mpeak, bins=bins)
        n_minfall, _ = np.histogram(sub_minfall, bins=bins)

        n_mvir, n_mpeak, = n_mvir/n_host, n_mpeak/n_host
        n_minfall = n_minfall/n_host

        n_mvir = np.cumsum(n_mvir[::-1])[::-1]
        n_mpeak = np.cumsum(n_mpeak[::-1])[::-1]
        n_minfall = np.cumsum(n_minfall[::-1])[::-1]

        # Everything beyond this point is plotting junk
        m_min = 300 * param["mp"]
        ok = (n_mvir > 0) & (edges[:-1] > m_min)
        ax.plot(edges[:-1][ok], n_mvir[ok], "--", c=colors[i_suite])
        ok = (n_mpeak > 0) & (edges[:-1] > m_min)
        ax.plot(edges[:-1][ok], n_mpeak[ok], "-", c=colors[i_suite])
        ok = (n_minfall > 0) & (edges[:-1] > m_min)
        ax.plot(edges[:-1][ok], n_minfall[ok], ":", c=colors[i_suite], lw=2)
        
    ax.set_ylabel(r"$N(>M_X)$")
    ax.set_xlabel(r"$M_X\,(M_\odot)$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    fig.savefig(path.join(out_base, "mass_func.png"))

if __name__ == "__main__": main()
