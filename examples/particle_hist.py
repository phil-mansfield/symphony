import matplotlib.pyplot as plt
import symlib
import matplotlib as mpl

try:
    import palette
    palette.configure(True)
except:
    pass

def main():
    fig, ax = plt.subplots(2, 2, figsize=(16, 16), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    base_dir = "/sdf/home/p/phil1/ZoomIns"
    sim_dir = symlib.get_host_directory(base_dir, "SymphonyLCluster", "Halo_042")

    rs, hist = symlib.read_rockstar(sim_dir)
    host, sub = rs[0], rs[1]
    lim = 1.5*host["rvir"][-1]

    for irow in range(2):
        for icol in range(2):
            if icol == 0: ax[irow,icol].set_ylabel(r"$Y\ ({\rm kpc})$")
            if irow == 1:
                ax[irow,icol].set_xlabel(r"$X\ ({\rm kpc})$")
                sub_ls, host_ls = "-", ":"
            else:
                sub_ls, host_ls = ":", "-"

            symlib.plot_circle(
                ax[irow,icol], host["x"][-1,0], host["x"][-1,1],
                host["rvir"][-1], c="white", ls=host_ls
            )
            symlib.plot_circle(
                ax[irow,icol], sub["x"][-1,0], sub["x"][-1,1],
                sub["rvir"][-1], c="white", ls=sub_ls
            )
            ax[irow,icol].set_xlim(-lim, lim)
            ax[irow,icol].set_ylim(-lim, lim)
            
    part = symlib.Particles(sim_dir)
    p = part.read(199)
    hp, sp = p[0], p[1]

    norm = mpl.colors.LogNorm(vmin=1, vmax=10000)
    kwargs = {"extent": [-lim, lim, -lim, lim],
              "norm": norm, "cmap": "inferno", "gridsize": 200}
    ax[0,0].hexbin(hp["x"][:,0], hp["x"][:,1], **kwargs)
    ax[1,0].hexbin(sp["x"][:,0], sp["x"][:,1], **kwargs)
    ax[0,1].hexbin(hp["x"][hp["smooth"],0], hp["x"][hp["smooth"],1], **kwargs)
    ax[1,1].hexbin(sp["x"][sp["smooth"],0], sp["x"][sp["smooth"],1], **kwargs)
    
    plt.savefig("particle_hist.png")

if __name__ == "__main__": main()
