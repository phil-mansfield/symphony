import symlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as random

try:
    import palette
    palette.configure(True)
except:
    pass

def main():
    sim_dir = symlib.get_host_directory(
        "/sdf/home/p/phil1/ZoomIns", "SymphonyLCluster", "Halo_042")
    i_sub = 15

    part = symlib.Particles(sim_dir)
    cores = part.core_indices(mode="smooth", halo=i_sub)

    x, ok = np.zeros((200, 32, 3)), np.zeros((200, 32), dtype=bool)
    
    for snap in range(200):
        p = part.read(snap, mode="smooth", halo=i_sub)
        x[snap], ok[snap] = p["x"][cores], p["ok"][cores]

    scale = symlib.scale_factors(sim_dir)
    def r(x): return np.sqrt(np.sum(x**2, axis=1))

    fig, ax = plt.subplots()

    rs, _ = symlib.read_rockstar(sim_dir)
    sy, _ = symlib.read_symfind(sim_dir)
    ok_rs, ok_sy = rs["ok"][i_sub], sy["ok"][i_sub]
    ax.plot(scale[ok_rs], r(rs["x"][i_sub,ok_rs]), lw=6, c="tab:red",
            label=r"${\rm Rockstar}$")
    ax.plot(scale[ok_sy], r(sy["x"][i_sub,ok_sy]), "--", lw=6, c="tab:blue",
            label=r"${\rm Symfind}$")
    ax.set_xlabel(r"$a(t)$")
    ax.set_ylabel(r"$r\ ({\rm physical}\ {\rm kpc})$")

    for i in range(32):
        ok_i = ok[:,i]
        ax.plot(scale[ok_i], r(x[ok_i,i]), lw=1, c="k", alpha=0.5)

    ax.legend(loc="lower left")

    fig.savefig("trajectory.png")

if __name__ == "__main__": main()
