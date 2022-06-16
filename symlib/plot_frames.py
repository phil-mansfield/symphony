import numpy as np
import matplotlib.pyplot as plt
import palette
from palette import pc
import subhalo_tracking as sh
import os.path as path
import lib
import util
import matplotlib.colors as mpl_colors
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import subfind
import os

#TARGETS = [3, 20, 21, 322, 325]
TARGETS = [20, 21, 98, 99, 100, 110, 138, 145, 257, 325]
MIN_SNAP = 0
BASE_PLOT_DIR = "/oak/stanford/orgs/kipac/users/phil1/project_data/symphony_movies/SymphonyMilkyWay/Halo023/"

def frame_dir(i_sub):
    return path.join(BASE_PLOT_DIR, "sub_%04d" % i_sub)

def frame_name(i_sub, i_frame):
    return path.join(frame_dir(i_sub), "frame_%03d.png" % i_frame)

def get_c_range(x, frac):
    x = np.sort(x.flatten())
    x_sum = np.cumsum(x)
    x_tot = x_sum[-1]
    i_min = np.searchsorted(x_sum, x_tot*(1-frac)/2)
    i_max = np.searchsorted(x_sum, x_tot*(1+frac)/2)
    return x[i_min], x[-1]

def plot_r_max(r_min, xs, rs, mult=1.2):
    r_maxes = np.asarray(rs) + np.sqrt(np.sum(np.asarray(xs)**2, axis=1))
    
    return max(r_min, max(r_maxes))*mult

def sph_density_2d_proj(x, m, test, dim_x=0, dim_y=1, k=64, return_tree=False):
    xx = np.zeros((len(x), 2))
    xx[:,0], xx[:,1] = x[:,dim_x], x[:,dim_y]
    if test.shape[1] == 2:
        yy = test
    else:
        yy = np.zeros((len(test), 2))
        yy[:,0], yy[:,1] = test[:,dim_x], test[:,dim_y]

    return subfind.sph_density(xx, m, yy, k=k)

def plotting_grid(r, center, pts):
    x_edge = np.linspace(center[0] - r, center[0] + r, pts+1)
    y_edge = np.linspace(center[1] - r, center[1] + r, pts+1)
    x = (x_edge[1:] + x_edge[:-1]) / 2
    y = (x_edge[1:] + x_edge[:-1]) / 2
    xy, yx = np.meshgrid(x, y)
    xy, yx = xy.flatten(), yx.flatten()
    return np.vstack((xy, yx)).T


def main():
    palette.configure(False)
    
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns"
    suite_name = "SymphonyMilkyWay"
    sim_dir = path.join(base_dir, suite_name, "Halo023")
    
    param = lib.parameter_table[suite_name]
    h, hist = lib.read_subhalos(param, sim_dir)
    h_cmov = np.copy(h)
    info = lib.ParticleInfo(sim_dir)

    scale = lib.scale_factors(sim_dir)
    cosmo = cosmology.setCosmology("", lib.colossus_parameters(param))
    h = util.set_units_halos(h_cmov, scale, param)

    targets = TARGETS
    
    #targets = np.where(h["ok"][:,234] & h["ok"][:,235])[0]
    #targets = targets[targets >= 3]

    n_core = 32
    n_grid_proj = 200
    tracks, infall_tracks = [None]*len(h), [None]*len(h)

    max_snap = len(scale) - 1
    starting_snap = np.maximum(hist["merger_snap"], MIN_SNAP)
    infall_cores = [None]*len(h)

    print("""# 0 - snap
# 1 - subhalo index
# 2-4 - x (infall tracking)
# 5-7 - x (tree tracking)
# 8 - R_tidal
# 9 - R_tidal,klypin
# 10 - R_half,bound
# 11 - M_tidal
# 12 - M_tidal,klypin
# 13 - M_bound
# 14 - M_rockstar""")
    fig, ax = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    for i_sub in targets:
        os.makedirs(frame_dir(i_sub), exist_ok=True)
    
    for snap in range(np.min(starting_snap[targets]), max_snap + 1):

        sd = sh.SnapshotData(info, sim_dir, snap, scale[snap], h_cmov, param)
        prof = sh.MassProfile(sd.param, snap, h, sd.x, sd.owner, sd.valid)
        
        for j in range(len(targets)):
            i_sub = targets[j]

            if snap < starting_snap[i_sub]: continue

            elif tracks[i_sub] is None:
                infall_cores[i_sub] = sh.n_most_bound(
                    h["x"][i_sub,snap], h["v"][i_sub,snap],
                    sd.x[i_sub], sd.v[i_sub], sd.ok[i_sub],
                    n_core, sd.param
                )

                tracks[i_sub] = sh.SubhaloTrack(
                    i_sub, sd, infall_cores[i_sub], param)
                infall_tracks[i_sub] = sh.SubhaloTrack(
                    i_sub, sd, infall_cores[i_sub], param)
            else:
                infall_tracks[i_sub].next_snap(sd, infall_cores[i_sub])
                tracks[i_sub].next_snap(sd, tracks[i_sub].cores[snap-1])
                
            xp, vp = sd.x[i_sub], sd.v[i_sub]
            mp = np.ones(len(xp))*sd.mp
        
            x_t1 = infall_tracks[i_sub].x[snap]
            x_t2 = tracks[i_sub].x[snap]
            x_rs = h["x"][i_sub,snap]
        
            v_t1 = infall_tracks[i_sub].v[snap]
            v_t2 = tracks[i_sub].v[snap]

            dxp_t1, dvp_t1 = sh.delta(xp, x_t1), sh.delta(vp, v_t1)
            dxp_t2, dvp_t2 = sh.delta(xp, x_t2), sh.delta(vp, v_t2)

            r_rs, m_rs = h["rvir"][i_sub,snap], h["mvir"][i_sub,snap]

            ok, valid = sd.ok[i_sub], sd.valid[i_sub]

            bound_only=True
            r_t1, m_t1 = prof.tidal_radius(
                x_t1, xp[valid], v_t1, vp[valid], bound_only=False)
            r_t2, m_t2 = prof.tidal_radius(
                x_t2, xp[valid], v_t2, vp[valid], bound_only=False)
            r_t2_klypin, m_t2_klypin = prof.tidal_radius(
                x_t2, xp[valid], v_t2, vp[valid],
                bound_only=False, method="klypin99")
            

            is_bound,_ = sh.is_bound_iter(5, sd.param, dxp_t1, dvp_t1, ok=valid)
            m_bound = np.sum(mp[is_bound])
            r = np.sqrt(np.sum(dxp_t2**2, axis=1))
            r_half_bound = np.median(r[sd.valid[i_sub]])

            print("%d %d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.3g %.3g %.3g %.3g" %
                  (snap, i_sub, x_t1[0], x_t1[1], x_t1[2],
                   x_t2[0], x_t2[1], x_t2[2], r_t2, r_t2_klypin,
                   r_half_bound,
                   m_t1, m_t2_klypin, m_bound, m_rs))

            host_rvir = h["rvir"][0,snap]
            r_max = plot_r_max(host_rvir, [x_rs, x_t1, x_t2],
                               [r_rs, r_t1, r_t2])

            dim_xs = [0, 2]
            dim_ys = [1, 1]
            x_labels = [r"$X\,({\rm kpc})$", r"$Z\,({\rm kpc})$"]

            grid = plotting_grid(r_max, [0, 0], n_grid_proj)
            ax[0].cla()
            ax[1].cla()
            
            for i in range(2):
                dim_x, dim_y = dim_xs[i], dim_ys[i]

                # Plot host
                util.plot_circle(ax[i], 0, 0, host_rvir,
                                 c="k", lw=3, ls="--")

                # Plot the rockstar location
                util.plot_circle(ax[i], h["x"][i_sub,snap,dim_x],
                                 h["x"][i_sub,snap,dim_y], r_rs,
                                 c=pc("b"), lw=2)
            
                # Plot infall tracking
                util.plot_circle(ax[i], x_t1[dim_x], x_t1[dim_y], r_t1,
                                 c="w", lw=2)
                ax[i].plot(xp[infall_cores[i_sub],dim_x],
                           xp[infall_cores[i_sub],dim_y],
                           ".", alpha=0.4, c=pc("o"))

                # Plot boostrapped tracking
                util.plot_circle(ax[i], x_t2[dim_x], x_t2[dim_y], r_t2,
                                 c="w", lw=2, ls="--")
                
                log_rho_proj = np.log10(sph_density_2d_proj(
                    xp[ok], mp[ok], grid, k=64,
                    dim_x=dim_xs[i], dim_y=dim_ys[i]
                ))
                log_rho_proj = np.reshape(
                    log_rho_proj, (n_grid_proj, n_grid_proj))
            
                vmin, vmax = 1.5, 6.5

                ax[i].imshow(log_rho_proj, vmin=vmin, vmax=vmax,
                             origin="lower", cmap="magma_r",
                             extent=[-r_max, +r_max, -r_max, +r_max])
                    
                ax[i].set_xlim(-r_max, +r_max)
                ax[i].set_ylim(-r_max, +r_max)
                ax[i].set_xlabel(x_labels[i])
            
            ax[0].set_ylabel(r"$Y\,({\rm kpc})$")
            ax[0].set_title(r"$a(z) = %.3f$" % scale[snap])

            i_frame = snap - starting_snap[i_sub]
            fig.savefig(frame_name(i_sub, i_frame))
            
    plt.show()
        
if __name__ == "__main__": main()
