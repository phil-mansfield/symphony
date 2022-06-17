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
import sys

MIN_SNAP = 125

def main():
    palette.configure(False)

    out_file = sys.argv[1] 
    
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

    targets = np.arange(1, len(h), dtype=int)
    
    n_core = 32
    tracks = [None]*len(h)

    max_snap = len(scale) - 1
    starting_snap = np.maximum(hist["merger_snap"], MIN_SNAP)
    infall_cores = [None]*len(h)

    with open(out_file, "a") as fp:
        print("""# 0 - snap
# 1 - subhalo index
# 2-4 - x (pkpc)
# 5-7 - v (km/s)
# 8 - R_tidal (pkpc)
# 9 - R_50,bound (pkpc)
# 10 - R_95,bound (pkpc)
# 11 - M_tidal (Msun)
# 12 - M_tidal,bound (Msun)
# 13 - M_bound (msun)""", file=fp)

    for snap in range(np.min(starting_snap[targets]), max_snap + 1):

        sd = sh.SnapshotData(info, sim_dir, snap, scale[snap], h_cmov, param)
        prof = sh.MassProfile(sd.param, snap, h, sd.x, sd.owner, sd.valid)
        
        for j in range(len(targets)):
            i_sub = targets[j]
            if snap < starting_snap[i_sub]: continue

            if tracks[i_sub] is None:
                infall_cores[i_sub] = sh.n_most_bound(
                    h["x"][i_sub,snap], h["v"][i_sub,snap],
                    sd.x[i_sub], sd.v[i_sub], sd.ok[i_sub],
                    n_core, sd.param
                )

                tracks[i_sub] = sh.SubhaloTrack(
                    i_sub, sd, infall_cores[i_sub], param)
            else:
                tracks[i_sub].next_snap(sd, infall_cores[i_sub])
                
            xp, vp = sd.x[i_sub], sd.v[i_sub]
            mp = np.ones(len(xp))*sd.mp
        
            xc = tracks[i_sub].x[snap]
            vc = tracks[i_sub].v[snap]

            dxp, dvp = sh.delta(xp, xc), sh.delta(vp, vc)
            ok, valid = sd.ok[i_sub], sd.valid[i_sub]

            bound_only=True
            r_tidal, m_tidal = prof.tidal_radius(
                xc, xp[valid], vc, vp[valid], bound_only=False)            

            is_bound,_ = sh.is_bound_iter(10, sd.param, dxp, dvp, ok=valid)
            m_bound = np.sum(mp[is_bound])
            
            r = np.sqrt(np.sum(dxp**2, axis=1))
            if np.sum(is_bound) > 2:
                r_50_bound = np.quantile(r[is_bound], 0.5)
                r_95_bound = np.quantile(r[is_bound], 0.95)
            else:
                r_50_bound, r_95_bound = 0, 0

            m_tidal_bound = np.sum(mp[(r < r_tidal) & is_bound])
            
            with open(out_file, "a") as fp:
                print(("%d %d "+"%.4f "*3+"%.4f "*3+"%.4f "*3+"%.4g "*3) %
                      (snap, i_sub, xc[0], xc[1], xc[2], vc[0], vc[1], vc[2],
                       r_tidal, r_50_bound, r_95_bound,
                       m_tidal, m_tidal_bound, m_bound), file=fp)

        
if __name__ == "__main__": main()
