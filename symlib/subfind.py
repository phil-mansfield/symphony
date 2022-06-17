import numpy as np
import scipy.spatial as spatial
import scipy.special as special

def spline(x):
    r1 = x <= 1
    r2 = (x > 1) & (x <= 2)
    
    out = np.zeros(len(x))
    out[r1] = 1 - (3/2)*x[r1]**2 + 0.75*x[r1]**3
    out[r2] = 0.25*(2 - x[r2])**3

    return out
    
def rho_sph(ri, mi, n_dim):
    r = np.max(ri)/2
    Vr = np.pi**(n_dim/2)*r**n_dim / special.gamma(n_dim/2 + 1)
    # No need to add the central particle's mass: It's already in the sample
    # if it's a real particel and its mass is zero if it's a test point.
    return np.sum(mi*spline(ri/r))/Vr

def mean_sph(ri, mi, val):
    r = np.max(ri)/2
    w = mi*spline(ri/r)
    return np.sum(val*w)/np.sum(w)


def sph_densest_neighbor(x, rho, test, k=64, tree=None):
    k = min(k, len(x))
    if tree is None:
        tree = spatial.cKDTree(x)

    i_max = np.zeros(rho.shape, dtype=int)
    for i in range(len(test)):
        _, idx = tree.query(test[i], k)
        i_max[i] = idx[np.argmax(rho[idx])]

    return i_max
    
def sph_density(x, m, test, k=64, return_tree=False):
    k = min(k, len(x))
    tree = spatial.cKDTree(x)
    
    rho = np.zeros(len(test))
    for i in range(len(test)):
        ri, idx = tree.query(test[i], k)
        rho[i] = rho_sph(ri, m[idx], len(x[0]))

    if return_tree:
        return rho, tree
    else:
        return rho

def sph_scalar_mean(x, val, m, test, k=64, tree=None):
    k = min(k, len(x))
    if tree is None: tree = spatial.cKDTree(x)

    out = np.zeros(len(test))
    for i in range(len(test)):
        ri, idx = tree.query(test[i], k)
        out[i] = mean_sph(ri, m[idx], val[idx])
    
    return out
        
def sph_vector_mean(x, vec, m, test, k=64, tree=None):
    k = min(k, len(x))
    if tree is None: tree = spatial.cKDTree(x)

    out = np.zeros((len(test), len(vec[0])))
    for dim in range(len(vec[0])):
        out[:,dim] = sph_scalar_mean(x, vec[:,dim], m, test, k=k, tree=tree)
    return out

def subfind(x, v, mp, k=64):
    if type(mp) in [float, np.float32, np.float64]:
        mp = np.ones(len(x))*mp
        
    rho, tree = sph_density(x, mp, x, k=k, return_tree=True)
    
    order = np.argsort(rho)[::-1]
    orig_idx = order
    rank = np.zeros(len(order), dtype=int)
    rank[order] = np.arange(len(rank), dtype=int)

    i_max = sph_densest_neighbor(x, rho, x, k=k, tree=tree)
    is_peak = np.arange(len(rho), dtype=int) == i_max
    peaks = np.where(is_peak)[0]
    peak_n = np.zeros(len(peaks), dtype=int)
    parent = np.arange(len(peaks))
    parent_densities = rho[peaks]
    
    owner = np.ones(len(x), dtype=int) * -1
    owner[is_peak] = np.arange(len(peaks), dtype=int)

    order = np.argsort(rho)[::-1]

    for i_sort in range(len(rho)):
        i = orig_idx[i_sort]

        if owner[i] != -1: continue
        
        k = min(k, len(x))
        _, n_idx = tree.query(x[i], k)
        n_parents = parent[owner[n_idx]]
        n_parents = np.unique(n_parents[n_parents != -1])

        assert(len(n_parents) > 0)
        if len(n_parents) == 1:
            owner[i] = n_parents[0]
        else: # saddle point
            densest_parent = n_parents[np.argmax(parent_densities[n_parents])]
            #biggest_parent = n_parents[np.argmax(peak_n[n_parents])]
            for j in range(len(n_parents)):
                parent[parent == n_parents[j]] == densest_parent
                
            owner[i] = densest_parent

        peak_n[owner[i]] += 1
        
    xp = x[peaks]
    vp = sph_vector_mean(x, v, mp, xp, k=k, tree=tree)
    rho_p = rho[peaks]
    
    return xp, vp, rho_p, owner, peak_n

def main():
    import matplotlib.pyplot as plt
    import palette
    from palette import pc
    import util
    import lib
    import matplotlib.colors as mpl_colors
    from colossus.cosmology import cosmology
    from colossus.halo import mass_so
    import util
    import os.path as path
    
    palette.configure(False)
    
    base_dir = "/home/phil/code/src/github.com/phil-mansfield/symphony_pipeline/tmp_data"
    suite_name = "SymphonyMilkyWay"
    sim_dir = path.join(base_dir, suite_name, "Halo023")
    param = lib.parameter_table[suite_name]

    h, hist = lib.read_subhalos(param, sim_dir)
    h_cmov = np.copy(h)
    info = lib.ParticleInfo(sim_dir)
    k = 64
    scale = lib.scale_factors(sim_dir)

    cosmo = cosmology.setCosmology("", lib.colossus_parameters(param))
    
    for snap in [235]:
        x_snap = lib.read_particles(info, sim_dir, snap, "x")
        v_snap = lib.read_particles(info, sim_dir, snap, "v")
        valid_snap = lib.read_particles(info, sim_dir, snap, "valid")
        owner_snap = lib.read_particles(info, sim_dir, snap, "ownership")

        a = scale[snap]
        for i in range(len(x_snap)):
            x_snap[i] = util.set_units_x(x_snap[i], h_cmov[0,snap], a, param)
            v_snap[i] = util.set_units_v(v_snap[i], h_cmov[0,snap], a, param)

        mp, eps = util.set_units_param(a, param)
        h = util.set_units_halos(h_cmov, a, param)
        r_max = h[0,-1]["rvir"]*1.25

        rho_m = cosmo.rho_m(1/a-1) * param["h100"]**2
        
        for sub_i in [13]:
            x, v = x_snap[sub_i], v_snap[sub_i]
            valid, owner = valid_snap[sub_i], owner_snap[sub_i]

            ok = valid & (owner == 0)

            for k in [20, 40, 80]:
                xp, vp, rho_p, p_owner, p_count = subfind(x[ok], v[ok], mp, k=k)
                real = (p_count > 10) & (rho_p/rho_m > 10)
                
                print(rho_p[real]/rho_m)
                print(p_count[real])
                print(p_count[real])
                
                order = np.argsort(p_count)[::-1]
                orig_idx = np.arange(len(order), dtype=int)[order]
                plot_color = np.ones(len(p_owner), dtype=int)*-1
                
                colors = [pc("r"), pc("o"), pc("g"), pc("b"), pc("p")]
                for i in range(len(colors)):
                    plot_color[orig_idx[i] == p_owner] = i
                    
                fig, ax = plt.subplots()
                ax.plot(x[ok][plot_color == -1][:,0],
                        x[ok][plot_color == -1][:,1], ".", c=pc("k"), alpha=0.2)
                for i in range(len(colors)):
                    ax.plot(x[ok][plot_color == i][:,0],
                            x[ok][plot_color == i][:,1],
                            ".", c=colors[i], alpha=0.2)

                ax.set_xlim(-r_max, +r_max)
                ax.set_ylim(-r_max, +r_max)

    plt.show()
    
if __name__ == "__main__": main()
