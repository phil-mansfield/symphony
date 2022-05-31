import numpy as np
import scipy.spatial as spatial
import scipy.special as special
import os.path as path
import star_tagging
import scipy.signal as signal
import scipy.interpolate as interpolate

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

def density(x, m, test, k=64, return_tree=False):
    tree = spatial.cKDTree(x)
    
    rho = np.zeros(len(test))
    for i in range(len(test)):
        ri, idx = tree.query(test[i], k)
        rho[i] = rho_sph(ri, m[idx], len(x[0]))

    if return_tree:
        return rho, tree
    else:
        return rho

def density(x, m, test, k=64, return_tree=False):
    tree = spatial.cKDTree(x)
    
    rho = np.zeros(len(test))
    for i in range(len(test)):
        ri, idx = tree.query(test[i], k)
        rho[i] = rho_sph(ri, m[idx], len(x[0]))

    if return_tree:
        return rho, tree
    else:
        return rho

def neighbor_density(x, rho, test, k=64, tree=None):
    if tree is None:
        tree = spatial.cKDTree(x)

    rho_neighbor = np.zeros(rho.shape)
    for i in range(len(test)):
        _, idx = tree.query(test[i], k)
        rho_neighbor[i] = np.max(rho[idx])

    return rho_neighbor
    
def density_2d_proj(x, m, test, dim_x=0, dim_y=1, k=64, return_tree=False):
    xx = np.zeros((len(x), 2))
    xx[:,0], xx[:,1] = x[:,dim_x], x[:,dim_y]
    if test.shape[1] == 2:
        yy = test
    else:
        yy = np.zeros((len(test), 2))
        yy[:,0], yy[:,1] = test[:,dim_x], test[:,dim_y]

    return density(xx, m, yy, k=k)

def dispersion_scale(x, m):
    xc = np.zeros(len(x[0]))
    for dim in range(len(x[0])):
        xc[dim] = np.sum(x[:,dim]*m, axis=0)/np.sum(m)
    dx = np.zeros(x.shape)
    for dim in range(len(x[0])):
        dx[:,dim] = x[:,dim] - xc[dim]
    return np.sqrt(np.mean(dx**2))

def density_phase_space(x, v, m, test_x, test_v, k=64, return_tree=False):
    dx, dv = dispersion_scale(x, m), dispersion_scale(v, m)
    xx = np.zeros((len(x), len(x[0]) + len(v[0])))
    yy = np.zeros((len(test_x), len(x[0]) + len(v[0])))
    xx[:,:len(x[0])] = x / dx
    xx[:,len(x[0]):] = v / dv
    yy[:,:len(x[0])] = test_x / dx
    yy[:,len(x[0]):] = test_v / dv
    return density(xx, m, yy, k=k, return_tree=return_tree)

def sph_scalar(x, val, m, test, k=64, tree=None):
    if tree is None: tree = spatial.cKDTree(x)

    out = np.zeros(len(test))
    for i in range(len(test)):
        ri, idx = tree.query(test[i], k)
        out[i] = mean_sph(ri, m[idx], val[idx])
    
    return out
        
def sph_vector(x, vec, m, test, k=64, tree=None):
    if tree is None: tree = spatial.cKDTree(x)

    out = np.zeros((len(test), len(vec[0])))
    for dim in range(len(vec[0])):
        out[:,dim] = sph_scalar(x, vec[:,dim], m, test, k=k, tree=tree)
    return out

def sph_vector_phase_space(x, v, vec, m, test_x, test_v, k=64, tree=None):
    dx, dv = dispersion_scale(x, m), dispersion_scale(v, m)
    yy = np.zeros((len(test_x), len(x[0]) + len(v[0])))
    yy[:,:len(x[0])] = np.asarray(test_x) / dx
    yy[:,len(x[0]):] = np.asarray(test_v) / dv
    xx = np.zeros((len(x), len(x[0]) + len(v[0])))
    xx[:,:len(x[0])] = x / dx
    xx[:,len(x[0]):] = v / dv
    
    return sph_vector(xx, vec, m, yy, k=k, tree=tree)

def fix_units_x(x, h, scale, param):
    dx = np.zeros(x.shape)
    for dim in range(3): dx[:,dim] = x[:,dim] - h["x"][dim]
    dx *= scale*1e3/param["h100"]
    return dx

def fix_units_v(v, h, scale, param):
    v *= np.sqrt(scale)
    dv = np.zeros(v.shape)
    for dim in range(3): dv[:,dim] = v[:,dim] - h["v"][dim]
    return dv

def fix_units_param(scale, param):
    return param["mp"]/param["h100"], param["eps"]*scale/param["h100"]

def fix_units_halos(h, scale, param):    
    for dim in range(3):
        x0 = np.copy(h["x"][0,:,dim])
        v0 = np.copy(h["v"][0,:,dim])
        
        for hi in range(len(h)):
            h["x"][hi,:,dim] -= x0
            h["v"][hi,:,dim] -= v0
                    
    h["x"] *= 1e3*scale/param["h100"]
    h["rvir"] *= 1e3*scale/param["h100"]
    
    for hi in range(len(h)):
        invalid = h[hi]["rvir"] < 0
        h[hi, invalid]["rvir"] = -1
        h[hi, invalid]["x"] = -1
        h[hi, invalid]["v"] = -1
        
    return h

def capped_min(x, cap):
    ok = x > cap
    return np.min(x[ok])

class MassProfile(object):
    # Positions should already be in physical kpc, no h100.
    def __init__(self, eps, mp, snap, h, x, own, valid):
        r_low = eps / 10
        r_high = h[0,snap]["rvir"]*10
        r_edges = 10**np.linspace(np.log10(r_low), np.log10(r_high), 100)

        h0 = h[0,snap]["x"]
        dm_hist = np.zeros(len(r_edges) - 1)
        
        for i in range(len(x)):
            dx = np.copy(x[i])
            for dim in range(3): dx[:,dim] -= h0[dim]
            r = np.sqrt(np.sum(dx**2, axis=1))

            ok = (own[i] == 0) & (valid[i])
            n, _, = np.histogram(r[ok], bins=r_edges)
            dm_hist += n*mp
            dm_hist[0] += np.sum(r < r_low)

        dlog_r = np.log10(r_edges[1]) - np.log10(r_edges[0])
        mr = np.cumsum(dm_hist)
        
        ok = (mr > 0) & (r_edges[1:] > eps)
        dlog_m_dlog_r = signal.savgol_filter(
            np.log10(mr[ok]), 25, 4, deriv=1,
            delta=dlog_r, mode="interp"
        )
        dlog_m_dlog_r[dlog_m_dlog_r < 0] = 0

        log_r = np.log10(r_edges[1:][ok])
        
        self._deriv_interp = interpolate.interp1d(
            log_r, dlog_m_dlog_r, kind="linear",
            fill_value=(dlog_m_dlog_r[0], dlog_m_dlog_r[-1]),
            bounds_error=False
        )

        self._m_interp = interpolate.interp1d(
            log_r, np.log10(mr[ok]), kind="linear",
            fill_value=(np.log10(mr[ok])[0], np.log10(mr[ok])[-1]),
            bounds_error=False
        )

        self.mp = mp
        self.eps = eps

    def m(self, r):
        return 10**self._m_interp(np.log10(r))

    def dln_m_dln_r(self, r):
        return self._deriv_interp(np.log10(r))

    def _tidal_radius_iter(self, m, dx, dv=None, method="centrifugal"):
        """ Helper function for tidal_radius

        methods are:
        "jacobi" -  Eq. 2 from vdB & Ogiya 2018
        "radial" - Eq. 3 from vdB & Ogiya 2018
        "centrifugal" - Eq. 4 from vdB & Ogia 2018
        "circular" - Eq. 5 from vdB & Ogiya 2018
        "klypin99" - Eq. 6 from vdB & Ogiya 2018

        This is not vectorized
        """
        
        if dx.shape != (3,):
            raise ValueError("dx must be single 3-vector.")

        R = np.sqrt(np.sum(dx**2))
        dln_M_dln_R = self.dln_m_dln_r(R)
        M = self.m(R)
        
        if method == "jacobi":
            factor = 3
        elif method == "radial":
            factor = 2 - dln_M_dln_R
        elif method == "centrifugal":
            if dv is None:
                raise ValueError("dv must be set if method = 'centrifugal'")
            # In units of  (km/s * kpc)/kpc^2 = km/kpc/s:
            Omega = np.sqrt(np.sum(np.cross(dv, dr)**2))/R**2
            # In units of 1/Gyr (funny how closw those two numbers are...)
            Omega *= 1.022
            Omega_circ = 0.750 * (M/1e12)**0.5 * (R/200)**-1.5
            factor = 2 + (Omega/Omega_circ)**2 - dln_M_dln_R
        elif method == "circular":
            factor = 3 - dln_M_dln_R
        elif method == "klypin99":
            factor = 1
        else:
            raise ValueError("Unrecognized tidal radius method, '%s'" % method)

        if factor <= 0:
            return np.inf
        
        return R*(m/(factor*M))**(1.0/3)

    def tidal_radius(dx_core, dx, dv_core=None, dv=None, bound_only=False):
        if bound_only:
            raise ValueError("bound_only not yet supported")
        
        if dx_core[0] == 0 and dx_core[1] == 0 and dix_core[2] == 0:
            return 0.0

        dr2 = np.zeros(len(dx))
        for dim in range(3):
            dr2 += (dx[:,dim] - dx_core[:,dim])**2
        dr = np.sqrt(dr2)

        m = len(dr)*self.mp
        conv_limit = 0.001
        
        while True:
            r_tidal = self._tidal_radius_iter(m, dx_core, )
            n_tidal = np.sum(dr < r_tidal)
            if n_tidal == 0: return 0.0
            
            if np.abs(m/m_tidal - 1) < conv_limit:
                return r_tidal
            
            m = m_tidal

        return m_tidal
            
def eval_grid(r, center, pts):
    x = np.linspace(center[0] - r, center[0] + r, pts)
    y = np.linspace(center[1] - r, center[1] + r, pts)
    xy, yx = np.meshgrid(x, y)
    xy, yx = xy.flatten(), yx.flatten()
    return np.vstack((xy, yx)).T

def find_core_phase_space(x_own, v_own, m_own, k=64):
    rho, tree = density_phase_space(
        x_own, v_own, m_own, x_own, v_own,
        k=k, return_tree=True)

    i_max = np.argmax(rho)
    xc_init, vc_init = x_own[i_max], v_own[i_max]

    xc = sph_vector_phase_space(x_own, v_own, x_own, m_own,
                                [xc_init], [vc_init], k=k, tree=tree)[0]
    vc = sph_vector_phase_space(x_own, v_own, v_own, m_own,
                                [xc_init], [vc_init], k=k, tree=tree)[0]

    return xc, vc

def find_core(x_own, v_own, m_own, k=64):
    rho, tree = density(x_own, m_own, x_own,
                        k=k, return_tree=True)

    i_max = np.argmax(rho)
    xc_init, vc_init = x_own[i_max], v_own[i_max]

    xc = sph_vector(x_own, x_own, m_own,
                    [xc_init], k=k, tree=tree)[0]
    vc = sph_vector(x_own, v_own, m_own,
                    [xc_init], k=k, tree=tree)[0]

    return xc, vc


def find_cores(x_own, v_own, m_own, rho_cut, k=64):
    rho, tree = density(x_own, m_own, x_own,
                        k=k, return_tree=True)
    rho_neighbor = neighbor_density(x_own, rho, x_own, k=k, tree=tree)

    core_idx = np.where((rho_neighbor == rho) & (rho > rho_cut))[0]
    if len(core_idx) == 0:
        core_idx = np.array([np.argmax(rho)], dtype=int)
    
    core_rho = rho[core_idx]
    core_idx = core_idx[np.argsort(core_rho)][::-1]
    
    xc_init = x_own[core_idx]

    xc = sph_vector(x_own, x_own, m_own,
                    xc_init, k=k, tree=tree)
    vc = sph_vector(x_own, v_own, m_own,
                    xc_init, k=k, tree=tree)

    return xc, vc

def main():
    import matplotlib.pyplot as plt
    import palette
    from palette import pc
    import util
    import lib
    import matplotlib.colors as mpl_colors
    from colossus.cosmology import cosmology
    from colossus.halo import mass_so

    palette.configure(False)
    
    base_dir = "/home/phil/code/src/github.com/phil-mansfield/symphony_pipeline/tmp_data"
    suite_name = "SymphonyMilkyWay"
    sim_dir = path.join(base_dir, suite_name, "Halo023")
    param = lib.parameter_table[suite_name]

    h, hist = lib.read_subhalos(param, sim_dir)
    info = lib.ParticleInfo(sim_dir)
    k = 64
    n_iter = 10
    scale = lib.scale_factors(sim_dir)

    cosmo = cosmology.setCosmology("", lib.colossus_parameters(param))
    
    #for h0 in range(len(h)):
    #    if h[h0,-1]["mvir"] < 0: continue
    #    print("%3d %.3f" % (h0, h[h0,-1]["mvir"]/np.max(h[h0]["mvir"])))

    fig_proj, ax_proj = plt.subplots(1, 2, figsize=(16, 8), sharey=True)
    fig_den, ax_den = plt.subplots()
    fig_prof, ax_prof
    
    for snap in [235]:
        x_snap = lib.read_particles(info, sim_dir, snap, "x")
        v_snap = lib.read_particles(info, sim_dir, snap, "v")
        valid_snap = lib.read_particles(info, sim_dir, snap, "valid")
        owner_snap = lib.read_particles(info, sim_dir, snap, "ownership")
        
        a = scale[snap]
        for i in range(len(x_snap)):
            x_snap[i] = fix_units_x(x_snap[i], h[0,snap], a, param)
            v_snap[i] = fix_units_v(v_snap[i], h[0,snap], a, param)

        mp, eps = fix_units_param(a, param)
        h = fix_units_halos(h, a, param)
        
        r_max = h[0,snap]["rvir"]*1.5
        
        pts = 201
        grid = eval_grid(r_max, [0, 0], pts)

        rho_c = cosmo.rho_c(1/a - 1) * param["h100"]**2
        rho_cut = mass_so.deltaVir(1/a - 1)*rho_c
        
        #for h0 in range(1, len(h)):
        #for h0 in [3, 80, 110, 253, 254, 257, 398]:
        for h0 in [21]:
            if h[h0,snap]["rvir"] <= 0: continue
            print(h0)
            
            x, ok, owner = x_snap[h0], valid_snap[h0], owner_snap[h0]
            v = v_snap[h0]
            own = ok & (owner == 0)
            mp = np.ones(len(x))*param["mp"]

            #xc, vc = find_core(x[own], v[own], mp[own], k=k)
            xcs, vcs = find_cores(x[own], v[own], mp[own], rho_cut, k=k)
            xc, vc = xcs[0], vcs[0]
            rho_ok = density(x[ok], mp[ok], x[ok], k=k)
            rho_own = density(x[own], mp[own], x[own], k=k)
            
            dxc, dvc = np.zeros(x.shape), np.zeros(v.shape)
            for dim in range(3):
                dxc[:,dim] = x[:,dim] - xc[dim]
                dvc[:,dim] = v[:,dim] - vc[dim]
            rc = np.sqrt(np.sum(dxc**2, axis=1))
            rc = np.sqrt(rc**2 + eps**2)
            ke = 0.5*np.sum(dvc**2, axis=1)

            rmax, vmax, phi, order = star_tagging.profile_info(param, dxc[ok])
            e = phi+ke[ok]/vmax**2

            bound = e < 0
            n_bound = np.sum(bound)

            for _ in range(n_iter):
                rmax, vmax, phi, order = star_tagging.profile_info(
                    param, dxc[ok], ok=bound)
                e = phi+ke[ok]/vmax**2
                bound = e < 0
                if (np.sum(bound) == 0 or
                    abs(np.sum(bound) - n_bound)/n_bound < 0.01): break
                n_bound = np.sum(bound)

            rho_proj_1 = density_2d_proj(x[ok], mp[ok], grid, k=k)
            rho_proj_2 = density_2d_proj(x[ok], mp[ok], grid, k=k, dim_x=2)
            
            ax = ax_den
            ax.cla()
            
            own_ok = own[ok]
            bound_own = bound[own_ok]
            
            ax.plot(np.log10(rc[ok][(~bound) & (~own_ok)]),
                    np.log10(rho_ok[(~bound) & (~own_ok)]),
                    ".", alpha=0.2, c=pc("b"))
            ax.plot(np.log10(rc[own][~bound_own]),
                    np.log10(rho_own[~bound_own]),
                    ".", alpha=0.2, c="k")
            ax.plot( np.log10(rc[own][bound_own]),
                     np.log10(rho_own[bound_own]),
                     ".", alpha=0.2, c=pc("r"))
            ax.plot(np.log10(rc[ok][(bound) & (~own_ok)]),
                    np.log10(rho_ok[(bound) & (~own_ok)]),
                    ".", alpha=0.2, c=pc("o"))
            
            rho_min = np.log10(np.min(rho_ok))
            rho_max = np.log10(np.max(rho_ok))
            r_low, r_high = ax.get_xlim()
            ax.set_xlim(r_low, r_high)
            
            ax.fill_between([r_low, r_high], [rho_min]*2, [np.log10(rho_cut)]*2,
                            alpha=0.2, color="k")
            
            for i in range(1, len(xcs)):
                rc = np.sqrt(np.sum((xcs[i] - xc)**2))
                plt.plot([np.log10(rc)]*2, [rho_min, rho_max], "--",
                         lw=1, c=pc("b"))

            log_rvir = np.log10(h[h0,snap]["rvir"])
            ax.set_ylim(rho_min, rho_max)
            
            ax.set_xlabel(r"$\log_{10}(r)\,({\rm kpc})$")
            ax.set_ylabel(r"$\rho_{\rm dm}\,(M_\odot\,{\rm kpc}^{-3})$")
            
            ax1, ax2 = ax_proj
            ax1.cla()
            ax2.cla()
            
            rho_proj_1 = np.log10(np.reshape(rho_proj_1, (pts, pts)))
            rho_proj_2 = np.log10(np.reshape(rho_proj_2, (pts, pts)))

            c_range = get_c_range(10**rho_proj_1, 0.99)
            vmin, vmax = np.log10(c_range[0]), np.log10(c_range[1])
            vmin, vmax = vmin, 6.5
            
            ax1.imshow(rho_proj_1, vmin=vmin, vmax=vmax, origin="lower",
                      cmap="afmhot", extent=[-r_max, r_max, -r_max, r_max])
            util.plot_circle(ax1, 0, 0, h[0,snap]["rvir"], c="w", lw=2.5)
            util.plot_circle(ax2, 0, 0, h[0,snap]["rvir"], c="w", lw=2.5)
            
            rc = r_max/20
            for i in range(len(xcs)):
                util.plot_circle(ax1, xcs[i,0], xcs[i,1], rc, lw=1.5, c="w")
                util.plot_circle(ax1, xcs[i,0], xcs[i,1], rc, lw=1.5,
                                 c="k", ls="--")
            util.plot_circle(ax1, xcs[0,0], xcs[0,1],
                             rc*1.25, c=pc("b"), lw=1.5)
            util.plot_circle(ax1, h[h0,snap]["x"][0], h[h0,snap]["x"][1],
                             rc*1.5, c=pc("p"), lw=1.5)
            ax1.set_xlim(-r_max, +r_max)
            ax1.set_ylim(-r_max, +r_max)
            
            ax1.set_xlabel(r"$X\,(h^{-1}{\rm Mpc})$")
            ax1.set_ylabel(r"$Y\,(h^{-1}{\rm Mpc})$")

            ax2.imshow(rho_proj_2, vmin=vmin, vmax=vmax, origin="lower",
                      cmap="afmhot", extent=[-r_max, r_max, -r_max, r_max])
            util.plot_circle(ax1, 0, 0, h[0,snap]["rvir"], c="w")
            
            for i in range(len(xcs)):
                util.plot_circle(ax2, xcs[i,2], xcs[i,1], rc, lw=1.5, c="w")
                util.plot_circle(ax2, xcs[i,2], xcs[i,1], rc, lw=1.5,
                                 c="k", ls="--")
            util.plot_circle(ax2, xcs[0,2], xcs[0,1],
                             rc*1.25, c=pc("b"), lw=1.5)
            util.plot_circle(ax2, h[h0,snap]["x"][2], h[h0,snap]["x"][1], 
                             rc*1.5, c=pc("p"), lw=1.5)
            
            ax2.set_xlabel(r"$Z\,(h^{-1}{\rm Mpc})$")
            ax2.set_xlim(-r_max, +r_max)
            ax2.set_ylim(-r_max, +r_max)
            
            fig_den.savefig(
                "../../symphony_pipeline/plots/" +
                "subhalo_tracking/h%03d_den.png" % h0
            )
            fig_proj.savefig(
                "../../symphony_pipeline/plots/" +
                "subhalo_tracking/h%03d_proj.png" % h0
            )
    plt.show()
    
def get_c_range_hist(x, y, mp, r_max, bins, frac):
    range = ((-r_max, r_max), (-r_max, r_max))
    H, _, _ = np.histogram2d(x, y, weights=mp, range=range, bins=bins)
    H = np.sort(H.flatten())

    H_sum = np.cumsum(H)
    H_tot = H_sum[-1]
    i_min = np.searchsorted(H_sum, H_tot*(1 - frac)/2)
    i_max = np.searchsorted(H_sum, H_tot*(1 + frac)/2)

    return H[i_min], H[i_max]

def get_c_range(x, frac):
    x = np.sort(x.flatten())
    x_sum = np.cumsum(x)
    x_tot = x_sum[-1]
    i_min = np.searchsorted(x_sum, x_tot*(1-frac)/2)
    i_max = np.searchsorted(x_sum, x_tot*(1+frac)/2)
    return x[i_min], x[-1]

def main2():
    import matplotlib.pyplot as plt
    import palette
    from palette import pc
    import util
    import lib
    import matplotlib.colors as mpl_colors
    from colossus.cosmology import cosmology
    from colossus.halo import mass_so

    palette.configure(False)
    
    base_dir = "/home/phil/code/src/github.com/phil-mansfield/symphony_pipeline/tmp_data"
    suite_name = "SymphonyMilkyWay"
    sim_dir = path.join(base_dir, suite_name, "Halo023")
    param = lib.parameter_table[suite_name]

    h, hist = lib.read_subhalos(param, sim_dir)
    info = lib.ParticleInfo(sim_dir)
    k = 64
    n_iter = 10
    scale = lib.scale_factors(sim_dir)

    cosmo = cosmology.setCosmology("", lib.colossus_parameters(param))

    for snap in [235]:
        x_snap = lib.read_particles(info, sim_dir, snap, "x")
        v_snap = lib.read_particles(info, sim_dir, snap, "v")
        valid_snap = lib.read_particles(info, sim_dir, snap, "valid")
        owner_snap = lib.read_particles(info, sim_dir, snap, "ownership")

        a = scale[snap]
        for i in range(len(x_snap)):
            x_snap[i] = fix_units_x(x_snap[i], h[0,snap], a, param)
            v_snap[i] = fix_units_v(v_snap[i], h[0,snap], a, param)

        mp, eps = fix_units_param(a, param)
        h = fix_units_halos(h, a, param)

        prof = MassProfile(eps, mp, snap, h, x_snap, owner_snap, valid_snap)

        r = 10**np.linspace(-1, 3, 200)
        dln_m_dln_r = prof.dln_m_dln_r(r)
        m  = prof.m(r)

        plt.figure()
        plt.plot(r, dln_m_dln_r, c=pc("r"))
        plt.xscale("log") 
 
        plt.figure()
        plt.plot(r, m, c=pc("r"))
        plt.yscale("log")
        plt.xscale("log")
        
        plt.xscale("log")
        
    plt.show()


if __name__ == "__main__":
    main()
    #main2()
