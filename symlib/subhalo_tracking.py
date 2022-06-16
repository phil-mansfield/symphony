import numpy as np
import matplotlib.pyplot as plt
import lib
import util
import star_tagging
import subfind
import matplotlib.colors as mpl_colors
from colossus.cosmology import cosmology
from colossus.halo import mass_so
import os.path as path
import palette
from palette import pc
import scipy.signal as signal
import scipy.interpolate as interpolate

def delta(x, x0):
    dx = np.zeros(x.shape)
    for dim in range(3):
        dx[:,dim] = x[:,dim] - x0[dim]
    return dx

def distance(x, x0):
    return np.sqrt(np.sum(delta(x, x0)**2, axis=1))

def n_most_bound(xc, vc, x, v, ok, n_core, param):
    dx, dv = delta(x, xc), delta(v, vc)

    ke = 0.5*np.sum(dv**2, axis=1)
    _, vmax, pe_scaled, _ = star_tagging.profile_info(param, dx, ok=ok)
    E = ke + pe_scaled*vmax**2
    
    order = np.argsort(E[ok])
    orig_idx = np.arange(len(x), dtype=int)[ok][order]
    
    assert(len(orig_idx) >= n_core)
    return orig_idx[:n_core]

def is_bound(param, dx, dv, ok=None, order=None):
    rmax, vmax, pe, order = star_tagging.profile_info(param, dx, ok, order)
    pe *= vmax**2
    ke = 0.5*np.sum(dv**2, axis=1)

    if ok is None:
        return pe + ke < 0, order
    else:
        return (pe + ke < 0) & ok, order

def is_bound_iter(n_iter, param, dx, dv, ok=None, order=None):
    ok, order = is_bound(param, dx, dv, ok, order)
    n_bound = np.sum(ok)

    for i in range(n_iter - 1):
        ok, order = is_bound(param, dx, dv, ok, order)
        n_bound_i = np.sum(ok)
        if n_bound_i == n_bound: break
        
    return ok, order

def rockstar_cores(snap_info, h, sub_idxs, n_core):
    core_idxs = np.ones((len(h), n_core), dtype=int)*-1
    param = snap_info.param

    snap = snap_info.snap
    
    for i_sub in sub_idxs:
        x, v = snap_info.x[i_sub], snap_info.v[i_sub]
        valid, owner =  snap_info.valid[i_sub], snap_info.owner[i_sub]
        ok = valid & (owner == 0)
        
        xc, vc = h["x"][i_sub,snap], h["v"][i_sub,snap]
        core_idxs[i_sub,:] = n_most_bound(xc, vc, x, v, ok, n_core, param)

        r = distance(x, xc)
        
    return core_idxs

class SnapshotData(object):
    def __init__(self, info, sim_dir, snap, a, h_cmov, param):
        self.x = lib.read_particles(info, sim_dir, snap, "x")
        self.v = lib.read_particles(info, sim_dir, snap, "v")
        self.valid = lib.read_particles(info, sim_dir, snap, "valid")
        self.owner = lib.read_particles(info, sim_dir, snap, "ownership")
        self.ok = [None]*len(self.x)
        
        for i in range(len(self.x)):
            self.x[i] = util.set_units_x(self.x[i], h_cmov[0,snap], a, param)
            self.v[i] = util.set_units_v(self.v[i], h_cmov[0,snap], a, param)
            self.ok[i] = self.valid[i] & (self.owner[i] == 0)
            
        self.mp, self.eps = util.set_units_param(a, param)
        self.snap = snap
        self.param = param

class SubhaloTrack(object):
    def __init__(self, i_sub, snap_data_init, core, param):
        x, v = snap_data_init.x[i_sub], snap_data_init.v[i_sub]
        
        self.i_sub = i_sub

        snap = snap_data_init.snap
        n_snap = param["n_snap"]
        n_core = len(core)
        self.n_core = n_core
        self.snaps = np.arange(n_snap, dtype=int)
        self.cores = np.ones((n_snap, n_core), dtype=int)*-1
        self.x = np.ones((n_snap, 3))*np.nan
        self.v = np.ones((n_snap, 3))*np.nan
        self.param = param

        self._set_xv(snap_data_init, core)
        
    def next_snap(self, snap_data, core):    
        snap = snap_data.snap
        self._set_xv(snap_data, core)

    def _set_xv(self, snap_data, prev_core):
        snap = snap_data.snap
        core = self.cores[snap]
    
        x, v = snap_data.x[self.i_sub], snap_data.v[self.i_sub]
        valid, owner = snap_data.valid[self.i_sub], snap_data.owner[self.i_sub]
        mp, eps = snap_data.mp, snap_data.eps
        
        ok = valid & (owner == 0)

        x_sf, v_sf, rho_sf, owner_sf, peak_n = subfind.subfind(x[ok], v[ok], mp)
        owner_all = np.ones(len(x), dtype=int)*-1
        owner_all[ok] = owner_sf
        
        owner_votes = owner_all[core]
        assert(-1 not in owner_votes)
        i_owner = np.argmax(np.bincount(owner_votes))
        xc_sf, vc_sf = x_sf[i_owner], v_sf[i_owner]

        core_sf = n_most_bound(xc_sf, vc_sf, x, v, ok, self.n_core, self.param)
        
        self.cores[snap] = core_sf
        self.x[snap] = xc_sf
        self.v[snap] = vc_sf   
        
class MassProfile(object):
    # Positions should already be in physical kpc, no h100.
    def __init__(self, param, snap, h, x, own, valid):
        r_low = param["eps"] / 10
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
            dm_hist += n*param["mp"]
            dm_hist[0] += np.sum(r < r_low)

        dlog_r = np.log10(r_edges[1]) - np.log10(r_edges[0])
        mr = np.cumsum(dm_hist)
        
        ok = (mr > 0) & (r_edges[1:] > param["eps"])
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

        self.param = param
        self.mp = param["mp"]
        self.eps = param["eps"]

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
            Omega = np.sqrt(np.sum(np.cross(dv, dx)**2))/R**2
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
        elif factor < 1:
            factor = 1.0
        
        return R*(m/(factor*M))**(1.0/3)

    def tidal_radius(self, dx_core, dx, dv_core=None, dv=None,
                     method="centrifugal", bound_only=True):
        if dx_core[0] == 0 and dx_core[1] == 0 and dx_core[2] == 0:
            return 0.0

        dx = delta(dx, dx_core)
        dv = delta(dv, dv_core)
        
        dr = np.sqrt(np.sum(dx**2, axis=1))
        order = np.argsort(dr)

        m = len(dr)*self.mp
        conv_limit = 0.001
        
        while True:
            r_tidal = self._tidal_radius_iter(m, dx_core, dv_core,
                                              method=method)
            ok = dr < r_tidal
            n_tidal = np.sum(ok)
            if bound_only:
                is_bound = is_bound_iter(5, self.param, dx, dv,
                                         ok=ok, order=order)
            else:
                is_bound = ok
                
            m_tidal = np.sum(is_bound)*self.param["mp"]

            if m_tidal == 0:
                return 0.0, 0.0
            if np.abs(m/m_tidal - 1) < conv_limit:
                return r_tidal, m_tidal

            
            m = m_tidal

        assert(0)

        
def main():
    palette.configure(False)
    
    base_dir = "/home/phil/code/src/github.com/phil-mansfield/symphony_pipeline/tmp_data"
    suite_name = "SymphonyMilkyWay"
    sim_dir = path.join(base_dir, suite_name, "Halo023")
    
    param = lib.parameter_table[suite_name]
    h, hist = lib.read_subhalos(param, sim_dir)
    h_cmov = np.copy(h)
    info = lib.ParticleInfo(sim_dir)

    scale = lib.scale_factors(sim_dir)
    cosmo = cosmology.setCosmology("", lib.colossus_parameters(param))
    h = util.set_units_halos(h_cmov, scale, param)
    
    targets = np.where(h["ok"][:,234] & h["ok"][:,235])[0]
    targets = targets[targets >= 3]
    
    snap234 = SnapshotData(info, sim_dir, 234, scale[234], h_cmov, param)
    snap235 = SnapshotData(info, sim_dir, 235, scale[235], h_cmov, param)

    n_core = 32
    
    rs_cores = rockstar_cores(snap234, h, targets, n_core)
    infall_cores = lib.read_particles(info, sim_dir, -1, "infall_core")
    infall_cores = infall_cores[:,:n_core]

    infall_tracks = [None]*len(targets)
    rs_tracks = [None]*len(targets)

    prof235 = MassProfile(snap235.param, 235, h, snap235.x,
                          snap235.owner, snap235.valid)
    
    for i in range(len(targets)):
        i_sub = targets[i]
        print("%3d" % i_sub, end=" ")
        infall_tracks[i] = SubhaloTrack(i_sub, snap234, infall_cores, param)
        rs_tracks[i] = SubhaloTrack(i_sub, snap234, rs_cores, param)

        x = h["x"][i_sub,235]
        print("%8.3f %8.3f %8.3f" % (x[0], x[1], x[2]), end=" ")
        
        infall_tracks[i].next_snap(snap235)
        x = infall_tracks[i].x[235]
        print("%8.3f %8.3f %8.3f" % (x[0], x[1], x[2]), end=" ")
        
        rs_tracks[i].next_snap(snap235)
        x = rs_tracks[i].x[235]
        print("%8.3f %8.3f %8.3f" % (x[0], x[1], x[2]), end=" ")

        v = rs_tracks[i].v[235]
        dx = delta(snap235.x[i_sub], x)
        dv = delta(snap235.v[i_sub], v)

        ok = snap235.ok[i_sub]
        rs_r_tidal, rs_m_tidal = prof235.tidal_radius(
            x, snap235.x[i_sub][ok], v, snap234.v[i_sub][ok], bound_only=False)
        
        rs_is_bound, _ = is_bound_iter(10, snap235.param, dx, dv, ok=ok)
        rs_dr = np.sqrt(np.sum(dx**2, axis=1))
        rs_m_bound = snap235.param["mp"]*np.sum(rs_is_bound)
        
        if rs_m_bound == 0:
            rs_r99 = 0.0
        else:
            rs_r99 = np.quantile(rs_dr[rs_is_bound], 0.99)
            
        print("%7.3f %7.3f %7.3f %.3g %.3g %.3g" % (
            (h["rvir"][i_sub,235], rs_r_tidal, rs_r99,
             h["mvir"][i_sub,235], rs_m_tidal, rs_m_bound)
        ))
        
if __name__ == "__main__": main()
