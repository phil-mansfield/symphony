import numpy as np
import heapq
import symlib

def tag_energy_cut(p, E, r50, n_min=8):
    """ tag_energy_cut tags particles according to an energy cut. Invalid
    partcles need to be removed and the particles already need to be centered
    on their subhalos.
    """
    assert(len(p) > n_min)
    
    order = np.argsort(E)

    r = np.sqrt(np.sum(p["x"][order]**2, axis=1))
    ok = p["smooth"]

    r_med = running_median(r[ok])

    too_small = r_med <= r50
    too_big   = r_med >r50

    candidates = np.where(too_small[:-1] & too_big[1:])[0]

    if len(candidates) == 0 or candidates[-1] < n_min:
        n_star = n_min
    else:
        n_star = candidates[-1]

    mp = np.zeros(len(E))
    mp[order[:n_star]] = 1/n_star

    return mp

def tag_fixed_n(E, n_star):
    order = np.argsort(E)
    mp = np.zeros(len(E))
    mp[order[:n_core]] = 1/n_star

    return mp

def tag_nimbus(p, ranks, r_enc, m_enc):
    ok = p["smooth"]
    
    M = profile_matrix(p[ok], ranks, r_enc)
        
    n = M.shape[1]
    
    dm = np.zeros(len(m_enc))
    dm[1:] = m_enc[1:] - m_enc[:-1]
    dm[0] = m_enc[0]

    res = optimize.lsq_linear(
        M, dm, bounds=(np.zeros(n), np.inf*np.ones(n))
    )

    mp_star_table = res.x

    mp_star = np.zeros(len(p))
    mp_star[ok] = mp_star_table[ranks[ok]]
    mp_star[ok][ranks[ok] == -1] = 0

    m_tot = np.sum(mp_star)
    correction_frac = m_enc[-1]/m_tot if m_tot > 0 else 1.0

    mp_star *= correction_frac
    mp_star_table *= correction_frac
        
    return mp_star, mp_star_table

def profile_table(p, ranks, r_bins):
    n_rank = 1 + np.max(ranks)
    n_bins = len(r_bins) - 1

    M = np.zeros((n_bins, n_rank))

    for i in range(n_rank):
        ri = r[ranks[idx] == i]
        N, _ = np.histogram(ri, bins=r_bins)
        M[:,i] = N
        
    return M

def running_median(x):
    if len(x) == 0: return np.array([], dtype=x.dtype)
    med = np.zeros(len(x))

    low, high = [], [x[0]]

    med[0] = x[0]

    # This can be sped up with more conditionals that replace separate
    # heappush and heappop calls with a combined heappushpop.
    for i in range(1, len(x)):
        # Grow the correct queue.
        if x[i] <= med[i-1]:
            heapq.heappush(low, -x[i])
        else:
            heapq.heappush(high, x[i])

        # If one queue is too big, equalise them.
        if len(low) == len(high) + 2:
            xx = heapq.heappop(low)
            heapq.heappush(high, -xx)
        elif len(high) == len(low) +2:
            xx = heapq.heappop(high)
            heapq.heappush(low, -xx)

        # Evaluate median.
        if len(low) == len(high):
            med[i] = (-low[0] + high[0]) / 2
        elif len(low) == len(high) + 1:
            med[i] = -low[0]
        elif len(high) == len(low) + 1:
            med[i] = high[0]

    return np.array(med)

def r50(p, E, mp):
    ok = p["ok"] & (E < 0)
    p, mp = p[ok], mp[ok]
    
    r = radius(p["x"])
    
    order = np.argsort(r)
    r, mp = r[order], mp[order]
    
    # Compute the mass CDF.
    m_enc = np.cumsum(mp)/np.sum(mp)
    
    # We need to add zeros to the start to avoid crashing on a degenerate
    # case where the first particle has more than half the mass. Although in
    # such a case, we're probably doomed anyway.
    r, m_enc = np.hstack([[0], r]), np.hstack([[0], m_enc])
    
    # interpolate
    i = np.searchsorted(m_enc, 0.5)
    return r[i-1] + (0.5-m_enc[i-1])*(r[i]-r[i-1])/(m_enc[i]-m_enc[i-1])

def energy(param, p):
    rmax, vmax, phi_vmax2, _ = symlib.profile_info(param, p["x"]) 
    return np.sum(p["v"]**2, axis=1)/2 + phi_vmax2*vmax**2

def radius(x):
    return np.sqrt(np.sum(x**2, axis=1))

def center(p, halo):
    out = np.copy(p)
    out["x"] -= halo["x"]
    out["v"] -= halo["v"]
    return out

def expand(x, ok):
    out = np.ones(len(ok))*np.nan
    out[ok] = x
    return out
    
def main():
    sim_dir = "/fs/ddn/sdf/group/kipac/g/cosmo/ki21/phil1/simulations/ZoomIns/MWest/Halo004"

    # Read stuff in
    rs, _    = symlib.read_rockstar(sim_dir)
    sf, hist = symlib.read_symfind(sim_dir)
    part = symlib.Particles(sim_dir)
    param = symlib.simulation_parameters(sim_dir)
    
    # Time range to look at
    snap_i = hist["first_infall_snap"]
    snap_f = 235

    # Only look at these halos
    targets = np.where(sf["ok"][:10,-1])[0]

    print("""# 0 - subhalo index
# 1 - r50/rvir,infall
# 2 - m/m_peak,pre
# 3 - r50,target (kpc)
# 4 - r50,i (kpc)
# 5 - r50,f (kpac)""")
    
    for i in targets:
        # You need all the particles if you compute the energy yourself.
        pi = part.read(snap_i[i], mode="all", halo=i)
        pf = part.read(snap_f,    mode="all", halo=i)
        pi, pf = center(pi, sf[i,snap_i[i]]), center(pf, sf[i,snap_f])

        # Compute energies
        ok = pi["ok"]
        Ei = expand(energy(param, pi[ok]), ok)
        Ef = expand(energy(param, pf[ok]), ok)

        m_ratio = sf["m"][i,snap_f]/hist["mpeak_pre"]
        
        rvir_infall = rs["rvir"][i, snap_i[i]]
        r50_mults = [1/200, 1/100, 1/50, 1/20]
        
        for j, r50_mult in enumerate(r50_mults):
            rvir_infall = rs["rvir"][i,snap_i[i]]
            r50_target = rvir_infall * r50_mult

            mp = expand(tag_energy_cut(pi[ok], Ei[ok], r50_target), ok)
            
            r50_i = r50(pi[ok], Ei[ok], mp[ok])
            r50_f = r50(pf[ok], Ef[ok], mp[ok])
            
            print("%2d %.3f %.3f %.2f %.2f %.2f" %
                  (i, r50_mult, m_ratio[i], r50_target, r50_i, r50_f))

        print()        
    
if __name__ == "__main__": main()
