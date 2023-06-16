import numpy as np
from . import util
import scipy.spatial as spatial

class BinnedKDTrees(object):
    def __init__(self, x, mpeak, ok, bin_factor=1.25):
        log_mpeak = np.log10(mpeak)
        log_m_min, log_m_max = np.min(log_mpeak), np.max(log_mpeak)
        d_log_m = np.log10(bin_factor)

        n_bins = int(np.ceil((log_m_max - log_m_min) / d_log_m))

        self.bins = np.linspace(log_m_min, log_m_min+(n_bins+1)*d_log_m, n_bins+1)
        self.trees = [None]*n_bins
        self.idxs = [None]*n_bins

        for i in range(n_bins):
            ok_i = ok & (log_mpeak >= self.bins[i]) & (log_mpeak < self.bins[i+1])
            if np.sum(ok_i) > 0:
                self.idxs[i] = np.where(ok_i)[0]
                self.trees[i] = spatial.KDTree(x[ok_i])

    def best_match(self, x, mpeak, search_factor=2):
        m_low = np.log10(mpeak/search_factor)
        m_high = np.log10(mpeak*search_factor)
        
        idx_high = min(np.searchsorted(self.bins, m_high), len(self.bins)-1)
        idx_low = max(np.searchsorted(self.bins, m_low), 0)

        match_idxs, match_rs = [], []

        for i in range(idx_low, idx_high):
            if self.trees[i] is None: continue
            r, j = self.trees[i].query(x)
            match_rs.append(r)
            match_idxs.append(self.idxs[i][j])

        if len(match_rs) == 0:
            return -1

        return match_idxs[np.argmin(match_rs)]
        

def match_subhalos(h_1, h_2, hist_1, hist_2, min_votes=4):
    """ match_subhalos matches subhalos from simulation 1 against simulation 2
    using their pre-infall trajectories. h_1 and h_2 are 2D ROCKSTAR_DTYPE
    arrays, as returned by read_rockstar(). Any similarly structured
    2D array will work as long as it has a position field called "x" and a
    field called "ok" that indicates whether the subhalo exists. hist_1
    and hist_2 are HISTORY_DYTPE arrays, as returned by read_rockstar(). Any
    similarly formatted strucutred array will work as long as it also has
    a pre-infall mpeak field called "mpeak_pre" and an infall snapshot field
    called "first_infall_snap".

    This can be done with either comoving positions or physical displacements
    from the central halo, but we've noticed that the matching is usually
    more stable to chaotic changes in host position if comoving units are used.

    Two arrays are returned, match_1 and match_2, each of length |hist_1| and
    |hist_2|, respectively. They give the index of the corresponding match in
    the other simulation. this index will be unique and can be -1 if no
    reliable match can be found.
    
    The optional argument min_votes can be used to tune the matching. Higher
    values lead to fewer false positives and lower values lead to fewer false
    positives. 4 works well for SyphonyMilkyWay, but if looking at a simulation
    with a different output cadence, you should experiment with this. The only
    good way to test this that we are aware of is to qualitatively look at
    MAHs and movies, such as through the use of plot_matched_subhalos().
    """
    n_snap = h_1.shape[1]
    mpeak_1 = hist_1["mpeak_pre"]
    mpeak_2 = hist_2["mpeak_pre"]

    votes_1, votes_2 = [], []

    for snap in range(n_snap):
        match_1 = np.ones(len(h_1), dtype=int)*-1
        match_2 = np.ones(len(h_2), dtype=int)*-1

        last_snap_1 = hist_1["first_infall_snap"]
        last_snap_2 = hist_2["first_infall_snap"]
        last_snap_1[last_snap_1 < 0] = n_snap+1
        last_snap_2[last_snap_2 < 0] = n_snap+1

        ok_1 = (snap < hist_1["first_infall_snap"]) & h_1["ok"][:,snap]
        ok_2 = (snap < hist_2["first_infall_snap"]) & h_2["ok"][:,snap]
        trees_1 = BinnedKDTrees(h_1["x"][:,snap,:], mpeak_1, ok_1)
        trees_2 = BinnedKDTrees(h_2["x"][:,snap,:], mpeak_2, ok_2)

        for i in np.where(ok_1)[0]:
            match_1[i] = trees_2.best_match(h_1["x"][i,snap], mpeak_1[i])
        for i in np.where(ok_2)[0]:
            match_2[i] = trees_1.best_match(h_2["x"][i,snap], mpeak_2[i])
                    
        for i in range(len(match_1)):
            if match_2[match_1[i]] == i:
                votes_1.append(i)
                votes_2.append(match_1[i])

    votes_1 = np.array(votes_1, dtype=int)
    votes_2 = np.array(votes_2, dtype=int)        

    vote_count = { }

    for i in range(len(votes_1)):
        pair = (votes_1[i], votes_2[i])
        vote_count[pair] = 1 + vote_count.get(pair, 0)

    idx_1, idx_2, n_votes = process_vote_count(vote_count)
    order = np.argsort(n_votes)[::-1]
    idx_1, idx_2, n_votes = idx_1[order], idx_2[order], n_votes[order]

    match_1 = np.ones(len(h_1), dtype=int)*-1
    match_2 = np.ones(len(h_2), dtype=int)*-1
    in_use_1, in_use_2 = {}, {}
    for i in range(len(n_votes)):
        if idx_1[i] in in_use_1 or idx_2[i] in in_use_2 or n_votes[i] < min_votes:
            continue
        i1, i2 = idx_1[i], idx_2[i]
        match_1[i1] = i2
        match_2[i2] = i1
        in_use_1[i1] = None
        in_use_2[i2] = None
    
    return match_1, match_2
    
def process_vote_count(vote_count):
    pairs = vote_count.keys()
    idx_1, idx_2 = list(zip(*pairs))
    idx_1, idx_2 = np.array(idx_1, dtype=int), np.array(idx_2, dtype=int)
    n_votes = np.array([vote_count[pair] for pair in pairs], dtype=int)
    return idx_1, idx_2, n_votes

def plot_matched_subhalos(ax, h_1, h_2, match_1, match_2, snap, r_max,
                          plot_set_1=None, plot_set_2=None):
    """ plot_matched_subhalos plots plots the connections between subhalos
    according to the results of match_subhalos(). ax is the pyplot axis that
    will be used, h_1 and h_2 are 2D ROCKSTAR_DTYPE arrays from the matched
    simulations, formatted as returned by read_rockstar(). match_1 and match_2
    are the results of match_subhalos(), snap is the snapshot being plotted,
    and r_max is the maximum radius to plot relative to the host center, given
    in physical kpc.

    Two optional argument, plot_set_1 and plot_set_2 can be provided. They give
    the indices of all the subhaloes to be plotted.

    Subhalos in sim 1 will be plotted with solid lines and subhalos in sim 2
    will be plotted as dashed lines. Mached subhalos will be plotted in blue
    with a line connecting them. Subhaloes which are matched to a subhalo that
    doesn't exist in the current snapshot will be red, and subhalos which are
    not matched to anything will be plotted in black.

    This is indended to be used to generate a single frame in a movie.
    """
    if plot_set_1 is None: plot_set_1 = np.arange(len(h_1), dtype=int)
    if plot_set_2 is None: plot_set_2 = np.arange(len(h_2), dtype=int)

    ls_1, ls_2 = "-", "--"
    main_c = "tab:gray"
    surv_c, dis_c, unmatch_c = "tab:blue", "tab:red", "black"

    for i in plot_set_1:
        if not h_1["ok"][i,snap]: continue

        m = match_1[i]

        if i == 0:
            c = main_c
        elif m == -1:
            c = unmatch_c
        elif not h_2[m,snap]["ok"]:
            c = dis_c
        else:
            c = surv_c
            
        util.plot_circle(ax, h_1["x"][i,snap,0],
                         h_1["x"][i,snap,1],
                         h_1["rvir"][i,snap],
                         c=c, ls=ls_1, lw=1.5)
        ax.plot([h_1["x"][i,snap,0]], [h_1["x"][i,snap,1]], ".", c=main_c)

        if c == surv_c:
            util.plot_circle(ax, h_2["x"][m,snap,0],
                             h_2["x"][m,snap,1],
                             h_2["rvir"][m,snap],
                             c=c, ls=ls_2, lw=1.5)
            ax.plot([h_2["x"][m,snap,0]], [h_2["x"][m,snap,1]], ".", c=main_c)
            ax.plot([h_2["x"][m,snap,0], h_1["x"][i,snap,0]],
                     [h_2["x"][m,snap,1], h_1["x"][i,snap,1]],
                     lw=1, c=main_c)

    for i in plot_set_2:
        if not h_2["ok"][i,snap]: continue

        m = match_2[i]

        if i == 0:
            c = main_c
        elif m == -1:
            c = unmatch_c
        elif not h_1[m,snap]["ok"]:
            c = dis_c            
        else:
            c = surv_c

        util.plot_circle(ax, h_2["x"][i,snap,0],
                           h_2["x"][i,snap,1],
                           h_2["rvir"][i,snap],
                           c=c, ls=ls_2, lw=1.5)
        ax.plot([h_2["x"][i,snap,0]], [h_2["x"][i,snap,1]], ".", c=main_c)

        if c == surv_c:
            util.plot_circle(ax, h_1["x"][m,snap,0],
                               h_1["x"][m,snap,1],
                               h_1["rvir"][m,snap],
                               c=c, ls=ls_1, lw=1.5)
            ax.plot([h_1["x"][m,snap,0]], [h_1["x"][m,snap,1]], ".", c=main_c)
            ax.plot([h_1["x"][m,snap,0], h_2["x"][i,snap,0]],
                    [h_1["x"][m,snap,1], h_2["x"][i,snap,1]],
                    lw=1, c=main_c)

    ax.set_xlim(-r_max, +r_max)
    ax.set_ylim(-r_max, +r_max)
    ax.set_xlabel(r"$X\ ({\rm kpc})$")
    ax.set_ylabel(r"$Y\ ({\rm kpc})$")
