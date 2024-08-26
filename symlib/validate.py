import numpy as np
from . import lib
from . import util
import pathlib
import datetime
import symlib
import matplotlib.pyplot as plt
import os
import numpy.random as random
import shutil

def validate_symfind(base_dir, suite, halo, base_out_dir,
                     suffix="fid4",
                     resolution_bins=[3e2, 3e3, 3e4, 3e5, 3e6, 3e7, 3e8],
                     major_merger_ratio=0.15,
                     examples_per_bin=10,
                     seed=None):

    if seed is not None: random.seed(seed)
    
    # Handle different ways of supplying suite(s)
    if type(suite) == str:
        suites = [suite]
    else:
        suites = suite

    resolution_bins = np.array(resolution_bins)

    out_dir = pathlib.Path(base_out_dir, suffix)
    text_fname = pathlib.Path(out_dir, "overview.txt")
    
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    
    f = open(text_fname, "w+")
    def fprint(*args): print(*args, file=f)
    
    # Header to remind you about test configuration in case this gets dumped
    # into a context-less text file.
    now = datetime.datetime.now()
    fprint("Running validate_symfind", now.strftime("%m/%d/%Y %H:%M:%S"))
    fprint("  writing to %s" % base_out_dir)
    fprint("  suites = %s" % suites)
    fprint("  suffix = '%s'" % suffix)
    fprint("  resolution_bins =", resolution_bins)
    fprint("  major_merger_ratio = %f" % major_merger_ratio)
    fprint("  examples_per_bin = %d" % examples_per_bin)
    if seed is None:
        fprint("  RNG seed = None")
    else:
        fprint("  RNG seed = %d" % seed)

    # Binned data that will be accumulated
    n_bin = len(resolution_bins) - 1
    n_sub = np.zeros(n_bin)
    n_infall_err = np.zeros(n_bin)
    n_mm = np.zeros(n_bin)
    
    n_rs_ok = np.zeros(n_bin)
    n_rs_err = np.zeros(n_bin)
    n_sf_ok = np.zeros(n_bin)
    n_sf_err = np.zeros(n_bin)

    n_rs_ok_mm = np.zeros(n_bin)
    n_rs_err_mm = np.zeros(n_bin)
    n_sf_ok_mm = np.zeros(n_bin)
    n_sf_err_mm = np.zeros(n_bin)

    # Arrays that halo data will be siloed off into and then concatenated.
    suite_all, host_all, sub_all = [], [], []
    hist_all, rs_all, sf_all = [], [], []
    npeak_all = []

    # Loop through all targeted suites and halos. Compute basic diagnostics
    # and join subhalos into big arrays that can be concatenated.
    for suite in suites:
        # Handle different ways of supplying halo(s)
        if type(halo) in [int, str]:
            if halo == -1:
                halos = np.arange(util.n_hosts(suite), dtype=int)
            else:
                halos = [halo]
        else:
            halos = halo
            
        for halo in halos:
            # Read stuff in
            sim_dir = util.get_host_directory(base_dir, suite, halo)
            print("Testing", sim_dir)
            
            param = lib.simulation_parameters(sim_dir)
            mp = param["mp"]/param["h100"]

            rs, hist = lib.read_rockstar(sim_dir)
            sf, _ = lib.read_symfind(sim_dir, suffix=suffix)
            sub_idx = np.arange(len(sf))
            host = rs[0]
            # Remove the host
            rs, hist, sf = rs[1:], hist[1:], sf[1:]

            is_sub = after_infall(rs, hist)
            
            is_mm = hist["merger_ratio"] > major_merger_ratio
            npeak = hist["mpeak_pre"]/mp

            n_sub += count_subhalos(npeak, resolution_bins)
            n_mm += count_subhalos(npeak, resolution_bins, is_mm)
            
            n_infall_err += count_subhalos(
                npeak, resolution_bins, is_infall_error(sf, hist))

            is_ok_rs = rs["ok"] & is_sub
            is_ok_sf = sf["ok"] & is_sub
            is_err_rs = ((sf["is_err_rs"] & is_ok_rs) |
                         (is_ok_sf & (~is_ok_rs)))
            is_err_sf = ((sf["is_err"] & is_ok_sf) |
                         (is_ok_rs & (~sf["is_err_rs"]) & (~is_ok_sf)))
            
            n_rs_ok += count_subhalo_grid(
                npeak[~is_mm], resolution_bins, is_ok_rs[~is_mm])
            n_rs_err += count_subhalo_grid(
                npeak[~is_mm], resolution_bins, is_err_rs[~is_mm])
            n_sf_ok += count_subhalo_grid(
                npeak[~is_mm], resolution_bins, is_ok_sf[~is_mm])
            n_sf_err += count_subhalo_grid(
                npeak[~is_mm], resolution_bins, is_err_sf[~is_mm])

            n_rs_ok_mm += count_subhalo_grid(
                npeak[is_mm], resolution_bins, is_ok_rs[is_mm])
            n_rs_err_mm += count_subhalo_grid(
                npeak[is_mm], resolution_bins, is_err_rs[is_mm])
            n_sf_ok_mm += count_subhalo_grid(
                npeak[is_mm], resolution_bins, is_ok_sf[is_mm])
            n_sf_err_mm += count_subhalo_grid(
                npeak[is_mm], resolution_bins, is_err_sf[is_mm])
            
            hist_all.append(hist)
            rs_all.append(rs)
            sf_all.append(sf)
            npeak_all.append(npeak)
            suite_all.append([suite for _ in range(len(sf))])
            host_all.append([halo for _ in range(len(sf))])
            sub_all.append(sub_idx)

    # Join halos together

    npeak = np.concatenate(npeak_all)
    hist = np.concatenate(hist_all)
    rs = np.concatenate(rs_all, axis=0)
    sf = np.concatenate(sf_all, axis=0)
    suite, host = np.hstack(suite_all), np.hstack(host_all)
    sub = np.hstack(sub_idx)

    # Print simple diagnostics out
    
    ok = n_sub > 3
    ok_mm = n_mm > 3
    fprint("populated mass bins")
    fprint("low: ", np.array(resolution_bins)[:-1][ok])
    fprint("high: ", np.array(resolution_bins)[1:][ok])
    fprint()
    fprint("n_sub (subhalo count per bin)")
    fprint(n_sub[ok])
    fprint()
    fprint("f_mm (major merger fraction)")
    fprint(n_mm[ok]/n_sub[ok])
    fprint()
    fprint("f_infall_error (fraction of subhalos that disappear after infall)")
    fprint(n_infall_err[ok]/n_sub[ok])
    fprint()
    fprint("Minor merger 'error rates'")
    fprint("i.e. subhalos found in one catalog but not the other")
    fprint("Rockstar:")
    fprint(n_rs_err[ok]/n_rs_ok[ok])
    fprint("Symfind:")
    fprint(n_sf_err[ok]/n_sf_ok[ok])
    fprint()
    fprint("Major merger 'error rates'")
    fprint("This is a very meaningful qauntitfy for minor mergers, but I'm")
    fprint("less convinced it means something for major mergers. These things")
    fprint("sink directly into the host center and most remain bound the")
    fprint("entire time, so the termination condition is the main issue here.")
    fprint("Rockstar:")
    fprint(n_rs_err_mm[ok_mm]/n_rs_ok_mm[ok_mm])
    fprint("Symfind:")
    fprint(n_sf_err_mm[ok_mm]/n_sf_ok_mm[ok_mm])

    f.close()
    
    # Plot individual halos
    fig, ax = plt.subplots(3, 2, figsize=(14, 21), sharex=True)

    is_example = select_examples(npeak, resolution_bins, examples_per_bin)
    prev_suite, prev_host = None, None
    for i in range(len(is_example)):
        if not is_example[i]: continue

        curr_suite, curr_host = suite[i], host[i]
        if curr_suite != prev_suite or curr_host != prev_host:
            sim_dir = symlib.get_host_directory(base_dir, suite[i], host[i])
        #    part = symlib.Particles(sim_dir)
        prev_suite, prev_host = curr_suite, curr_host
        
        clear_axes(ax)

        scale = symlib.scale_factors(sim_dir)

        # I don't want to set rcParams because I don't know what environment
        # this is being run in. So this will be a bit annoying.

        r_rs = np.sqrt(np.sum(rs["x"][i]**2, axis=1))
        r_sf = np.sqrt(np.sum(rs["x"][i]**2, axis=1))
        
        ok = rs["ok"][i]
        ax[0,0].plot(scale[ok], rs["m"][i,ok], lw=3, c="tab:orange",
                     label=r"${\rm Rockstar\ (errors)}$")
        ax[1,0].plot(scale[ok], r_rs[ok], lw=3, c="tab:orange")
        
        ok = sf["ok_rs"][i]
        ax[0,0].plot(scale[ok], rs["m"][i,ok], lw=3, c="tab:red",
                     label=r"${\rm Rockstar\ (no\ errors)}$")
        ax[1,0].plot(scale[ok], r_rs[ok], lw=3, c="tab:red")
        
        ok = sf["ok"][i]
        ax[0,0].plot(scale[ok], sf["m"][i,ok], lw=3, c="tab:blue",
                     label=r"${\rm Symfind}$")
        ax[1,0].plot(scale[ok], r_sf[ok], lw=2, c="tab:blue")

        ax[0,0].legend(loc="lower left", fontsize=16)
        ax[0,0].set_ylabel(r"$m\ (M_\odot)$", fontsize=18)
        ax[0,0].set_yscale("log")

        ax[1,0].set_ylabel(r"$r\ ({\rm kpc})$", fontsize=18)
        ax[1,0].set_yscale("log")
        
        ax[2,0].set_xlabel(r"$a(t)$", fontsize=18)
        ax[2,0].set_xlim(0, 1)

        set_tick_fontsizes(ax, 18)
        
        fname = pathlib.Path(out_dir, "evo_%s_h%d_s%d.png" %
                             (suite[i], host[i], sub[i]))
        fig.savefig(fname)
        
def clear_axes(ax):
    for i in range(ax.shape[0]):
        for j in range(ax.shape[1]):
            ax[i,j].cla()
            
def set_tick_fontsizes(ax, fontsize):
    for i in range(ax.shape[0]):
        for j in range(ax.shape[1]):
            ax[i,j].tick_params(
                direction="in",
                labelsize=fontsize
            )            
    
            
def after_infall(rs, hist):
    snap = np.arange(rs.shape[1], dtype=int)
    out = np.zeros(rs.shape, dtype=bool)
    for i in range(len(rs)):
        out[i,snap >= hist["first_infall_snap"][i]] = True
    return out
    
def select_examples(npeak, resolution_bins, examples_per_bin):
    flags = np.zeros(len(npeak), dtype=bool)
    idx = np.arange(len(flags), dtype=int)
    for i in range(len(resolution_bins) - 1):
        ok = (npeak > resolution_bins[i]) & (npeak < resolution_bins[i+1])
        if np.sum(ok) == 0: continue

        n_example = min(np.sum(ok), examples_per_bin)
        idx_i = random.choice(idx[ok], n_example, replace=False)
        flags[idx_i] = True

    return flags
        
def count_subhalos(npeak, res_bins, flag=None):
    if flag is None: flag = np.ones(len(npeak), dtype=bool)
    return np.histogram(npeak[flag], res_bins)[0]

def count_subhalo_grid(npeak, res_bins, flag):    
    out = np.zeros(len(res_bins) - 1)
    for i in range(len(res_bins) - 1):
        ok = (res_bins[i] < npeak) & (npeak < res_bins[i+1])
        out[i] = np.sum(flag[ok])
    return out
        
def is_infall_error(sf, hist):
    out = np.zeros(len(hist), dtype=bool)
    for i in range(len(out)):
        snap = hist[i]["first_infall_snap"]+1
        snap = np.minimum(snap, sf.shape[1]-1)
        out[i] = not sf["ok"][i,snap]
    return out
