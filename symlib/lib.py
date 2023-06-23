import array
import struct
import numpy as np
import os
import os.path as path
import scipy.interpolate as interpolate
import numpy.random as random
import scipy.optimize as optimize
import glob
from . import util

""" SUBHALO_DTYPE is the numpy datatype used by the main return value of
read_subhalos().
"""
SUBHALO_DTYPE = [("id", "i4"), ("mvir", "f4"), ("vmax", "f4"), ("rvmax", "f4"),
                 ("x", "f4", (3,)), ("v", "f4", (3,)), ("ok", "?"),
                 ("rvir", "f4"), ("cvir", "f4")]


""" ROCKSTAR_DTYPE is the main numpy datatype returned by read_rockstar().
 - id: unique identifier for the halo.
 - m (Msun): mass of the subhalo.
 - vmax (km/s): maximum value of the subhalo's rotation curve.
 - rvmax (physical kpc): radius of vmax.
 - x (physical kpc): displacement from the host's centre to the subhalo's
   centre.
 - v (km/s) velocity of the subhalo relative to the host's centre.
 - rvir (physical kpc): radius of the subhalo according to Bryan & Norman 1998.
   This is rarely a meaningful quantity for subhaloes as their tidal radii are
   usually smaller than rvir.
 - cvir: the NFW concentration parameter, r_-2/rvir for the subhalo, as
   estimated from vmax/vvir. This is rarely a meaninfiul quantity for
   subhaloes, because their profiles are very different from NFW..
 - ok: True if the subhalo exists and False otherwise. Note that if Rockstar
   has any errors but still believes this subhalo is real, this will be true.
   This can only be calculated for haloes which symfind has been run on
   You must use the HISTORY_DTYPE or SYMFIND_DTYPE to identify
"""
ROCKSTAR_DTYPE = [("id", "i4"), ("m", "f4"), ("vmax", "f4"), ("rvmax", "f4"),
                 ("x", "f4", (3,)), ("v", "f4", (3,)), ("ok", "?"),
                 ("rvir", "f4"), ("cvir", "f4")]

""" HISTORY_DTYPE is a numpy datatype representing data about a subhalo which 
is indepdent of time, such as the maximum mass that it takes on, its merger
snapshot, and its location in the full merger tree file.

 - The main branch of a given halo within the depth-first merger tree can be 
   found with x[start:end].
 - mpeak: (Msun) the peak value of Mvir across the entire branch
 - mpeak_pre: (Msun) the peak value of Mvir before becoming a subhalo
 - vpeak: (km/s) the peak value of vmax across the entire branch
 - merger_snap: the snapshot that the subhalo first became a subhalo of the
   main host halo.
 - first_infall_snap: the snapshot that the subhalo first became a subhalo of
   any other object, excluding glitches during major mergers.
 - conv_snap_*: the snaphot at which the subhalo is predicted to be
   unconverged, according to various criteria. See mansfield et al. 2023 for
   discussion
 - disrupt_snap: the snapshot at which the subhalo disrupts in the symfind
   catalogue
 - disrupt_snap_rs: the snapshot at which the subhalo disrupts in Rockstar + 
   consistent-trees. If the Rockstar has erroneous subhalo
   identifications at the end of its branch, this flags the last non-error
   snapshot.
 - several other non-informative fields are included for consistency with
   BRANCH_DTYPE.
"""
HISTORY_DTYPE = [("mpeak", "f4"), ("vpeak", "f4"), ("merger_snap", "i4"),
                 ("merger_ratio", "f4"),
                 ("start", "i4"), ("end", "i4"), ("is_real", "?"),
                 ("is_disappear", "?"), ("is_main_sub", "?"),
                 ("preprocess", "i4"), ("first_infall_snap", "i4"),
                 ("branch_idx", "i4"), ("false_selection", "?"),
                 ("mpeak_pre", "f4"), ("conv_snap_discrete", "i4"),
                 ("conv_snap_eps", "i4"), ("conv_snap_relax", "i4"),
                 ("conv_snap_relax_hydro", "i4"), ("disrupt_snap", "i4"),
                 ("disrupt_snap_rs", "i4")]


""" BRANCH_DTYPE is a numpy datatype representing the main branch of of halo's
merger tree. This forms a component of HISTORY_DTYPE, which 
 - The main branch of a given halo within the depth-first merger tree can be 
   found with x[start:end].
 - is_real: false if a halo is definitely a numerical artefact (e.g. if it's
   already a subhalo in its first snapshot).
 - is_disappear: true if the subhalo disappears without merging. This is also
   bad and probably means the thing is a numerical artefact.
 - is_main_sub: true is the halo was ever a subhalo of the zoom-in box's main 
   halo (includes splashback subhaloes).
 - preprocess: the index of the branch of the largest halo that hosted this 
   halo before it became a subhalo of the main halo. If the halo was never a
   subhalo of the main halo, this is just the index of largest halo to have
   every hosted this halo. If no other halo every hosted this halo before
   it entered Rvir of the main halo, this is -1. Includes splashback subhaloes
   and doesn't include collisions between subhaloes once they've already been
   accreted.
 - branch_idx: the index of this halo in the full merger tree.
"""
BRANCH_DTYPE = [("start", "i4"), ("end", "i4"), ("is_real", "?"),
                 ("is_disappear", "?"), ("is_main_sub", "?"),
                 ("preprocess", "i4"), ("first_infall_snap", "i4")]


""" CORE_DTYPE is a numpy datatype representing the properties of the
particle-tracked subhalo.
TODO: UPDATE THIS COMMENT
 - x: position in pkpc
 - v: velocity in km/s
 - r_tidal: tidal radius in pkpc. Accounts for angular momentum and the
   mass profile of the central halo.
 - r50_bound: radius enclosing 50% of bound particles
 - r50_bound_rs: radius enclosing 50% of bound particles
 - r95_bound: radius enclosing 95% of bound particles
 - m_tidal: total mass within tidal radius
 - m_tidal_bound: total bound mass within tidal radius
 - m_bound: total bound mass
 - ok: true if the core is being tracked and false otherwise
 - ok_rs: true if a rockstar halo exists and at least one core particle is
   still inside Rockstar's r_50 (almost all cases where this fails are cases
   where Rockstar has made an error)
 - is_err: Rockstar and core-tracking disagree, but more core particles are
   assoicated with Rockstar
 - is_err_rs: Rockstar and core-tracking disagree, but more core particles are
   with the particle-tracking subhalo.
 - interp: true if the core in this snapshot is interpolated betwen two
   nearby snapshots.
"""
CORE_DTYPE = [("x", "f4", (3,)), ("v", "f4", (3,)), ("r_tidal", "f4"),
              ("r50_bound", "f4"), ("r50_bound_rs", "f4"), ("m_tidal", "f4"),
              ("m_tidal_bound", "f4"), ("m_bound", "f4"), ("vmax", "f4"),
              ("f_core", "f4"), ("f_core_rs", "f4"), ("d_core_mbp", "f4"),
              ("ok", "?"), ("ok_rs", "?"), ("is_err", "?"), ("is_err_rs", "?"),
              ("interp", "?")]

""" SYMFIND_DTYPE is the main numpy datatype returned by read_symfind().
 - x (physical kpc): displacement from the host halo to the subhalo.
 - v (km/s): velocity relative to the host.
 - m (Msun): bound mass of the subhalo. Post-porcessing has been done to reduce
   noise by enfocing monotonicity (i.e. m(z) = max(m_raw(>= z)))
 - m_raw (Msun): bound mass of the suhalo with no post-processing.
 - r_half (physical kpc) half-mass radius of the subhalo.
 - r_tidal (physical kpc): tidal radius of the subhalo. Accounts for angular
   momentum and the mass profile of the central halo. Tidal radii can be
   poor approximations if the size of the subhalo is comparable to the
   distance to the host centre or if the subhalo is at pericentre. There is no
   strict criteria for identifying this.
 - m_tidal (Msun): total mass within tidal radius.
 - m_tidal_bound (Msun): total bound mass within tidal radius, using only
   particles within the tidal radius. This value can be a bit fragile at
   pericentre.
 - ok: true if the subhalo exists and is not an error, and false otherwise.
   Interpolated subhaloes will be marked as true.
 - interp: true if the subhalo in this snapshot is interpolated betwen two
   nearby snapshots.

Some values for Rockstar haloes are also calculated here:
 - r_half_rs (physical kpc): half-mass radius of tracked particles using
   Rockstar's position and velocity. This is the r_half used for error
   detection, not the r_half in the Rockstar catalogues.
 - ok_rs: true if a rockstar halo exists and at least one core particle is
   still inside Rockstar's r_half. A slightly more conservative way to
   identify Rockstar errors would be to check if ok_rs & (~is_err_rs).

Some diagnostic values are included that will not be useful for a typical user:
 - is_err: Rockstar and symfind disagree, but more core particles are
   assoicated with Rockstar. If this is true, the subhalo is either
   interpolated or flagged as disrupted.
 - is_err_rs: Rockstar and symfind disagree, but more core particles are
   associated with symfind.
 - f_core: the fraction of core particles associated with Rockstar
 - f_core_rs: the fraction of core particles associated with symfind.
 - d_core_mbp (physical kpc): the distance between the subhalo and the
   most-bound particle. Sometimes the most-bound particle can numerically
   diffuse out of the subhalo.
"""
SYMFIND_DTYPE = [("x", "f4", (3,)), ("v", "f4", (3,)), ("r_tidal", "f4"),
                 ("r_half", "f4"), ("m_tidal", "f4"),
                 ("m_tidal_bound", "f4"), ("m", "f4"), ("vmax", "f4"),
                 ("m_raw", "f4"),
                 ("ok", "?"), ("interp", "?"),
                 ("f_core", "f4"), ("f_core_rs", "f4"), ("d_core_mbp", "f4"),
                 ("r_half_rs", "f4"), ("ok_rs", "?"),
                 ("is_err", "?"), ("is_err_rs", "?")]

""" UM_DTYPE is a numpy datatype representing UniverseMachine predictions for
galaxy properties.
- m_star: the stellar mass of the galaxy in Msun
- m_in_situ: the stellar mass of a galaxy that formed inside the galaxy (as
             opposed to coming in via mergers) in Msun
- m_icl: the stellar mass of the galaxy's stellar halo
- sfr: the starformaiton rate in Msun/yr
- mvir: predicted virial mass of the halo in msun
- vmax: predicted V_max of the halo in km/s
- x: position of the galaxy relative to the host halo in pkpc
- v: velocity of the galaxy relative to the host halo in km/s
- rank: the galaxy's rank in delta vmax
- is_orphan: true if the the galaxy is un-tracked and replaced with an orphan
             model during this timestep and false otherwise
- ok: true if UniverseMachine has an entry for this galaxy at this snapshot and
      false otherwise
"""
UM_DTYPE = [("m_star", "f4"), ("m_in_situ", "f4"), ("m_icl", "f4"),
            ("sfr", "f4"), ("x", "f4", (3,)), ("v", "f4", (3,)),
            ("rank", "f4"), ("mvir", "f4"), ("vmax", "f4"),
            ("is_orphan", "?"), ("ok", "?")]

""" PARTICLE_DTYPE is a numpy datatype representing the properties of tracked 
particles within a subhalo. Positions (x) and velocities (v) are given in
physical units.
- x: position in pkpc
- v: velocity in km/s
- snap: integer value for the desired snapshot
 - ok: true if the particle is being tracked and false otherwise
 - smooth: true if the particle smoothly accreted onto the subhalo and false otherwise.
"""
PARTICLE_DTYPE = [
    ("x", "f4", (3,)), ("v", "f4", (3,)),
    ("snap", "i4"), ("id", "i4"),
    ("ok", "?"), ("smooth", "?")
]

""" TREE_COL_NAMES is the mapping of variable names to columns in the
consistent-trees file. These are the variable names you need to pass to
read_tree().
"""
TREE_COL_NAMES = {
    "dfid": 28,
    "id": 1,
    "desc_id": 3,
    "upid": 6,
    "phantom": 8,
    "snap": 31,
    "next_co_prog": 32,
    "mvir": 10,
    "rs": 12,
    "vmax": 16,
    "m200b": 39,
    "m200c": 40,
    "m500c": 41,
    "xoff": 43,
    "spin_bullock": 45,
    "b_to_a": 46,
    "c_to_a": 47,
    "t_to_u": 56,
    "r_vmax": 60,
    "x": 17,
    "v": 20,
    "j": 23,
    "a": 48,
}
_chinchilla_cosmology = { 'flat': True, 'H0': 70.0, 'Om0': 0.286,
                          'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95 }
_banerjee_cosmology = { 'flat': True, 'H0': 70.0, 'Om0': 0.3,
                        'Ob0': 0.049, 'sigma8': 0.85, 'ns': 0.95 }
_carmen_cosmology = { 'flat': True, 'H0': 70.0, 'Om0': 0.25,
                      'Ob0': 0.049, 'sigma8': 0.8, 'ns': 1 }

""" parameter_table contains simulation-defining parameters. This includes a 
table of cosmological parameters ("flat", "H0", "Om0", "Ob0", "sigma8", and
"ns"), the particle mass in Msun/h, "mp", the Plummer-equivalent force
softening scale in comoving kpc/h, "eps", and for convience, "h100". "n_snap" is
the number of snapshots in the suite.
"""
parameter_table = {
    "SymphonyLMC": _chinchilla_cosmology.copy(),
    "SymphonyMilkyWay": _chinchilla_cosmology.copy(),
    "SymphonyMilkyWayDisk": _chinchilla_cosmology.copy(),
    "SymphonyMilkyWayDiskDMO": _chinchilla_cosmology.copy(),
    "SymphonyMilkyWayLR": _chinchilla_cosmology.copy(),
    "SymphonyMilkyWayHR": _chinchilla_cosmology.copy(),
    "SymphonyGroup": _chinchilla_cosmology.copy(),
    "SymphonyLCluster": _banerjee_cosmology.copy(),
    "SymphonyCluster": _carmen_cosmology.copy(),
    "MWest": _chinchilla_cosmology.copy(),
}

parameter_table["SymphonyLMC"]["eps"] = 0.080
parameter_table["SymphonyMilkyWay"]["eps"] = 0.170
parameter_table["SymphonyMilkyWayDisk"]["eps"] = 0.170*2
parameter_table["SymphonyMilkyWayDiskDMO"]["eps"] = 0.170*2
parameter_table["SymphonyMilkyWayLR"]["eps"] = 0.170
parameter_table["SymphonyMilkyWayHR"]["eps"] = 0.170/2
parameter_table["SymphonyGroup"]["eps"] = 0.360
parameter_table["SymphonyLCluster"]["eps"] = 1.200
parameter_table["SymphonyCluster"]["eps"] = 3.250
parameter_table["MWest"]["eps"] = 0.170

parameter_table["SymphonyLMC"]["mp"] = 3.52476e4
parameter_table["SymphonyMilkyWay"]["mp"] = 2.81981e5
parameter_table["SymphonyMilkyWayDisk"]["mp"] = 2.81981e5*8
parameter_table["SymphonyMilkyWayDiskDMO"]["mp"] = 2.81981e5*8
parameter_table["SymphonyMilkyWayLR"]["mp"] = 2.81981e5
parameter_table["SymphonyMilkyWayHR"]["mp"] = 2.81981e5/8
parameter_table["SymphonyGroup"]["mp"] = 2.25585e6
parameter_table["SymphonyLCluster"]["mp"] = 1.51441632e8
parameter_table["SymphonyCluster"]["mp"] = 1.26201360e8
parameter_table["MWest"]["mp"] = 2.81981e5

parameter_table["SymphonyLMC"]["n_snap"] = 236
parameter_table["SymphonyMilkyWay"]["n_snap"] = 236
parameter_table["SymphonyMilkyWayDisk"]["n_snap"] = 236
parameter_table["SymphonyMilkyWayDiskDMO"]["n_snap"] = 236
parameter_table["SymphonyMilkyWayLR"]["n_snap"] = 236
parameter_table["SymphonyMilkyWayHR"]["n_snap"] = 236
parameter_table["SymphonyGroup"]["n_snap"] = 236
parameter_table["SymphonyLCluster"]["n_snap"] = 200
parameter_table["SymphonyCluster"]["n_snap"] = 200
parameter_table["MWest"]["n_snap"] = 236

for sim in parameter_table:
    param = parameter_table[sim]
    param["h100"] = param["H0"]/100

parameter_table["ExampleSuite"] = parameter_table["MWest"]

def suite_names():
    """ suite_names returns the names of all the simulation suites.
    """
    return util.SUITE_NAMES

def simulation_parameters(suite_name):
    """ parameter_table contains simulation-defining parameters. This includes a
    table of cosmological parameters ("flat", "H0", "Om0", "Ob0", "sigma8", and
    "ns"), the particle mass in Msun/h, "mp", the Plummer-equivalent force
    softening scale in comoving kpc/h, "eps", and for convience, "h100".
    "n_snap" is the number of snapshots in the suite.
    """
    if suite_name in parameter_table:
        return parameter_table[suite_name]
    else:
        sim_dir = suite_name
        suite_name = halo_dir_to_suite_name(sim_dir)
        if suite_name not in parameter_table:
            raise ValueError("'%s' is neither a recognized suite name nor a simulation direoctry containing a recognized suite name. Recognized suite names are: " + str(suite_names()))
        return parameter_table[suite_name]

def colossus_parameters(param):
    """ colossus_parameters converts a parameter dictionary into a colossus-
    comparitble parameter dictionary.
    """
    var_names = ["Om0", "flat", "H0", "Ob0", "sigma8", "ns"]
    return { name: param[name] for name in var_names }

P_FILE_FMT = "%s/particles/part_%03d.%d"

def halo_dir_to_suite_name(dir_name):
    """ halo_dir_to_suite_name returns the name of the suite that a halo in the
    given directory belongs to.
    """
    if dir_name[-1] == "/": dir_name = dir_name[:-1]
    suite_dir, halo_name = os.path.split(dir_name)
    base_dir, suite_name = os.path.split(suite_dir)
    return suite_name

_SCALE_FACTOR_CACHE = { }

def scale_factors(dir_name):
    """ scale_factors returns the scale factors of each snapshot for a halo in
    the given directory.
    """
    if dir_name in _SCALE_FACTOR_CACHE:
        return _SCALE_FACTOR_CACHE[dir_name]

    # TODO: individual halo-by-halo scale factors
    suite_name = halo_dir_to_suite_name(dir_name)
    if (suite_name in ["SymphonyLMC", "SymphonyGroup", "SymphonyMilkyWayDisk",
                       "SymphonyMilkyWayDiskDMO",
                      "SymphonyMilkyWay", "MWest", "SymphonyMilkyWayLR",
                       "SymphonyMilkyWayHR"] or
        dir_name in ["SymphonyLMC", "SymphonyGroup", "SymphonyMilkyWayDisk",
                       "SymphonyMilkyWayDiskDMO",
                     "SymphonyMilkyWay", "MWest", "SymphonyMilkyWayLR",
                     "SymphonyMilkyWayHR"]):
        default = 10**np.linspace(np.log10(0.05), np.log10(1), 236)
    elif (suite_name in ["SymphonyLCluster",  "SymphonyCluster"] or
          dir_name in ["SymphonyLCluster",  "SymphonyCluster"]):
        default = 10**np.linspace(np.log10(0.075), np.log10(1), 200)
    else:
        raise ValueError(("The halo in %s does not belong to a " + 
                          "recognized suite.") % dir_name)

    if dir_name in ["SymphonyLMC", "SymphonyGroup",
                    "SymphonyMilkyWay", "MWest", "SymphonyMilkyWayLR",
                    "SymphonyMilkyWayHR", 
                    "SymphonyLCluster",  "SymphonyCluster"]:
        return default

    file_name = path.join(dir_name, "halos", "snap_scale.dat")
    if not path.exists(file_name):
        _SCALE_FACTOR_CACHE[dir_name] = default
        return default

    f = open(file_name, "rb")
    n_snap = struct.unpack("q", f.read(8))[0]
    scale = np.fromfile(f, np.float64, n_snap)
    f.close()

    _SCALE_FACTOR_CACHE[dir_name] = scale
    return scale


def mvir_to_rvir(mvir, a, omega_M):
    """ mvir_to_rvir converts a Bryan & Norman virial mass in Msun/h to a virial
    radius in comoving Mpc/h at a given scale factor a, and omega_M.
    """
    
    if type(mvir) == np.ndarray:
        mvir = np.maximum(np.zeros(len(mvir)), mvir)
    else:
        mvir = max(0, mvir)

    omega_L = 1 - omega_M
    Ez = np.sqrt(omega_M/a**3 + omega_L)
    rho_crit = 2.77519737e11*Ez**2
    omega_Mz = (omega_M/a**3)/Ez**2

    x = omega_Mz - 1
    delta_vir = 18*np.pi**2 + 82*x - 39.0*x**2
    rho_vir = rho_crit*delta_vir

    r_phys = (mvir/(rho_vir * (4*np.pi / 3)))**(1.0/3)
    r_cmov = r_phys/a

    return r_cmov

def flatten(arrays):
    """ flatten takes a list of numpy arrays and flattens it into a single
    array. arrays[i].shape[0] can be any value, but all other components of the
    shape vectors must be the same. THis is needed because hstack doesn't work
    on tensor-arrays.
    """

    N = sum(arr.shape[0] for arr in arrays)

    shape, dtype = arrays[0].shape, arrays[0].dtype
    if len(shape) == 1:
        out = np.zeros(N, dtype=dtype)
    else:
        out = np.zeros((N,) + shape[1:], dtype=dtype)

    start, end = 0, 0
    for i in range(len(arrays)):
        end += arrays[i].shape[0]
        out[start: end] = arrays[i]
        start = end

    return out

def merger_snap(h, x_sub, snap_sub):
    """ merger_snap returns the snapshot where a subhalo with a given set of
    positions and snapshots first falls within h, a halo (i.e. a return value
    of read_mergers()).
    """
    dist = np.sqrt(np.sum((h["x"][snap_sub] - x_sub)**2, axis=1))
    within = dist < h["rvir"][snap_sub]
    merger = h["ok"][snap_sub] & within

    if np.sum(merger) == 0:
        return -1
    else:
        return np.min(snap_sub[merger])

def pristine_merger_indices(b):
    return np.where(b["is_real"] & (~b["is_disappear"]) &
                    b["is_main_sub"] & (b["preprocess"] == -1))[0]

def merger_stats(b, m, x, mvir, snap):
    """ merger_stats returns several useful values: mpeak of the subhalo, the
    scale factor of the merger, the snapshot of the merger, and the mass ratio
    between the host and subhalo at the merger snapshot.
    """
    sub_idx = pristine_merger_indices(b)
    mw = m[0]

    ratio = np.zeros(len(sub_idx))
    m_snap = np.zeros(len(sub_idx), dtype=int)
    mpeak = np.zeros(len(sub_idx))
    
    mw_mass = np.max(mw["mvir"])
    for j, i in enumerate(sub_idx):
        mvir_i = mvir[b["start"][i]: b["end"][i]]
        snap_i = snap[b["start"][i]: b["end"][i]]
        x_i = x[b["start"][i]: b["end"][i]]
                
        m_snap_i = merger_snap(mw, x_i, snap_i)

        m_snap_sub = np.searchsorted(snap_i[::-1], m_snap_i)
        mpeak[j] = np.max(mvir_i)
        m_snap[j] = m_snap_i
        ratio[j] = mvir_i[::-1][m_snap_sub]/mw["mvir"][m_snap_i]

    return mpeak, m_snap, ratio, sub_idx

ROCKSTAR_DTYPE = [("id", "i4"), ("m", "f4"), ("vmax", "f4"), ("rvmax", "f4"),
                 ("x", "f4", (3,)), ("v", "f4", (3,)), ("ok", "?"),
                 ("rvir", "f4"), ("cvir", "f4")]

def read_rockstar(dir_name, comoving=False):
    h, hist = read_subhalos(dir_name, comoving=comoving)
    r = np.zeros(h.shape, dtype=ROCKSTAR_DTYPE)
    r["id"], r["m"], r["vmax"] = h["id"], h["mvir"], h["vmax"]
    r["rvmax"], r["x"], r["v"] = h["rvmax"], h["x"], h["v"]
    r["ok"], r["rvir"], r["cvir"] = h["ok"], h["rvir"], h["cvir"]
    return r, hist

def read_subhalos(dir_name, comoving=False, include_false_selections=False):
    """ read_subhalos reads major merger data from the halo directory dir_name.
    It returns two arrays. The first, m_idx, is the indices of the major
    mergers within the branches arrays. Index 0 is the main halo, index 1 is the
    biggest merger (by Mpeak), index 2 is the second biggest, etc. As a
    convenience, data from the main branches of those haloes are returned as
    the second argument, m. m[halo_idx, snap] gives the  properties of the
    halo at the index halo_idx in m_idx and the snapshot snap. The fields are
    given by SUBHALO_DTYPE (see the header of this file). If a halo doesn't
    exist at a given snapshot, values are set -1 ok will be set to false.

    If comoving is True, the subhaloes will be returned in Rockstar's default,
    comoving units. Otherwise it's returned in symlib's units. If
    include_false_slections is true, haloes which would be included with a
    straight M_peak > 300*mp cutoff are included. If False (the default)
    """
    params = simulation_parameters(dir_name)

    fname = path.join(dir_name, "halos", "subhalos.dat")
    f = open(fname, "rb")

    n_snap = struct.unpack("i", f.read(4))[0]
    n_merger = struct.unpack("i", f.read(4))[0]

    idx = np.fromfile(f, np.int32, n_merger)
    out = np.zeros((n_merger, n_snap), dtype=SUBHALO_DTYPE)
    a = scale_factors(dir_name)

    for i in range(n_merger):
        out["mvir"][i,:] = np.fromfile(f, np.float32, n_snap)
        out["ok"][i,:] = out["mvir"][i,:] > 0
        out["rvir"][i,:] = mvir_to_rvir(out["mvir"][i,:], a, params["Om0"])

    out["rvir"][out["rvir"] == 0] = -1
        
    for i in range(n_merger):
        out["vmax"][i,:] = np.fromfile(f, np.float32, n_snap)
        ok = out["rvir"][i,:] > 0
        out["cvir"][i,ok] = vmax_to_cvir_nfw(
            out["vmax"][i,ok], out["mvir"][i,ok], out["rvir"][i,ok]*a[ok]
        )

    for i in range(n_merger):
        out["rvmax"][i,:] = np.fromfile(f, np.float32, n_snap)
    for i in range(n_merger):
        out["id"][i,:] = np.fromfile(f, np.int32, n_snap)
    for i in range(n_merger):
        out["x"][i,:,:] = np.fromfile(f, (np.float32, (3,)), n_snap)
    for i in range(n_merger):
        out["v"][i,:,:] = np.fromfile(f, (np.float32, (3,)), n_snap)
        
    f.close()

    histories = get_subhalo_histories(out, idx, dir_name)

    file_name = path.join(dir_name, "halos", "core_pos.dat")
    if path.exists(file_name):
        with open(file_name, "rb") as f:
            n_halo, n_snap = struct.unpack("qq", f.read(16))
            x = np.fromfile(f, np.float32, 3*n_halo*n_snap)
            v = np.fromfile(f, np.float32, 3*n_halo*n_snap)
            x = x.reshape((n_halo, n_snap, 3))
            v = v.reshape((n_halo, n_snap, 3))
        
    if not comoving:
        scale = scale_factors(dir_name)
        out = util.set_units_halos(out, scale, params)
        histories = util.set_units_histories(histories, scale, params)

    if not include_false_selections:
        out = out[~histories["false_selection"]]
        histories = histories[~histories["false_selection"]]

    return out, histories

def get_subhalo_histories(s, idx, dir_name):
    if s["x"][0,-1,0] == 0:
        raise ValueError("The input subhalo array has already been " + 
                         "converted to symlib's default units.")

    param = simulation_parameters(dir_name)
    b = read_branches(dir_name)
    h = np.zeros(len(s), dtype=HISTORY_DTYPE)

    # These can just be copied over
    h["start"], h["end"] = b[idx]["start"], b[idx]["end"]
    h["is_real"], h["is_disappear"] = b[idx]["is_real"], b[idx]["is_disappear"]
    h["is_main_sub"] = b[idx]["is_main_sub"]
    # TODO: point into subhalos, not branches
    h["preprocess"] = b[idx]["preprocess"]
    h["first_infall_snap"] = b[idx]["first_infall_snap"]
    h["branch_idx"] = idx

    central = s[0]

    snap = np.arange(len(s[0]))
    
    for i in range(len(s)):
        mvir, vmax, x = s[i]["mvir"], s[i]["vmax"], s[i]["x"]
        ok = mvir > 0
        
        h["mpeak"][i], h["vpeak"][i] = np.max(mvir[ok]), np.max(vmax[ok])
        if i == 0: continue

        m_snap = merger_snap(central, x[ok], snap[ok])
        m_ratio = mvir[m_snap] / central["mvir"][m_snap]

        target = (snap >= m_snap) & s["ok"][i]

        if h["first_infall_snap"][i] == -1: h["first_infall_snap"][i] = m_snap
        h[i]["merger_snap"], h[i]["merger_ratio"] = m_snap, m_ratio
        h["first_infall_snap"][i] = min(h["merger_snap"][i],
                                        h["first_infall_snap"][i])

        infall_snap = h["first_infall_snap"][i]
        mpeak_infall = np.max(s[i,:infall_snap+1]["mvir"])
        h["mpeak_pre"][i] = mpeak_infall
        h["false_selection"][i] = mpeak_infall < 300*param["mp"]

    h["mpeak_pre"][0] = h["mpeak"][0]

    (snap_discrete, snap_eps,
     snap_relax, snap_relax_hydro,
     disrupt_snap, disrupt_snap_rs) = read_convergence_limits(dir_name)
    if snap_discrete is not None:
        ok = ~h["false_selection"]
        h["conv_snap_discrete"][ok] = snap_discrete
        h["conv_snap_eps"][ok] = snap_eps
        h["conv_snap_relax"][ok] = snap_relax
        h["conv_snap_relax_hydro"][ok] = snap_relax_hydro
        h["disrupt_snap"][ok] = disrupt_snap
        h["disrupt_snap_rs"][ok] = disrupt_snap_rs

    return h

def read_convergence_limits(sim_dir):
    file_name = path.join(sim_dir, "convergence_limits.dat")
    if not path.exists(file_name): return None, None, None, None, None, None

    try:
        with open(file_name, "rb") as fp:
            n_halo, = struct.unpack("q", fp.read(8))
            snap_discrete = np.fromfile(fp, np.int32, n_halo)
            snap_eps = np.fromfile(fp, np.int32, n_halo)
            snap_relax = np.fromfile(fp, np.int32, n_halo)
            snap_relax_hydro = np.fromfile(fp, np.int32, n_halo)
            disrupt_snap = np.fromfile(fp, np.int32, n_halo)
            disrupt_snap_rs = np.fromfile(fp, np.int32, n_halo)

            return snap_discrete, snap_eps, snap_relax, snap_relax_hydro, disrupt_snap, disrupt_snap_rs
    except:
        return None, None, None, None, None, None

def read_symfind(dir_name, suffix="fid3"):
    """ read_symfind returns a pair of arrays representing each symfind
    subhalo. dir_name is the simulation's directory. The first is a 2D array
    of type SYMFIND_DTYPE where halos[i,j] is subhalo i at snapshot j.
    Symfind only tracks subhaloes, so there are no entries prioer to first
    infall. The second return value is an array of type HISTORY_DTYPE with
    one element for each subhalo. Subhaloes are ordered by decreasing
    pre-infall Mpeak, and the first element corresponds to the host halo.
    (Because the host is not a subhalo by definition, symfind does not have
    entries for this object).
    """
    h, hist = read_subhalos(dir_name)
    c = read_cores(dir_name)

    s = np.zeros(c.shape, dtype=SYMFIND_DTYPE)
    s["x"], s["v"] = c["x"], c["v"]
    s["r_tidal"], s["r_half"] = c["r_tidal"], c["r50_bound"] 
    s["m_tidal_bound"], s["m"] = c["m_tidal_bound"], c["m_bound"]
    s["m_raw"] = c["m_bound"]
    s["vmax"], s["ok"], s["interp"] = c["vmax"], c["ok"], c["interp"]
    s["f_core"], s["f_core_rs"] = c["f_core"], c["f_core"]
    s["d_core_mbp"] = c["d_core_mbp"]
    s["r_half_rs"], s["ok_rs"] = c["r50_bound_rs"], c["ok_rs"]
    s["is_err"], s["is_err_rs"] = c["is_err"], c["is_err_rs"]

    # Enforce monotonicity to reduce noise in the commonly used "m"
    for i in range(len(s)):
        s["m"][i] = monotonic_m(s["m_raw"][i], s["ok"][i])

    return s, hist

def monotonic_m(m_in, ok):
    m_out = np.zeros(len(m_in))
    idx = np.where(ok)[0]

    prev_max = 0
    for i in idx[::-1]:
        m_out[i] = max(m_in[i], prev_max)
        prev_max = m_out[i]

    return m_out

def read_cores(dir_name, include_false_selections=False, suffix="fid3"):
    """ read_cores read the particle-tracked halo cores of the halo in the
    given directory. The returned array is a structured array of type
    CORE_DTYPE with shape (n_halos, n_snaps).
    """
    if suffix == "":
        file_name = path.join(dir_name, "halos", "cores.dat")
    else:
        file_name = path.join(dir_name, "halos", "cores_%s.dat" % suffix)

    param = simulation_parameters(dir_name)
    mp = param["mp"]/param["h100"]
    h, hist = read_subhalos(dir_name, include_false_selections=True)
    with open(file_name, "rb") as fp:
        n_halo, n_snap = struct.unpack("qq", fp.read(16))
        n = n_halo*n_snap
        out = np.zeros(n, dtype=CORE_DTYPE)

        x = np.fromfile(fp, np.float32, 3*n)
        v = np.fromfile(fp, np.float32, 3*n)

        out["x"] = x.reshape((n,3))
        out["v"] = v.reshape((n,3))
        out["r_tidal"] = np.fromfile(fp, np.float32, n)
        out["r50_bound"] = np.fromfile(fp, np.float32, n)
        out["r50_bound_rs"] = np.fromfile(fp, np.float32, n)
        out["m_tidal"] = np.fromfile(fp, np.float32, n)
        out["m_tidal_bound"] = np.fromfile(fp, np.float32, n)
        out["m_bound"] = np.fromfile(fp, np.float32, n)
        out["vmax"] = np.fromfile(fp, np.float32, n)
        out["f_core"] = np.fromfile(fp, np.float32, n)
        out["f_core_rs"] = np.fromfile(fp, np.float32, n)
        out["d_core_mbp"] = np.fromfile(fp, np.float32, n)
        out["ok"] = np.fromfile(fp, np.bool, n)
        out["ok_rs"] = np.fromfile(fp, np.bool, n)
        out["is_err"] = np.fromfile(fp, np.bool, n)
        out["is_err_rs"] = np.fromfile(fp, np.bool, n)
        out["interp"] = np.fromfile(fp, np.bool, n)
       

    out = out.reshape((n_halo, n_snap))

    if not include_false_selections:
        h, hist = read_subhalos(dir_name, include_false_selections=True)
        out = out[~hist["false_selection"]]

    return out

def read_um(dir_name):
    fname = path.join(dir_name, "um", "um.dat")
    h, _ = read_subhalos(dir_name)
    n_halo, n_snap = h.shape
    n = n_snap*n_halo

    h, hist = read_subhalos(dir_name)

    out = np.zeros(h.shape[0]*h.shape[1], UM_DTYPE)

    with open(fname) as fp:
        out["m_star"] = np.fromfile(fp, np.float32, n)
        out["m_in_situ"] = np.fromfile(fp, np.float32, n)
        x = np.fromfile(fp, np.float32, 3*n)
        v = np.fromfile(fp, np.float32, 3*n)
        out["x"] = x.reshape((n,3))
        out["v"] = v.reshape((n,3))
        out["sfr"] = np.fromfile(fp, np.float32, n)
        out["rank"] = np.fromfile(fp, np.float32, n)
        out["mvir"] = np.fromfile(fp, np.float32, n)
        out["vmax"] = np.fromfile(fp, np.float32, n)
        out["is_orphan"] = np.fromfile(fp, np.bool, n)
        out["ok"] = np.fromfile(fp, np.bool, n)
        
    out = out.reshape((n_halo, n_snap))
    out["is_orphan"] = out["ok"] & (~h["ok"])
    return out

def read_branches(dir_name):
    """ read_branches reads main branch data from the halo directory dir_name.
    It returns an array with length n_branches where each element has type
    BRANCH_DTYPE.
    """
    fname = path.join(dir_name, "halos", "tree_header.dat")
    f = open(fname, "rb")
    
    n = struct.unpack("i", f.read(4))[0]
    central_idx = struct.unpack("i", f.read(4))[0]
    out = np.zeros(n, dtype=BRANCH_DTYPE)

    edges = np.fromfile(f, np.int32, n+1)
    out["start"] = edges[:-1]
    out["end"] = edges[1:]
    out["is_real"] = np.fromfile(f, bool, n)
    out["is_disappear"] = np.fromfile(f, bool, n)
    out["is_main_sub"] = np.fromfile(f, bool, n)
    out["preprocess"] = np.fromfile(f, np.int32, n)
    out["first_infall_snap"] = np.fromfile(f, np.int32, n)

    f.close()
    return out    

def read_tree(dir_name, var_names):
    """ read_tree reads variables from the halo directory halo_dir in
    depth-first order. var_names is a list of variables to be read (see
    TREE_COL_NAMES). A list of arrays is returned. Use the branches and merger
    files to identify main branches and important haloes, respectively.
    """
    tree_files = sorted(glob.glob(path.join(dir_name, "halos", "tree_*.dat")))
    tree_files = [fname for fname in tree_files if
                  path.basename(fname) != "tree_header.dat"]

    out = []
    for i in range(len(var_names)):
        var = []
        for j in range(len(tree_files)):
            hd = read_tree_header(tree_files[j])
            offset = tree_var_offset(hd, var_names[i])
            f = open(tree_files[j], "rb")
            col = tree_var_col(hd, var_names[i])
            f.seek(offset)
        
            if is_int(col, hd):
                var.append(np.fromfile(f, np.int32, hd.n))                
            elif is_float(col, hd):
                var.append(np.fromfile(f, np.float32, hd.n))
            else:
                var.append(np.fromfile(f, (np.float32, (3,)), hd.n))

        out.append(flatten(var))
    return out

def propagate_parent_indices(idx):
    """ propagate_parent_indices performs a union-find on the array idx so that
    if subhalo A was preprocessed by B which was preprocessed by C, A's index
    points to C.

    This is currently brute-force and very slow. There's a well-known union-find
    algorithm which is faster. Can switch to this if it ever matters, but it's
    also ~130 lines in my old Julia implementation, so I'm not sure that it's 
    worht it.
    """
    n_changed = -1
    while n_changed != 0:
        n_changed = 0
        for i in range(len(idx)):
            if idx[i] == i: idx[i] = -1
            if idx[i] != -1 and idx[[idx[i]]] != -1:
                n_changed += 1
                idx[i] = idx[idx[i]]

def read_merger_idxs(dir_name):
    halo_name = dir_name.split("/")[-1]
    parent_dir = "/".join(dir_name.split("/")[:-1])
    merger_idxs_name = path.join(parent_dir, "merger_idxs.txt")
    lmc_idx_cosmo, lmc_idx_zoom, gse_idx = np.loadtxt(
        merger_idxs_name, usecols=(1, 4, 7), dtype=int).T
    halo_names = np.loadtxt(merger_idxs_name, usecols=(0,), dtype=str)
    for i in range(len(halo_names)):
        if halo_name == halo_names[i]:
            return lmc_idx_cosmo[i], lmc_idx_zoom[i], gse_idx[i]
    raise ValueError("Could not find %s in %s" % (halo_name, merger_idxs_name))

# Everything else in this file is an internal helper function

def is_int(col, hd): return col < hd.n_int
def is_float(col, hd): return col >= hd.n_int and col < hd.n_float
     
def read_tree_header(fname):
    f = open(fname, "rb")
    class Header(object): pass
    hd = Header()
    hd.n, hd.n_int, hd.n_float, hd.n_vec = struct.unpack("iiii", f.read(16))
    cols = array.array("i")
    cols.fromfile(f, hd.n)
    hd.cols = np.asarray(cols, dtype=np.int32)
    return hd

def tree_var_col(hd, var_name):
    if var_name not in TREE_COL_NAMES:
        raise ValueError("Variable name '%s' unrecognized" % var_name)
    
    target_col = TREE_COL_NAMES[var_name]
    for i in range(len(hd.cols)):
        if hd.cols[i] == target_col: break

    return i

def tree_var_offset(hd, var_name):
    i = tree_var_col(hd, var_name)
    hd_size = 4*4 + 4*(hd.n_int + hd.n_float + hd.n_vec)
    return hd.n*4*(i + 2*max(0, i - hd.n_int - hd.n_float)) + hd_size    

def merger_lookup_table(b, dfid):
    """ merger_lookup_table creates a look-up table for finding the branches of
    mergers. It's used along side find_merger_branch and
    find_all_merger_branches.

    b is an array of branches (BRANCH_DTYPE) and dfid is an array of
    depth-first IDs with merger tree ordering (can be read by calling 
    read_tree on "dfid")
    """
    return dfid[b["start"]]

def find_merger_branch(lookup_table, co_prog):
    """ find_merger_branch returns the index of the branch corresponding to
    a given co-progenitor ID.
    """
    return np.searchsorted(lookup_table, co_prog)

def find_all_merger_branches(b, lookup_table, co_prog, i):
    """ find_all_merger_branches returns an array with the indexes of all the
    branches that merge with a halo in the given snapshot (i.e. this is the
    last resolved snapshot of those branches). b is an array of branches,
    (BRANCH_DTYPE), lookup_table is a table returned by merger_lookup_table,
    co_prog is an array of co-progenitor IDs (can be read by calling read_tree
    on "next_co_prog"), and i is the index of the halo whose mergers you're
    trying to find within the merger tree.
    """
    branches = []
    while co_prog[i] != -1:
        bi = find_merger_branch(lookup_table, co_prog[i])
        branches.append(bi)
        i = b["start"][bi]

    return np.array(branches, dtype=int)
        
def read_particle_header(base_dir):
    n_snap = struct.unpack("i", f.read(4))[0]
    n_merger = struct.unpack("i", f.read(4))[0]

    idx = np.fromfile(f, np.int32, n_merger)    
    
class ParticleHeader(object):
    def __init__(self, base_dir):
        file_name = path.join(base_dir, "particles", "particle_header.dat")
        try:
            f = open(file_name, "rb")
        except FileNotFoundError:
            raise FileNotFoundError("The particle data for this halo does not exist.") from None
        # TODO: if not working, try os.path.exists() and os.path.isfile() after that

        self.n_file = struct.unpack("i", f.read(4))[0]
        self.n_halo = struct.unpack("i", f.read(4))[0]
        self.n_particle = struct.unpack("i", f.read(4))[0]
        
        self.file_lengths = np.fromfile(f, np.int32, self.n_file)
        self.offsets = np.fromfile(f, np.int32, self.n_halo)
        self.sizes = np.fromfile(f, np.int32, self.n_halo)
        self.file_idxs = np.fromfile(f, np.int32, self.n_halo)
        self.n0 = np.fromfile(f, np.int32, self.n_halo)

        self.base_dir = base_dir
        
        f.close()

class ParticleTags(object):
    def __init__(self, base_dir, part_hd):
        self.id = [None]*part_hd.n_halo
        self.snap = [None]*part_hd.n_halo
        self.flag = [None]*part_hd.n_halo

        for i_file in range(part_hd.n_file):
            file_name = path.join(base_dir, "particles",
                                  "tags.%d.dat" % i_file)
            f = open(file_name, "rb")
            
            to_read = np.where(i_file == part_hd.file_idxs)[0]
            for i in to_read:
                self.id[i] = np.fromfile(f, np.int32, part_hd.sizes[i])
                self.snap[i] = np.fromfile(f, np.int16, part_hd.sizes[i])
                self.flag[i] = np.fromfile(f, np.uint8, part_hd.sizes[i])

            f.close()

        self.n0 = part_hd.n0
        self.n_file = part_hd.n_file
        self.n_halo =  part_hd.n_halo
        self.file_idxs = part_hd.file_idxs

        self._n0_snaps = [self.snap[i][self.flag[i] == 0]
                          for i in range(part_hd.n_halo)]
        self._offset_lookup = {}

    def halo_offset(self, snap, i_halo):
        if snap in self._offset_lookup:
            return self._offset_lookup[snap][i_halo]
        else:
            offsets = np.zeros(self.n_halo, dtype=int)
            for i_file in range(self.n_file):
                to_read = np.where(i_file == self.file_idxs)[0]
                curr_offset = 4
                for i in to_read:
                    offsets[i] = curr_offset

                    if len(self._n0_snaps[i]) == 0: continue
                    n = np.sum(self._n0_snaps[i] <= snap)
                    if n == 0: continue

                    curr_offset += 6*n + 8 + 12 + 12
            self._offset_lookup[snap] = offsets
            return self._offset_lookup[snap][i_halo]

class TagLookup(object):
    def __init__(self, base_dir):
        file_name = path.join(base_dir, "particles", "tags.lookup_table.dat")
        f = open(file_name, "rb")
        Np = struct.unpack("i", f.read(4))[0]
        self.halo = np.fromfile(f, np.int16, Np)
        self.index = np.fromfile(f, np.int32, Np)
        
class ParticleInfo(object):
    def __init__(self, base_dir):
        self.part_hd = ParticleHeader(base_dir)
        self.tags = ParticleTags(base_dir, self.part_hd)
        self.lookup = TagLookup(base_dir)

        global_offset = np.zeros(self.part_hd.n_halo, dtype=int)
        global_offset[1:] = np.cumsum(self.part_hd.n0[:-1])
        self.global_offset = global_offset
        self.global_index = global_offset[self.lookup.halo] + self.lookup.index
        
        # Used to find false selections
        self.h_false, self.hist_false = read_subhalos(
            base_dir, include_false_selections=True)

def is_real_confirmed(part_info, h, i_sub):
    first_snap = np.where(h[i_sub]["ok"])[0][0]
    ok = read_particles(part_info, None, first_snap, "valid", owner=i_sub)
    return np.sum(ok) > 0

class Particles(object):
    """ A wrapper for accessing particle data associated with the desired halo suite."""
    def __init__(self, sim_dir):
        self.sim_dir = sim_dir
        self.part_info = ParticleInfo(sim_dir)
        self.params = simulation_parameters(sim_dir)
        self.scale = scale_factors(sim_dir)
        self.h_cmov, _ = read_subhalos(sim_dir, comoving=True)

    def read(self, snap, halo=-1, mode="current", comoving=False):
        """ read returns a list of 1D arrays containing frequently used variables 
        stored in the particle data for a given subhalo. 

        Parameters
        ----------
        snap: int, snapshot.
        halo: int, index for subhalo of interest, default returns all subhalos.
        comoving: boolean, default is False to return physical units.
        mode: string 
            "Current" is the default setting and returns variables for valid particles.
            "All" returns variables for all particles, valid or invalid.
            "Smooth" returns variables for smoothly accreted particles. 

        Returns
        -------
        p: list containing 1D arrays of the following variables from stored particle data:
            id, x, v, snap, ok, and smooth. See PARTICLE_DTYPE for more detail.

        """

        sim_dir = self.sim_dir
        part_info = self.part_info
        h0 = self.h_cmov[0,snap]
        a = self.scale[snap]

        # These modes need to read in all particles simultaneously: cannot
        # just read one halo.
        if mode == "all" or mode == "current":
            idp = read_particles(part_info, sim_dir, snap, "id")
            x = read_particles(part_info, sim_dir, snap, "x")
            v = read_particles(part_info, sim_dir, snap, "v")
            if not comoving:
                for s in range(len(x)):
                    x[s] = util.set_units_x(x[s], h0, a, self.params)
                    v[s] = util.set_units_v(v[s], h0, a, self.params)
            
            snaps = read_particles(part_info, sim_dir, snap, "snap")
            valid = read_particles(part_info, sim_dir, snap, "valid")
            smooth = read_particles(part_info, sim_dir, snap, "ownership")
            
            # Remove invalid particles.
            if mode == "current":
                for i in range(len(x)):
                    ok = valid[i]
                    x[i], v[i], idp[i] = x[i][ok], v[i][ok], idp[i][ok]
                    snaps[i], smooth[i] = snaps[i][ok], smooth[i][ok]
                    valid[i] = valid[i][ok]

        if mode == "smooth":
            idp = [None]*len(self.h_cmov)
            x = [None]*len(self.h_cmov)
            v = [None]*len(self.h_cmov)
            snaps = [None]*len(self.h_cmov)
            valid = [None]*len(self.h_cmov)
            smooth = [None]*len(self.h_cmov)
            for sh in range(len(self.h_cmov)):
                # A single halo read is fast for smooth particles, so we only
                # have to read the specified halo in this mode
                if halo != -1 and sh != halo: continue
                idp[sh] = read_particles(part_info, sim_dir, snap, "id", owner=sh)
                x[sh] = read_particles(part_info, sim_dir, snap, "x", owner=sh)
                v[sh] = read_particles(part_info, sim_dir, snap, "v", owner=sh)
                if not comoving:
                    x[sh] = util.set_units_x(x[sh], h0, a, self.params)
                    v[sh] = util.set_units_v(v[sh], h0, a, self.params)
            
                snaps[sh] = read_particles(part_info, sim_dir, snap,
                                           "snap", owner=sh)
                valid[sh] = read_particles(part_info, sim_dir, snap,
                                           "valid", owner=sh)
                smooth[sh] = read_particles(part_info, sim_dir, snap,
                                            "ownership", owner=sh)
        p = [None]*len(x)
        for i in range(len(x)):
            if halo == -1 or i == halo:
                p[i] = np.zeros(len(x[i]), dtype=PARTICLE_DTYPE)
                p[i]["id"] = idp[i]
                p[i]["x"] = x[i]
                p[i]["v"] = v[i]
                p[i]["snap"] = snaps[i]
                p[i]["ok"] = valid[i]
                p[i]["smooth"] = smooth[i] == 0

        if halo == -1:
            return p 
        else:
            return p[halo]

    def core_particles(self, snap, halo=-1, comoving=False):
        """ TODO: write docstring
        """
        pass

def read_particles(part_info, base_dir, snap, var_name,
                   owner=None, include_false_selections=False):
    hd, tags = part_info.part_hd, part_info.tags
    
    h_idx = np.arange(len(part_info.hist_false), dtype=int)
    if include_false_selections:
        false_remapping = h_idx
    else:
        false_remapping = h_idx[~part_info.hist_false["false_selection"]]
    if owner is not None and not include_false_selections:
        owner = false_remapping[owner]

    if var_name == "id":
        if owner is None:
            return [tags.id[i] for i in false_remapping]
        else:
            return tags.id[owner][tags.flag[owner] == 0]
    elif var_name == "snap":
        if owner is None:
            return [tags.snap[i] for i in false_remapping]
        else:
            return tags.snap[owner][tags.flag[owner] == 0]
    elif var_name == "ownership":
        if owner is None:
            return [tags.flag[i] for i in false_remapping]
        else:
            return tags.flag[owner][tags.flag[owner] == 0]
    elif var_name == "valid":
        if owner is None:
            valid = [None]*len(tags.snap)
            for i in range(len(valid)):
                valid[i] = tags.snap[i] <= snap
            return [valid[i] for i in false_remapping]
        else:
            valid = (tags.snap[owner] <= snap)[tags.flag[owner] == 0]
            return valid
    elif var_name in ["x", "v"]:
        snap_name = "snap_%03d" % snap
        if owner is not None:
            first_snap = np.min(np.where(part_info.h_false["ok"][owner,:])[0])
            if snap < first_snap:
                num = int(np.sum(tags.flag[owner] == 0))
                return np.ones((num, 3)) * np.nan

            i_file = hd.file_idxs[owner]
            file_name = "%s.%d.dat" % (var_name, i_file)

            path = os.path.join(hd.base_dir, "particles", snap_name, file_name)
            f = open(path, "rb")

            code = struct.unpack("i", f.read(4))[0]
            if code != 0: raise ValueError("%s is not a vector file." % path)

            offset = tags.halo_offset(snap, owner)
            f.seek(offset, 0)
            
            size = struct.unpack("q", f.read(8))[0]
            min = np.array(struct.unpack("fff", f.read(12)))
            max = np.array(struct.unpack("fff", f.read(12)))
            
            qx_i = np.fromfile(f, np.uint16, size*3)
            x_i = _dequantize_vector(qx_i, min, max)
            x_i = _expand_vector(tags, x_i, owner, snap)

            f.close()
            return x_i

        x_full = np.zeros((hd.n_particle, 3))

        for i_file in range(hd.n_file):
            file_name = "%s.%d.dat" % (var_name, i_file)
            path = os.path.join(hd.base_dir, "particles", snap_name, file_name)
            f = open(path, "rb")
            
            code = struct.unpack("i", f.read(4))[0]
            if code != 0: raise ValueError("%s is not a vector file." % path)

            to_read = np.where(hd.file_idxs == i_file)[0]
            for i_halo in to_read:
                if hd.n0[i_halo] == 0: continue
                ok = tags.snap[i_halo][tags.flag[i_halo] == 0] <= snap
                if np.sum(ok) == 0: continue

                curr_offset = f.seek(0, 1)

                size = struct.unpack("q", f.read(8))[0]
                min = np.array(struct.unpack("fff", f.read(12)))
                max = np.array(struct.unpack("fff", f.read(12)))
                
                qx_i = np.fromfile(f, np.uint16, size*3)
                x_i = _dequantize_vector(qx_i, min, max)

                x_i = _expand_vector(tags, x_i, i_halo, snap)
                offset = part_info.global_offset
                x_full[offset[i_halo]: offset[i_halo]+len(x_i)] = x_i

            f.close()
            
        out = [None]*hd.n_halo
        for i_halo in range(hd.n_halo):
            idx = part_info.global_index[tags.id[i_halo] - 1]
            out[i_halo] = x_full[idx]

        return [out[i] for i in false_remapping]
    elif var_name == "infall_core":
        # Setting owner doesn't do anything for infall_core.

        file_name = os.path.join(base_dir, "halos", "infall_cores.dat")
    
        with open(file_name, "rb") as fp:
            n_halo, n_core = struct.unpack("qq", fp.read(16))
            idxs = np.fromfile(fp, dtype=np.int32, count=n_halo*n_core)
            idxs = idxs.reshape((n_halo, n_core))
        return [idxs[i] for i in false_remapping]
    else:
        raise ValueError("Unknown property name, '%s'" % var_name)    
    
def _dequantize_vector(qx, min, max):
    dx = max - min
    n = len(qx)//3
    
    qx_vec = qx.reshape((n, 3))
    uint16_max = float(np.iinfo(np.uint16).max)
    out = np.zeros((n, 3))
    for dim in range(3):
        out[:,dim] = dx[dim]/uint16_max*(qx_vec[:,dim]+ random.random(n)) + min[dim]

    return out
         
def _expand_vector(tags, x, i, snap):
    ok = tags.snap[i][tags.flag[i] == 0] <= snap
    out = np.ones((len(ok), 3)) * np.nan
    if np.sum(ok) != len(x):
        print("   ", np.sum(ok), len(x))
    assert(np.sum(ok) == len(x))
    out[ok] = x
    return out
    
# NFW math
def _two_way_interpolate(x, y):
    """ two_way_interpolate returns interpolation functions that map from x -> y
    and y -> x. Both input arrays must monotonic.
    """
    if x[0] < x[1]:
        x_to_y = interpolate.interp1d(x, y)
    else:
        x_to_y = interpolate.interp1d(x[::-1], y[::-1])

    if y[0] < y[1]:
        y_to_x = interpolate.interp1d(y, x)
    else:
        y_to_x = interpolate.interp1d(y[::-1], x[::-1])

    return x_to_y, y_to_x

def _f(x): return np.log(1+x) - x/(1+x)
def _x_max_nfw(): return 2.1626
def _v_vmax_nfw(x): return 2.1506 * np.sqrt(_f(x) / x)
def _m_enc_nfw(x): return _f(x)
def _alpha_nfw(x): return -1 - 2*x/(1 + x)

def _half_mass_nfw(x, mass_fraction):
    def f_mass_ratio(xx):
        return _m_enc_nfw(xx) / _m_enc_nfw(x) - mass_fraction
    sol = optimize.root_scalar(f_mass_ratio, bracket=[1e-4*x, x])
    return sol.root

_c = 10**np.linspace(0, 3, 1000)
_cv = np.sqrt(_f(_x_max_nfw()) / _x_max_nfw() * _c/_f(_c))
# Make sure we're solving the correct part of the cV - cvir relation
_cv[_c < _x_max_nfw()] = 1.0 
c_to_cv_nfw, cv_to_c_nfw = _two_way_interpolate(_c, _cv)

def vmax_to_cvir_nfw(vmax, mvir, rvir):
    vvir = 655.8 * (mvir/1e14)**0.5 * (rvir/1.0)**-0.5
    cv = vmax/vvir
    cv[cv < 1] = 1
    return cv_to_c_nfw(cv)    

def main():
    base_dir = "/oak/stanford/orgs/kipac/users/phil1/simulations/ZoomIns"
    suite = "SymphonyMilkyWay"
    i_halo = 0
    sim_dir = get_host_directory(base_dir, suite, i_halo)

    param = simulation_parameters(suite)
    h, hist = read_subhalos(param, sim_dir)
    
    info = ParticleInfo(sim_dir)

    x = read_particles(info, sim_dir, 235, "x")
    v = read_particles(info, sim_dir, 235, "v")
    ok = read_particles(info, sim_dir, 235, "valid")
            
if __name__ == "__main__": main()
