import numpy as np
import matplotlib.pyplot as plt
import numpy.random as random
import scipy.stats as stats
import abc
import scipy.optimize as optimize
import symlib.lib as lib

# cosmology information
from colossus.cosmology import cosmology
from colossus.halo import mass_so

""" 
Variables names:

'mpeak' - Peak Brayn & Norman virial mass (Msun)
'mvir' - Bryan & Norman virial mass (Msun)
'rvir' - Bryan & Norman virial radius (pkpc)
'cvir' - NFW concentration relative to the Bryan & Norman virial radius
'z' redshift
"""

NIL_RANK = -1

DEFAULT_QUANTILE_EDGES = np.zeros(13)
DEFAULT_QUANTILE_EDGES[1:] = 10**np.linspace(-4, 0, 12)

DEFAULT_R_BINS = np.zeros(102)
DEFAULT_R_BINS[1:] = 10**np.linspace(-2.5, 0, 101)

MIN_PARTICLES_PER_QUANTILE = 10

DEFAULT_CORE_PARTICLES = 30

#########################
# Abstract base classes #
#########################

class ProfileModel(abc.ABC):
    """ AbstractProfile is an abstract base class for galaxy profile models. It
    allows you to convert a half-light radius into an enclosed stellar mass
    profile.
    """
    
    @abc.abstractmethod
    def m_enc(self, m_star, r_half, r):
        """ m_enc returns the enclosed mass profile as a function of 3D radius,
        r. m_star is the asymptotic stellar mass of the galaxy, r_half is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units as m_star.
        """
        pass

    @abc.abstractmethod
    def density(self, m_star, r_half, r):
        """ density returns the local density as a function of 3D radius, r.
        m_star is the asymptotic stellas mass of the galaxy, r_hald is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units of m_star.
        """
        pass

class RHalfModel(abc.ABC):
    """ RHalfModel is an abstract base class for models of galaxy half-mass
    radii. 
    """
    
    @abc.abstractmethod
    def r_half(self, **kwargs):
        """ r_half returns the half-mass radius of a a galaxy in physical kpc.
        """
        pass

    @abc.abstractmethod
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        pass

    def trim_kwargs(self, kwargs):
        out = {}
        for key in self.var_names():
            out[key] = kwargs[key]
        return out
    
class MStarModel(abc.ABC):
    """ MStarModel is an abstract base class for models of the Mhalo-Mstar
    relation.
    """

    @abc.abstractmethod
    def m_star(self, **kwargs):
        """ m_star returns the stellar mass of a galaxy in Msun.
        """
        pass

    @abc.abstractmethod
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        pass

    def trim_kwargs(self, kwargs):
        out = {}
        for key in self.var_names():
            out[key] = kwargs[key]
        return out
    
class AbstractRanking(abc.ABC):
    """ AbstractRanking is an abstract base class for various models that 
    rank particles.
    """

    def __init__(self, ranks, core_idx):
        """ ranks is an array of ranks for every particle in the halo. If that
        particle has been accreted by the time of tagging, it should have an
        non-negative integer rank. If the particle hasn't been accreted yet,
        it should be set to NIL_RANK

        core_idx is the index of the core particles of the halo according to
        this ranking scheme.

        quantile_edges are the upper and lower quantile bounds for the each rank
        (i.e. the quantile, q, of a particle in rank i should be in the range
        quantile_edges[i] < q < quantile_edges[i+1]).
        """
        self.ranks = ranks
        self.n_max = np.max(ranks)
        self.core_idx = core_idx
        self.mp_star = np.zeros(len(ranks))
        
        self.idx = None
        self.x = None
        self.v = None
        self.xc = None
        self.vc = None
        
    def load_particles(self, x, v, idx):
        """ load_particles loads paritcle properties into the ranking. Must
        be called before mp_star() or either of the core functions. v may be
        None unless you call a later method which works with velocities (i.e.
        core_v)
        """
        self.x = np.copy(x)
        if v is not None:
            self.v = np.copy(v)
        else:
            self.v = None
        self.idx = idx

        self.xc = self.core_x()
        for dim in range(3): self.x[:,dim] -= self.xc[dim]
        
        if self.v is not None:
            self.vc = self.core_v()
            for dim in range(3): self.v[:,dim] -= self.vc[dim]

    def core_x(self):
        """ core_x returns the position of the halo core. You must have called
        load_particles() first.
        """
        if self.x is None: raise ValueError("x not set")
        core = self.x[np.searchsorted(self.idx, self.core_idx)]
        return np.median(core, axis=0)
        
    def core_v(self):
        """ core_v returns the velocity of the halo core. You must have called
        load_particles() first with a non-None argument for v.
        """
        if self.v is None: raise ValueError("v not set")
        
        core = self.v[np.searchsorted(self.idx, self.core_idx)]
        return np.median(core, axis=0)
        
    def set_mp_star(self, rvir, profile_model, r_half, m_star,
                r_bins=DEFAULT_R_BINS):
        r = np.sqrt(np.sum(self.x**2, axis=1))/rvir
        M = ranked_np_profile_matrix(self.ranks, self.idx, r, r_bins)
        n = M.shape[1]

        m_star_enc_target = profile_model.m_enc(
            m_star, r_half, r_bins[1:]*rvir)
        
        dm_star_enc_target = np.zeros(len(m_star_enc_target))
        dm_star_enc_target[1:] = m_star_enc_target[1:] - m_star_enc_target[:-1]
        dm_star_enc_target[0] = m_star_enc_target[0]
        
        res = optimize.lsq_linear(
            M, dm_star_enc_target, bounds=(np.zeros(n), np.inf*np.ones(n))
        )
        self.mp_star_table = res.x

        self.mp_star[self.idx] = self.mp_star_table[self.ranks[self.idx]]
        is_nil = self.ranks[self.idx] == NIL_RANK
        self.mp_star[self.idx[is_nil]] = 0

        correction_frac = m_star/np.sum(self.mp_star)
        self.mp_star *= correction_frac
        self.mp_star_table *= correction_frac
        
        return self.mp_star

    def ranked_halfmass_radius(self):
        """ ranked_halfmass_radii returns an array of the current half-mass
        radii of each rank bin in pkpc.
        """
        r = np.sqrt(np.sum(self.x**2, axis=1))

        rhalf = np.zeros(self.n_max+1)
        for i in range(self.n_max+1):
            rhalf[i] = np.median(r[self.ranks[self.idx] == i])
            
        return rhalf

    def relaxation_time(self, mp, eps):
        """ relaxation_time returns the relaxation time in Gyr of each particle
        in the halo, assuming that it is on a circular orbit. mp and eps are in
        Msun and pkpc, respectively.
        """
        
        r = np.sqrt(np.sum(self.x**2, axis=1))
        
        order = np.argsort(r)
        Nr = np.zeros(len(r))
        Nr_sort = np.arange(len(Nr)) + 1
        Nr[order] = Nr_sort

        t_orb = t_orbit(Nr*mp, r)
        t_relax = t_relax_t_orbit(Nr, r, eps)*t_orb

        return t_relax
    
    def ranked_relaxation_time(self, mp, eps):
        """ ranked_relaxation_times returns an array of the relaxation times of
        each rank bin in Gyr. Ranks begin expanding when 0.18*t_relax has
        passed.
        """

        t_relax = self.relaxation_time(mp, eps)
        mean_t_relax = np.zeros(self.n_max+1)
        for i in range(len(mean_t_relax)):
            mean_t_relax[i] = np.mean(t_relax[self.ranks[self.idx] == i])
        
        return mean_t_relax
    
##################################
# Specific Model Implementations #
##################################

    
class PlummerProfile(ProfileModel):
    """ PlummerProfile models a galaxy's mass distribution as a Plummer sphere.
    """
    def m_enc(self, m_star, r_half, r):
        """ m_enc returns the enclosed mass profile as a function of 3D radius,
        r. m_star is the asymptotic stellar mass of the galaxy, r_half is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units as m_star.
        """
        a = r_half
        return m_star*r**3 /(r**2 + a**2)**1.5
    
    def density(self, m_star, r_half, r):
        """ density returns the local density as a function of 3D radius, r.
        m_star is the asymptotic stellas mass of the galaxy, r_hald is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units of m_star.
        """
        a = r_half
        return 3*m_star/(4*np.pi*a**3) * (1 + r**2/a**2)**(-5/2)
    
class Nadler2020RHalf(RHalfModel):
    """ Nadler20202RHalf models galaxies according to the z=0 size-mass relation
    in Nadler et al. 2020. (https://arxiv.org/abs/1912.03303)

    This is calibrated by fitting a galaxy-halo model to Milky Way satellites.
    """
    
    def __init__(self, A=27e-6, n=1.07, R0=10e-3, sigma_log_R=0.63):
        """ The constructor for Nadler2020RHalf allows you to change the
        parameters of the fit. Please consult the paper for discussion on what
        these parameters mean. This allows you to sample the posterior
        distribution for these parameters.
        """
        self.A = A
        self.n = n
        self.R0 = R0
        self.sigma_log_R = sigma_log_R

    def r_half(self, **kwargs):
        """ r_half returns the half-mass radius of a a galaxy in physical kpc.
        Required keyword arguments:
         - rvir
        """
        # inputs and outputs are in pMpc (no h).
        log_R = np.log10(self.A * (rvir/self.R0)**self.n)
        log_scatter = self.sigma_log_R*random.normal(0, 1, size=np.shape(rvir))
        return 10**(log_R + log_scatter)

    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["rvir"]

class Jiang2019RHalf(RHalfModel):
    """ Jiang2019Rhalf models galaxies according to the z-dependent size-mass
    relation in Jiang et al. 2019 (https://arxiv.org/abs/1804.07306; Appendix
    D). This model is fit to   galaxies in NIHAO/VELA and has been rescaled
    to  match the noramlization of the z-dependence seen in GAMA+CANDLES when
    combined with abundance matching (Somerville et al. 2018;
    https://arxiv.org/abs/1701.03526; M* > 1e8 Msun).
    """
    def __init__(self, sigma_log_R=0.2):
        """ The constructor for Jiang2019RHalf allows you to change the
        intrinsic scatter in the relation. I eyeballed sigma_log_R at 0.2 from
        their plots, so this value in particular is worth modifying.
        """
        # I eyeballed the 0.2 from their plots.
        self.sigma_log_R = sigma_log_R
    
    def r_half(self, rvir=None, cvir=None, z=None, no_scatter=False):
        """ r_half returns the half-mass radius of a a galaxy in physical kpc.
        Required keyword arguments:
         - rvir
         - cvir
         - z
        """
        if rvir is None: raise ValueError("rvir not supplied")
        if cvir is None: raise ValueError("cvir not supplied")
        if z is None: raise ValueError("z not supplied")
        
        # From Appendix D
        fz = 0.02*(1 + z)**-0.2
        R = fz * (cvir/10)**-0.7 * rvir
        log_scatter = self.sigma_log_R*random.normal(0, 1, size=np.shape(rvir))
        if not no_scatter:
            return 10**(np.log10(R) + log_scatter)
        else:
            return R

    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["rvir", "cvir", "z"]
    
class UniverseMachineMStar(MStarModel):
    def m_star(self, mpeak=None, z=None, no_scatter=False):
        """
         Required keyword arguments:
         - rvir
         - cvir
         - z
        """

        if mpeak is None: raise ValueError("mpeak not supplied")
        if z is None: raise ValueError("z not supplied")
        
        mpeak = mpeak
        
        a = 1/(1 + z)

        e0 = -1.435
        al_lna = -1.732

        ea = 1.831
        alz = 0.178

        e_lna = 1.368
        b0 = 0.482

        ez = -0.217
        ba = -0.841

        m0 = 12.035
        bz = -0.471

        ma = 4.556
        d0 = 0.411

        m_lna = 4.417
        g0 = -1.034

        mz = -0.731
        ga = -3.100

        al0 = 1.963
        gz = -1.055

        ala = -2.316
        
        log10_M1_Msun = m0 + ma*(a-1) - m_lna*np.log(a) + mz*z
        e = e0 + ea*(a - 1) - e_lna*np.log(a) + ez*z
        al = al0 + ala*(a - a) - al_lna*np.log(a) + alz*z
        b = b0 + ba*(a - 1) + bz*z
        d = d0
        g = 10**(g0 + ga*(a - a) + gz*z)

        x = np.log10(mpeak/10**log10_M1_Msun)

        al = 2
        
        log10_Ms_M1 = (e - np.log10(10**(-al*x) + 10**(-b*x)) +
                       g*np.exp(-0.5*(x/d)**2))
                       
        log10_Ms_Msun = log10_Ms_M1 + log10_M1_Msun

        if not no_scatter:
            log_scatter = 0.2*random.normal(0, 1, size=np.shape(mpeak))
            log10_Ms_Msun += log_scatter
        
        Ms = 10**log10_Ms_Msun
        
        return Ms
    
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["mpeak", "z"]
    

class RadialEnergyRanking(AbstractRanking):
    def __init__(self, params, x, v, idx, n_max,
                 core_particles=DEFAULT_CORE_PARTICLES,
                 quantile_edges=DEFAULT_QUANTILE_EDGES):
        """ Takes mp (Msun; may be a scalar), x (pkpc), v (pkm/s), idx (the
        index into the full particle array), and n_max (the number of particles
        in the full particle arrays).
        
        only_tag may be set so that only a subset of particles may be tagged
        (e.g. early-accreted particles).
        """        
        self.rmax, self.vmax, self.pe_vmax2, self.order = profile_info(params,x)

        self.ke_vmax2 = 0.5*np.sum(v**2, axis=1) / self.vmax**2
        self.E_vmax2 = self.ke_vmax2 + self.pe_vmax2
        
        core_q = min(core_particles / len(self.E_vmax2), 1)
        core_E_vmax2 = np.quantile(self.E_vmax2, core_q)
        
        core_idx = idx[self.E_vmax2 <= core_E_vmax2]

        E_edges = np.linspace(-10, 0, 41)
        quantile_edges = [np.sum(self.E_vmax2< E_edges[i])/
                          len(self.E_vmax2) for i in range(len(E_edges))]
        
        ranks, self.E_edges = rank_by_quantile(
            quantile_edges, self.E_vmax2, idx, n_max)
                
        super(RadialEnergyRanking, self).__init__(ranks, core_idx)
        
#################################
# General classes and functions #
#################################

def tag_stars(base_dir, params, galaxy_halo_model, mergers, halo_idx, tag_snap,
              E_snap=None):
    """ tag_stars returns an array of star particle masses and an object
    that implements AbstractRanking (for core-tracking purposes) of the
    specified subhalo.

    base_dir is the base directory of the simulation output, galaxy_halo_model
    is a GalaxyHaloModel instance, mergers is the second result of
    lib.read_mergers(), halo_idx is the index of the halo within mergers, and
    tag_snap is the snapshot where star particles will be tagged. If you'd like
    to force default_tag to measure energies at a specific snapshot, you can do
    that with E_snap.
    """
    halo = mergers[halo_idx]
    scale = lib.scale_factors(base_dir)
    
    star_tag_snap = tag_snap
    if E_snap is None:
        E_tag_snap = look_back_orbital_time(
            params, scale, star_tag_snap, 0.125, halo, 0.5)
    else:
        E_tag_snap = E_snap
        
    x_star = lib.read_particles(base_dir, star_tag_snap, halo_idx, ["x"])
    x_E, v_E = lib.read_particles(base_dir, E_tag_snap, halo_idx, ["x", "v"])
    n_max = len(x_E)

    x_E, v_E, idx_E = clean_particles(
        params, x_E, v_E, mergers[halo_idx,E_tag_snap],
        scale[E_tag_snap]
    )
    x_star, _, idx_star = clean_particles(
        params, x_star, None, mergers[halo_idx,star_tag_snap],
        scale[star_tag_snap]
    )
    
    ranks = RadialEnergyRanking(params, x_E, v_E, idx_E, n_max)
    ranks.load_particles(x_star, None, idx_star)

    kwargs = galaxy_halo_model.get_kwargs(
        params, scale, mergers[halo_idx], star_tag_snap)
    mp_star = galaxy_halo_model.set_mp_star(ranks, kwargs)
    
    return mp_star, ranks
    
def look_back_orbital_time(params, scale, snap, dt_orbit, halo, min_mass_frac):
    """ look_back_orbital_time returns the snapshot of the time which has
    experienced dt_orbit orbital times before the given snapshot, snap. To
    protect against looking back too far, you may specify a maximum amount that
    the halo can grow between the look-back snapshot. halo is an object from
    the subhalos.dat file.
    """
    cosmo = cosmology.setCosmology("", lib.colossus_parameters(params))
    
    z = 1/scale - 1
    age = cosmo.age(z)
    T_orbit = mass_so.dynamicalTime(z, "vir", "orbit")
    
    dT = (age - age[snap])
    dT_orbit = dT/T_orbit

    orbit_start = np.searchsorted(dT_orbit, -dt_orbit)

    for snap_start in range(snap, orbit_start, -1):
        # This also catches cases where mvir = -1
        if halo["mvir"][snap_start-1] < min_mass_frac*halo["mvir"][snap]:
            break

    return snap_start
            

class GalaxyHaloModel(object):
    def __init__(self, m_star_model, r_half_model, profile_model,
                 no_scatter=False):
        """ GalaxyHaloModel requires a model for the M*-Mhalo relation,
        m_star_model, a model for how the projected half-mass radius and Mhalo
        are related, and a model for the halo profile, profile_model. These
        should be types that inherit from AbstractMstarModel,
        AbstractRHalfModel, and AbstractProfileModel, respectively.

        If you'd like to remove scatter from your model, set no_scatter=True.
        """
        self.m_star_model = m_star_model
        self.r_half_model = r_half_model
        self.profile_model = profile_model
        self.no_scatter = no_scatter
        
    def set_mp_star(self, ranks, kwargs, r_half=None, m_star=None):
        """ set_mp_star sets the stellar masses of a halo's dark matter
        particles given their positions relative to the halo center, and
        ranking, ranks (type: inherits from AbstractParticleRanking). This
        function accepts the same keyword arguments as its m_star_model and
        r_half_model arguments. You may also fix the half-mass radius and
        stellar mass to whatever you want with r_half (units: pkpc) and m_star
        (units: Msun).
        """
        if m_star is None:
            check_var_names(kwargs, self.m_star_model)
            m_star = self.m_star_model.m_star(
                no_scatter=self.no_scatter,
                **self.m_star_model.trim_kwargs(kwargs))
        if r_half is None:
            check_var_names(kwargs, self.r_half_model)
            r_half = self.r_half_model.r_half(
                no_scatter=self.no_scatter,
                **self.r_half_model.trim_kwargs(kwargs))

        return ranks.set_mp_star(
            kwargs["rvir"], self.profile_model, r_half, m_star)
    
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        r_half_names = self.r_half_model.var_names()
        m_star_names = self.m_star_model.var_names()
        
        out_dict = { }
        for i in range(len(r_half_names)):
            out_dict[r_half_names[i]] = None
        for i in range(len(m_star_names)):
            out_dict[m_star_names[i]] = None

        if "rvir" not in out_dict: out_dict["rvir"] = None
            
        return sorted(list(out_dict.keys()))

    def get_kwargs(self, params, scale, halo, snap):
        """ get_kwargs creates the kwargs that that are needed for a call to
        set_mp_star.
        """
        h100 = params["h100"]
        
        z = 1/scale - 1
        
        var_names = self.var_names()

        kwargs = { }
        
        for i in range(len(var_names)):
            if var_names[i] == "z":
                kwargs["z"] = z
                continue

            if var_names[i] == "mpeak":
                x = np.max(halo["mvir"])
            else:
                x = halo[snap][var_names[i]]
            
            kwargs[var_names[i]] = x

        return kwargs
    
def check_var_names(kwargs, model):
    """ check_var_names checks whether kwargs contains all the arguments that
    model requires.
    """
    model_names = model.var_names()
    for i in range(len(model_names)):
        if model_names[i] not in kwargs:
            raise ValueError(("The variable '%s' is required, but was not " +
                              "given as an argument to set_mp_star.") %
                             model_names[i])
    
        
def clean_particles(params, x, v, s, scale, ok=None):
    """ clean_particles performs various normalizing and centering operations
    on particle positions and velocities. x and v are the positions and
    velocities of the particles in Gadget code units (cMpc/h, ckm/s). h is a
    halo from a merger.dat file. h100 is H0 / (100 km/s/Mpc) and scale is the
    scale factor. vp can be None and will be ignored if so.

    Returns x_clean, v_clean, idx. Particles that have not yet been accreted
    or which are False in the ok array will be removed. x_clean are centered
    on the halo and in units of pkpc. v_clean are also centered on the halo and
    are in units of pkm/s. idx is the indices of x_clean and v_clean in the
    original x and v arrays.
    """

    h100 = params["h100"]
    
    if ok is None:
        ok = x[:,0] > 0
    else:
        ok = ok & (x[:,0] > 0)

    x = np.copy(x)
    if v is not None: v = np.copy(v)
        
    for dim in range(3):
        if v is not None:
            v[:,dim] *= np.sqrt(scale)
            v[:,dim] -= s["v"][dim]
        x[:,dim] -= s["x"][dim]
        x[:,dim] *= scale
        
    idx = np.arange(len(x))
    x, idx = x[ok], idx[ok]
    if v is not None:
        v = v[ok]

    return 1e3*x/h100, v, idx

def v_circ(m, r):
    """ v_circ returns the circular velocity (in km/s) for a given mass (Msun)
    and radius (pkpc)
    """
    return 655.8 * (m/1e14)**0.5 * (r/1e3)**-0.5

def profile_info(params, x, ok=None, order=None):
    """ profile_info returns basic information about the spherically averaged
    profile of a halo. x is the the position of particles (pkpc) relative to
    the halo center, mp is the particle mass[es] (Msun), and ok flags 
    which particles should be used ignored (e.g. if you're iteratively removing
    particles after unbinding). Set to None if not needed.

    returns r_max, v_max, PE/Vmax, and the order of particles according to
    their radii.
    """

    mp, h100 = params["mp"], params["h100"]
    mp /= h100
    
    r = np.sqrt(np.sum(x**2, axis=1))

    if ok is not None:
        r[np.isnan(r)] = np.max(r[~np.isnan(r)])
    
    if order is None: order = np.argsort(r)
    if ok is not None and np.sum(ok) == 0:
        return 0.0, 0.0, np.zeros(len(x)), order
    
    r_sort = np.sqrt(r[order]**2 + (params["eps"])**2)
    dm = np.ones(len(r_sort))*mp
    if ok is not None: dm[~ok] = 0
    m_enc = np.cumsum(dm[order])

    np.sum(np.isnan)
    v_rot = v_circ(m_enc, r_sort)
    i_max = np.argmax(v_rot)
    rmax, vmax = r_sort[i_max], v_rot[i_max]
    
    dW = np.zeros(len(m_enc))

    dr = r_sort[1:] - r_sort[:-1]
    dr[dr == 0] = np.finfo(dr.dtype).eps
    r_scaled = (r_sort[1:]*r_sort[:-1]) / dr
    dW[:-1] = v_circ(m_enc[:-1], r_scaled)**2

    # I don't think this is true when there's a tidal radius. Come back to this
    # later.
    vesc_lim = v_circ(m_enc[-1], r_sort[-1]) * np.sqrt(2)

    W = (np.cumsum(dW[::-1])[::-1] + 2*vesc_lim**2)/vmax**2
    out = np.zeros(len(W))
    out[order] = W
    
    return rmax, vmax, -out, order

def rank_by_quantile(quantiles, x, idx, n_max):
    """ rank_by_quantile breaks the paritcles represented by x (any value) and
    idx (index within the full set of particles, of length n_max). quantiles
    gives the upper and lower quantile edges for each rank.
    """
    q = np.quantile(x, quantiles)
    valid_edges = np.zeros(len(q), dtype=bool)
    ranks = np.ones(n_max, dtype=int)*NIL_RANK

    n_valid_edges, curr_rank = 0, 0
    for i in range(len(q) - 1):
        ok = (ranks[idx] == -1) & (x < q[i+1])
        if np.sum(ok) < MIN_PARTICLES_PER_QUANTILE: continue
        ranks[idx[ok]] = curr_rank
        curr_rank += 1
        
        valid_edges[i+1] = True
        if n_valid_edges == 0: valid_edges[i] = True
        n_valid_edges += 1
        
    return ranks, q[valid_edges]

def accreted_early(curr_snap, accrete_snap, f_accrete):
    accreted_so_far = accrete_snap <= curr_snap
    mid_snap = np.quantile(accrete_snap[accreted_so_far], f_accrete)
    return accrete_snap <= mid_snap

def t_relax_t_orbit(Nr, r, eps):
    return Nr/4 * (np.log(r**2/eps**2 + 1) +
                   eps**2-2*r**2/(3*(eps**2+r**2)) -
                   np.log(3/2))**-1

def t_orbit(Mr, r):
    return 7.5 * (r/40)**1.5 * (1e10/Mr)**0.5

def ranked_np_profile_matrix(ranks, idx, r, bins):
    n_rank = 1 + np.max(ranks)
    n_bins = len(bins)

    M = np.zeros((n_bins-1, n_rank))

    for i in range(n_rank):
        ri = r[ranks[idx] == i]
        N, _ = np.histogram(ri, bins=bins)
        #M[:,i] =  np.cumsum(N)
        M[:,i] = N
        
    return M
