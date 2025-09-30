import numpy as np
import matplotlib.pyplot as plt
import numpy.random as random
import scipy.stats as stats
import abc
import scipy.optimize as optimize
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.special as special
from . import lib
from . import util
from . import sampling
import time

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
'mstar' - The stellar mass assinged to this halo (Msun)
"""

NIL_RANK = -1

DEFAULT_QUANTILE_EDGES = np.zeros(13)
DEFAULT_QUANTILE_EDGES[1:] = 10**np.linspace(-4, 0, 12)
#DEFAULT_QUANTILE_EDGES = np.zeros(25)
#DEFAULT_QUANTILE_EDGES[1:] = 10**np.linspace(-4, 0, 24)

DEFAULT_R_BINS = np.zeros(102)
DEFAULT_R_BINS[1:] = 10**np.linspace(-2.5, 0, 101)

MIN_PARTICLES_PER_QUANTILE = 10

DEFAULT_CORE_PARTICLES = 30

#########################
# Abstract base classes #
#########################

class ProfileShapeModel(abc.ABC):
    """ AbstractProfile is an abstract base class for galaxy profile models. It
    allows you to convert a half-light radius into an enclosed stellar mass
    profile.
    """ 
    
    @abc.abstractmethod
    def m_enc(self, m_star, r_half, r, **kwargs):
        """ m_enc returns the enclosed mass profile as a function of 3D radius,
        r. m_star is the asymptotic stellar mass of the galaxy, r_half is 3D
        half-light radius of the galaxy. is_2d is True if the input radius is
        a projected 2d radius and false otherwise. Returned masses will be in
        the same units as m_star.
        """
        pass


    @abc.abstractmethod
    def set_r_half_is_2d(self, r_half_is_2d):
        """ set_rhalf_is_2d sets the variable r_half_is_2d, which determines
        how the r_half value in the m_enc function is interpreted. This
        variable should be set to True.
        """
        pass

    @abc.abstractmethod
    def r2d_r3d(self):
        """ returns the ratio of r2d/r3d for whatever internally generated
        parameters were used int he previous call to m_enc().
        """
        pass
    
    @abc.abstractmethod
    def density(self, m_star, r_half, r, **kwargs):
        """  DEPRECATED

        density returns the local density as a function of 3D radius, r.
        m_star is the asymptotic stellas mass of the galaxy, r_hald is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units of m_star. Only the r input will be vecotrized.
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
    def r_half_is_2d(self):
        """ r_half_is_2d returns True if the model returns 2D radii or 3D
        radii.
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

class SFHModel(abc.ABC):
    """ SFHModel is an abstract base class for models that compute star
    formation histories.
    """
    def sfh(self, **kwargs):
        """ sfh returns a (2, n_snap) array. sfh[0,:] = time, sf[1,:] = M*.
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


class FeHMeanModel(abc.ABC):
    """ MeanFeHModel is an abstract base class for models of the
    mean [Fe/H] values of galaxies
    """
    @abc.abstractmethod
    def Fe_H(self, **kwargs):
        """ Fe_H returns the mean [Fe/H] of the galaxy
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
    
class FeHMDFModel(abc.ABC):
    """ FeHMDFModel is an abstract base class for models of the metallicitiy
    distirbution funciton of [Fe/H]. 
    """
    @abc.abstractmethod
    def mdf(self, Fe_H_mean, **kwargs):
        """ mdf returns an unnormalized  1d funtion defined on the range
        [-10, 5] which gives the pdf of metallicities in the galaxy, given 
        a target mean metallicity.
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

class FeHProfileModel(abc.ABC):
    """ FeHProfileModel is an abstract base class for models of the metallicitiy
    distirbution funciton of [Fe/H]. 
    """
    @abc.abstractmethod
    def Fe_H_profile(self, r_half, FeH, ranks, **kwargs):
        """ FeH_profile returns the indices that a given set of particle
        metallicities would need to have in order to obey the target
        metallicity profile. Also returns the target metallicity gradient.
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

class MetalCorrelationModel(abc.ABC):
    """ MetalCorrelationModel computes stellar ages based on [Fe/H] values
    """

    def a_form(self, Fe_H, sfh, **kwargs):
        """ ages returns the ages of star particles
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

############################
# Particle ranking classes #
############################
    
class AbstractRanking(abc.ABC):
    """ AbstractRanking is an abstract base class for various models that 
    rank particles.
    """

    def __init__(self, ranks, core_idx, rvir, r_bins=DEFAULT_R_BINS):
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
        self.core_idx = core_idx
        if len(ranks) == 0:
            self.n_max = 0
            self.mp_star = np.empty(0)
        else:
            self.n_max = np.max(ranks)
            self.mp_star = np.zeros(len(ranks))
        
        self.idx = None
        self.x = None
        self.v = None
        self.xc = None
        self.vc = None
        self.r_bins = r_bins
        self.rvir=rvir
        self.M = None
        self.r = None
        
    def load_particles(self, x, v, idx):
        """ load_particles loads paritcle properties into the ranking. Must
        be called before mp_star() or either of the core functions. v may be
        None unless you call a later method which works with velocities (i.e.
        core_v)
        """
        
        # Initially, this did a shallow copy of the arrays. 
        self.x = np.array(x)
        if v is not None:
            self.v = np.array(v)
        else:
            self.v = None
        self.idx = idx

        self.xc = self.core_x()
        for dim in range(3): self.x[:,dim] -= self.xc[dim]
        
        if self.v is not None:
            self.vc = self.core_v()
            for dim in range(3): self.v[:,dim] -= self.vc[dim]

        self.r = np.sqrt(np.sum(self.x**2, axis=1))/self.rvir
        self.M = ranked_np_profile_matrix(self.ranks, self.idx,
                                          self.r, self.r_bins)

    def core_x(self):
        """ core_x returns the position of the halo core. You must have called
        load_particles() first.
        """
        if self.x is None: raise ValueError("x not set")
 
        core = self.x[np.searchsorted(self.idx, self.core_idx)]
        # I'm extremely confused by something that I think is a bug in
        # np.median here???
        return np.median(core, axis=0)
        
    def core_v(self):
        """ core_v returns the velocity of the halo core. You must have called
        load_particles() first with a non-None argument for v.
        """
        if self.v is None: raise ValueError("v not set")
        
        core = self.v[np.searchsorted(self.idx, self.core_idx)]
        return np.median(core, axis=0)
        
    def set_mp_star(self, kwargs, profile_model, r_half, m_star):
        M = self.M

        # This is some sort of crazy error caused by Rockstar problems.
        if M is None:
            self.mp_star[:] = m_star/len(self.mp_star)
            return self.mp_star
        
        n = M.shape[1]

        profile_kwargs = profile_model.trim_kwargs(kwargs)
        m_star_enc_target = profile_model.m_enc(
            m_star, r_half, self.r_bins[1:]*self.rvir,
            **profile_kwargs)
        
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

        mp_star_tot = np.sum(self.mp_star)
        if mp_star_tot == 0:
            correction_frac = 1.0
        else:
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

class EnergyRanking(AbstractRanking):
    def __init__(self, p, E, rvir, vmax,
                 core_particles=DEFAULT_CORE_PARTICLES,
                 quantile_edges=DEFAULT_QUANTILE_EDGES):
        """ Takes a set of particles read in "smooth" mode, their corresponding
        energies, and ok, a boolean flag indicating which particles have valid
        energy values. Regardless of what ok is set to, 
        """

        if len(p) == 0:
            self.E_edges = np.empty(0)
            super(EnergyRanking, self).__init__(np.empty(0), np.empty(0), 0)
            return

        ok = p["ok"]
        idx = np.arange(len(ok), dtype=int)[ok]

        # This can happen if there are rockstar errors which screw up the
        # particle tagging. It's rare but happens ever few-thousand halos.
        if len(E) <= core_particles:
            core_idx = np.arange(len(E), dtype=int)
            ranks = np.zeros(len(E), dtype=int)
            super(EnergyRanking, self).__init__(ranks, core_idx, rvir)
            return
        
        core_q = core_particles / len(E)
        core_E = np.quantile(E, core_q)
        core_idx = np.where(E <= core_E)[0]

        E_edges = np.linspace(-15, 0, 61)
        quantile_edges = [np.sum(E/vmax**2 < E_edges[i])/ len(E)
                          for i in range(len(E_edges))]

        ranks, self.E_edges = rank_by_quantile(
            quantile_edges, E[idx]/vmax**2, idx, len(p))
        
        super(EnergyRanking, self).__init__(ranks, core_idx, rvir)
    
##################################
# Specific Model Implementations #
##################################

class PlummerProfile(ProfileShapeModel):
    """ PlummerProfile models a galaxy's mass distribution as a Plummer sphere.
    """

    def __init__(self, r_half_is_2d=False):
        self.r_half_is_2d = r_half_is_2d

    def set_r_half_is_2d(self, r_half_is_2d):
        self.r_half_is_2d = r_half_is_2d
        
    def m_enc(self, m_star, r_half, r, **kwargs):
        """ m_enc returns the enclosed mass profile as a function of 3D radius,
        r. m_star is the asymptotic stellar mass of the galaxy, r_half is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units as m_star.
        """

        if not self.r_half_is_2d:
            r_half *= self.r2d_r3d()

        a = r_half
        return m_star*r**3 /(r**2 + a**2)**1.5

    def r2d_r3d(self):
        return (2**(2/3) - 1)**0.5
    
    def density(self, m_star, r_half, r, **kwargs):
        """ density returns the local density as a function of 3D radius, r.
        m_star is the asymptotic stellas mass of the galaxy, r_hald is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units of m_star.
        """
        a = r_half
        return 3*m_star/(4*np.pi*a**3) * (1 + r**2/a**2)**(-5/2)
    
    def var_names(self):
        return []
    
class HernquistProfile(ProfileShapeModel):
    """ HernquistProfile models a galaxy's mass distribution as a Hernquist sphere.
    """
    def __init__(self, r_half_is_2d=False):
        self.r_half_is_2d = r_half_is_2d

    def set_r_half_is_2d(self, r_half_is_2d):
        self.r_half_is_2d = r_half_is_2d
    
    def m_enc(self, m_star, r_half, r, **kwargs):
        """ m_enc returns the enclosed mass profile as a function of 3D radius,
        r. m_star is the asymptotic stellar mass of the galaxy, r_half is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units as m_star.
        """
        if not self.r_half_is_2d:
            r_half *= self.r2d_r3d()
        
        a = r_half/1.8153
        return m_star*r**2 /(r+a)**2.
    
    def density(self, m_star, r_half, r, **kwargs):
        """ density returns the local density as a function of 3D radius, r.
        m_star is the asymptotic stellas mass of the galaxy, r_hald is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units of m_star.
        """
        a = r_half
        return m_star/(2*np.pi*a**3) * (a**4/(r*(r+a)**3))

    def r2d_r3d(self):
        return 1/1.33#1/(1 + np.sqrt(2))
    
    def var_names(self):
        return []

class EinastoProfile(ProfileShapeModel):
    """ EinastoProfile is is a profile following an Einasto profile. This
    can be a somewhat reasonable approximation of a de-projected Sersic profile
    and more generally is just a clean way to carve out a range of profile
    shapes. (See Cardone et al. 2005 for a bunch of analytical Einasto math,
    although note that Cardone used gamma as the index)

    rho = rho_s * exp(-2/alpha * ((r/rs)^alpha - 1))

    M(r) = Mtot*Gamma(3/alpha, 2*(r/rs)^alpha/alpha)

    (here, Gamma is the normalized lower incomplete gamma function and rs
    is the radius where the log slope = -2.)
    """
    def __init__(self, alpha):
        self.alpha = alpha
        def f(x):
            return special.gammainc(3/alpha, 2*x**alpha/alpha) - 0.5
        self.r_half_rs = optimize.root_scalar(f, x0=1, x1=2).root

    def set_r_half_is_2d(self, r_half_is_2d):
        self.r_half_is_2d = r_half_is_2d
        
    def m_enc(self, m_star, r_half, r, **kwargs):
        """ m_enc returns the enclosed mass profile as a function of 3D radius,
        r. m_star is the asymptotic stellar mass of the galaxy, r_half is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units as m_star.
        """

        if self.r_half_is_2d:
            r_half /= self.r2d_r3d()
        
        rs = r_half/self.r_half_rs
        x = r/rs
        return special.gammainc(
            3/self.alpha, 2*x**self.alpha/self.alpha)*m_star

    def r2d_r3d(self):
        return (special.gammaincinv(3/self.alpha, 0.5)*self.alpha/2)**(1/self.alpha)
    
    def density(self, m_star, r_half, r, **kwargs):
        """ density returns the local density as a function of 3D radius, r.
        m_star is the asymptotic stellas mass of the galaxy, r_hald is the 2D
        half-light radius of the galaxy. Returned masses will be in the same
        units of m_star.
        """
        rs = r_half/self.r_half_rs
        rho_s = (m_star*(2/self.alpha)**(3/self.alpha) *
                 (self.alpha/(4*np.pi*rs**3)))
        return rho_s * np.exp(-2/self.alpha * ((r/rs)**self.alpha - 1))
        
    def var_names(self):
        return []

class DeprojectedSersicProfile(ProfileShapeModel):
    def __init__(self, n_sersic="mansfield25", r_half_is_2d=False):
        """ n_sersic can be one of two different values:
        - float: this will be the n_sersic value used
        - string: the name of an sersic index relation. Curently, only 
        "mansfield25" is supported.
        - r_half_is_2d determines whether the r_half input to
        """
        self.n_sersic = n_sersic
        self.previous_n_sersic_sample = 2.0 # just to prevent crashes
        self.r_half_is_2d = r_half_is_2d

    def set_r_half_is_2d(self, r_half_is_2d):
        self.r_half_is_2d = r_half_is_2d
        
    def density(self, m_star, r_half_2d, r, **kwargs):
        raise ValueError("Not yet implemented")
            
    def m_enc(self, m_star, r_half_2d, r, **kwargs):
        if type(self.n_sersic) == str:
            if self.n_sersic == "mansfield25":
                n_sersic = self._mansfield25_sample(m_star)
            else:
                raise ValueError(
                    "Unrecognized sersic model, '%s'" % self.n_sersic)
        else:
            n_sersic = self.n_sersic

        self.previous_n_sersic_sample = n_sersic

        if not self.r_half_is_2d:
            r_half_2d = self.r2d_r3d()*r_half_2d
            
        # Use Lima Neto+ (1999) fit to n.
        # According to Vitral & Mamon (2020), this runs into problems for
        # n ~< 2 /and/ r ~< Re/10. This is fine for us. (I would have used
        # Vitral & Mamon's fit, but it's much harder to use.)
        pn = 1 - 0.6097/n_sersic + 0.05463/n_sersic**2
        
        a = 2*n_sersic
        bn = special.gammaincinv((3-pn)*n_sersic, 0.5)
        
        r_eff_a = np.exp((0.6950 - np.log(1/n_sersic))/n_sersic - 0.1789)
        a = r_half_2d/r_eff_a
        r3d_a = r_eff_a / self.r2d_r3d() * 10

        r_half_3d = r_half_2d / self.r2d_r3d()
        return special.gammainc((3-pn)*n_sersic, bn*(r/r_half_3d)**(1/n_sersic))*m_star
    
    def r2d_r3d(self):
        # From Lima & Neto
        nu = 1/self.previous_n_sersic_sample
        return 1.356 - 0.0293*nu + 0.0023*nu**2
    
    
    def _mansfield25_sample(self, mstar):
        # Fits to the 16th, 50th, and 84th quantiles
        x = np.log10(mstar)
        q16 = -0.18 + 0.35*(1 + special.erf(1.65*(x - 10.95)))
        q50 = -0.07 + 0.35*(1 + special.erf(1.44*(x - 10.56)))
        q84 = +0.08 + 0.35*(1 + special.erf(1.10*(x - 10.15)))
        
        # 68% skew parameter
        s68_signed = np.abs((q16 + q84 - 2*q50)/(q84 - q16))
        s68 = np.abs(s68_signed)

        # fit to ln(kappa)
        z = -np.log(1 - s68)
        a = 0.6332
        b = 0.3734
        c = -0.008641
        d = 0.003886
        log_k = b*np.log(z/a) + c*np.log(z/a)**2 + d*z/a
        k = 10**log_k * np.sign(s68_signed)
        # fit to alpha/(Q_84 - Q_16)
        a = 0.5029
        b = 1.986
        c = 0.3789
        d = 1.289
        alpha_qq = a*(1 - s68**(b + c*s68))**2
        alpha = alpha_qq * (q84 - q16)

        mu = q50
        q = random.random()
        # inverse transfer sampling
        log_n = mu + (alpha/k)*(np.exp(k*np.sqrt(2)*special.erfinv(2*q - 1)) - 1)
        return 10**log_n
    
    def var_names(self):
        return []
    
class Nadler2020RHalf(RHalfModel):
    """ Nadler20202RHalf models galaxies according to the z=0 size-mass relation
    in Nadler et al. 2020. (https://arxiv.org/abs/1912.03303)

    This is calibrated by fitting a galaxy-halo model to Milky Way satellites.
    """
    
    def __init__(self, A=27e-6, n=1.07, R0=10e-3, scatter=0.63):
        """ The constructor for Nadler2020RHalf allows you to change the
        parameters of the fit. Please consult the paper for discussion on what
        these parameters mean. This allows you to sample the posterior
        distribution for these parameters. Scatter is logarithmic in base-10.
        """
        self.A = A
        self.n = n
        self.R0 = R0
        self.sigma_log_R = scatter

    def r_half(self, **kwargs):
        """ r_half returns the half-mass radius of a a galaxy in physical kpc.
        Required keyword arguments:
         - rvir
        """
        # inputs and outputs are in pMpc (no h).
        log_R = np.log10(self.A * (rvir/self.R0)**self.n)
        if scatter <= 0:
            return 10**log_R
        else:
            log_scatter = self.sigma_log_R*random.normal(
                0, 1, size=np.shape(rvir))
            return 10**(log_R + log_scatter)

    def r_half_is_2d(self):
        return True
        
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["rvir"]


class FixedRHalf(RHalfModel):
    def __init__(self, ratio=0.015, scatter=0.0):
        self.ratio = ratio
        self.sigma_log_R = scatter

    def r_half(self, rvir=None):
        if self.sigma_log_R <= 0.0:
            return rvir*self.ratio
        else:
            log_scatter = self.sigma_log_R*random.normal(
                0, 1, size=np.shape(rvir))
            return rvir*self.ratio * 10**(log_scatter)

    def r_half_is_2d(self):
        return False
        
    def var_names(self):
        return ["rvir"]

class Jiang2019RHalf(RHalfModel):
    """ Jiang2019Rhalf models galaxies according to the z-dependent size-mass
    relation in Jiang et al. 2019 (https://arxiv.org/abs/1804.07306; Appendix
    D). This model is fit to   galaxies in NIHAO/VELA and has been rescaled
    to  match the noramlization of the z-dependence seen in GAMA+CANDLES when
    combined with abundance matching (Somerville et al. 2018;
    https://arxiv.org/abs/1701.03526; M* > 1e8 Msun).
    """
    def __init__(self, scatter=0.2):
        """ The constructor for Jiang2019RHalf allows you to change the
        intrinsic scatter in the relation. I eyeballed sigma_log_R at 0.2 from
        their plots, so this value in particular is worth modifying.
        """
        # I eyeballed the 0.2 from their plots.
        self.sigma_log_R = scatter
    
    def r_half(self, rvir=None, cvir=None, z=None):
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
        if self.sigma_log_R > 0.0:
            log_scatter = self.sigma_log_R*random.normal(
                0, 1, size=np.shape(rvir))
            return 10**(np.log10(R) + log_scatter)
        else:
            return R

    def r_half_is_2d(self):
        return False
        
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["rvir", "cvir", "z"]

class Carlsten2021RHalf(RHalfModel):
    """ Carlsten2022RHalf models galaxies according to a z=0 size-mass relation
    calibrated off of observations in the local volume.
    (https://arxiv.org/pdf/2105.03435.pdf) (10^5.5 Msun < M* < 10^8.5 Msun).
    """
    def __init__(self, scatter=0.181, a=1.071, b=0.247):
        """ The constructor for Carlsten2021RHalf is just a power law
        """
        self.sigma_log_R = scatter
        self.a = a
        self.b = b
    
    def r_half(self, rvir=None, cvir=None, z=None):
        """ r_half returns the half-mass radius of a a galaxy in physical kpc.
        Required keyword arguments:
         - rvir
         - cvir
         - z
        """
        R = 10**(self.a + self.b*np.log10(0.247))
        if self.sigma_log_R > 0:
            log_scatter = self.sigma_log_R*random.normal(
                0, 1, size=np.shape(rvir))
            return 10**(np.log10(R) + log_scatter)
        else:
            return R

    def r_half_is_2d(self):
        return True
        
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["rvir", "cvir", "z"]

class UniverseMachineSFH(SFHModel):
    """ SFHModel is an abstract base class for models that compute star
    formation histories.
    """
    def sfh(self, um_sfh=None):
        return um_sfh

    def var_names(self):
        return ["um_sfh"]

class DarkMatterSFH(SFHModel):
    """ SFHModel is an abstract base class for models that compute star
    formation histories.
    """
    def sfh(self, dm_mah=None):
        return dm_mah

    def var_names(self):
        return ["dm_mah"]

class Kirby2013Metallicity(FeHMeanModel):
    """ Kirby2013Metallicity models galaxy metallicity according to the
    z-agnostic fit in Kirby et al. 2013 (https://arxiv.org/pdf/1310.0814.pdf;
    Equation 4).
    """
    def __init__(self, scatter=0.17):
        """ The constructor for Kirby2013Metallicity allows you to change the
        intrinsic scatter in the relation.

        That 0.17 comes from scraping the data points in Fig 9 with WPD.
        """
        self.sigma_Fe_H = scatter
    
    def Fe_H(self, mstar=None):
        """ Fe_H returns the metallicity of a given galaxy.
        Required keyword arguments:
         - mstar
        """
        if mstar is None: raise ValueError("mstar not supplied")
        
        Fe_H_mean = -1.69 + 0.30*np.log10(mstar/1e6)
        if self.sigma_Fe_H > 0.0:
            Fe_H_mean = random.normal(Fe_H_mean, self.sigma_Fe_H)
        return Fe_H_mean

        #return random.normal(Fe_H_mean, self.Fe_H_dist_width, size=n_part)


    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["mstar"]

class Kirby2013MetallicityVariable(FeHMeanModel):
    """ 
    """
    def __init__(self, offset = -1.69, sigma_Fe_H=0.17, mass_slope=0.3, z_slope=-0.11):
        """

        """
        self.sigma_Fe_H = sigma_Fe_H
        self.offset = offset 
        self.slope = z_slope 
        self.mass_slope = mass_slope
        
    def Fe_H(self, mstar=None, z=None, no_scatter=False):
        """ Fe_H returns the metallicity of a given galaxy.
        Required keyword arguments:
         - mstar
        """
        if mstar is None: raise ValueError("mstar not supplied")
        if z is None: raise ValueError("z not supplied")

        Fe_H_mean = self.offset + self.mass_slope*np.log10(mstar/1e6) + self.slope*z 
        if not no_scatter:
            Fe_H_mean = random.normal(Fe_H_mean, self.sigma_Fe_H)
        return Fe_H_mean
       
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["mstar", "z"]
    
class Kirby2013MDF(FeHMDFModel):
    def __init__(self, model_type="gaussian"):
        """ accepted model_types are:
          - "gaussian"
          - "leaky box"

        Maybe I'll add pre-enriched and accretion later if people care.

        Note: the leaky box model is nice because it's a single parameter
        model, but it doesn't fit higher mass galaxies as well as more
        complicated models.
        """
        if model_type not in ["leaky box", "gaussian"]:
            raise ValueError("Unknown model_type, %s" % model_type)
        self.model_type = model_type

        if model_type == "gaussian":
            sigmas = [0.29, 0.28, 0.44, 0.36, 0.38, 0.34, 0.35, 0.39, 0.45,
                      0.6, 0.59, 0.57, 0.6, 0.39, 0.33, 0.47, 0.32, 0.55, 0.54,
                      0.36, 0.31, 0.47]
            self.mdf_sigma_mean = np.mean(sigmas)
            self.mdf_sigma_sigma = np.std(sigmas)
        elif model_type == "leaky box":
            def leaky_mean(p_eff):
                pdf = lambda Fe_H: 10**Fe_H/p_eff*np.exp(-10**Fe_H/p_eff)
                return (integrate.quad(lambda Fe_H: Fe_H*pdf(Fe_H), -10, 5)[0] /
                        integrate.quad(lambda Fe_H: pdf(Fe_H), -10, 5)[0])
            p_eff = 10**np.linspace(-4, 2, 200)
            mean_Fe_H = [leaky_mean(p_eff[i]) for i in range(len(p_eff))]
            self.f_leaky_mean = interpolate.interp1d(mean_Fe_H, np.log10(p_eff))
        
    def mdf(self, Fe_H_mean): 
        if self.model_type == "gaussian":
            sigma = random.randn()*self.mdf_sigma_sigma + self.mdf_sigma_mean
            return sampling.PDF(lambda x: np.exp(-(x-Fe_H_mean)**2 / 
                                                 (2*sigma**2)), -10, 5)
        elif self.model_type == "leaky box":
            log_p_eff = self.f_leaky_mean(Fe_H_mean)
            p_eff = 10**log_p_eff
            return sampling.PDF(lambda Fe_H: 10**Fe_H/p_eff*np.exp(-10**Fe_H/p_eff), -10, 5)
        else:
            assert(0)

    def var_names(self):
        return []

class FlatFeHProfile(FeHProfileModel):
    def Fe_H_profile(self, r_half, FeH, ranks):
        return FeH, 0.0

    def var_names(self):
        return []

class Taibi2022FeHProfile(FeHProfileModel):
    """ Requires that you're using EnergyRanking
    """
    def __init__(self, r_r50_max=6):
        FeHProfileModel.__init__(self)
        delta_Fe_H = [
            -0.101, -0.23, -0.09, -0.14, -0.21, -0.04, -0.32, -0.15, -0.15,
            -0.07,   0.00,  0.00,  0.00,  0.00, -0.20,  0.00,  0.00, -0.39,
            -0.10,  -0.46, -0.06,  0.00, -0.12, -0.16, -0.29,  0.00, -0.10,
            -0.08,  -0.35, -0.07
        ]
        self.mean_delta_Fe_H = np.mean(delta_Fe_H)
        self.std_delta_Fe_H = np.std(delta_Fe_H)
        self.r_r50_max = r_r50_max

    def Fe_H_profile(self, r_half, Fe_H, ranks):
        delta_Fe_H = random.randn()*self.std_delta_Fe_H + self.mean_delta_Fe_H

        def f_Fe_H(r):
            out = delta_Fe_H*(r/r_half)
            min_Fe_H = np.abs(delta_Fe_H*self.r_r50_max)
            out[out < -min_Fe_H] = -min_Fe_H
            out[out > min_Fe_H] = min_Fe_H
            return out

        """
        # You need to be unusually careful here since you're using both
        # information from the energy and star snapshots. I've appended E or
        # star to each variable based on the length/indexing of the underlyng
        # array.
        
        idx_star, idx_E = ranks.idx, ranks.E_idx
        
        x_star = ranks.x
        x_all = np.zeros((len(ranks.ranks),3))
        x_all[idx_star] = x_star
        x_E = x_all[idx_E]
        
        r_E = np.sqrt(np.sum(x_E**2, axis=1))
        order_E = ranks.order
        mp_E = ranks.mp_star[idx_E]
        is_bound_E = (ranks.ranks != NIL_RANK)[idx_E]
        
        # f_Fe_H is the metallicity gradient of the galaxy with an arbitary
        # offset. [Fe/H] is kept constant outside r_r50_max.t
        
        r_E_sort = r_E[order_E]
        mp_E_sort = ranks.mp_star[ranks.order]
        
        # The idea here is that you want to evaluate the metallicity at a
        # radius corresponding to a random particle with a similar energy
        # to remove biases where particles are pericenter are given higher
        # metallicities. You don't want this because correlations with orbital
        # phase lead to immediate time evolution after tagging.

        window = 10 # This number can be whatever, it just can't be large
        idx = np.arange(len(r_E_sort), dtype=int)
        min_idx = np.maximum(idx - window, 0)
        max_idx = np.minimum(idx + window, len(r_E_sort)-1)
        idx = random.randint(min_idx, max_idx)
        r_eval = r_E_sort[idx]

        Fe_H_i = f_Fe_H(r_eval) # intial Fe/H before transforms/scatter
        input_std = np.std(Fe_H)
        w_std = weighted_std(Fe_H_i, mp_E_sort)

        target_std = np.sqrt(max(input_std**2 - w_std**2, 0))
        # This std/5 comes from basically nowhere: I don't think the
        # literature has much guidance on what this should be, but I don't
        # want it to be a delta function.
        scatter = min(input_std/5, target_std)
        Fe_H_i += random.randn(len(Fe_H_i)) * scatter

        avg = np.average(Fe_H_i, weights=mp_E_sort)
        Fe_H_i = Fe_H_i - avg + np.mean(Fe_H)

        Fe_H_i = weighted_abundance_match(Fe_H, Fe_H_i, mp_E_sort)
        
        # Unsort Fe/H back into the orignal order
        idx_sort = np.arange(len(Fe_H_i), dtype=int)[ranks.order]
        Fe_H_0 = np.zeros(len(Fe_H_i))
        Fe_H_0[idx_sort] = Fe_H_i

        out = f_Fe_H()

        out = np.zeros(len(ranks.ranks))
        out[idx_E] = Fe_H_0
        """        

        out = np.zeros(len(ranks.ranks))
        r = np.sqrt(np.sum(ranks.x**2, axis=1))
        out[ranks.idx] = f_Fe_H(r)
        
        return out, delta_Fe_H
    
    def var_names(self):
        return []

def weighted_std(x, w):
    avg = np.average(x, weights=w)
    var = np.average((x - avg)**2, weights=w)
    return np.sqrt(var)

def weighted_abundance_match(x, y, wy):
    x = np.sort(x)

    order = np.argsort(y)
    orig_idx = np.arange(len(y))[order]
    sy, swy = y[order], wy[order]

    sy = np.asarray(sy, dtype=np.float64)
    swy = np.asarray(swy, dtype=np.float64)

    P = np.cumsum(swy) / np.sum(swy)
    f_idx = P*(len(x) - 1)

    prev_idx = np.array(np.floor(f_idx), dtype=int)
    frac = f_idx - prev_idx

    # This stuff here is needed due to floating point numbers not being real
    # numbers
    prev_idx[prev_idx >= len(x)] = len(x) - 1
    prev_idx[prev_idx < 0] = 0
    frac[frac > 1] = 1
    
    next_idx = prev_idx + 1
    next_idx[next_idx >= len(x)] = len(x) - 1

    dx = x[next_idx] - x[prev_idx]
    zero_dx = dx <= 0
    dx[zero_dx], frac[zero_dx] = 1, 0
    out = np.zeros(len(y))
    out[orig_idx] = (x[prev_idx] + dx*frac)
    return out

class GaussianCoupalaCorrelation(MetalCorrelationModel):
    """ GaussianCoupalaCorrelation connects ages to metallicities by assuming 
    that thjey can be represented by a Gaussian coupala.
    """

    def __init__(self, rho=0.75):
        self.rho = rho

    def a_form(self, Fe_H, sfh):
        return sampling.gaussian_coupala_sample(Fe_H, sfh, self.rho)
        
    def var_names(self):
        return []

    def trim_kwargs(self, kwargs):
        out = {}
        for key in self.var_names():
            out[key] = kwargs[key]
        return out

class UniverseMachineMStarFit(MStarModel):
    def __init__(self, scatter=0.2):
        self.scatter = scatter
    
    def m_star(self, mpeak=None, z=None):
        """
         Required keyword arguments:
         - mpeak
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
        al = al0 + ala*(a - 1) - al_lna*np.log(a) + alz*z
        b = b0 + ba*(a - 1) + bz*z
        d = d0
        g = 10**(g0 + ga*(a - 1) + gz*z)

        x = np.log10(mpeak/10**log10_M1_Msun)
      
        log10_Ms_M1 = (e - np.log10(10**(-al*x) + 10**(-b*x)) +
                       g*np.exp(-0.5*(x/d)**2))
                       
        log10_Ms_Msun = log10_Ms_M1 + log10_M1_Msun

        if self.scatter > 0.0:
            log_scatter = self.scatter*random.normal(
                0, 1, size=np.shape(mpeak))
            log10_Ms_Msun += log_scatter
        
        Ms = 10**log10_Ms_Msun
        
        return Ms
    
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["mpeak", "z"]
    

class UniverseMachineMScatterGrowing(MStarModel):
    
    def __init__(self, slope=1.963, mpivot=10, gamma=0.0, scatter=0.2):
        self.gamma = gamma #rate of growth 
        self.scatter = scatter #high mass scatter
        self.slope = slope #z=0 slope, rescaled at higher redshifts by UM relation
        self.mpivot = mpivot #pivot mass below which scatter grows 

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
        al = al0 + ala*(a - 1) - al_lna*np.log(a) + alz*z
        
        al = al*self.slope/al0 


        b = b0 + ba*(a - 1) + bz*z
        d = d0
        g = 10**(g0 + ga*(a - 1) + gz*z)

        x = np.log10(mpeak/10**log10_M1_Msun)
      
        log10_Ms_M1 = (e - np.log10(10**(-al*x) + 10**(-b*x)) +
                       g*np.exp(-0.5*(x/d)**2))
                       
        log10_Ms_Msun = log10_Ms_M1 + log10_M1_Msun
        
        if(np.log10(mpeak))>self.mpivot:
            growing_scatter = self.scatter
        else:
            growing_scatter = self.scatter + self.gamma*(np.log10(mpeak)-self.mpivot)
        
        log_scatter = growing_scatter*random.normal(0, 1, size=np.shape(mpeak))
        log10_Ms_Msun += log_scatter
        
        Ms = 10**log10_Ms_Msun
        
        return Ms
    
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """
        return ["mpeak", "z"]
    
class UniverseMachineMStar(MStarModel):
    def m_star(self, um_mstar=None):
        return np.maximum(1, um_mstar)

    def var_names(self):
        return ["um_mstar"]
        
#################################
# General classes and functions #
#################################

def tag_stars(sim_dir, galaxy_halo_model, star_snap=None, E_snap=None,
              target_subs=None, seed=None, energy_method="E_sph"):
    # energy_method can be E_sph, E, or smooth
    if seed is not None:
        random.seed(seed)
        
    # Basic simulation information
    param = lib.simulation_parameters(sim_dir)
    h, hist = lib.read_subhalos(sim_dir)
    scale = lib.scale_factors(sim_dir)

    # Not everyone is going to run UM on their zoom-ins, so don't force this
    # file to be read.
    if "um_mstar" in galaxy_halo_model.var_names():
        um = lib.read_um(sim_dir)
    else:
        um = [None]*len(h)

    if target_subs is None:
        # Don't read in the host: wouldn't make sense.
        target_subs = np.arange(len(h), dtype=int)
    target_subs = np.asarray(target_subs, dtype=int)

    # Set the optional arguments to their default values
    if star_snap is None:
        # Match profiles at the snapshot when the subhalo first crosses the
        # virial radius.
        star_snap = hist["first_infall_snap"]

    idx0 = np.where(target_subs == 0)
    if len(idx0) != 0:
        star_snap[idx0] = h.shape[1]-1
    

    if E_snap is None:
        E_snap = np.zeros(len(h), dtype=int)
        for i in target_subs:
            # Look back an eighth of an orbital time, or to the half-mass
            # scale, whichever is later.
            if star_snap[i] == 0:
                E_snap[i] = 0
            else:
                E_snap[i] = look_back_orbital_time(
                    param, scale, star_snap[i], 0.125, h[i,:], 0.5)
                


    n_smooth = np.zeros(len(h), dtype=int)
    # particles in "all" mode at the energy and tagging snapshots
    p_E, p_star = [None]*len(h), [None]*len(h)

    if energy_method in ["E", "E_sph"]:
        part = lib.Particles(sim_dir, include=[energy_method])
    elif energy_method == "smooth":
        part = lib.Particles(sim_dir)
    else:
        raise ValueError("Unrecognized E_method, '%s'" % E_method)    
    
    for snap in range(h.shape[1]):
        # If nobody needs this snapshot, don't read it in.
        if (snap not in E_snap[target_subs] and
            snap not in star_snap[target_subs]): continue
        
        star_ok = np.zeros(len(h), dtype=bool)
        E_ok = np.zeros(len(h), dtype=bool)
        star_ok[target_subs], E_ok[target_subs] = True, True
        i_star = np.where(star_ok & (star_snap == snap))[0]
        i_E = np.where(E_ok & (E_snap == snap))[0]
        
        # Read in data.

        # Add information at the snapshot stellar masses are assigned
        for i in i_star:
            p = part.read(snap, mode="smooth", halo=i)
            n_smooth[i] = np.sum(p["smooth"])
            
            sx = h[i,snap]["x"]
            sv = h[i,snap]["v"]
            p["x"] -= sx
            p["v"] -= sv
            p_star[i] = p

        # Add information at the snapshot energies are evaluated
        for i in i_E:
            p = part.read(snap, mode="smooth", halo=i)
            
            sx = h[i,snap]["x"]
            sv = h[i,snap]["v"]
            p["x"] -= sx
            p["v"] -= sv
            p_E[i] = p
            
    # Initialize various arrays.
    mp_star, ranks = [None]*len(h), [None]*len(h)
    r_half, m_star = np.ones(len(h))*-1, np.ones(len(h))*-1
    Fe_H = [None]*len(h)
    gal_hists = np.zeros(len(h), dtype=lib.GALAXY_HISTORY_DTYPE)
    
    stars = [None]*len(h)
    for i in target_subs:
        stars[i] = np.zeros(n_smooth[i], dtype=lib.STAR_DTYPE)
    
    # Finally, assign stellar masses and create ranking objects for each halo.
    for i in target_subs:
        if energy_method == "smooth":
            ok = p_E[i]["ok"]
            _, vmax, pe_vmax2 , _ = profile_info(param, p_E[i]["x"], ok=ok)
            pe = pe_vmax2 * vmax**2
            ke = np.sum(p_E[i]["v"]**2, axis=1)/2
            E = pe + ke
        elif energy_method in ["E", "E_sph"]:
            E = p_E[i][energy_method]
            vmax = h["vmax"][i,E_snap[i]]
        else:
            assert(0) # Already tested for.

        rvir = h["rvir"][i,E_snap[i]]
        E[~p_E[i]["ok"]] = np.inf
        ranks[i] = EnergyRanking(p_E[i], E, rvir, vmax)
        
        if len(p_E[i]) == 0: continue

        idx = np.where(p_E[i]["ok"])[0]
        if len(idx) == 0 or len(idx) <= DEFAULT_CORE_PARTICLES: continue
        
        ranks[i].load_particles(p_star[i]["x"][idx], p_star[i]["v"][idx], idx)
        
        # It's some tiny subhalo that Rockstar made a mistake on. Doesn't have
        # any bins with more than 10 particles in them.
        if np.max(ranks[i].ranks) == -1:
            continue

        kwargs = galaxy_halo_model.get_kwargs(
            param, scale, h[i], um[i], star_snap[i])

        stars[i], gal_hists[i] = galaxy_halo_model.star_properties(
            ranks[i], kwargs)

    return stars, gal_hists, ranks

class RetagStarsState(object):
    def __init__(self, sim_dir, galaxy_halo_model):
        self.param = lib.simulation_parameters(sim_dir)
        self.h, self.hist = lib.read_subhalos(sim_dir)
        self.h_cmov, _ = lib.read_subhalos(sim_dir, comoving=True)
        self.scale = lib.scale_factors(sim_dir)

        if "um_mstar" in galaxy_halo_model.var_names():
            self.um = lib.read_um(sim_dir)
        else:
            self.um = [None]*len(self.h)
    
    def get_all(self):
        return self.param, self.h, self.hist, self.h_cmov, self.scale, self.um

def retag_stars(sim_dir, galaxy_halo_model, ranks,
                state=None, star_snap=None,
                target_subs=None, seed=None):
    if seed is not None:
        random.seed(seed)
    # Basic simulation information

    if state is None:
        state = RetagStarsState(sim_dir, galaxy_halo_model)
    param, h, hist, h_cmov, scale, um = state.get_all()
    if star_snap is None:
        star_snap = hist["first_infall_snap"]
    if target_subs is None:
        # Don't read in the host: wouldn't make sense.
        target_subs = np.arange(1, len(h), dtype=int)

    if 0 in target_subs:
        raise ValueError("Cannot do star-tagging on the central halo. " + 
                         "Remove 0 from the target_subs array.")

    gal_hists = np.zeros(len(h), dtype=lib.GALAXY_HISTORY_DTYPE)
    stars = [None]*len(h)
    for i in target_subs:
        stars[i] = np.zeros(len(ranks[i].ranks), dtype=lib.STAR_DTYPE)

    for i in target_subs:
        # It's some tiny subhalo that Rockstar made a mistake on. Doesn't have
        # any bins with more than 10 particles in them.
        if len(ranks[i].ranks) == 0 or  np.max(ranks[i].ranks) == -1: continue

        kwargs = galaxy_halo_model.get_kwargs(
            param, scale, h[i], um[i], star_snap[i])

        stars[i], gal_hists[i] = galaxy_halo_model.star_properties(
            ranks[i], kwargs)

    return stars, gal_hists, state            
    
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

    if snap == orbit_start: return snap

    for snap_start in range(snap, orbit_start, -1):
        if halo["mvir"][snap_start-1] < min_mass_frac*halo["mvir"][snap]:
            break

    return snap_start

class StellarMassModel(object):
    def __init__(self, m_star_model, sfh_model):
        self.m_star_model = m_star_model
        self.sfh_model = sfh_model

class ProfileModel(object):
    def __init__(self, r_half_model, profile_shape_model):
        self.r_half_model = r_half_model
        self.profile_shape_model = profile_shape_model

class MetalModel(object):
    def __init__(self, Fe_H_model, Fe_H_mdf_model, Fe_H_profile_model,
                 correlation_model):
        self.Fe_H_model = Fe_H_model
        self.Fe_H_mdf_model = Fe_H_mdf_model
        self.Fe_H_profile_model = Fe_H_profile_model
        self.correlation_model = correlation_model
        
class GalaxyHaloModel(object):
    def __init__(self, stellar_mass_model, profile_model, metal_model):
        """ GalaxyHaloModel requires a model for the M*-Mhalo relation,
        m_star_model, a model for how the projected half-mass radius and Mhalo
        are related, and a model for the halo profile, profile_model. These
        should be types that inherit from AbstractMstarModel,
        AbstractRHalfModel, AbstractProfileModel, and AbstractMetallicityModel
        respectively.
        """
        self.m_star_model = stellar_mass_model.m_star_model
        self.sfh_model = stellar_mass_model.sfh_model
        self.r_half_model = profile_model.r_half_model
        self.profile_shape_model = profile_model.profile_shape_model
        self.Fe_H_model = metal_model.Fe_H_model
        self.Fe_H_mdf_model = metal_model.Fe_H_mdf_model
        self.Fe_H_profile_model = metal_model.Fe_H_profile_model
        self.correlation_model = metal_model.correlation_model
        
    def star_properties(
            self, ranks, kwargs, r_half=None, m_star=None, Fe_H=None):
        """ star_properties sets the stellar masses of a halo's dark matter
        particles given their positions relative to the halo center, and
        ranking, ranks (type: inherits from AbstractParticleRanking). This
        function accepts the same keyword arguments as its m_star_model and
        r_half_model arguments. You may also fix the half-mass radius and
        stellar mass to whatever you want with r_half (units: pkpc) and m_star
        (units: Msun).
        
        This function also returns the M_star, r_half, and [Fe/H] values that
        it assigned to the halo.
        """
        if m_star is None:
            check_var_names(kwargs, self.m_star_model)
            m_star = self.m_star_model.m_star(
                **self.m_star_model.trim_kwargs(kwargs))
            
        kwargs["mstar"] = m_star

        if r_half is None:
            check_var_names(kwargs, self.r_half_model)
            r_half = self.r_half_model.r_half(
                **self.r_half_model.trim_kwargs(kwargs))

        self.profile_shape_model.set_r_half_is_2d(
            self.r_half_model.r_half_is_2d()
        )
            
        mp_star = ranks.set_mp_star(kwargs, self.profile_shape_model,
                                    r_half, m_star)
        gal_hist = np.zeros(1, dtype=lib.GALAXY_HISTORY_DTYPE)[0]

        if Fe_H is None:
            check_var_names(kwargs, self.Fe_H_model)
            Fe_H_mean = self.Fe_H_model.Fe_H(
                **self.Fe_H_model.trim_kwargs(kwargs))

        mdf = self.Fe_H_mdf_model.mdf(Fe_H_mean,
            **self.Fe_H_mdf_model.trim_kwargs(kwargs))
        sfh = self.sfh_model.sfh(**self.sfh_model.trim_kwargs(kwargs))
        Fe_H = mdf.sample(len(mp_star))
        Fe_H, delta_Fe_H = self.Fe_H_profile_model.Fe_H_profile(
            r_half, Fe_H, ranks,
            **self.Fe_H_profile_model.trim_kwargs(kwargs))
        a_form = self.correlation_model.a_form(Fe_H, sfh,
            **self.correlation_model.trim_kwargs(kwargs))

        gal_hist["sigma_Fe_H_i"] = np.std(Fe_H)
        gal_hist["Fe_H_i"] = Fe_H_mean
        gal_hist["delta_Fe_H_i"] = delta_Fe_H
        gal_hist["a50"] = sfh[0,np.searchsorted(sfh[1,:], sfh[1,-1]*0.5)]
        gal_hist["a90"] = sfh[0,np.searchsorted(sfh[1,:], sfh[1,-1]*0.9)]
        
        stars = np.zeros(len(mp_star), dtype=lib.STAR_DTYPE)
        stars["mp"] = mp_star
        stars["Fe_H"] = Fe_H
        stars["a_form"] = a_form

        gal_hist["m_star_i"] = m_star

        r2d_r3d = self.profile_shape_model.r2d_r3d()
        if self.r_half_model.r_half_is_2d():
            gal_hist["r_half_2d_i"] = r_half
            gal_hist["r_half_3d_i"] = r_half/r2d_r3d
        else:
            gal_hist["r_half_2d_i"] = r_half*r2d_r3d
            gal_hist["r_half_3d_i"] = r_half
        
        return stars, gal_hist
    
    def var_names(self):
        """ var_names returns the names of the variables this model requires.
        """

        models = [self.m_star_model,
                  self.sfh_model,
                  self.r_half_model,
                  self.profile_shape_model,
                  self.Fe_H_model,
                  self.Fe_H_mdf_model,
                  self.Fe_H_profile_model,
                  self.correlation_model]
        name_set = [m.var_names() for m in models]

        out_dict = { }
        for n in name_set:
            for i in range(len(n)):
                out_dict[n[i]] = None
            
        if "rvir" not in out_dict: out_dict["rvir"] = None
        return sorted(list(out_dict.keys()))

    def get_kwargs(self, params, scale, halo, um, snap):
        """ get_kwargs creates the kwargs that that are needed for a call to
        set_mp_star.
        """
        h100 = params["h100"]
        
        z = 1/scale[snap] - 1
        
        var_names = self.var_names()

        kwargs = { }
        
        for i in range(len(var_names)):
            if var_names[i] == "z":
                kwargs["z"] = z
                continue

            if var_names[i] == "mpeak":
                x = np.max(halo["mvir"][:snap+1])
            elif var_names[i] == "mstar":
                continue
            elif var_names[i] == "um_mstar":
                x = um[snap]["m_star"]
            elif var_names[i] == "dm_mah":
                x = np.zeros((2, len(scale)))
                x[0,:] = scale
                x[1,:] = running_mpeak(halo["mvir"])
            elif var_names[i] == "um_sfh":
                col_param = lib.colossus_parameters(params)
                cosmo = cosmology.setCosmology("",
                    lib.colossus_parameters(params))
                age = cosmo.age(1/scale - 1)
                dt = age[1:] - age[:-1]
                avg_dt = np.zeros(age.shape)
                avg_dt[0], avg_dt[-1] = dt[0], dt[-1]
                avg_dt[1:-1] = (dt[1:]+dt[:-1]) / 2

                dm = um["sfr"]*avg_dt
                
                x = np.zeros((2, len(scale)))
                x[0,:] = scale
                x[1,:] = np.cumsum(dm) / np.sum(dm) * np.max(um["m_star"])
            else:
                x = halo[snap][var_names[i]]
            
            kwargs[var_names[i]] = x

        return kwargs
    
def running_mpeak(m):
    out = np.zeros(m.shape)
    out[0] = m[0]
    for i in range(1, len(m)):
        out[i] = max(out[i-1], m[i])
    return out

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

def rmax_vmax(params, x, ok=None, bins=300):
    eps = params["eps"]/params["h100"]
    mp, h100 = params["mp"], params["h100"]
    mp /= h100
    
    if ok is None:
        r = np.sqrt(np.sum(x**2, axis=1))
    else:
        r = np.sqrt(np.sum(x[ok]**2, axis=1))

    log_r_min = np.log10(params["eps"]/params["h100"])
    r[r < eps/10] = eps/10
    log_r_max = np.log10(np.max(r))
    d_log_r = (log_r_max - log_r_min)/bins

    log_r = np.log10(r)
    ri = np.asarray(np.floor(((log_r - log_r_min)/d_log_r)), dtype=int)
    ri[ri < -1] = -1
    ri += 1
    
    n = np.bincount(ri)
    mass = np.cumsum(n*mp)
    r = 10**(np.arange(len(n))*d_log_r + log_r_min)
    Vr = v_circ(mass, r)
    
    i_max = np.argmax(Vr) 
    return r[i_max], Vr[i_max]


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
    """ rank_by_quantile breaks the particles represented by x (any value) and
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


DWARF_GALAXY_HALO_MODEL = GalaxyHaloModel(
    StellarMassModel(
        UniverseMachineMStar(),
        UniverseMachineSFH()
    ),
    ProfileModel(
        Jiang2019RHalf(),
        DeprojectedSersicProfile()
    ),
    MetalModel(
        Kirby2013Metallicity(),
        Kirby2013MDF(model_type="gaussian"),
        #Kirby2013MDF(model_type="leaky box"),
        #Taibi2022FeHProfile(),
        FlatFeHProfile(),
        GaussianCoupalaCorrelation()
    )
)

DWARF_GALAXY_HALO_MODEL_NO_UM = GalaxyHaloModel(
    StellarMassModel(
        UniverseMachineMStarFit(),
        DarkMatterSFH()
    ),
    ProfileModel(
        Jiang2019RHalf(),
        DeprojectedSersicProfile()
    ),
    MetalModel(
        Kirby2013Metallicity(),
        Kirby2013MDF(model_type="gaussian"),
        #Kirby2013MDF(model_type="leaky box"),
        #Taibi2022FeHProfile(),
        FlatFeHProfile(),
        GaussianCoupalaCorrelation()
    )
)
