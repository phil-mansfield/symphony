Symlib Documentation
====================
			 
			 
Datatypes
---------

``symlib`` halo data is generally returned as a numpy `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_, which allows for fields to be accessesed and subselected easily. See the :doc:`Getting Started <getting_started>` page for usage examples.

.. note::
   TODO: I don't like the formatting of the fields. The alignment is ugly. Any ideas for fixes?

   Also, the vertical spacing between different entries is super small and I think it looks ugly. But Sphinx doesn't let you just add newlines in, so I don't know how to fix that either.

.. data:: symlib.SUBHALO_DTYPE
		   
    Time-dependent subhalo information (e.g., position) from the Rockstar halo finder. You can get it for all the host's subhalos by calling :func:`symlib.read_subhalos` or :func:`symlib.set_units_halos`. Different instances of this type can use different units, so look at the function that created the array for this information.
	
    **DATA FIELDS:**
	
    * ``"id"`` (*numpy.int32*) - A unique integer identifying each object. Changes from snapshot to snapshot.
    * ``"mvir"``  (*numpy.float32*) - The mass of the halo, :math:`M_{\rm vir}`. When isolated, this an overdensity mass from the Bryan & Norman (1998) definition of the virial overdensity. When deep in a host halo, this is the bound mass. The transition between these two definitions is fuzzy.
    * ``"rvir"`` (*numpy.float32*) - The overdensity radius of the halo, :math:`R_{\rm vir}`.
    * ``"vmax"`` (*numpy.float32*) - The maximum value of the halo's circular rotation curve, :math:`V_{\rm max} = {\rm max}\left\{V_{\rm rot}(r) = \sqrt{G M(<r)/r}\right\}`.
    * ``"rvmax"`` (*numpy.float32*) - The radius where this maximum rotation velocity occurs.
    * ``"cvir"`` (*numpy.float32*) - An estimate of how centrally concentrated the subhalo's mass is, :math:`c_{\rm vir}=R_s/R_{\rm vir}`., :math:`R_s` is the transition radius between shallow inner density slopes (:math:`d \ln(\rho)/d \ln(r)` > -2) and steep outer slopes (i.e. :math:`d \ln(\rho)/d \ln(r)` < -2). *cvir* is estimated  by measuring :math:`V_{\rm max}/V_{\rm rot}(R_{\rm vir})`, assuming an NFW profile, and solving for :math:`R_s`. Because of this, the *value* of :math:`c_{\rm vir}` is only meaningful for halos where the assumption of NFW profiles is reasonable (non-subhalos). However, the *relative ordering* of concentrations will be correct regardless.
    * ``"x"`` (*numpy.float32*) - :math:`x`, the position of the subhalo.
    * ``"v"`` (*numpy.float32*) - :math:`v`, The velocity of the subhalo.
    * ``"ok"`` (*bool*) - True if the subhalo exists during the specified snapshot and False otherwise.
		
.. data:: symlib.HISTORY_DTYPE

    Time-independent subhalo information about the subhalo's entire history in the simulation (e.g. when it first fell into the host halo). You can get it for all the host's subhalos by calling :func:`symlib.read_subhalos`.

	**DATA FIELDS**
	
    * ``"mpeak"`` (*numpy.float32*) - :math:`M_{\rm peak}`, the largest :math:`M_{\rm vir}` that the subhalo ever had. This quantity is often useful for reasoning about subhalo disruption or as a component in models of galaxy mass.
    * ``"vpeak"`` (*numpy.float32*) - :math:`V_{\rm peak}`, the largest :math:`V_{\rm max}` that the subhalo ever had. This is useful in the same places that :math:`M_{\rm peak}` is.
    * ``"merger_snap"`` (*numpy.int32*) - The snapshot where the subhalo first fell within the virial radius of the host halo.
    * ``"merger_ratio"`` (*numpy.float32*) - The ratio of masses at this infall snapshot.
    * :data:`symlib.HISTORY_DTYPE` also contains all the fields in :func:`symlib.BRANCH_DTYPE`. Note, however, that subhalos where ``is_disappear`` is True or ``is_real`` is False have already been removed, so there is no need to make cuts on this.

.. data:: symlib.BRANCH_DTYPE

    Information about the main branch of a subhalo in the full consistent-trees merger tree. You probably will not need this unless you walk through the full consistent-trees merger tree, which is an advanced action. You can get it by calling :func:`symlib.read_branches`.
	
	**DATA FIELDS**
	
    * ``"start"`` (*numpy.int32*) - The index of the last halo along this main branch. It is labeled "start" because the tree is ordered from later times to earlier times. See the documentation on :func:`read_tree` for more details on tree structure.
    * ``"end"`` (*numpy.int32*) - The index after the first halo in the branch. This means that the full main branch can be accessed by using index slicing: ``branch = tree[start: end]``.
    * ``"is_real"`` (*bool*) - False if the first tracked halo of this branch is a subhalo and True otherwise. Branches where this is False are virtually always tree-linking errors.
    * ``"is_disappear"`` (*bool*) True if the last tracked halo of this branch disrupts without merging with any other halos and True otherwise. Branches where this is True are virtually always barely-resolved object fluctuating in-and-out of existence near the resolution barrier.
    * ``"is_main_sub"`` (*bool*) - True if any halo in the branch was every a subhalo of the main host.
    * ``"preprocess"`` (*numpy.in32*) - A non-negative integer if the branch was ever the subhalo of a larger halo prior to becoming a subhalo of the host and -1 otherwise. If the first case is true, this variable is the index of the largest branch that this branch was a subhalo of. There's some non-trivial bookkeeping required to deal with tree errors caused by major mergers, which will be described in a future paper. For now, suffice to say that it is a generalized version of Section 2.3.1 of Mansfiled & Kravtsov (2020).
    * ``"fist_infall_snap"`` (*numpy.int32*) - If ``"preprocess"`` is non-negative, the snapshot when this branch first fell into a halo of the branch pointed to by ``"preprocess"``.

		  
General Functions
-----------------					

.. function:: symlib.n_hosts(suite_name)

    Returns the number of hosts in a simulation suite. Can be used with :func:`symlib.get_host_directory` to loop over all the host halos in a suite.

    :param str suite_name: The name of the simulation suite.
    :rtype: int

.. function:: symlib.get_host_directory(base_dir, suite_name, halo_name)

    Returns the name of a simulation directory given the base directory that all the suites are stored in, the suite, and the halo name. The halo name can either be the literal halo name (e.g., ``"Halo023"``) or a number in the range [0, *N_hosts*). This can be combined wiht :func:`symlib.n_hosts` to loop over all the hosts in the suite.

    :param str base_dir: Base directory containing all suites.
    :param str suite_name: Name of the simulation suite.
    :param halo_name: Name or index of the halo.
    :type halo_name: str or int
    :rtype: str, the name of the host's simulation directory.
    
.. function:: symlib.scale_factors(sim_dir)

    Returns an array of the scale factors, :math:`a(z)`, of each of snapshot. Sorted from earliest to latest.

    The scale factor arrays of two simulations in different suites may be very different from one another. The scale factor arrays of two simulations in the same suite may be slightly different from one another, depending on whether simulations needed to be restarted midway through.

	:param str sim_dir: The directory of the target host halo.

.. function:: symlib.simulation_parameters(dim_dir)

    Returns a dictionary containing parameters of the simulation suite. These parameters are returned as a dictionary which maps the string names of variables to their values.

    * ``"eps"`` - :math:`\epsilon`, the effective radius of dark matter particles in comoving :math:`h^{-1}{\rm kpc}` (i.e. the "Plummer-equivalent force softening scale").
    * ``"mp"`` - :math:`m_p`, the mass of dark matter particles in :math:`h^{-1}M_\odot`.
    * ``"n_snap"`` - :math:`N_{\rm snap}`, the number of snapshots in the simulation.
    * ``"h100"`` - :math:`h_{100} = H_0 / (100\ {\rm km/s/Mpc})`, the scaled Hubble parameter.

    It also contains `colossus <https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html>`_-compatible cosmology parameters. These are not the same between all suites.
	
    * ``"flat"`` - True if the universe is flat and False otherwise.
    * ``"H0"`` - :math:`H_0`, the Hubble constant in units of km/s/Mpc.
    * ``"Om0"`` - :math:`\Omega_{m,0}`, the total matter density relative to the citical density at :math:`z=0`.
    * ``"Ob0"`` - :math:`\Omega_{m,0}` baryon density relative to the critical density at :math:`z=0`.
    * ``"sigma8"`` - :math:`\sigma_8` the amplitude of the power spectrum at :math:`8\ h^{-1}{\rm Mpc}`.
    * ``"ns"`` - :math:`n_s`, the spectral tilt of the power spectrum.
    
    :param sim_dir: The directory of the target host halo. You may also just pass it the name of the simulation suite (e.g. ``"SymphonyMilkyWay"``)
    :rtype: dict
	

.. function:: symlib.set_units_parameters(scale, param)

.. function:: symlib.set_units_halos(h, scale, param)


Halo Functions
--------------
				  
.. function:: symlib.read_subhalos(params, sim_dir)

    Reads the subhalo data for a single host halo. Two arrays are returned.

    The first return value is a 2D :data:`symlib.SUBHALO_DTYPE` array representing the time-dependent behavior of each subhalo (e.g. postions). The array first indexes over subhaloes in order of their peak :math:`M_{\rm vir}` value and then indexes over snapshots from first to last. The host halo is at the first index. The second argument is a 1D :data:`symlib.SUBHALO_DTYPE` array which represents time-independent information about each subhalo (e.g. merger time). It has the same ordering as the first index of the :data:`symlib.SUBHALO_DTYPE` array.
	
    Subhalos are determined by the Rockstar halo finder and consistent-trees merger tree code. All objects which have ever been within :math:`R_{\rm vir,host}` of the host halo are included, meaning that disrupted, merged, and "splashback" subhalos are included.

    The output arrays use Rockstar's unit conventions by default: all masses, positions, and distances have :math:`h_{100}`-scalings: masses have units of :math:`h^{-1}M_\odot`, positions comoving :math:`h^{-1}{\rm Mpc}`, and radii comoving :math:`h^{-1}{\rm kpc}`. Positions are centered on the zero-point of the box. Almost all users will want to call ``symlib.set_units_halos`` on the returned :data:`symlib.SUBHALO_DTYPE` array to convert to more convenient conventions.
	
    :param dict params: Simulation parameters, as returned by :func:`symlib.simulation_parameters`
    :param str sim_dir: The directory of the target host halo.
    :rtype: (``h``, ``hist``): ``h`` is a ``symlib.SUBHALO_DTYPE`` ``np.array`` with shape (:math:`N_{\rm subhalos}`, :math:`N_{\rm snaps}`), ``hist`` is is a ``symlib.HISTORY_DTYPE`` ``np.array`` with length :math:`N_{\rm subhalos}`.
		
.. function:: symlib.read_tree(sim_dir)

.. function:: symlib.read_branches(sim_dir)


Particle Functions
------------------

.. note::
   Coming with a future paper release
				  
Star Tracking
-------------

.. note::
   Coming with a future paper release


Halo Core Tracking
------------------

.. note::
   Coming with a future paper release


Utility Functions
-----------------

.. function:: symlib.colossus_parameters(param)
				  
.. function:: symlib.suite_names()

.. function:: symlib.plot_circle(ax, x, y, r, **kwargs)
