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
    * ``"mvir"``  (*numpy.float32*) - The mass of the halo. When isolated, this an overdensity mass from the Bryan & Norman (1998) definition of the virial overdensity. When deep in a host halo, this is the bound mass. The transition between these two definitions is fuzzy.
    * ``"rvir"`` (*numpy.float32*) - The overdensity radius of the halo.
    * ``"vmax"`` (*numpy.float32*) - The largest value that the halo's rotation curve, *V(<r) = sqrt(G M(<r)/r)* takes on.
    * ``"rvmax"`` (*numpy.float32*) - The radius where this maximum rotation velocity occurs.
    * ``"cvir"`` (*numpy.float32*) - An estimate of how centrally concentrated the subhalo's mass is. Specifcially, this is *Rvir/Rs*, where Rs is the transition radius between shallow inner density slopes (*d ln rho/d ln r > -2*) and steep outer slopes (*d ln rho/d ln r < -2*). *cvir* is estimated  by measuring *Vmax/V(<Rvir)*, assuming an NFW profile, and solving for *Rs*. Because of this, the value of *cvir* is only meaningful for halos where the assumption of NFW profiles is reasonable (non-subhalos). However, the *relative ordering* of concentrations will be correct regardless.
    * ``"x"`` (*numpy.float32*) - The position of the subhalo.
    * ``"v"`` (*numpy.float32*) - The velocity of the subhalo.
    * ``"ok"`` (*bool*) - True if the subhalo exists during the specified snapshot and False otherwise.
		
.. data:: symlib.HISTORY_DTYPE

    Time-independent subhalo information about the subhalo's entire history in the simulation (e.g. when it first fell into the host halo). You can get it for all the host's subhalos by calling :func:`symlib.read_subhalos`.

	**DATA FIELDS**
	
    * ``"mpeak"`` (*numpy.float32*) - The largest mass that the subhalo ever had.
    * ``"vpeak"`` (*numpy.float32*) - The highest Vmax that the subhalo ever had.
    * ``"merger_snap"`` (*numpy.int32*) - The snapshot where the subhalo first fell within the virial radius of the host halo.
    * ``"merger_ratio"`` (*numpy.float32*) - The ratio of masses at this infall snapshot.
	* :data:`symlib.HISTORY_DTYPE` also contains all the fields in :func:`symlib.BRANCH_DTYPE`. Note, however, that subhalos where ``is_disappear`` is True or ``is_real`` is False have already been removed, so there is no need to make cuts on this.

.. data:: symlib.BRANCH_DTYPE

    Information about the main branch of a subhalo in the full consistent-trees merger tree. You probably will not need this unless you walk through the full consistent-trees merger tree, which is an advanced action. You can get it by calling :func:`symlib.read_branches`.
	
	**DATA FIELDS**
	
    * ``"start"`` (*numpy.int32*)
    * ``"end"`` (*numpy.int32*)
    * ``"is_real"`` (*bool*)
    * ``"is_disappear"`` (*bool*)
    * ``"is_main_sub"`` (*bool*)
    * ``"preprocess"`` (*numpy.in32*)
    * ``"fist_infall_snap"`` (*numpy.int32*)

		  
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

    Returns an array of the scale factors, *a(z)*, of each of snapshot. Sorted from earliest to latest.

    The scale factor arrays of two simulations in different suites may be very different from one another. The scale factor arrays of two simulations in the same suite may be slightly different from one another, depending on whether simulations needed to be restarted midway through.

	:param str sim_dir: The directory of the target host halo.

.. function:: symlib.simulation_parameters(dim_dir)

    Returns a dictionary containing parameters of the simulation suite. These parameters are returned as a dictionary which maps the string names of variables to their values.

    * ``"eps"`` - The effective radius of dark matter particles in comoving kpc/h (i.e. the "Plummer-equivalent force softening scale").
    * ``"mp"`` - The mass of dark matter particles in Msun/h.
    * ``"n_snap"`` - The number of snapshots in the simulation.
    * ``"h100"`` - H0 / (100 km/s/Mpc).

    It also contains `colossus <https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html>`_-compatible cosmology parameters
	
    * ``"flat"`` - True if the universe is flat and False otherwise.
    * ``"H0"`` - Hubble constant in units of km/s/Mpc.
    * ``"Om0"`` - Total matter density relative to the citical density at z=0.
    * ``"Ob0"`` - Baryon density relative to the critical density at z=0.
    * ``"sigma8"`` - Amplitude of the power spectrum at 8 Mpc/h.
    * ``"ns"`` - Spectral tilt of the power spectrum.
    
    :param sim_dir: The directory of the target host halo. You may also just pass it the name of the simulation suite (e.g. ``"SymphonyMilkyWay"``)
    :rtype: dict
	

.. function:: symlib.set_units_parameters(scale, param)

.. function:: symlib.set_units_halos(h, scale, param)


Halo Functions
--------------
				  
.. function:: symlib.read_subhalos(params, sim_dir)

    Reads the subhalo data for a single host halo. Two arrays are returned.

    The first return value is a 2D :data:`symlib.SUBHALO_DTYPE` array representing the time-dependent behavior of each subhalo (e.g. postions). The array first indexes over subhaloes in order of their peak *M_vir* value and then indexes over snapshots from first to last. The host halo is at the first index. The second argument is a 1D :data:`symlib.SUBHALO_DTYPE` array which represents time-independent information about each subhalo (e.g. merger time). It has the same ordering as the first index of the :data:`symlib.SUBHALO_DTYPE` array.
	
    Subhalos are determined by the Rockstar halo finder and consistent-trees merger tree code. All objects which have ever been within *R_vir* of the host halo are included, meaning that disrupted, merged, and "splashback" subhalos are included.

    The output arrays use Rockstar's unit conventions by default: all masses, positions, and distances have *h_100*-scalings: masses have units of Msun/h, positions comoving Mpc/h, and radii comoving kpc/h. Positions are centered on the zero-point of the box. Almost all users will want to call ``symlib.set_units_halos`` on the returned :data:`symlib.SUBHALO_DTYPE` array to convert to more convenient conventions.
	
    :param dict params: Simulation parameters, as returned by :func:`symlib.simulation_parameters`
    :param str sim_dir: The directory of the target host halo.
    :rtype: (``h``, ``hist``): ``h`` is a ``symlib.SUBHALO_DTYPE`` ``np.array`` with shape (*N_subhalos*, *N_snaps*), ``hist`` is is a ``symlib.HISTORY_DTYPE`` ``np.array`` with length *N_subhalos*.
		
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
