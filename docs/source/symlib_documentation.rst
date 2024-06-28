``symlib`` Documentation
=========================

This page is meant to be a full technical reference to ``symlib``. First-time users will probably find it easiest to start with the tutorial pages on :doc:`Data Access <data_access>`, :doc:`Working with Subhalos <working_with_subhalos>`, :doc:`Working with Particles <working_with_particles>`, and :doc:`Working with Merger Trees <working_with_trees>` rather than trying to read through this page from start to finish.

.. _units_ref:

Units
-----

``symlib`` functions generally expect variables to be in the following units:

- Masses: :math:`M_\odot`
- Distances/radii: physical :math:`{\rm kpc}`
- Positions: physical :math:`{\rm kpc}`, centered on the host halo when appropriate
- Velocities: physical :math:`{\rm km/s}`, centered on the host halo when appropriate

Much of the data that ``symlib`` reads in is a processed form of data from another code, which may use different conventions. To convert to ``symlib``'s conventions, use the ``set_units_*`` family of functions: :func:`symlib.set_units_halos`, :func:`symlib.set_units_parameters`, :func:`symlib.set_units_branches` as needed.
			 
Datatypes
---------

``symlib`` halo data is generally returned as a numpy `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_, which allows for fields to be accessesed and subselected easily. See the :doc:`Working with Subhalos <working_with_subhalos>` page for usage examples.

.. data:: symlib.SUBHALO_DTYPE
		   
    Time-dependent subhalo information (e.g., position) from the Rockstar halo finder. You can get this information for all a host's subhalos by calling :func:`symlib.read_subhalos`.
	
    **DATA FIELDS:**
	
    * ``"id"`` (*numpy.int32*) - A unique integer identifying each subhalo. Changes from snapshot to snapshot.

    * ``"mvir"``  (*numpy.float32*) - The mass of the halo, :math:`M_{\rm vir}`. When isolated, this is an overdensity mass from the Bryan & Norman (1998) definition of the virial overdensity. When deep in a host halo, this is the bound mass. The transition between these two definitions is ill-defined.

    * ``"rvir"`` (*numpy.float32*) - The overdensity radius of the halo, :math:`R_{\rm vir}`.

    * ``"vmax"`` (*numpy.float32*) - The maximum value of the halo's circular rotation curve, :math:`V_{\rm max} = {\rm max}\left\{V_{\rm rot}(r) = \sqrt{G M(<r)/r}\right\}`.

    * ``"rvmax"`` (*numpy.float32*) - The radius where this maximum rotation velocity occurs.

    * ``"cvir"`` (*numpy.float32*) - An estimate of how centrally concentrated the subhalo's mass is, :math:`c_{\rm vir}=R_s/R_{\rm vir}`. :math:`R_s` is the transition radius between shallow inner density slopes (:math:`d \ln(\rho)/d \ln(r)` > -2) and steep outer slopes (i.e. :math:`d \ln(\rho)/d \ln(r)` < -2). :math:`c_{\rm vir}` is estimated  by measuring :math:`V_{\rm max}/V_{\rm rot}(R_{\rm vir})`, assuming an NFW profile, and solving for :math:`R_s`. Because of this, the *value* of :math:`c_{\rm vir}` is only meaningful for halos where the assumption of NFW profiles is reasonable (non-subhalos). However, the *relative ordering* of concentrations will be correct regardless.

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

    * ``"branch_idx"`` (*numpy.int32*) - The index of this halo's branch in the full merger tree. This allows you to switch back and forther between the two data structures as needed.

    * :data:`symlib.HISTORY_DTYPE` also contains all the fields in :func:`symlib.BRANCH_DTYPE`. Note, however, that subhalos where ``is_disappear`` is True or ``is_real`` is False have already been removed, so there is no need to make cuts on this.

    * ``false_selection`` (*bool*) - True if the branch has :math:`M_{\rm peak} \geq 300\cdot m_p` after infall, but :math:`M_{\rm peak} < 300\cdot m_p`.
    

.. data:: symlib.BRANCH_DTYPE

    Information about the main branch of a subhalo in the full consistent-trees merger tree. You probably will not need this unless you walk through the full merger tree, which is an advanced action. You can get it by calling :func:`symlib.read_branches`.
	
	**DATA FIELDS**
	
    * ``"start"`` (*numpy.int32*) - The index of the last halo along this main branch. It is labeled "start" because the tree is ordered from later times to earlier times. See the documentation on :func:`read_tree` for more details on tree structure.

    * ``"end"`` (*numpy.int32*) - The index after the first halo in the branch. This means that the full main branch can be accessed by using index slicing: ``branch = tree[start: end]``.

    * ``"is_real"`` (*bool*) - False if the first tracked halo of this branch is a subhalo and True otherwise. Branches where this is False are virtually always tree-linking errors.

    * ``"is_disappear"`` (*bool*) True if the last tracked halo of this branch disrupts without merging with any other halos and False otherwise. Branches where this is True are virtually always barely-resolved object fluctuating in-and-out of existence near the resolution barrier.

    * ``"is_main_sub"`` (*bool*) - True if any halo in the branch was ever a subhalo of the main host.

    * ``"preprocess"`` (*numpy.int32*) - A non-negative integer if the branch was ever the subhalo of a larger halo prior to becoming a subhalo of the host and -1 otherwise. If the first case is true, this variable is the index of the largest branch that this branch was a subhalo of. There's some non-trivial bookkeeping required to deal with tree errors caused by major mergers, which will be described in a future paper. For now, suffice to say that it is a generalized version of Section 2.3.1 of Mansfield & Kravtsov (2020).

    * ``"first_infall_snap"`` (*numpy.int32*) - If ``"preprocess"`` is non-negative, the snapshot when this branch first fell into a halo of the branch pointed to by ``"preprocess"``.

.. data:: symlib.PARTICLE_DTYPE
    
    Information about the particles associated with a given subhalo. 

    **DATA FIELDS:**
    
    * ``"id"`` (*numpy.int32*) - A unique integer identifying each particle. Changes from snapshot to snapshot.

    * ``"x"`` (*numpy.float32*) - :math:`x`, the position of the particle.

    * ``"v"`` (*numpy.float32*) - :math:`v`, the velocity of the particle.

    * ``"snap"`` (*numpy.int32*) - The particle's snapshot.

    * ``"ok"`` (*bool*) - True if the particle exists during the specified snapshot and False otherwise.

    * ``"smooth"`` (*bool*) - True if the particle smoothly accreted onto the subhalo and False otherwise.
      

.. _merger_tree_variables:
      
Merger Tree Variables
---------------------

The following variables can be read in from merger trees with the :func:`symlib.read_tree` function. These variables are taken directly from the consistent-trees output files and still retain its units and ID conventions.

* ``"dfid"`` - The depth-first ID of the halo.

* ``"id"`` - The ID of the halo.

* ``"desc_id`` - The ID (``id``, not ``dfid``) of the descendant. -1 if the halo has no descendants.

* ``"upid"`` - The UpID of a halo. This is -1 if the halo is not within a larger halo's virial radius, otherwise it is the ID (``id``, not ``dfid``) of that larger halo.

* ``"phantom"`` - A flag indicating whether consistent-trees was able to track the object during this snapshot. 1 if so, and 0 otherwise. If 0, this halo's properties were interpolated during this snapshot.

* ``"snap"`` -  This halo's snapshot.

* ``"next_co_prog"`` - The depth-first ID (``dfid``, not ``id``) of this halo's co-progenitor, if it exists. If this halo doesn't have a co-progenitor, this variable is -1. See :doc:`Intro to Merger Trees <intro_to_merger_trees>` for a description of what this is.

* ``"mvir"`` -  The mass of the halo, :math:`M_{\rm vir}`. When isolated, this an overdensity mass from the Bryan & Norman (1998) definition of the virial overdensity. When deep in a host halo, this is the bound mass. The transition between these two definitions is ill-defined.

* ``"rs"`` - The NFW scale radius of the halo, :math:`R_s`. Units are comoving :math:`h^{-1}{\rm kpc}`.

* ``"vmax"`` -  The maximum value of the halo's circular rotation curve, :math:`V_{\rm max} = {\rm max}\left\{V_{\rm rot}(r) = \sqrt{G M(<r)/r}\right\}`. Units are physical km/s.

* ``"m200b"`` - The overdensity mass, :math:`M_{\rm 200b}`, corresponding to :math:`200\times \rho_m`.

* ``"m200c"`` - The overdensity mass, :math:`M_{\rm 200b}`, corresponding to :math:`200\times \rho_c`.

* ``"m500c"`` - The overdensity mass, :math:`M_{\rm 200b}`, corresponding to :math:`500\times \rho_c`.

* ``"xoff"`` - The distance between the center of mass and the densest part fo the halo. units are comoving :math:`h^{-1}{\rm kpc}`.

* ``"spin_bullock"`` - Unitless paramater that tracks the specific anular momentum of the halo. :math:`|\vec{J}|/(\sqrt{2}\,M_{\rm vir}\,V_{\rm vir}\,R_{\rm vir})`.

* ``"c_to_a"`` - The unitless minor-to-major axis ratio of the halo.

* ``"b_to_a"`` - The unitless intermediate-to-major axis ratio of the halo.

* ``"t_to_u"`` - The virial ratio, :math:`T/|U|`.

* ``"r_vmax"`` - The radius, :math:`R_{\rm vmax}`, at which :math:`V_{\rm max}` occurs.

* ``"x"`` - A 3-vector, :math:`\vec{x}` giving the position of the halo in comivng :math:`h^{-1}{\rm Mpc}`.

* ``"v"`` - A 3-vector, :math:`\vec{v}`, giving the velocity of the halo in physical km/s.

* ``"j"`` - A 3-vector, :math:`\vec{J}`, giving the angular momentum of the halo in physical :math:`h^{-2}M_\odot\cdot{\rm Mpc}\cdot{\rm km/s}`.

* ``"a"`` - A 3-vector, :math:`\vec{A}`, pointing in the direction of the halo's major axis with length equal to that major axis. Units are comoving :math:`h^{-1}{\rm kpc}`.

  
General Functions
-----------------					

.. function:: symlib.n_hosts(suite_name)

    Returns the number of zoom-in simulations, each of which is associated with one "target" host halo, in a simulation suite (the Symphony suites are: LMC, Milky Way, Group, L-Cluster, and Cluster). Can be used with :func:`symlib.get_host_directory` to loop over all target host halos in a suite.

    :param str suite_name: The name of the simulation suite.
    :rtype: int


.. function:: symlib.get_host_directory(base_dir, suite_name, halo_name)

    Returns the name of a simulation directory given the base directory that all the suites are stored in, the suite, and the halo name. The halo name can either be the literal halo name (e.g., ``"Halo023"``) or a number in the range :math:`[0,\,N_{\rm host})`. This can be combined with :func:`symlib.n_hosts` to loop over all the hosts in a suite.

    :param str base_dir: Base directory containing all suites.
    :param str suite_name: Name of the simulation suite.
    :param halo_name: Name or index of the target host halo.
    :type halo_name: str or int
    :rtype: str, the name of the host's simulation directory.
    

.. function:: symlib.scale_factors(sim_dir)

    Returns an array of the scale factors, :math:`a(z)`, of each of snapshot. Sorted from earliest to latest.

    The scale factor arrays of two simulations in different suites may be different from one another. The scale factor arrays of two simulations in the same suite sometimes also slightly differ, depending on whether simulations needed to be restarted midway through.

    :param str sim_dir: The directory of the target host halo.
    :rtype: ``np.array`` containing the scale factors of each snapshot in the simulation.


.. function:: symlib.simulation_parameters(dim_dir)

    Returns a dictionary containing parameters of the simulation suite. These parameters are returned as a dictionary which maps the string names of variables to their values.

    * ``"eps"`` - :math:`\epsilon`, the effective radius of dark matter particles in comoving :math:`h^{-1}{\rm kpc}` (i.e. the "Plummer-equivalent force softening scale").

    * ``"mp"`` - :math:`m_p`, the mass of dark matter particles in :math:`h^{-1}M_\odot`.

    * ``"n_snap"`` - :math:`N_{\rm snap}`, the number of snapshots in the simulation.

    * ``"h100"`` - :math:`h_{100} = H_0 / (100\ {\rm km/s/Mpc})`, the scaled Hubble parameter.

    It also contains `colossus <https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html>`_-compatible cosmology parameters. Note that these are not the same between all suites.
	
    * ``"flat"`` - True if the universe is flat and False otherwise.

    * ``"H0"`` - :math:`H_0`, the Hubble constant in units of km/s/Mpc.

    * ``"Om0"`` - :math:`\Omega_{m,0}`, the total matter density relative to the citical density at :math:`z=0`.

    * ``"Ob0"`` - :math:`\Omega_{m,0}` baryon density relative to the critical density at :math:`z=0`.

    * ``"sigma8"`` - :math:`\sigma_8` the amplitude of the power spectrum at :math:`8\ h^{-1}{\rm Mpc}`.
    
    * ``"ns"`` - :math:`n_s`, the spectral tilt of the power spectrum.
    
    :param sim_dir: The directory of the target host halo. You may also just pass it the name of the simulation suite (e.g. ``"SymphonyMilkyWay"``)
    :rtype: dict
	

.. function:: symlib.set_units_parameters(scale, param)
	      
   Converts the particle mass (:math:`m_p`, ``"mp"``) and particle size (:math:`\epsilon`, ``"eps"``) to ``symlibs``'s default units.

   :param float mp: particle mass in :math:`M_\odot`
   :param float eps: Plummer-equivalent force softening scale in physical :math:`{\rm kpc}`.
	  
	      
.. function:: symlib.set_units_halos(h, scale, param)
	      
   Converts the units of a 2D ``np.array`` with type :data:`symlib.SUBHALO_DTYPE` to ``symlib``'s default units. All masses will be in units of :math:`M_\odot`, all positions and radii will be units of physical :math:`{\rm kpc}`. Positions will be centered on the first halo in the array at the given snapshot. Velocities will be in physical :math:`{\rm km/s}` and similarly centered on the velocity of the first halo at each snapshot.

   This function only needs to be called if `comoving=True` in :func:`symlib.read_subhalos`. This is not true by default
   
   :param symlib.SUBHALO_DTYPE np.array h: A 2D array of subhalos, with the first halo indexing over halos and the second over snapshots. (see :func:`symlib.read_subhalos`).
   :param np.array scale: An array of the scale factors of each snapshot (see :func:`symlib.scale_factors`)
   :param dict param: The simulation parameters (see :func:`symlib.simulation_parameters`)


.. function:: symlib.set_units_histories(hist, scale, param)
	      
   Converts the units of an ``np.array`` with type :data:`symlib.HISTORY_DTYPE` to ``symlib``'s default units. All masses will be in units of :math:`M_\odot`, all positions and radii will be units of physical :math:`{\rm kpc}`. Positions will be centered on the first halo in the array at the given snapshot. Velocities will be in physical :math:`{\rm km/s}` and similarly centered on the velocity of the first halo at each snapshot.

   This function only needs to be called if `comoving=True` in :func:`symlib.read_subhalos`. This is true by default.
   
   :param symlib.HISTORY_DTYPE np.array h: Array of subhalo histories (see :func:`symlib.read_subhalos`).
   :param np.array scale: An array of the scale factors of each snapshot (see :func:`symlib.scale_factors`)
   :param dict param: The simulation parameters (see :func:`symlib.simulation_parameters`)


Halo Functions
--------------
				  
.. function:: symlib.read_subhalos(sim_dir, comoving=False, include_false_selections=False)

    Reads the subhalo data for a single host halo. Two arrays are returned.

    The first return value is a 2D :data:`symlib.SUBHALO_DTYPE` array representing the time-dependent behavior of each subhalo (e.g. positions). The array first indexes over subhaloes in order of their peak :math:`M_{\rm vir}` value and then indexes over snapshots from first to last. The host halo is at the first index. The second argument is a 1D :data:`symlib.SUBHALO_DTYPE` array which represents time-independent information about each subhalo (e.g. merger time). It has the same ordering as the first index of the :data:`symlib.SUBHALO_DTYPE` array.
	
    Subhalos are determined by the Rockstar halo finder and consistent-trees merger tree code. All objects that have ever been within :math:`R_{\rm vir,host}` of the host halo are included, meaning that disrupted, merged, and "splashback" subhalos are included.

    If ``comoving=False``, ``symlib``'s default units are used. Positions and velocities are centered on the host halo. Otherwise, the output arrays use Rockstar's unit conventions by default: all masses, positions, and distances have :math:`h_{100}`-scalings: masses have units of :math:`h^{-1}M_\odot`, positions comoving :math:`h^{-1}{\rm Mpc}`, and radii comoving :math:`h^{-1}{\rm kpc}`. In this case positions will be centered on the zero-point of the box.

    By default, subhalos which have :math:`M_{\rm peak}` above the 300-particle cutoff, but were below the cutoff when they first became a subhalo are considered numerical artifacts and are _not_ included. They can be reintroduced to the catalog by setting ``include_false_selections=True``
    
    :param str sim_dir: The directory of the target host halo.
    :param bool comoving=False: Controls whether the resturn values are in default Rockstar/consistent-trees units (``False``) or default symlib units (``True``).
    :param bool include_false_selections=False: Controls whether subhalos which only have :math:`M_{\rm peak}` above the catalog cutoff due toa consistent-trees error are included (``True``) or excluded (``False``).
    :rtype: (``h``, ``hist``): ``h`` is a :data:`symlib.SUBHALO_DTYPE` ``np.array`` with shape (:math:`N_{\rm subhalos}`, :math:`N_{\rm snaps}`), ``hist`` is is a :data:`symlib.HISTORY_DTYPE` ``np.array`` with length :math:`N_{\rm subhalos}`.

	
.. function:: symlib.read_tree(sim_dir, var_names)

   Reads the time-dependent properties of every halo in the simulation, not just the subhalos of the target host in a "depth-first merger tree" format.

   The user supplies a list of variable names and a single, 1D array is returned for each variable. Each element of each array is a halo at a specific snapshot, and these arrays are ordered in a way that encodes which halos evolve and merge into which other halos. To decode this structure, you will need to use the results of :func:`symlib.read_branches`, which breaks the tree into smaller structures, or "branches."

   The full strucutre of this merger tree is too large of a topic to be covered here. A writeup can be found on the :doc:`Intro to Merger Trees <intro_to_merger_trees>` page.
	      
   :param str sim_dir: The directory of the target host halo.
   :param str list var_names: The names of variables.
   :rtype: tuple of ``np.array``, one for each element in ``var_names``.

	      
.. function:: symlib.read_branches(sim_dir)
	      	      
   Reads information about the time-independent properties of every halo in the simulation, not just the subhalos of target host. Each element corresonds to a single branch in the tree (i.e. the evolution of a single halo over time) and gives information on the properties and location of the branch.

   The full strucutre of this merger tree is too large of a topic to be covered here. A writeup can be found on the :doc:`Intro to Merger Trees <intro_to_merger_trees>` page.
   
   :param str sim_dir: The directory of the target host halo.
   :rtype: :data:`symlib.BRANCH_DTYPE` ``np.array`` 


.. function:: symlib.merger_lookup_table(b, dfid)

   Creates a lookup table to aid with finding the branches of merging halos. The details of this table are not important and may be changed at any time to improve performance.

   :param b:
   :type b: :data:`symlib.BRANCH_DTYPE` np.array
   :param int np.array dfid:
   :rtype: int np.array

   
.. function:: symlib.find_merger_branch(lookup_table, co_prog)

   Searches for the index of the branch corresponding of a given merging subhalo. The subhalo is identified by a "co-progenitor" ID. See the writeup in :doc:`Intro to Merger Trees <intro_to_merger_trees>` for more discussion on what this means.

   In practice, most users will want to use :func:`symlib.find_all_merger_branches`.

   :param int np.array lookup_table: A look up table, as created by :func:`symlib.merger_lookup_table`.
   :param int co_prog: a single "co-progenitor depth-first ID" (``"next_co_prog"`` in calls to :func:`read_tree`).
   :rtype: int

		       
.. function:: symlib.find_all_merger_branches(b, lookup_table, co_prog, i)

   Returns the indices of all the branches that merge with a given halo. (i.e. branches that exist in the current snapshot but disrupt in the next snapshot).

   :param b: The branch information for the merger tree.
   :type b: :data:`symlib.BRANCH_DTYPE` np.array
   :param int np.array lookup_table: A look up table, as created by :func:`symlib.merger_lookup_table`.
   :param int np.array co_prog: A tree-ordered array of co-progenitor IDs (``"next_co_prog"`` in calls to :func:`read_tree`).
   :param int i: The index of the halo in the tree that you are interested in. 
   :rtype: int np.array


Particle Functions
------------------

.. class:: symlib.Particles(sim_dir)
    
    Wraps around the particle data files generated by Symphony for easier user access.

    :param str sim_dir: The directory of the target host halo.

.. function:: read(snap, halo, mode, comoving)
    
    Reads the particle data for all subhalos associated with a host halo. A list containing 1D arrays corresponding to six different particle data variables is returned.

    :param int snap: The desired snapshot.
    :param int halo: The index for the desired subhalo. Defaults to -1 to return all subhalos. 
    :param bool comoving=False: Controls whether the return values are in default Rockstar/consistent-trees units (``False``) or default symlib units (``True``).
    :param str mode: Takes 3 different string values corresponding to the desired output:
        "Current" is the default setting and returns variables for only valid particles.
        "All" returns variables for all particles, whether they are valid or invalid.
        "Smooth" returns variables for only smoothly accreted particles. For advanced users, 
        this is the fastest option for getting particle data.

    :rtype: list np.array
        Returns a list containing the following variables:[``id``, ``x``, ``v``, ``snap``, ``ok``, 
        ``smooth``]. Each variable within the returned list is a 1D ``np.array``.
    

Utility Functions
-----------------

.. function:: symlib.colossus_parameters(param)
	      
   Converts a ``symlib`` parameter dictionary to a parameter dictionary that can be passed to a call to `colossus.cosmology.cosmology.setCosmology <https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html#cosmology.cosmology.setCosmology>`_. This will allow you to calculate cosmological quantities (e.g. the mass-concentration relation) using the colossus library.

   :param dict param: A ``symlib`` parameter dictionary returned by :func:`symlib.simulation_parameters`.
   :rtype: A ``colossus`` parameter dictionary.
	      
				  
.. function:: symlib.suite_names()
	      
   Returns a list of all the valid suite names.

   :rtype: string list 

	      
.. function:: symlib.plot_circle(ax, x, y, r, **kwargs)

   Plots the a circle to a given `matplotlib.pyplot.Axes <https://matplotlib.org/stable/api/axes_api.html#the-axes-class>`_. This is a convenience function that helps with example code in the tutorial.

   All keyword arguments accepted by `matplotlib.pyplot.plot <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html>`_ are accepted as keywords arguments by this function.

   :param matplotlib.pyplot.Axes ax: The axis to plot the circle on.
   :param float x: The :math:`x` coordinate of the circle.
   :param float y: The :math:`y` coordinate of the circle.
   :param float r: The radius of the circle.


File Management
---------------

.. function:: symlib.download_files(user, password, suite, halo_name, base_out_dir, target="halos", logging=True)

   Downloads data associated with a set of halos/suites. See :doc:`Data Access <data_access>` for usage examples.

   This download has two stages. First, all the data is downloaded in "packed" ``tar`` files. Once this finishes, all the ``tar`` files are expanded into data directories and deleted. This first step is handled with :func:`symlib.download_packed_files()` and the second with :func:`symlib.unpack_files()`. If you are running a large download job that stops halfway and don't want to repeat work when you restart it, you can use these two functions do do it.

   :param user: The username you would like to use to perform the download. See instructions on the :doc:`Data Access <data_access>` page for obtaining a username/password.
   :param password: The password associated with your username. See instructions on the :doc:`Data Access <data_access>` page for obtaining a username/password.
   :param suite_name: The suite to download a halo from. This may either be the full name of a symlib suite or None. If None, :func:`symlib.download_files()` will be applied to every simulation suite with the given value of ``halo_name``.
   :type suite_name: str or None
   :param halo_name: The halo to download. This can either be an int giving the index of the halo in the suite, a string giving the name of the halo, or None. If None, all the halos in the given suite[s] will be downloaded.
   :type halo_name: str, int, or None
   :param base_out_dir: The directory where data is stored.
   :param target="halos": What type of data to download. Possible options are ``"halos"`` and ``"trees"``.
   :param logging=True: True if you would like output printed telling the user what stage in the download they are at and False if you would like to turn off as much printing as possible. 

   
.. function:: symlib.download_packed_files(user, password, suite, halo_name, base_out_dir, target="halos", logging=True)
	      
   Downloads "packed" ``tar`` files containing the requested data for a given set of halos/suites. This function represents half of the :func:`symlib.download_files()` command and may be useful to users whose download stops halfway through and would like to restart. Note that in such a case, the _last_ downloaded ``tar`` file is likely an incomplete download and is probably corrupted. It should be repeated.

   :param user: The username you would like to use to perform the download. See instructions on the :doc:`Data Access <data_access>` page for obtaining a username/password.
   :param password: The password associated with your username. See instructions on the :doc:`Data Access <data_access>` page for obtaining a username/password.
   :param suite_name: The suite to download a halo from. This may either be the full name of a symlib suite or None. If None, :func:`symlib.download_packed_files()` will be applied to every simulation suite with the given value of ``halo_name``.
   :type suite_name: str or None
   :param halo_name: The halo to download. This can either be an int giving the index of the halo in the suite, a string giving the name of the halo, or None. If None, all the halos in the given suite[s] will be downloaded.
   :type halo_name: str, int, or None
   :param base_out_dir: The directory where data is stored.
   :param target="halos": What type of data to download. Possible options are ``"halos"`` and ``"trees"``.
   :param logging=True: True if you would like output printed telling the user what stage in the download they are at and False if you would like to turn off as much printing as possible. 

   
.. function:: unpack_files(suite, halo_name, base_out_dir, target="halos", logging=True)
	      
   Opens "packed" ``tar`` files containing the requested data for a given set of halos/suites. This function represents the second half of the :func:`symlib.download_files()` command and may be useful to users whose download stops halfway through and would like to restart. Note that in such a case, the _last_ downloaded ``tar`` file is likely an incomplete download and is probably corrupted. It should be repeated.
	      
   :param suite_name: The suite to download a halo from. This may either be the full name of a symlib suite or None. If None, :func:`symlib.unpack_files()` will be applied to every simulation suite with the given value of ``halo_name``.
   :type suite_name: str or None
   :param halo_name: The halo to download. This can either be an int giving the index of the halo in the suite, a string giving the name of the halo, or None. If None, all the halos in the given suite[s] will be downloaded.
   :type halo_name: str, int, or None
   :param base_out_dir: The directory where data is stored.
   :param target="halos": What type of data to download. Possible options are ``"halos"`` and ``"trees"``.
   :param logging=True: True if you would like output printed telling the user what stage in the download they are at and False if you would like to turn off as much printing as possible. 
