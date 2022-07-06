Symlib Documentation
====================
			 
			 
Datatypes
---------

``symlib`` halo data is generally returned as a numpy `structured array <https://numpy.org/doc/stable/user/basics.rec.html>`_, which allows for fields to be accessesed and subselected easily. See the :doc:`Getting Started <getting_started>` for usage examples.

**SUBHALO_DTYPE** contains information on the instantaneous properties of a subhalo over time according to the Rockstar halo finder and is returned by ``read_subhalos()``.

**HISTORY_DTYPE** contains information about the full history of a halo across the entire simulation.

**BRANCH_DTYPE**


	
General Functions
-----------------					

symlib.n_hosts

symlib.get_host_directory
    

symlib.scale_factors

symlib.simulation_parameters

symlib.set_units_parameters

symlib.set_units_halos


Halo Functions
--------------
				  
.. function:: symlib.read_subhalos(params, sim_dir)

    Reads the subhalo data for a single host halo. Two arrays are returned.

    The first return value is a 2D ``symlib.SUBHALO_DTYPE`` array representing the time-dependent behavior of each subhalo (e.g. postions). The array first indexes over subhaloes in order of their peak *M_vir* value and then indexes over snapshots from first to last. The host halo is at the first index.

    The second argument is a 1D ``symlib.HISTORY_DTYPE`` array which represents time-independent information about each subhalo (e.g. merger time). It has the same ordering as the first index of the ``symlib.SUBHALO_DTYPE`` array.

	Subhalos are determined by the Rockstar halo finder and consistent-trees merger tree code. All objects which have ever been within *Rvir* of the host halo are included, meaning that disrupted, merged, and "splashback" subhalos are included.

	:param params: Simulation parameters, as returned by ``symlib.simulation_parameters``
    :type params: dict
    :type sim_dir: The directory of the target host halo.
    :type sim_dir: str
    :rtype: (``h``, ``hist``): ``h`` is a ``symlib.SUBHALO_DTYPE`` ``np.array`` with shape (``n_subhalos``, ``n_snaps``), ``hist`` is is a ``symlib.HISTORY_DTYPE`` ``np.array`` with length ``n_subhalos``.

symlib.read_branches

symlib.read_tree


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

symlib.colossus_parameters
				  
symlib.suite_names

symlib.plot_circle
