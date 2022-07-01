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

.. autofunction:: symlib.scale_factors

.. autofunction:: symlib.simulation_parameters

Halo Functions
--------------
				  
.. autofunction:: symlib.read_subhalos

.. autofunction:: symlib.read_branches

.. autofunction:: symlib.read_tree


Particle Functions
------------------
				  
Star Tracking
-------------

Utility Functions
-----------------

.. autofunction:: symlib.set_units_parameters

.. autofunction:: symlib.set_units_halos

.. autofunction:: symlib.colossus_parameters
				  
.. autofunction:: symlib.n_hosts

.. autofunction:: symlib.get_host_directory

.. autofunction:: symlib.suite_names
