Symlib Documentation
====================
			 
			 
Datatypes
---------

.. py:data:: symlib.SUBHALO_DTYPE = [("id", "i4"), ("mvir", "f4"), ("vmax", "f4"), ("rvmax", "f4"), ("x", "f4", (3,)), ("v", "f4", (3,)), ("ok", "?"), ("rvir", "f4"), ("cvir", "f4")]

.. py:data:: symlib.HISTORY_DTYPE = [("mpeak", "f4"), ("vpeak", "f4"), ("merger_snap", "i4"), ("merger_ratio", "f4"), ("start", "i4"), ("end", "i4"), ("is_real", "?"), ("is_disappear", "?"), ("is_main_sub", "?"), ("preprocess", "i4"), ("first_infall_snap", "i4")]

.. py:data:: symlib.BRANCH_DTYPE = [("start", "i4"), ("end", "i4"), ("is_real", "?"), ("is_disappear", "?"), ("is_main_sub", "?"), ("preprocess", "i4"), ("first_infall_snap", "i4")]

General Functions
-----------------					

.. autofunction:: symlib.scale_factors

.. autofunction:: symlib.read_subhalos

.. autofunction:: symlib.read_branches

.. autofunction:: symlib.read_tree

.. autofunction:: symlib.read_particles

.. autofunction:: symlib.read_particle_tags
				  
Star Tracking
-------------

Utility Functions
-----------------

.. autodata:: symlib.parameter_table

.. autofunction:: symlib.scale_factors

.. autofunction:: symlib.pristine_merger_indices

.. autofunction:: symlib.propagate_parent_indices

.. autofunction:: symlib.colossus_parameters

.. autofunciton:: symlib.clean_particles
