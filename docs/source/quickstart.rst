Quickstart Instructions
=======================

This page is designed for experienced users who have worked with halo data before and don't mind a bit of jargon. It also contains FAQs outlining several convenience functions. You can find a more detailed explanation in the :doc:`Data Access <data_access>` page and in the the full tutorials on data analysis, the :doc:`Subalo Tutorial <working_with_subhalos>`, :doc:`Particle Tutorial <working_with_particles>`, and :doc:`Merger Tree Tutorial <working_with_trees>`.

First, install/update the Python 3 analysis library, ``symlib``, using pip:

.. code-block:: console

        pip install symlib -U

Next, fill out `this form <https://docs.google.com/forms/d/e/1FAIpQLSdud6b4i51AP13glVibkzyLAtT9b2ctVx516_hvy5nm76uq1Q/viewform?usp=sf_link>`__ to get emailed a user name and password. To download all the Milky Way mass hosts, call the following within Python:

.. code-block:: python

	import symlib
	symlib.download_files(
	    "my_user_name", "my_password",
	    "SymphonyMilkyWay", None, "my/local/directory"
	)

The files will be downloaded and unpacked in ``my/local/directory/SymphonyMilkyWay/HaloXXX``, with a different ``HaloXXX`` directory for each zoom-in host. To load a host's data into Python, run

.. code-block:: python

	r, hist = symlib.read_rockstar("my/local/directory/SymphonyMilkyWay/HaloXXX")
	s, hist = symlib.read_symfind("my/local/directory/SymphonyMilkyWay/HaloXXX")

The first line reads in subhalo data accoridng to the Rockstar halo finder and the second according to the Symfind halo finder. Both functions return a pair of similarly formatted arrays: ``r`` and ``s`` are 2D arrays which track the evolution of subhalo properties along the main branches of each of the host's subhaloes, including disrupted and splashback subhaloes. ``hist`` contains useful annotations on each main branch (e.g. the snapshot of first infall, :math:`V_{\rm peak}`, etc.) Some subhalos that are clearly artifacts have been removed.
	
``r`` and ``s`` are 2D structured numpy arrays which means that multiple values are stored at each element and can be accessed withstrings. The first index goes over subhaloes and the second over snapshots. The first object is the host halo and the remaining are subhlaos sorted by decreasing :math:`M_{\rm peak}`, down to :math:`N_{\rm peak} > 300`. For example, ``r["mvir"][3,200]`` gives the mass of of the third largest subhalo during snapshot 200. ``hist`` is a 1D structured array, so ``hist["mpeak"][3]`` gives :math:`M_{\rm peak}` for the same subhalo.

A full list of the values in ``h`` can be found :data:`here <symlib.SUBHALO_DTYPE>`, and a similar list for ``hist`` can be found by combining :data:`this list <symlib.HISTORY_DTYPE>` and :data:`this list <symlib.BRANCH_DTYPE>`. Units can be found :ref:`here <units_ref>`.

One ``h`` variable requires special note, ``h["ok"]``, which is true when a subhalo exists and false when it doesn't exist. For example, if you want to analyze all the subhalos in snapshot 200, you should only analyze halos where ``h["ok"][:,200]`` is true. Symfind does not track subhalos prior to infall.

Quickstart FAQs
---------------

.. dropdown:: Where can I find example code?
	:animate: fade-in
			  
			  
	The tutorial has example code that plots :ref:`subhalo locations <halo_position_example>`, looks at the :ref:`mass accretion history <mah_example>` of subhalos, and constructs a :ref:`subhalo mass function <shmf_example>` from all the hosts in a suite.

.. dropdown:: How do I loop over all the hosts in a suite?
	:animate: fade-in
			  
	.. code-block:: python

		symlib.get_host_directory("my/base/directory", "SymphonyMilkyWay", 3)

	returns the directory of host 3 in the Milky Way-mass suite. Use a for loop ranging from 0 to ``symlib.n_hosts("SymphonyMilkyWay")`` to access all the directories.

.. dropdown:: How do I get scale factors?
	:animate: fade-in
		  
	.. code-block:: python

		symlib.scale_factors("SymphonyMilkyWay")

.. dropdown:: How do I get simulation parameters?
	:animate: fade-in

	.. code-block:: python

		param = symlib.simulation_parameters("path/to/HaloXXX")

	``params`` is a dictionary with various cosmological and numerical parameters

	.. code-block:: python

		{'flat': True, 'H0': 70.0, 'Om0': 0.286, 'Ob0': 0.049,
		 'sigma8': 0.82, 'ns': 0.95, 'eps': 0.17, 'mp': 281981.0,
		 'h100': 0.7}

	Note that ``eps`` is in comoving :math:`h^{-1}\,{\rm kpc}` and ``mp`` is in :math:`h^{-1}M_\odot`.

.. dropdown:: How do I get halo properties in comoving units?
    :animate: fade-in
			  
	.. code-block:: python

		h, hist = symlib.read_subhalos("path/to/HaloXXX", comoving=True)

.. dropdown:: How do I get halos/properties not included in the "halos" dataset?
    :animate: fade-in
			  
	The default "halos" dataset (i.e. the data read in by :func:`symlib.read_subhalos`) contains the main branches of every object that has ever been a subhalo of the host as long as the three following conditions are met:

	- :math:`N_{\rm peak} > 300`, where :math:`N_{\rm peak}` is measured prior to the subhalo's first infall. First infall includes halos other than the host and does not include temporary Rockstar errors caused by major mergers.
	- The halo is not a subhalo during its first snapshot.
	- If the halo disrupts, consistent-trees merges it with any other halo.

	If you want other objects, you will need to analyze the full merger tree. This must be :doc:`downloaded separately <data_access>`. Symphony's merger trees use a different format than consistent-trees, so it would be best to read through the :doc:`full tutorial <working_with_trees>`. The full merger tree also contains `additional variables <merger_tree_variables>` not included in the standard halo dataset.
