Data Access
===========

Available data
--------------

As discussed in the :doc:`Simulations <simulations>` page, Symphony is organized into several suites of zoom-in simulations that respectively correspond to host halos containing the following solar masses (:math:`M_\odot`): the LMC (:math:`10^{11}`), Milky Way (:math:`10^{12}`), Group (:math:`10^{13}`), L-Cluster (:math:`4\times 10^{14}`), and Cluster (:math:`10^{15}`) suites.

For each halo in these suites, Symphony contains information on the evolution of halos and particles over time. Because of the size of these data sets, a limited number of them are currently available for download.

**Halos:** This dataset contains information on the subhalos of the simulated host halos, including the evolution of those subhalos over time. Disrupted subhalos and "splashback" subhalos are also included. The subhalos, their properties, and their evolution over time were calculated with a combination of Rockstar, consistent-trees, and custom routines designed to remove numerical artifacts. A custom format is used for these halos, as descibed in the :doc:`Data Access Tutorial <getting_started>` and the :doc:`Data Access Reference <symlib_documentation>`. *This dataset is light-weight enough to analyze on a personal computer, is relatively easy to use, and will be sufficient for the majority of halo-based research tasks.*

**Trees:** This dataset contains the full "merger trees" for all the high-resolution objects in the box, not just the subhalos near the central host. The trees also have the full merger history of each object. These trees were generated with Rockstar and consistent-trees, but have been converted to a more efficient and easy-to-use format and heavily annotated with several custom routines, as described in the :doc:`Data Access Tutorial <getting_started>` and the :doc:`Data Access Reference <symlib_documentation>`. *This dataset is large enough that it will be difficult to analyze on a personal computer, and it will be unnecessary for most halo-based analysis tasks, even those that require looking at the evolution of objects over time.*

The table below shows the sizes of the different different datasets across the different simulation suites.

.. list-table::
	:header-rows: 1
		
	* - Suite Name
	  - :math:`N_{\rm hosts}`
	  - subhalos
	  - trees
	* - SymphonyLMC
	  - 39
	  - 0.16 GB
	  - 18 GB
	* - SymphonyMilkyWay
	  - 45
	  - 0.77 GB
	  - 62 GB
	* - SymphonyGroup
	  - 49
	  - 0.83 GB
	  - 227 GB
	* - SymphonyLCluster
	  - 33
	  - 0.18 GB
	  - 21 GB
	* - SymphonyCluster
	  - 96
	  - 1.1 GB
	  - 160 GB

Downloading Data
----------------

You can download this data with Symphony's data analysis library, ``symlib``. A full reference for this library can be found :doc:`here <symlib_documentation>`.

``symlib`` can be installed with pip. It depends on the standard scientific python libraries, `numpy <https://numpy.org/install/>`__, `scipy <https://scipy.org/install/>`__, and `matplotlib <https://matplotlib.org/stable/users/installing/index.html>`__. It also depends on the cosmology library `colossus <https://bdiemer.bitbucket.io/colossus/installation.html>`__. ``symlib`` will automatically install these libraries if they are not already on your machine, but if you encounter issues, you can find installation instructions for all four libraries on their respective web pages.

If you don't have pip installed, follow the instructions `here <https://pip.pypa.io/en/stable/installation/>`__. Then, from the command line, run the following command.

.. code-block:: console

	pip install symlib -U

You can use the :func:`symlib.download_files` function to download whichever datasets you want to use. The entire dataset will be stored in a single base directory. Each suite has is its own sub-directory within the base, and each zoom-in simulation has a subdirectory within its suite. 

The following Python code shows examples of how to use this function.

.. code-block:: python

	import symlib

	# The base directory where data will be downloaded to.
	# All the suites and halos will be properly ordered
	# within this, so each halo's location will end up being
	# {data_dir}/{suite_name}/{halo_name}.
	data_dir = "path/to/storage/location"

	# The type of data you want to download. "halos" is
	# the basic halo information associated with the
	# central host, and "trees" is that plus the full
	# merger tree of the simulation.
	target = "halos"

	# Exmaple 1
	# Download the first host halo in the Milky Way-mass suite.
	symlib.download_files("SymphonyMilkyWay", 0,
		data_dir, target=target)

	# Example 2
	# Download all the host halos in the Milky Way-mass suite.
	symlib.download_files("SymphonyMilkyWay", None,
		data_dir, target=target)

	# Example 3
	# Download all the host halos across all the suites.
	symlib.download_files(None, None,
		data_dir, target=target)

	# Example 4
	# Download a specific halo that you know the name of.
	symlib.download_files("SymphonyMilkyWay", "Halo023",
		data_dir, target=target)

You can also get a list of suite names with :func:`symlib.suite_names()` and host counts for a given suite with :func:`symlib.n_hosts()` so you can use a fine-tuned for loop instead of ``None``.

If you are running tests on a machine where you don't have much storage space, the smallest host is Halo933 in SymphonyLMC, with a ``"halos"`` size of 2.3 MB and ``"trees"`` size of 227 MB.

symlib also offers the functions :func:`symlib.download_packed_files()` and :func:`symlib.unpack_files()`, which might be helpful if you are running a long download request that gets interrupted midway through.
