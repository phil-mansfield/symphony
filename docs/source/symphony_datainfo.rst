Symphony Data
================

Available data
--------------

As discussed in the :doc:`Simulations <symphony_overview>` page, Symphony is organized into several suites of zoom-in simulations that respectively correspond to host halos containing the following solar masses (:math:`M_\odot`): the LMC (:math:`10^{11}`), Milky Way (:math:`10^{12}`), Group (:math:`10^{13}`), LCluster (:math:`4\times 10^{14}`), and Cluster (:math:`10^{15}`) suites.

Each halo in these suites has several processed datasets describing its subhalos and particles. (More detailed descriptions of each dataset for experts can be found below.)

* ``halos`` - Lightweight dataset giving the evolution of subhalos over time according to the Rockstar and Symgind halo finders. *This is the main dataset that most users will be interested in. It is light-weight enough to use on a personal computer.*
* ``trees`` - Full Rockstar + consistent-trees merger treees for the entire high-resolution region around the halo.
* ``particles`` - Particles of the halo and all its subhalos, tracked over time.
* ``full_snapshots`` - Full Gadget-2 snapshots for the entire zoom-in simulation

.. list-table::
	:header-rows: 1
		
	* - Suite Name
	  - :math:`N_{\rm hosts}`
	  - ``halos``
	  - ``trees``
	  - ``particles``
	  - ``full_snapshots``
	  - Total
	* - SymphonyLMC
	  - 39
	  - 0.16 GB
	  - 18 GB
	  - 121 GB
	  - 183 GB
	  - 332 GB
	* - SymphonyMilkyWay
	  - 45
	  - 0.77 GB
	  - 62 GB
	  - 206 GB
	  - 387 GB
	  - 654 GB
	* - SymphonyGroup
	  - 49
	  - 0.83 GB
	  - 227 GB
	  - 338 GB
	  - 1.1 TB
	  - 1.6 TB
	* - SymphonyLCluster
	  - 33
	  - 0.18 GB
	  - 21 GB
	  - 71 GB
	  - 162 GB
	  - 253 GB
	* - SymphonyCluster
	  - 96
	  - 1.1 GB
	  - 160 GB
	  - N/A
	  - 806 GB
	  - 965 GB
	* - Total
	  - 262
	  - 3.8 GB
	  - 482 GB
	  - 734 GB
	  - 2.6 TB
	  - 3.8 TB

No particle-tracking routines have been run on SymphonyCluster. This means that there are no Symfind catalogues for this suite and no ``particles`` dataset.

The ``trees`` and ``particles`` datasets rely on the ``halos`` dataset also being present.

**Dataset details**:

- ``halos``: This dataset contains information host halos and subhalos, including the full main-branch evolution of these objects. Disrupted subhalos and "splashback" subhalos are also included. Separate catalogues are provided, constructed with (a) a combination of Rockstar (`Behroozi et al. 2013a <https://ui.adsabs.harvard.edu/abs/2013ApJ...762..109B/abstract>`__) and consistent-trees (`Behroozi et al. 2013b <https://ui.adsabs.harvard.edu/abs/2013ApJ...763...18B/abstract>`__) and (b) Symfind (Mansfield et al., in prep). Numerous numerical artifacts in the Rockstar catalogs have been fixed. A custom format is used for these halos, as descibed in the :doc:`Working with Subhalos <working_with_subhalos>` and the :doc:`Library Documentation <symlib_documentation>` pages. *This dataset is light-weight enough to analyze on a personal computer, is relatively easy to use, and will be sufficient for the majority of halo-based research tasks.*
- ``trees``: This dataset contains merger trees for all objects in the high-resolution region of the box, not just the subhalos near the central host. These trees were generated with Rockstar (`Behroozi et al. 2013a <https://ui.adsabs.harvard.edu/abs/2013ApJ...762..109B/abstract>`__) and consistent-trees (`Behroozi et al. 2013b <https://ui.adsabs.harvard.edu/abs/2013ApJ...763...18B/abstract>`__), but have been converted to a more efficient and easy-to-use format and have been heavily annotated with several custom routines, as described in the :doc:`Working with Merger Trees <working_with_trees>` and the :doc:`Library Documentation <symlib_documentation>` pages. *This dataset is large enough that it will be difficult to analyze on a personal computer, and it will be unnecessary for most halo-based analysis tasks, even those that require looking at the evolution of objects over time.* You should use this dataset only if you want to look at poorly-resolved subhalos, objects far from the host, or completed mergers that occured between future subhalos prior to their infall into the host.
- ``particles``: This dataset contains all the particles associated with the host and its subhalos. Particles are organized according to which subhalo they orbited prior to becoming a subhalo (see Mansfield et al., 2023), making it easy to follow the evolution of the mass around subhalos over time. This sata is stored in a custom, compressed format and will need to be read with library functions (see :doc:`Working with Particles <working_with_particles>` and the :doc:`Library Documentation <symlib_documentation>`).
- ``full_snapshots``: This dataset contains full Gadget-2 snapshots for the entire zoom-in simulation at :math:`a(t)=0.2,\ 0.33,\ 0.5\ 0.67, 1.00`. These snapshots can be read with any software that can read Gadget-2 snapshots. IDs are 32-bit integers, positions and velocities use 32-bit floats. High-reoslution particles are stored at type-1 particles, and lower resolution particles are stored as successively higher particle types.


Accessing Data
--------------

Data is stored as a series of ``.tar`` files at the password-protected `s3df.slac.stanford.edu/data/kipac/symphony <s3df.slac.stanford.edu/data/kipac/symphony>`__.

.. button-ref:: data_access
   :ref-type: doc
   :color: secondary

   Get data access


If you're comfortable working with halo data already, you'll want to see the :doc:`Quickstart <quickstart>` page! Otherwise, check out the :doc:`Library Info pages <symlib>` to learn how to use the data.
   
