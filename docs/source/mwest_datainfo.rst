Milky Way-est Data
==================================

Suite Information
---------------------------
For information about the Milky Way hosts, LMC analogs, and GSE analogs in the Milky Way-est suite, please see Table 2 of Buch et al. (2024).

We use the following cosmological parameters in our simulations: :math:`h=0.7, \Omega_{\rm m} = 0.286, \Omega_{\Lambda} = 0.714, \sigma_8 = 0.82, n_s=0.96`.

For the definitions of these parameters, see the ``symlib.simulation_parameters`` description in :doc:`symlib_documentation`.

Initial conditions for the zoom-in simulations were generated using MUSIC (`Hahn and Abel (2011) <https://academic.oup.com/mnras/article/415/3/2101/1045260>`_), the simulations themselves were run using GADGET-2 code (`Springel (2005) <https://academic.oup.com/mnras/article/364/4/1105/1042826>`_), and we generated halo catalogs and merger trees using ROCKSTAR (`Behroozi et al. (2013a) <https://iopscience.iop.org/article/10.1088/0004-637X/762/2/109>`_) and CONSISTENT-TREES (`Behroozi et al. (2013b) <https://iopscience.iop.org/article/10.1088/0004-637X/763/1/18>`_).

Data Access
-----------------

There are 2 approaches to using Milky Way-est data: 1) load ``symlib``-compatible files stored at ``s3df.slac.stanford.edu/data/kipac/symphony``, or 2) load simulation analysis output like the Milky Way-est authors with files stored at ``s3df.slac.stanford.edu/data/kipac/symphony/mwest``.

You can load whichever data format suits your needs!

Both have the same access point:

.. button-ref:: data_access
   :ref-type: doc
   :color: secondary
           
   Get data access

Analysis Options
-----------------

How would you like to make use of the Milky Way-est suite?

.. grid:: 2 2 2 2
   
   .. grid-item-card::

      .. dropdown:: Use ``symlib`` functionality with Milky Way-est

         This is a great option for folks who already use ``symlib`` tools for their existing data analysis pipeline, want to make use of the functions and documentation that ``symlib`` offers, or otherwise want a more guided/abstracted approach to working with these halos. See the :doc:`symlib` page for more information.

   .. grid-item-card::

      .. dropdown:: Load Milky Way-est halo catalogs directly

         This might be a good option for you if, upon gaining access to the data files, you'd like to get started making plots and iterating directly from the halo catalogs themselves. This is how the authors approached their analysis, and we've put together a pipeline you can replicate to directly load and analyze the data. See :doc:`mwest_analysis` for more!
