Zoom-in Simulations
===================

"Zoom-in" simulations are designed to study individual halos at high resolutions. To do this, researchers run a large but low-resolution simulation, identify an interesting halo, and rerun the simulation with additional resolution added to the region which will eventually become that halo.

This site currently hosts two sets of zoom-in simulations. The :doc:`Symphony <symphony_overview>` simulation suite is a collection of five different smaller zoom-in simulations suites which cover masses ranging from LMC-mass dawrf galaxies to very high mass galaxy clusters. The :doc:`Milky Way-est <mwest_overview>` suite is a set of Milky-Way-mass dark matter halos which have been selected to match large mergers which were similar to those experienced by the actual Milky Way.

Instructions for downloading data can be found at the :doc:`Data Access Instructions <data_access>` page and a high level description of all the available datasets for each simulation can be found at the :doc:`Data Specifications <data_specifications>` page.

Most of the available data must be read through the lightweight ``symlib`` library, which allows you to access our specialized, space-efficient files and also structures the data in a way which is easier to use than most traditional data formats, particularly when tracking the evolution of subhalos and particles over time. Tutorials and documentations can be found on the :doc:`Symlib page <symlib>`.

.. toctree::
   :hidden:

   symphony_overview
   mwest_overview
   data_access
   data_specifications
   symlib
