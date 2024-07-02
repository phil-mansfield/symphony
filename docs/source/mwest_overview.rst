Milky Way-est Simulation Suite
==================================

The "Milky Way-est" Suite is a publicly-accessible collection of 20 high-resolution zoom-in simulations of dark matter halos with Milky Way (MW)-like properties, including having a Large Margellanic Cloud (LMC)-like subhalo around present day and a Gaia-Sausage-Enceladus (GSE)-like subhalo that merged with it roughly 10 billion years ago. The LMC and GSE both mark important accretion events that shape the Milky Way's formation, and this suite allows us to study the Milky Way's formation history and subhalo population in this context.

The grid below shows visualizations of the 20 Milky Way-est hosts, arranged from left to right in order of increasing host halo mass and from bottom to top in order of increasing host halo concentration. Visualizations were created using these phase-space tessellation method described in `Kaehler 2017 <https://www.sciencedirect.com/science/article/abs/pii/S2213133716301536?via%3Dihub>`_, `2018 <https://library.imaging.org/ei/articles/30/1/art00005>`_.

.. image:: mwest_grid.pdf
   :width: 800

Please see the other pages to learn how to use the Milky Way-est suite!
           
Please refer to `D. Buch et al. (2024) <https://arxiv.org/abs/2404.08043>`_ for complete technical details, including numerical and cosmological parameters, additional properties of MW hosts, LMC analogs, and GSE analogs, and analyses of Milky Way-est subhalo populations.

FAQs
------
.. dropdown:: Why "Milky Way-est"?
   :animate: fade-in

   The name of the suite, "Milky Way-est," was born out of a need to differentiate the halos in this suite from a generic halo that resembles the Milky Way. Many researchers refer to objects that have properties similar to the Milky Way as "Milky Way-like," but this can be used for halos with Milky Way-mass, with LMC-like analogs, or with other Milky Way-like properties. So, we opt for a term that conveys the inclusion of many important Milky Way-like criteria (and one that's shorter than, for example, "much more Milky Way-like than your average Milky Way-like object").

.. dropdown:: Is this suite part of Symphony?
   :animate: fade-in

   Nope! But, they are somewhat related. Milky Way-est halos and Symphony MW-mass halos are run using the same cosmological parameters, initial conditions code, simulation code, and post-processing code. This makes them great candidates for comparative studies, and the Milky Way-est paper goes into detail about analyses done on Milky Way-est compared to Symphony MW. However, halos in Symphony MW have a different mass and concentration distribution than Milky Way-est halos, and unlike Milky Way-est, they are not explicitly selected to include LMC and GSE analogs.

.. dropdown:: If I need to abbreviate Milky Way-est while referring to it, what should I do?
   :animate: fade-in

   We usually refer to Milky Way-est by its full name for clarity, but if you really really need to abbreviate (e.g., if you are making a very small plot with very little space for labels), you can use "MW-est"! (It is also valid to use "MWest" but we caution that readers unfamiliar with the suite might misread this as "M West" instead of "MW est"!)
   
.. Milky Way-est Subpages
.. -----------------------
.. toctree::
   :hidden:
	  
   mwest_datainfo
   mwest_analysis
   mwest_credits
   
