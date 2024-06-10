The Milky Way-est Simulation Suite
==================================

The "Milky Way-est" Suite is a publicly-accessible collection of 20 high-resolution zoom-in simulations of dark matter halos with Milky Way (MW)-like properties, including having a Large Margellanic Cloud (LMC)-like subhalo around present day and a Gaia-Sausage-Enceladus (GSE)-like subhalo that merged with it roughly 10 billion years ago. The LMC and GSE are both important accretion events that have shaped the Milky Way's formation, and this suite allows us to study the Milky Way's formation history and subhalo population.

.. image:: mwest_grid.pdf
   :width: 800

..
   can do the following for an iframe: [newline] .. raw:: html [leave one empty line] <iframe...>

.. grid:: 2
          
   .. grid-item-card::
      
      .. raw:: html

         <iframe id="igraph" scrolling="no" style="border:none;" seamless="seamless" src="https://plotly.com/~chris/1638.embed" height="525" width="100%"></iframe>
         
   .. grid-item-card:: hello

FAQs
------
.. dropdown:: Why "Milky Way-est"?
   :animate: fade-in

   The name of the suite, "Milky Way-est," was born out of a need to differentiate the halos in this suite from a generic halo that resembles the Milky Way. Many researchers refer to objects that have properties similar to the Milky Way as "Milky Way-like," but this can be used for halos with Milky Way-mass, with LMC-like analogs, or with other Milky Way-like properties. So, we opt for a term that conveys the inclusion of many important Milky Way-like criteria (and one that's shorter than, for example, "much more Milky Way-like than your average Milky Way-like object").

.. dropdown:: Is this suite part of Symphony?
   :animate: fade-in

   Nope! But, they are somewhat related. Milky Way-est halos and Symphony MW-mass halos are run using the same cosmological parameters, initial conditions code, simulation code, and post-processing code. This makes them great candidates for comparative studies, and the Milky Way-est paper goes into detail about analyses done on Milky Way-est compared to Symphony MW. However, halos in Symphony MW have a different mass and concentration distribution than Milky Way-est halos, and unlike Milky Way-est, they are not explicitly selected to include LMC and GSE analogs.
