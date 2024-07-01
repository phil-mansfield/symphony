Background Reading
==================

Simulations can be a complicated subject and may be intimidating to students and other first-time users. This page contains links to papers, lectures, and other resources that we have found helpful for learning about simulations, dark matter, and cosmology.

Cosmology Background Reading
----------------------------

When first starting to learn about cosmology, it's nice to get a bird's eye view before delving into technical specifics. The Wikipedia pages on `dark matter <https://en.wikipedia.org/wiki/Dark_matter>`__ and the cosmological model `Lambda-CDM <https://en.wikipedia.org/wiki/Lambda-CDM_model>`__ are good overviews. (Be warned that as you get to the ends of these articles, you are more likely to be reading about the personal interests of individual editors rather than the consensus of experts. This is especially true for many of the topics listed as "challenges" on the Lambda-CDM page.)

These review notes by `Olive & Peacock <https://pdg.lbl.gov/2006/reviews/bigbangrpp.pdf>`__ (up to section 19.2.3) serve as a good technical introduction to the types of calculations that one does in cosmology, but may be overkill if you only want to use these simulations. The minimum technical knowledge cosmology that you will need to work with these simulations is an understanding of *comoving distances,* *proper/physical distances,* and *scale factors*. This page by `Ned Wright <https://www.astro.ucla.edu/~wright/cosmo_02.htm>`__ gives an overview of these topics. If you would like to look at alternative explanations, sections IV and V from Anthony Lewis's `lecture notes <https://cosmologist.info/teaching/Cosmology/Cosmology.pdf>`__ are also quite good.

Depending on what you want to do with these simulations, you may encounter values that have an :math:`h` factor in their units (e.g. :math:`h^{-1}M_\odot`). In this case, you should read `Croton (2013) <https://arxiv.org/pdf/1308.4150.pdf>`__ (deligtfully titled "Damn You, Little :math:`h`!").

Lastly, you may be interested in reading a curated reivew of potential problems for the Lambda-CDM cosmological model in `Bullock & Boylan-Kolchin (2019) <https://arxiv.org/pdf/1707.04256.pdf>`__.  Be aware that many of these topics are as controversial as they are exciting. If you are interested in any of them, it would be best to read several papers on each "side" of the issue rather than a single author's summary.

Dark Matter Halo Background Reading
-----------------------------------

Dark matter halos are an important component of all the simulations hosted here. Dark matter halos are massive, egg-shaped structures made out of dark matter that galaxies form within. They are about fifty times wider than the galaxies inside of them and can can be tens to thousands of times more massive.

These graduate-level `lecture notes <https://campuspress.yale.edu/vdbosch/>`__ from Frank van den Bosch discuss the properties of dark matter halos, including their densities, shapes, spins, and the swarms of smaller subhalos that orbit inside them. Section 2.4 in `Mansfield & Aveztruz (2021) <https://arxiv.org/pdf/2008.08591.pdf>`__ gives a technical description of how many halo properties are computed by the population halo-finding tool Rockstar, and most other tools use similar techniques.

We encourage readers interested in understanding the numerical limitations of density profiles to read `Ludlow et al. (2019) <https://arxiv.org/pdf/1812.05777.pdf>`__, readers interested in the limitations of coarse-grained halo catalogs to read `Mansfield & Aveztruz (2021) <https://arxiv.org/pdf/2008.08591.pdf>`__, and readers interested in the limitations of orbiting subhalos to read `Mansfield et al. (2024) <https://arxiv.org/pdf/2308.10926.pdf>`__.

Simulation Background Reading
-----------------------------

There are a number of review papers which give an overview of general topics related to dakr matter simulations. All of these papers have different goals, different target audiences, and different sets of topics. We have found that `Kuhlen et al. (2012) <https://arxiv.org/pdf/1209.5745.pdf>`__ is a good introduction for starting students.

Frank van den Bosch also has high-quality `lecture notes <http://www.astro.yale.edu/vdbosch/astro610_lecture20.pdf>`__ which will be useful to readers interested in understanding how simulations work and what their limitations are.

Software Reading
----------------

Many software libraries are needed to run and analyze a dark matter simulation.

First, a program is needed to generate the starting locations of millions of particles in the early universe. Doing this involves working out the statistics of early fluctuations in the density of the unvierse given our chosen cosmological model and then picking out a random set of particle locations that end up following these statistics. In our case, we're trying to essentially rerun an existing simulation, but with more particles stuffed into a specific interesting region, so there are additional constraints on where particles can be placed. The zoom-in simulations hosted on this site used `MUSIC <https://www-n.oca.eu/ohahn/MUSIC/>`__ (`Hahn & Abel 2011 <https://arxiv.org/pdf/1103.6031.pdf>`__) to do this. Other simulations linked here use other schemes that share many basic properties, but are typically designed for ue on larger scales. Details can be found in the repsective release papers.

The next step is simulating where those particles move over time. The positions and velocities of particles are controlled by a massive system of ordinary differential equations, and running the simulation reduces to integrating these equations. This is known as the "N-body problem." The simulations hosted and linked to here make use of the most popular code for solving the N-body problem at the scales we are interested in: `Gadget-2 <https://wwwmpa.mpa-garching.mpg.de/gadget/>`__ (`Springel 2005 <https://arxiv.org/pdf/astro-ph/0505010.pdf>`__).

As the simulation runs, it will regularly stop and then write out the positions and velocities of particles to a series of files before continuing. Each set of files like this is called a "snapshot." Once the simulation finishes, we run a program called a "halo finder" on each snapshot. Halo finders try to identify the dense clumps of particles that correspond to the dark matter halos that we want to analyze. This is a relatively easy task when a halo is sitting on its own far away from other structures, but halo finders live or die based on how well they can identify small "subhalos" orbiting around larger ones. Many of our simulations use the halo finder `Rockstar <https://bitbucket.org/gfcstanford/rockstar/src/main/>`__ (`Behroozi et al. 2013a <https://arxiv.org/pdf/1110.4372.pdf>`__). Unfortunately, Rockstar and most other popular halo finding tools have substantial but poorly-characterized numerical biases. Some of our simulations also make use of another subhalo finder, Symfind, which is less biased and better characterized `Mansfield et al. (2024) <https://arxiv.org/pdf/2308.10926.pdf>`__.

Once halos have been identified, one could be interested in how a given halo evolves and changes over time. To do this, a program is needed to match subhalos up across snapshots. Halos can smash into one another and merge over time, but halos essentially never split apart into smaller halos after forming. This gives these connections a tree-like structure, with a large span of small objects at earlier times that merge together into massive objects at later times. The resulting data structures are called "merger trees." We use the merger tree code `consistent-trees <https://bitbucket.org/pbehroozi/consistent-trees/src/main/>`__ (`Behroozi et al. 2013b <https://arxiv.org/abs/1110.4370>`__), which is designed to reconstruct parts of the tree that might have been missed during the halo finder stage. Similar to halo finders, many merger trees have signficant numerical problems. Symfind also performs merger tree construction which avoids many of these problems.

In the real universe, we mostly observe galaxies, not dark matter halos. To make predictions that can be directly compared with observations, you will often need a "galaxy--halo model." A Galaxy--halo model is any method that places galaxies in dark matter halos and decides on their properties. These can range from complicated methods, like full, expensive simulations that track the evolution of gas and stars, to simple methods, like applyng a fitting formula. You can see an overview of these methods in section 2 of the reivew paper `Wechsler & Tinker (2018) <https://arxiv.org/pdf/1804.03097.pdf>`__. Many of our simulations  have had a method called Universe Machine (`Behroozi et al. 2018 <https://arxiv.org/pdf/1806.07893.pdf>`__) applied to them. Universe Machine starts with a wide range of galaxy obsevations and experiments with different ways that the growth of galaxies and the growth of halos could be connected. It checks these different "experiements" against observations. Then it uses the experiment that matches observations to assign stars to dark matter halos.
