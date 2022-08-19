Background Reading
==================

Simulations can be a complicated subject and may be intimidating to students and other first-time users. This page contains links to papers, lectures, and other resources that we have found helpful for learning about simulations, dark matter, and cosmology.

Cosmology Background Reading
----------------------------

When first starting to learn about cosmology, it's nice to get a bird's eye view before delving into technical specifics. The Wikipedia pages on `dark matter <https://en.wikipedia.org/wiki/Dark_matter>`__ and the cosmological model `Lambda-CDM <https://en.wikipedia.org/wiki/Lambda-CDM_model>`__ are good overviews, although be warned that as you get to the ends of these articles, you are more likely to be reading about the personal interests of individual editors rather than the consensus of experts. This is especially true for many of the topics listed as "challenges" on the Lambda-CDM page. 

These review notes by `Olive & Peacock <https://pdg.lbl.gov/2006/reviews/bigbangrpp.pdf>`__ (up to section 19.2.3) serve as a good technical introduction to the types of calculations that one does in cosmology, but may be overkill if you only want to use these simulations. The minimum technical knowledge cosmology that you will need to work with these simulations is an understanding of *comoving distances,* *proper/physical distances,* and *scale factors*. This page by `Ned Wright <https://www.astro.ucla.edu/~wright/cosmo_02.htm>`__ gives an overview of these topics. If you would like to look at alternative explanations, sections IV and V from Anthony Lewis's `lecture notes <https://cosmologist.info/teaching/Cosmology/Cosmology.pdf>`__ are also quite good.

Depending on what you want to do with these simulations, you may encounter values that have an :math:`h` factor in their units (e.g. :math:`h^{-1}M_\odot`). In this case, you should read `Croton (2013) <https://arxiv.org/pdf/1308.4150.pdf>`__ (deligtfully titled "Damn You, Little :math:`h`!").

Lastly, you may be interested in reading a curated reivew of potential problems for the Lambda-CDM cosmological model in `Bullock & Boylan-Kolchin (2019) <https://arxiv.org/pdf/1707.04256.pdf>`__.  Be aware that many of these topics are as controversial as they are exciting. If you are interested in any of them, it would be best to read several papers on each "side" of the issue rather than a single author's summary.

Dark Matter Halo Background Reading
-----------------------------------

The Symphony simulation suite focuses on dark matter halos. Dark matter halos are massive, egg-shaped structures made out of dark matter that galaxies form within. They are about fifty times wider than the galaxies inside of them and can can be tens to thousands of times more massive.

These graduate-level `lecture notes <https://campuspress.yale.edu/vdbosch/>`__ from Frank van den Bosch discuss the properties of dark matter halos, including their densities, shapes, spins, and the swarms of smaller subhalos that orbit inside them. Section 2.4 in `Mansfield & Aveztruz (2021) <https://arxiv.org/pdf/2008.08591.pdf>`__ gives a technical description of all the halo properties contained in Symphony's data sets.


Simulation Background Reading
-----------------------------

There are a number of review papers which give an overview of general topics related to dakr matter simulations. All of these papers have different goals, different target audiences, and different sets of topics. We have found that `Kuhlen et al. (2012) <https://arxiv.org/pdf/1209.5745.pdf>`__ is a good introduction for starting students.

Frank van den Bosch also has high-quality `lecture notes <http://www.astro.yale.edu/vdbosch/astro610_lecture20.pdf>`__ which will be useful to readers interested in understanding how simulations work and what their limitations are.

Software Reading
----------------

Beyond its own intenral libraries and analysis pipelines, Symphony makes use of many existing software projects. 

First, a program is needed to generate the starting locations of millions of particles in the early universe. Doing this involves working out the statistics of early fluctuations in the density of the unvierse givem our chosen cosmological model and then picking out a random set of particle locations that end up following these statistics. In our case, we're trying to essentially rerun an existing simulation, but with more particles stuffed into a specific interesting region, so there are additional constraints on where particles can be placed. We used the code `MUSIC <https://www-n.oca.eu/ohahn/MUSIC/>`__ (`Hahn & Abel 2011 <https://arxiv.org/pdf/1103.6031.pdf>`__) to do this.

The next step is simulating where those particles move over time. The positions and velocities of particles are controled by a massive system of ordinary differential equations, and running the simulation reduces to integrating these equations. This is known as the "N-body problem."  We make use of the most popular code for solving the N-body problem at the scales we are interested in: `Gadget-2 <https://wwwmpa.mpa-garching.mpg.de/gadget/>`__ (`Springel 2005 <https://arxiv.org/pdf/astro-ph/0505010.pdf>`__).

As the simulation runs, it will regularly stop and then write out the positions and velocities of particles to a series of files before continuing. Each set of files like this is called a "snapshot." Once the simulation finishes, we run a program called a "halo finder" on each snapshot. Halo finders try to identify the dense clumps of particles that correspond to the dark matter halos that we want to analyze. This is a relatively easy task when a halo is sitting on its own far away from other structures, but halo finders live or die based on how well they can identify small "subhalos" orbiting around larger ones. We use the halo finder `Rockstar <https://bitbucket.org/gfcstanford/rockstar/src/main/>`__ (`Behroozi et al. 2013a <https://arxiv.org/pdf/1110.4372.pdf>`__). Rockstar usually performs as well or better than other halo finders in common use (`Knebe et al. 2011 <https://arxiv.org/pdf/1104.0949.pdf>`__), but there is controversy over how well halo finders work in general. One of the goals of the Symphony suite is to study this issue in more depth.

Once halos have been identified, one could be interested in how a given halo evolves and changes over time. To do this, a program is needed to match subhalos up across snapshots. Halos can smash into one another and merge over time, but haloes essentially never split apart into smaller halos after forming. This gives these connections a tree-like structure, with a large fan small objects at earlier times that merge together into massive objects at later times. The resulting data structures are called "merger trees." We use the merger tree finder `consistent-trees <https://bitbucket.org/pbehroozi/consistent-trees/src/main/>`__ (`Behroozi et al. 2013b <https://arxiv.org/abs/1110.4370>`__), which is designed to reconstruct parts of the tree that might have been missed during the halo finder stage.

In the real universe, we mostly observe galaxies, not dark matter halos. To make predictions that can be directly compared with observations, you will often need a "galaxy--halo model." A Galaxy--halo model is any method that places galaxies in dark matter halos and decides on their properties. These can range from complicated methods, like full, expensive simulations that track the evolution of gas and stars, to simple methods, like applyng a fittign formula. You can see an overview of these methods in section 2 of the reivew paper `Wechsler & Tinker (2018) <https://arxiv.org/pdf/1804.03097.pdf>`__. In Symphony, we use a method somewhere in the middle in terms of complexity called Universe Machine (`Behroozi et al. 2018 <https://arxiv.org/pdf/1806.07893.pdf>`__). Universe Machine match starts with a wide range of galaxy obsevations and experiments with different ways that the growth of galaxies and the growth of halos could be connected. It checks these different experiements against observations Then it uses the experiment that matches observations to assign stars to dark matter halos.
