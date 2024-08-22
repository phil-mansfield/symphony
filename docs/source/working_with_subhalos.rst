Subhalo Tutorial
================

Symphony data is stored in cutsom data formats and must be read with its analysis library, ``symlib``. This page contains tutorials and example code showing how to use this library. A complete, detailed reference for this library can be found :doc:`here <symlib_documentation>`.

Installation
------------

As disucssed in the :doc:`Data Access <data_access>` page, ``symlib``  can be downloaded with the pip tool.

.. code-block:: console

	$ pip install symlib -U


Reading in Subhalo Data
-----------------------

The first step to working with simulation data is loading the simulation
parameters. This is a dictionary containing cosmological and numerical
information. These can be looked up with :func:`symlib.simulation_parameters` by giving the directory of the halo as an argument.

.. code-block:: python

	import symlib

	# sim_dir is the location of a single halo. We'll talk about
	# auto-generating this later, but if you downloaded "Halo933"
	# in the suite "SymphonyLMC" to the directory "data/", this
	# would be "data/SymphonyLMC/Halo933"
	sim_dir = "path/to/ExampleHalo"
	
	params = symlib.simulation_parameters(sim_dir)
	print(params)

.. code-block:: python
				
    {'flat': True, 'H0': 70.0, 'Om0': 0.286, 'Ob0': 0.049,
     'sigma8': 0.82, 'ns': 0.95, 'eps': 0.17, 'mp': 281981.0,
     'h100': 0.7}

The first six values are `colossus <https://bdiemer.bitbucket.io/colossus/>`__-compatible cosmological parameters, the next two are numerical parameters (``"eps"`` is the radius of each particle in comoving :math:`h^{-1}{\rm kpc}`, and ``mp`` is the mass of each particle in :math:`h^{-1}M_\odot`). The last value is :math:`h_{100} = H_0/(100\ {\rm km/s})`, which is often convenient to have.

Next, we read in the subhalos with two functions

.. code-block::

   s, hist = symlib.read_symfind(sim_dir)
   r, _ = symlib.read_rockstar(sim_dir)

The first function returns information from our Symfind subhalo finder, and the second returns information from the widely-used Rockstar halo finder.

Both functions have two return values: the first is information a 2D array representing the evolution of the host halo and subhaloes snapshot-by-snapshot, and the second is a 1D array summarizing various statistics about the overall evolution of each object.

In ``s`` and ``r``, the first index accesses different halos and
second different times. It contains information like mass and
position. ``r[0,235]`` is the host halo at snapshot
235, the last snapshot in the simulation. ``r[3, 100]`` is the third
largest subhalo, including disrupted subhalos, at snapshot 100. Subhalos are ordered according to the largest mass they achieved prior to becoming a subhalo, :math:`M_{\rm peak}`. Halos stay at the same index across their lifetimes.

``hist`` contains summary information about a halo's full history, including :math:`M_{\rm peak}` and when that subhalo fell into the host. Its length and ordering are the same as the first index of ``s`` and ``r``. 

``s`` and ``r`` are numpy structured arrays has type :data:`symlib.SYMLIB_DTYPE` and :data:`symlib.ROCKSTAR_DTYPE`, and ``hist`` is a structured array with type :data:`symlib.HISTORY_DTYPE`. Structured arrays are arrays that have different fields which can be accessed with strings. For example, ``r[3,100]["m"]`` and ``r["m"][3,100]`` both return the mass, :math:`M_{\rm vir}` of the third most massive halo at snapshot 100. The non-string indices obey normal numpy indexing rules, so you can use slicing, boolean indexing, axis removal and whatever other tricks you use with normal numpy arrays.

The full set of fields for ``s`` and ``r``  are described in the :data:`symlib.SYMFIND_DTYPE` and :data:`symlib.ROCKSTAR_DTYPE` documentation. In this tutorial we will only use:

* ``"x"`` - three-dimensional position vector (x, y, z)
* ``"v"`` - three-dimensional velocity vector (v_x, v_y, v_z)
* ``"m"`` - Mass (virial mass for non-subhaloes and bound mass for subhaloes)
* ``"rvir"`` (Rockstar) ``"r_half"`` (Symfind) - Radius (virial radius and half-mass radius, respectively)
* ``"ok"`` - ``True`` if the halo was tracked by the halo finder the given snapshot, ``False`` if the halo was not tracked by the halo finder at the given snapshot.

Fields in ``hist`` will be explained as needed, but can be found in full in the :data:`symlib.HISTORY_DTYPE` documentation.

Example Subhalo Analysis: Plotting Postions
-------------------------------------------
   
Our first step with analyzing any simulation data will be to look at it
qualitatively. We'll start by looking at the positions of the major ROkcstar subhalos around our central halo at the last snapshot of the simulation. We will plot the central halo in one color and the subhalos in another. We'll also need to avoid plotting any of the subhalos that were destroyed before the end of the simulation.

We'll use a utility function, :func:`symlib.plot_circle` to make circles that will represent each halo.

.. _halo_position_example:

.. code-block:: python

    import symlib
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    
    sim_dir = "path/to/ExampleHalo"
	r, hist = symlib.read_rockstar(sim_dir)
    
    host = r[0,-1] # First halo, last snapshot.
    symlib.plot_circle(ax, host["x"][0], host["x"][1],
                       host["rvir"], c="tab:red")
		       
    for i in range(1, len(r)):
        sub = r[i,-1] # i-th halo, last snapshot.
        if not sub["ok"]: continue
        symlib.plot_circle(
            ax, sub["x"][0], sub["x"][1],
            sub["rvir"], c="tab:blue"
        )
    
With a little bit of additional pyplot work that we've ellided here, this gives us the following.

.. image:: positions.png
   :width: 500

From this, we can see that our host halo is surrounded by a swarm of subhalos. Bigger subhalos are rarer and generally closer to the center of the host. Some subhalos are outside the radius of the host. These "splashback subhalos" had been inside the host in the past but have temporarily orbited outside of it. They are included in the symlink catalogs by default.
	   
Let's review the concepts that went into creating this image:

* We read in simulation parameters and halo information with :func:`symlib.read_rockstar`.
* We got the host halo at the last snapshot with ``r[0,-1]`` and the subhalos with ``r[i,-1]``.
* We got a vector representing the postion of the host by accessing ``host["x"]`` and the radius with ``host["rvir"]`` and were able to get similar quantities for subhalos.
* We needed to check ``sub["ok"]`` to make sure that the halo still existed at the snapshot we were interested in.

Here, the central halo at index 0 is red and all is subhalos are blue.
We used a built-in utility function called ``plot_circle`` and
needed to skip over some subhalos which disrupted before the final snapshot.

**Example exercise**

In the ``histories`` array, there is a field called ``merger_snap`` that gives the snapshot when a subhalo first fell into the host. Try coloring subhalos that fell in from the left side of the halo (:math:`x_{\rm infall} < 0`) differently from ones that fell in from the right.

Example Analysis: Mass Growth
-----------------------------

Now, we'll try analysis that's a bit more quantitative. We'll look at the growth of subhalos over time: looking at the growth of the host halo and its five most massive subhalos over time. To do this, we'll need to get the scale factors, :math:`a(z)`, for each snapshot with :func:`symlib.scale_factors`. We'll also use one of the fields in ``histories``, ``"merger_snap"`` which is the snapshot when the subhalo first fell into the host. We'll use it to plot times before infall as dashed lines and times afterwards as solid lines.

.. _mah_example:

.. code-block:: python
		
    sim_dir = "path/to/ExampleHalo"

    scale = symlib.scale_factors(sim_dir)
    s, hist = symlib.read_symfind(sim_dir)
    r, _ = symlib.read_rockstar(sim_dir)

    fig, ax = plt.subplots()
    colors = ["tab:red", "tab:orange", "tab:green",
              "tab:blue", "tab:purple"]

    # First, plot the host
    ok = r[0,:]["ok"]
    ax.plot(scale[ok], r[0,ok]["m"], c="k")

    # For now, let's only plot minor mergers which dirupt before
    # the end of the simulation in the Rockstar catalogue.
    is_target = (hist["merger_ratio"] < 0.1) & (~r["ok"][:,-1])
    targets = np.where(is_target)[0][:5]
    for i_color, i in enumerate(targets):
        # Plot the mass history of the rockstar subhalo
        ok = r[i,:]["ok"]
        ax.plot(scale[ok], r[i,ok]["m"], c=colors[i_color])
        # Plot the mass history of the symfind subhalo
        ok = s[i,:]["ok"]
        ax.plot(scale[ok], s[i,ok]["m"], c=colors[i_color], lw=1.5)
	

With a little bit of additional pyplot work, this gives us the following. The full script used to create this image, including the omitted pyplot code is shown in `examples/mah.py <https://github.com/phil-mansfield/symphony/blob/main/examples/mah.py>`__.

.. image:: mah.png
   :width: 500

Here we see that our subhalos spend most of their time in the simulation building up mass prior to falling in. The earlier-infalling halos shown here don't last for very long: they disrupt in a few snapshots! Others, like the green subhalo survive much longer.

Let's review the concepts that went into creating this image:

* We needed to read in scale factors with :func:`symlib.scale_factors` to figure out when each snapshot occured.
* We were able to figure out the snapshot when a subhalo fell into the host with ``histories``'s ``"merger_snap"`` field.
* The indices of structured arrays work just like normal numpy arrays, so we were able to select parts of them with the boolean arrays ``ok`` and ``is_sub``.

**Example exercise**

Try remaking this

You might have noticed that subhalos start losing mass before they actually start falling into the host (look at the transition from a dashed to solid line on the green curve in particular). Create a histogram showing :math:`R_{\rm peak}`/ :math:`R_{\rm virial}`, where :math:`R_{\rm peak}` is the distance between the subhalo and the host halo and :math:`R_{\rm virial}` is the virial radius of the host halo, both calculated at the time the subhalo reaches its peak mass.

Example Analysis: The Subhalo Mass Functions
--------------------------------------------

Lastly, let's try some more rigorous statistical analysis. So far we’ve been looking at a population of subhalos surrounding one host halo. Now, we’re going to measure the subhalo mass function for all of the host halos in the Milky Way suite. The subhalo mass function is a statistic that counts the number of subhalos orbiting a host halo as a function of the subhalo’s mass. It is essentially a cumulative histogram of subhalo mass. We'll need to look at :math:`N(>M_{\rm peak})`, the average number of subhalos per host halo whose maximum mass was larger than :math:`M_{\rm peak}`. 

In the previous exercise, we did analysis on the time when a subhalo reached its maximum mass, or :math:`M_{\rm peak}`. We can calculate that value ourselves or use the ``"mpeak"`` field of the ``histories`` array.

More importantly, to get good statistics we'll need to loop over all the host halos in the Milky Way suite, ``SymphonyMilkyWay``. One way to do this would be to manually store the names of all the halo directories, but instead we'll use library functions to do it. First, we'll count the number of halos in the Milky Way-mass suite with :func:`symlib.n_hosts`. Then, we can get directory names :func:`symlib.get_host_directory`, which takes the base directory, suite name, and the index of the halo you want to read. Together this lets you loop over halo directories.

Constructing a mass function has a bit more code overhead than the earlier examples: the important part is how the loop over files works.

.. _shmf_example:

.. code-block:: python

    base_dir = "path/to/base/dir"
    suite_name = "SymphonyMilkyWay"
    
    # Mass function bins and empty histogram.
    log_m_min, log_m_max, n_bin = 8, 12, 200
    bins = np.logspace(log_m_min, log_m_max, n_bin+1)
    N_vir = np.zeros(n_bin)

    n_hosts = symlib.n_hosts(suite_name)
    for i_host in range(n_hosts):
        sim_dir = symlib.get_host_directory(base_dir, suite_name, i_host)
	h, hist = symlib.read_subhalos(sim_dir)

	# Only count objects within R_vir
        host_rvir = h[0,-1]["rvir"]
        sub_x = h[:,-1]["x"]
        r = np.sqrt(np.sum(sub_x**2, axis=1))
        ok = h["ok"][:,-1] & (r < host_rvir)

        # Put in bins and add to cumulative histogram
        n_vir, _ = np.histogram(hist["mpeak"][ok][1:], bins=bins)
	N_vir += np.cumsum(n_vir[::-1])[::-1]/n_hosts

    plt.plot(bins[:-1], N_vir, "k")

With a little bit of additional pyplot work, this gives us the following. The full script used to create this image, including the omitted pyplot code is shown in `examples/mass_func.py <https://github.com/phil-mansfield/symphony/blob/main/examples/mass_func.py>`__.

.. image:: mass_func.png
   :width: 500

Here, we can see the classic form of the subhalo mass function. At smaller subhalo masses, decreasing the subhalo mass increases the number of subhalos and there’s an exponential cutoff as the subhalos approach the mass of the host halo.
   
Let's review the concepts that went into creating this image: 

* We needed to use :func:`symlib.n_hosts` to find the number of host halos in our target suite
* We needed to use :func:`symlib.get_host_directory` to find the names of the directories in the host halo.
* We needed the ``"mpeak"`` field of ``histories``
* We needed to do a little bit of array magic with numpy arrays, although this could also have been done in a less concise way.

**Example exercise**

You might notice that the plot above only includes subhalos with positions within the virial radius of the host halo. Try adding a curve for the mass function of surviving “splashback” subhalos, subhalos which have temporarily orbited outside of the host halo's virial radius, to this plot.

# TODO: add another tutorial for working with particle data and an example exercise.

