Getting Started
===============

Most users will be interested in ``symlib``, the Python library for reading and
analyzing Symphony's data. This page will walk you through using this library
for the first time. It will be enough for you to do real and meaningful
analysis, but will not try to be a complete reference.

If you are new to analysing simulations, you may want to consult
the terminology page, which will explain some terms and jargon.

Installation
------------

``symlib`` can be installed with pip.

.. code-block:: console

	$ pip install -i https://test.pypi.org/simple/ symlib

.. note::
   Currently hosted on the test PyPi server. Will move to the real one later.

Downloading Data
----------------

Symphony data is organized by halo. You can download an example halo,
``HaloExample.tar.gz`` for the purposes of this tutorial from my
`Google Drive <https://drive.google.com/file/d/1pnSqMtXDT_cE8lD3tMsA1NIxZEfQNS2z/view?usp=sharing>`_. Expand this file into a directory with the following
command. You can delete the tar file afterwards.

.. code-block:: console

	$ tar xzf HaloExample.tar.gz 

(This halo is a portion of Halo169 from the MWest suite.)
	
.. note::
   Obviously we'll do something different later. We'll use a non-MWest halo and
   won't use my Drive account for the example file.
   
Working With Subhaloes
----------------------

The first step to working with simulation data is loading the simulation
parameters. This is a dictionary containing cosmological and numerical
information. These can be looked up in ``parameter_table`` by supplying the
name of the simulation suite.

.. code-block:: python

	import symlib

	params = symlib.parameter_table["ExampleSuite"]
	print(params)

.. code-block:: python
				
	>>> {'flat': True, 'H0': 70.0, 'Om0': 0.286, 'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95, 'eps': 0.17, 'mp': 281981.0, 'h100': 0.7}

The first six values are `colossus <https://bdiemer.bitbucket.io/colossus/>`_-compatible cosmological parameters, and the last three are the Plummer equivalent force softening scale in comoving kpc/h, the dark matter particle mass in Msun/h, and H0/(100 km/s/Mpc).

Read in the largest subhalos of the example host halo with the
``read_subhaloes(parameters, directory)`` function.

.. code-block::

   subhalos, histories = symlib.read_subhalos(
       params, "path/to/HaloExample")

There are two return return values, ``subhalos`` and ``histories``. In your
code, you'll probably want to abbreviate these as ``s`` and ``h``.

``subhalos`` is a 2D array where the first index accesses different halos and
second different times. It contains time-dependent information, like mass and
position. ``subhalos[0,235]`` is the central subhalo at snapshot
235, the last snapshot in the simulation. ``subhalos[3, 100]`` is the third
largest subhalo, surviving or disrupted, at snapshot 100. Subhalos are ordered
according to the largest mass they ever had, Mpeak. ``histories`` contains
time-independent information on each halo, like Mpeak and and the snapshot
when the subhalo first falls into the central.

``subhalos`` is a numpy structured array has type ``SUBHALO_DTYPE`` and has
various fields that can be accessed by their names:

* ``"x"`` - Position in comoving Mpc/h
* ``"v"`` - Velocity in km/s
* ``"mvir"`` - Mass in Msun/h
* ``"rvir"`` - Radius in comoving Mpc/h
* ``"ok"`` - ``True`` if the halo exists at the given snapshot, ``False``
  otherwise.

``subhalos`` contains other value, too. These are described in the full symlib
documentation.

Our first step with analyzing any simulation data should be to look at it
qualitatively. We'll start by looking at the positions of the major subhalos
around our central halo at the last snapshot of the simulation (z=0).

.. code-block:: python

	import matplotlib.pyplot as plt
				
	for i in range(len(subhalos)):
		subhalo = subhalos[i]
		if not subhalo[-1]["ok"]: continue

		if i == 0:
			# Color the central halo differently.
			color = "tab:red"
		else:
			color = "tab:blue"

		x = subhalo["x"][-1,0]
		y = subhalo["x"][-1,1]
		r = subhalo["rvir"][-1]
		symlib.plot_circle(x, y, r, c=color, lw=3)

	plt.xlim(-0.3, +0.3)
	plt.ylim(-0.3, +0.3)
	plt.xlabel(r"$X\ (h^{-1}{\rm Mpc})$")
	plt.ylabel(r"$Y\ (h^{-1}{\rm Mpc})$")


Here, the central halo at index 0 is red and all is subhalos are blue.
We used a built-in utility function called ``plot_circle`` and
needed to skip over some subhalos which disrupted before the final snapshot.

.. note::
   Currently only the ten largest subhalos are stored in this file. That will be
   changed later.

We can also plot the growth of the largest subhalos over time. We will track
time through the scale factor, a(z), which we can get from the library
function ``scale_factors``. We'll check against the maximum mass value tabulated
in ``histories``, which are labeled at ``"mpeak"``.

.. code-block:: python

	colors = ["k", "tab:red",
	          "tab:orange", "tab:green",
			  "tab:blue", "tab:purple"]

	scales = symlib.scale_factors()
			  
	for i in range(6):
		subhalo = subhalos[i]

		# Plot growth history
		ok = subhalo["ok"]
		plt.plot(scales[ok], subhalo["mvir"][ok],
		    color=colors[i], lw=3)

		# Plot M_peak.
		mpeak = histories[i]["mpeak"]
		plt.plot([1/50, 1], [mpeak, mpeak],
		    "--", lw=1.5, color=colors[i])
		
	plt.xscale("log")
	plt.yscale("log")
	plt.xlabel(r"$a(z)$")
	plt.ylabel(r"$M_{\rm vir}$")

Here the central halo is in black and its subhalos are in color. We can see that
the maximum masses tracked by the ``histories`` look correct and that a number
of significant subhalos disrupted within our halo long ago. The red curve is an
analog for our Milky Way's LMC, and the orange curve is an analog for the
Gaia-Enceladus event.

.. note::
  Will put other insights here.


Working With Particles
----------------------

Particles are split up by subhalo. A subhalo owns all the particles that were
part of it before it became a subhalo, except the ones that belong to its own
subhalos or that already belonged to a bigger halo.

.. note::
   Need to put in an option that allows people to ignore that first constraint.

Let's analyze the particles of a subhalo at the time that it first falls into
out central halo. We'll focus on subhalo 4 and use the function
``read_particles``.

.. code-block:: python

	merger_snap = histories[4]["merger_snap"]
	x, v = symlib.read_particles(
	    "path/to/HaloExample", 4, merger_snap, ["x", "v"])

The function requires that you tell it the subhalo you want, the snapshot you
want, and the names of the variables you're reading. We'll make two plots: the
first will be the positions of the subhalo's particles around it and the second
will be a "phase diagram" showing their radii and radial velocities. We'll use
the function ``clean_particles``, which handles transforming everything out
of the simulation's code units, transforming into the inertial frame of the
subhalo, and removing particles that don't belong to the subhalo yet.

.. code-block:: python

	import matplotlib.colors as mpl_colors
				
	fig, ax = plt.subplots(2, figsize=(16, 8))

	scale = scales[merger_snap]
	x, v, idx = symlib.clean_particles(
	    params, x, v, subhalos[4], scale)

	r_max = 100 # kpc
	v_max = 100 # km/s

	r = np.sqrt(np.sum(x**2, axis=1))
	r_hat = x / r
	vr = np.dot(r_hat, v)
	
	ax[0].hist2d(
		x[:,0], x[:,1], bins=100,
		range=((-r_max, r_max), (-r_max, r_max)),
		norm=mpl_colors.LogNorm(1, 1000)
	)

	ax[1].hist2d(
		r, v, bins=100,
		range=((-r_max, r_max), (-r_max, r_max)),
		norm=mpl_colors.LogNorm(1, 1000)
	)

Put text discussing this here. (``"merger_snap"``, ``idx``, what's in the plots,
etc.)

.. note::
   Note that the units changed here. Need to figure out how to handle this...
