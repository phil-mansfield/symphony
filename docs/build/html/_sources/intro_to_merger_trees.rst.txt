Introduction to Merger Trees
============================

What Are Merger Trees?
----------------------

For some analysis, the set of subhaloes described above won't be enough. In some cases you might want to know about the extended history of subhalos, including objects that merged with those subhalos before they fell into the host. You would need to use a merger tree.

A merger tree is an array that contains all the halos and subhalos in a simulation across all times. It also has additional structure and information which allows one to figure out which halos at an earlier snapshot evolve into which halos at a later snapshot. This includes "mergers," events that occur when subhalos disrupt and contribute most of their mass to a larger host. This is a bit of a different definition than what we qualitatively think of as a merger: tree-mergers can happen many orbits after a subhalo falls into a host.

The merger trees in Symphony (generated with the merger tree code `consistent-trees <https://bitbucket.org/pbehroozi/consistent-trees>`__) are 1D arrays. These arrays are made of of separate sequences called "branches." A branch contains a single halo as it evolves over time. In Symphony's trees, halos in a branch are stored together, starting with the the halo's *last* and ending with the *first*. Every snapshot between the first and the last is included. Below is an example of the merger tree in a simulation with four snapshots, which has one halo that survived through all of these snapshots. The numbers in each circle show the index of that halo's data in the array.

.. image:: tree_1.png
    :height: 400
	     
Real simulations will have many halos, and these halos might not exist at all snapshots. When one halo's branch ends, the next halo's branch will begin at the next element. The image below shows an example of what a simulation with several halos might look like.
	     
.. image:: tree_2.png
    :height: 400

The tree also contains information on what happens to a halo after it disrupts. Some halos that disappear do so because they were orbiting a larger halo and were shredded apart. In these cases, the merger tree marks the event as a merger. The arrows below show an example set of mergers.
	     
.. image:: tree_3.png
    :height: 400

Note that more than one merger can occur within a halo in a single snapshot. This image also illustrates the ordering of branches within the tree ("depth-first ordering").

Merger information is stored in the tree the snapshot before the merger occurs. Each halo keeps track of its "co-progenitor", the last halo in the next branch that merges in that snapshot. It's easiest to understand visually:
	     
.. image:: tree_4.png
    :height: 400

In practice, this means that finding all the mergers for a host halo invovles hopping around from halo to halo.

Mergr tree analysis is an advanced technique and may not be necessary for many users. If one is only interested in resolved substructure of the host halo, virtually everything can be done with the subhalo arrays above. The tree essentially only allows one to analyze subresolution subhalos and objects far away from the host halo. Both tasks must be done with caution.

Using Merger Trees with Symlib
------------------------------

Full merger trees are more computationally intense than subhalo catalogs, meaning that they are a little more complicated to use than the host's subhalos. Instead of using structured arrays, individual variables are read from disk as 1D arrays that have the tree's ordering. This allows you to only load variables you need. Additionally, snapshots where a halo does not exist are not included in the tree. This makes indexing more complicated, but saves memory.

One last caveat is that tree variables are stored in consistent-trees's native units. These are inhomogenous and are listed in the full ``symlib`` documentation page.

As a first example, we'll do some analysis that doesn't care about the connections between different branches. We'll calculate :math:`M_{\rm peak}` (the maximum :math:`M_{\rm vir}` that the halo ever had). Because tree reading is a more advanced and less useful procedure than reading the host's subhalo information, this will be done purely by example. The full symlib documentation contains more detailed information on tree-reading and tree-manipulating functions.

Omitting some standard preamble and most of the plotting code, the following code block will:

* Read mass data in from the tree
* Loop over branches and calculate :math:`M_{\rm peak}`
* Split those branches into "real" branches and branches that are probably artifacts.
* Create mass functions for those groups

.. code-block:: python

    # Read in tree data
    b = symlib.read_branches(sim_dir)
    param = symlib.simulation_parameters(sim_dir)
    # Tree variables are always returned as a list, so if you
    # only specify one, unpack it as a length-1 tuple.
    mvir, = symlib.read_tree(sim_dir, ["mvir"])
    # Convert units
    mvir, mp = mvir/param["h100"], param["mp"]/param["h100"]

    # Flag halo branches which are probably not artifacts.
    ok = b["is_real"] & (~b["is_disappear"])

    mpeak = np.zeros(len(b))

    # Loop over all branches and calculate M_peak.
    for i in range(len(mpeak)):
        if not ok[i]: continue
        start, end = b[i]["start"], b[i]["end"]
        mpeak[i] = np.max(mvir[start: end])
    
    # Find host subhalos
    mpeak_host = mpeak[b["is_main_sub"]]

    # Calculate the mass function of both groups of halos
    bins = np.logspace(np.log10(param["mp"]/param["h100"]), 13, 200)
    n_host, _ = np.histogram(mpeak_host, bins=bins)
    n_all, _ = np.histogram(mpeak, bins=bins)
    N_host = np.cumsum(n_host[::-1])[::-1]
    N_all = np.cumsum(n_all[::-1])[::-1]

    # Plot
    fig, ax = plt.subplots()
    left_bins = bins[:-1]
    plt.plot(left_bins, N_host, c="tab:red", label=r"${\rm Host\ subhalos}$")
    plt.plot(left_bins, N_all, c="tab:blue", label=r"${\rm All}$")

With a bit more ellided plotting code, this results in the following plot
    
.. image:: tree_mass_func.png
    :width: 500

The dashed vertical line has been added to show a rough resolution limit, 300 particles. As this plot shows, most of the contents of the tree file are outside the target host halo and most of the corresponding objects are poorly resolved. The objects returned by :func:`symlib.read_subhalos` only consist of the red curve down to the black dashed line.

As a second example, we will navigate through the merger tree to find the number of mergers in each snapshot. These will be split into real mergers and artifacts. The latter includes a few criteria, but mostly it's objects whose first snapshots occured already inside the host halo. The time resolution of these simulations is high enough that any instances where this happens are either statistical noise or a halo which the merger tree had previously lost track of (i.e. a halo that merges twice). These objects have already been removed in the standard subhalo arrays.

This code will do the following:

 * Read in tree data, including snapshot and connectivity information.
 * Create a merger lookup table
 * Loop through the host's branch and use the lookup table to find the mergers in each snapshot.
 * Classify those mergers based on whether or not they are artifacts.
 * Convert snapshot information into cosmological times.

As before, some plotting code and standard setup code that reads in parameters and subhalo information has been omitted.

.. code-block:: python

    sim_dir = "path/to/ExampleHalo"
    h, hist = symlib.read_subhalos(sim_dir)
		
    # Read in tree data
    b = symlib.read_branches(sim_dir)
    dfid, next_co_prog, snap = symlib.read_tree(
        sim_dir, ["dfid", "next_co_prog", "snap"]
    )

    host_branch = b[hist["branch_idx"][0]]
    host_start = host_branch["start"]
    host_end = host_branch["end"]
    
    # Flag halo branches which are probably not artifacts.
    ok = b["is_real"] & (~b["is_disappear"])
    
    # Counting the number of mergers. Requires a lookup
    # table, which we construct from the branch
    # information and the depth-first IDs ("dfid")
    table = symlib.merger_lookup_table(b, dfid)
    n_mergers = np.zeros(host_end - host_start, dtype=int)
    n_artifacts = np.zeros(host_end - host_start, dtype=int)
    for i in range(host_start, host_end):
        branch_idx = symlib.find_all_merger_branches(
            b, table, next_co_prog, i)
        n_mergers[i - host_start] = np.sum(ok[branch_idx])
        n_artifacts[i - host_start] = np.sum(~ok[branch_idx])

    # Getting the scale factor of each snapshot.
    host_snap = snap[host_start: host_end]
    host_scale = scale[host_snap]

    # Plot
    fig, ax = plt.subplots()
    ax.plot(host_scale, n_mergers, "tab:blue",
            label=r"$N_{\rm merger}$")
    ax.plot(host_scale, n_mergers, "tab:blue",
            label=r"$N_{\rm artifacts}$")
		
		
		
.. image:: tree_mergers.png
   :width: 500
	    
