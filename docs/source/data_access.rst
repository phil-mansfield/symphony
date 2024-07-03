Data Access
===========

This page consists of data access instructions for the :doc:`Symphony <symphony_overview>` and :doc:`Milky Way-est <mwest_overview>` suites.

Accessing Data
--------------

Data is stored as a series of ``.tar`` files at the password-protected `s3df.slac.stanford.edu/data/kipac/symphony <s3df.slac.stanford.edu/data/kipac/symphony>`__. If you fill out `this form <https://docs.google.com/forms/d/e/1FAIpQLSdud6b4i51AP13glVibkzyLAtT9b2ctVx516_hvy5nm76uq1Q/viewform?usp=sf_link>`__, a username and password will be emailed to you. Once you have access, any method for downloading these ``.tar`` files will work, although we have several suggested methods below.

Each halo has a single ``.tar`` file per dataset, formatted as ``{suite}_{halo anme}_{dataset}`` so ``SymphonyLCluster_Halo_001_halos.tar`` contains the ``halos`` dataset for ``Halo_001`` in the LCluster suite. If you download all your ``.tar`` files into the same directory, they will automatically unpack into the correct directory structure expected by our analysis libraries by manually opening the ``.tar`` files or by running the following Python script in the data directory:

.. code-block:: python

   import tarfile, glob, os.path

   input_directory = "."
   output_directory = "."

   for file_name in glob.glob(os.path.join(input_directory, "*.tar")):
       f = tarfile.open(file_name)
       f.extractall(output_directory)
       f.close()

Download Methods
----------------

In principle, one could just manually download files from `s3df.slac.stanford.edu/data/kipac/symphony <s3df.slac.stanford.edu/data/kipac/symphony>`__ after getting password access, but we provide two suggested automated methods for downloading data. The first is to use some simple downloading library functions in our Symlib analysis library, and the second is to use an external tool named `rclone <https://rclone.org/>`__. rclone is a far better and more reliable approach than our library, especially when downloading the larger datasets, but it requires more setup.

**Downloading with Symlib**

First, install/update Symlib.

.. code-block:: console

	pip install symlib -U

You can use the :func:`symlib.download_files` function to download whichever datasets you want to use. The entire dataset will be stored in a single base directory. Each suite has is its own sub-directory within the base, and each zoom-in simulation has a subdirectory within its suite. 

The following Python code shows examples of how to use this function.

.. code-block:: python

	import symlib
	user = "my_user_name"
	password = "my_password"
	
	# The base directory where data will be downloaded to.
	data_dir = "path/to/storage/location"

	# The dataset you want to download.
	target = "halos"

	# Download the first host halo in the Milky Way-mass suite.
	symlib.download_files(user, password, "SymphonyMilkyWay", 0,
		data_dir, target=target)

	# Download all the host halos in the Milky Way-mass suite.
	symlib.download_files(user, password, "SymphonyMilkyWay", None,
		data_dir, target=target)

	# Download all the host halos across all the suites.
	symlib.download_files(user, password, None, None,
		data_dir, target=target)

You can also get a list of suite names with :func:`symlib.suite_names()` and host counts for a given suite with :func:`symlib.n_hosts()` so you can use a for loop instead of ``None``.

:func:`symlib.download_files()` is certainly capable of downloading all our datasets, but the chances of it encountering a network error that it cannot fix and recover from is quite large when attempting to download the three larger datasets. We recommend only relying on it if you want to download ``halos`` datasets or if you are unable to set up rclone on your computer.

**Downloading with rclone**

`rclone <https://rclone.org/>`__ is a cross-platform command-line tool which allows for easy file transfer between systems. Follow the instructions `here <https://rclone.org/install/>`__ to install. Once you have installed rclone, you will be asked to run

.. code-block:: console

	rclone config

to register the location of the simulation datasets and save your username and password. Run this command, then set up a "New remote" (``n``), call this remote ``symphony``, set the storage type of this remote to ``http``, and the url of this remote to ``https://xxxxx:yyyyyy@s3df.slac.stanford.edu/data/kipac/symphony``, except that ``xxxxx`` should be your username and ``yyyyy`` should be your password.

To test that you have set up rclone correctly, run

.. code-block:: console

	rclone ls symphony:

If successful, you will see a list of all the simulation ``.tar`` files printed to the screen. You can then download files using the command

.. code-block:: console
				
    rclone copy symphony: --include "target_file_name.tar" . --verbose

Here, ``"target_file_name.tar"`` can either be a single file or it can be a string with wildcards in it which will download each file that it matches. For example ``"SymphonyMilkyWay_*_halos.tar"`` would download the ``halos`` dataset for all the SymphonyMilkyWay hosts, and ``"*_*_particles.tar"`` would download the particle-tracking data for all suites and all halos. Remember to put quotes around the file name.

