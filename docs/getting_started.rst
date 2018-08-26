===============
Getting Started
===============

This source code will let you identify all the eddies in the ocean,
but also It can be adapted for any other normal function in a 2D surface.

To get the code:
================
1.- Make a new directory where you want the repository.

2.- Clone the TrackEddy repository from Github. In the command prompt, type::

        git clone https://github.com/Josue-Martinez-Moreno/trackeddy.git
        cd trackeddy

or set up SSH keys for github::

        git clone git@github.com:Josue-Martinez-Moreno/trackeddy.git

3.- Install the package globally::

        pip install -e .

or ::

        pip install --force-reinstall -e.

or for a local installation::

        pip install --user -e .

Package structure
=================

Suggested structure for the use of this package after cloning it from github.::

    trackeddy
    ├── LICENSE
    ├── README.md
    ├── data.input    --> Simbolic link to inputs.
    ├── data.output   --> Simbolic link to output.
    ├── docs          --> Documentation.
    │   ├── README.md
    │   ├── about.rst
    │   ├── conf.py
    │   ├── getting_started.rst
    │   ├── images
    │   ├── index.rst
    │   ├── pages
    │   │   ├── Diagnostics.rst
    │   │   └── Tests.rst
    │   ├── references.rst
    │   ├── related_projects.rst
    │   └── using_trackeddy.rst
    ├── examples      --> Notebooks implementing some of the functions in
    |                     the package.
    │   ├── Eddies_geostrophic_velocity_field.ipynb
    │   ├── Eddies_ssha_satellite.ipynb
    │   ├── Eddies_velocity_field.ipynb
    │   ├── Eddies_vertical_profiles.ipynb
    │   ├── MULTIPLE_STEPS_eddies_southern_ocean.ipynb
    │   ├── ONE_STEP_eddies_southern_ocean.ipynb
    │   ├── Potential_Vorticity.ipynb
    │   ├── eddy_V1_bk.ipynb
    │   ├── eddyidentification_fitness.ipynb
    │   ├── eddyidentification_fitness_south_africa.ipynb
    │   ├── eddyidentification_fitness_specific_eddy.ipynb
    │   ├── multiprocess_bk.ipynb
    │   ├── potential_vorticity_bk.ipynb
    │   ├── save_eddy_data.ipynb
    │   ├── ssh_mean.ipynb
    │   ├── ssh_mean_global_bk.ipynb
    │   ├── test_geodesics_bk.ipynb
    │   ├── track_eddy_bk.ipynb
    │   ├── track_eddy_southern_ocean_bk.ipynb
    │   └── vorticity_tracking_eddy.ipynb
    ├── output        --> Figure output or small files
    ├── setup.py
    ├── src           --> Work in progress: The core of the trackeddy algoritm
    |                     will be coded in Fortran or C.
    ├── tests         --> Folder full of tests used to check the proper
    |                     extraction and analysis of eddies.
    │   ├── Centroid_eddy.ipynb
    │   ├── Synthetic_fields.ipynb
    │   ├── gaussian_fitting_multiple_eta_level.ipynb
    │   ├── gaussian_fitting_one_eta_level.ipynb
    │   ├── improving_time.ipynb
    │   ├── join_files_func.ipynb
    │   ├── time_datastruct.ipynb
    │   ├── trackeddy_okubo.ipynb
    │   ├── trackeddy_ssh.ipynb
    │   └── trackeddy_ssh_bk.ipynb
    └── trackeddy     --> Functions included in the package.
    ├── __init__.py
    ├── datastruct.py
    ├── geometryfunc.py
    ├── init.py
    ├── physics.py
    ├── plotfunc.py
    ├── printfunc.py
    ├── savedata.py
    └── tracking.py

Test the code
=============
The source code have been compiled and tested into the `Travis CI <https://travis-ci.org/Josue-Martinez-Moreno/trackeddy>`_ environment
(Check the build status on the `trackeddy GitHub <https://github.com/Josue-Martinez-Moreno/trackeddy>`_ ).

1.- Move to the test directory::

    cd /path2trackeddy/test/

2.- Run any of the scripts located in that folder::

    # Example:
    python test_2d_gaussian_one_level.py

.. note::
    If you want to display the diagnostics for each test, just replace:
    "diagnostics=False" by "diagnostics=True" at the beginning of the test file.
..

.. warning::
    The testing code it's in a early version, so please submit all the Issues
    to `trackeddy GitHub <https://github.com/Josue-Martinez-Moreno/trackeddy>`_.
..
