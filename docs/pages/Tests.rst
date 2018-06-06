==========
Test Suite
==========

Some tests are provided to check the stability of the code and it's reproducibility.

Inside the folder ``tests`` are two available tests:

- ``*.py`` which correspond to the tests used in Travis CI and you can run them using::

        pytest -m testme

- ``*.ipynb`` which correspond to the same tests used in Travis CI, but displaying more information and diagnostics. Inside the trackeddy directory run::
        
        jupyther notebook

then, move into the ``tests`` folder and select the test you want to check.

Tests description
=================
.. warning::
  This section continues under development!
..

Now, a detailed description of each test is provided.

test_2d_gaussian_one_level
--------------------------

test_2d_gaussian_multiple_level
-------------------------------

test_ellipsoid_fitting
----------------------

test_timestep_tracking
----------------------

test_deta_tracking
------------------


