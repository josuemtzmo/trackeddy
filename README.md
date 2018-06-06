# TrackEddy (Eddy Identification Algorithm)


| Travis CI (Python 3.6) | Read the Docs | Code Coverage |
|:----------------------:|:-------------:|:-------------:|
| [![Build Status](https://travis-ci.org/Josue-Martinez-Moreno/trackeddy.svg?branch=develop)](https://travis-ci.org/Josue-Martinez-Moreno/trackeddy) | [![Documentation Status](https://readthedocs.org/projects/trackeddy/badge/?version=latest)](http://trackeddy.readthedocs.io/en/latest/?badge=latest) | [![codecov](https://codecov.io/gh/Josue-Martinez-Moreno/trackeddy/branch/develop/graph/badge.svg)](https://codecov.io/gh/Josue-Martinez-Moreno/trackeddy) |


**Supports Python 2.7 and Python 3**
This code will let you identify and track any normal shape in a 2D space. The principal objective of this development is finding and tracking all the eddies in the ocean. 

![Alt Text](https://github.com/Josue-Martinez-Moreno/trackeddy/blob/master/output/eddyn_13.gif "Eddy trajectory in the Souther Ocean")

# Beta Fase Release:
Now the source code supports the extraction of eddies each time step and moving in Z. Some bugs need to be fixed but right now it can be implemented in some data analysis.

![Alt Text](https://github.com/Josue-Martinez-Moreno/trackeddy/blob/develop/output/eke.png "Decomposition of the Kinetic energy in the Southern Ocean [Data provided by Adele Morrison].")

## To get the code

1. Make a new directory where you want the repository.
1. Clone the TrackEddy repository from Github. In the command prompt, type:
`git clone https://github.com/Josue-Martinez-Moreno/trackeddy.git`
1. Install the package globally:
`pip install -e .`
This make the package an editable install so that it can be updated with future additions to TrackEddy. To instead install the package locally:
`pip install --user .`


## To update the code later

1. Move into your TrackEddy directory.
1. Update your GitHub repository.
`git pull`
1. Edit your install of TrackEddy.
`pip install -e .` 
or
`pip install --force-reinstall -e .`
or, for local installation: 
`pip install --ignore-installed --user .`


## To test the code

Work in progress!

## To learn more about the module TrackEddy

Work in progress!

### References:
* Faghmous, J. H., Frenger, I., Yao, Y., Warmka, R., Lindell, A., & Kumar, V. (2015). A daily global mesoscale ocean eddy dataset from satellite altimetry. Scientific data, 2.
* Chang, Y. L., & Oey, L. Y. (2014). Analysis of STCC eddies using the Okubo–Weiss parameter on model and satellite data. Ocean Dynamics, 64(2), 259-271.
* Chelton, D. B., Schlax, M. G., Samelson, R. M., & de Szoeke, R. A. (2007). Global observations of large oceanic eddies. Geophysical Research Letters, 34(15).
Work in progress!
