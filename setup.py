from setuptools import setup, find_packages

setup(
    name='trackeddy',
    version='0.1',
    description='Tracking of eddies in the ocean or any other gaussian like function in a 3D [f(t,x,y)] space.',
    url='https://github.com/Josue-Martinez-Moreno/eddy_identification',
    author='Josue Martinez Moreno',
    author_email='josue.martinezmoreno@anu.edu.au',
    license='MIT License',
    packages=find_packages(),
    install_requires=['gsw','seawater'],
    packages=['trackeddy'],
    ext_package='trackeddy', 
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'],
 )
