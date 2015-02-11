commah
=======
[![Build Status](https://travis-ci.org/astroduff/commah.svg?branch=master)](https://travis-ci.org/astroduff/commah)

[commah](https://github.com/astroduff/commah) is an analytic formalism
based on Correa et al. 2015a,b,c) to calculate a range of physical
properties at any redshift for dark matter haloes of any mass using
cosmological parameters such as mass and dark energy densities.

Written in python this uses cosmological calculations based on cosmology
and routines in numpy and scipy to create a structured dataset with
mass accretion rates, halo mass history, NFW concentrations,
mass 'sigma' variance, fluctuation parameter and formation redshifts
for any halo mass at any redshift.

Note that implicitly this assumes halo virial mass is 200 times critical overdensity.

### Getting started

If you have `pip' then commah can be run using:
```
$ pip install commah
$ python
>>> import commah
>>> cosmo = 'WMAP9'
>>> output = commah.run(cosmo, zi=0., Mi=1e12, z=[0.,1.0,2.0], filename=cosmo+'_commah.txt')
```
### Examples

Several plot and command line usage examples are given in examples.py
and can be called using
```
$ python
>>> from commah import examples
>>> examples.runcommand()
>>> examples.plotcommand()
```
### Support and Contact

If you have trouble with commah or you have feature requests/suggestions please
send an email to (mail@alanrduffy.com) or twitter (astroduff)
