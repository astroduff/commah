##!/usr/bin/env ipython
# -*- coding: utf-8 -*-

"""Routine for creating Mass Accretion Histories and NFW profiles."""

__author__ = 'Camila Correa and Alan Duffy'
__email__ = 'mail@alanrduffy.com'
__version__ = '0.2.0'

from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import RectBivariateSpline

import numpy as np
import cosmolopy as cp
import cosmology_list as cg
import itertools
import matplotlib.pylab as plt

try:
  import cPickle as pkl
except:
  import pickle as pkl

def getcosmo(cosmology):
  """ Find the cosmological parameters for user provided named cosmology using cosmology.py list """

  defaultcosmologies = {'dragons' : cg.DRAGONS(), 'wmap1' : cg.WMAP1_2dF_mean(), 
  'wmap3' : cg.WMAP3_ML(), 'wmap5' : cg.WMAP5_ML(), 'wmap7' : cg.WMAP7_ML(), 
  'wmap9' : cg.WMAP9_ML(), 'planck' : cg.Planck_2013()}
  if cosmology.lower() in defaultcosmologies.keys():
    cosmo = defaultcosmologies[cosmology.lower()]
    A_scaling = getAscaling(cosmology.lower())
    cosmo.update({'A_scaling':A_scaling})
  elif isinstance(cosmology,dict):
    cosmo = cosmology
    if ('a_scaling').lower() not in cosmology.keys():
      ## Just assume WMAP5
      A_scaling = getAscaling(cosmology, newcosmo=True)
  else:
    print "You haven't passed a dict of cosmological parameters OR a recognised cosmology, you gave ",cosmology
  cosmo = cp.distance.set_omega_k_0(cosmo) ## No idea why this has to be done by hand but should be O_k = 0
  ## Use the cosmology as **cosmo passed to cosmolopy routines
  return cosmo

def cosmotodict(cosmo=None):
  """ Convert astropy cosmology to cosmolopy cosmology dict, still to do! Currently directly using cosmology.py list """

#def inv_h(z, **cosmo):
#  """ Inverse Hubble factor to be integrated """
#  return (1.+z)/(cosmo['omega_M_0']*(1.+z)**3.+cosmo['omega_lambda_0'])**(1.5)

def delta_sigma(**cosmo):
  """ Calculate delta_sigma (propto formation time) to convert best-fit parameter of rho_crit - rho_2 relation between cosmologies, WMAP5 standard (Correa et al 2014a) """

  """
  +
   NAME:
         delta_sigma
  
   PURPOSE:
         Perturb the rho_crit - rho_2 getAscaling result from the WMAP5 default to the user supplied cosmology

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy

   CALLING SEQUENCE:
      from commah import *
      Result = delta_sigma(**cosmo)
  
   INPUTS:
         cosmo: dictionary of cosmology
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS:

   OUTPUTS:
          A float 

   MODIFICATION HISTORY (by Alan Duffy):
          
          1/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """


  M8_cosmo = cp.perturbation.radius_to_mass(8.,**cosmo)
  return (0.796/cosmo['sigma_8'])*(M8_cosmo/2.5e14)**( (cosmo['n']-0.963)/6. )

def getAscaling(cosmology, newcosmo=False):
  """ Returns the normalisation constant between Rho_-2 and Rho_mean(z_formation) for a given cosmology """

  """
  +
   NAME:
         getAscaling
  
   PURPOSE:
         This function takes in a recognised string cosmology (such as DRAGONS, WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck) 
         or a cosmo dict that requires the sigma_8, n_s 
   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy

   CALLING SEQUENCE:
      from comma import *
      Result = getAscaling(cosmology, [newcosmo=False] )
  
   INPUTS:
         cosmology: Either a name for a cosmology, default WMAP7 (aka DRAGONS), such as DRAGONS, WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck
                or a dictionary like: 
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS (set to False):
         newcosmo: If false then cosmology is a string for a default cosmology & use Correa14b estimate, if true then recalculate

   OUTPUTS:
          A float that can scale between rho_mean at the formation redshift to rho_2 for NFW halo

   MODIFICATION HISTORY (by Alan Duffy):
          
          1/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """



  if newcosmo == False:
    defaultcosmologies = {'dragons' : 850., 'wmap1' : 787.01, 'wmap3' : 850.37, 
    'wmap5' : 903.75, 'wmap9' : 820.37, 'planck' : 798.82} # Values from Correa 14a

    if cosmology.lower() in defaultcosmologies.keys():
      A_scaling = defaultcosmologies[cosmology.lower()]    
    else:
      print "Error, don't recognise your cosmology for A_scaling ", cosmology

  else:
    # Scale from default WMAP5 cosmology using Correa et al 14b eqn C1 
    A_scaling = defaultcosmologies['wmap5'] * delta_sigma(**cosmology)

  return A_scaling

def int_growth(z, **cosmo):
  """ Returns the integral of the linear growth factor from z=200 to z=z """
  zmax = 200.

  if hasattr(z, "__len__"):
    for zval in z:
      assert(zval < zmax);
  else:
    assert(z < zmax);

  inv_h3 = lambda z: (1.+z)/(cosmo['omega_M_0']*(1.+z)**3.+cosmo['omega_lambda_0'])**(1.5)
  y, yerr = quad(inv_h3, z, zmax)  
  return y

def deriv_growth(z, **cosmo):
  """ Returns the derivative of the linear growth factor for redshift and cosmo (0m, 0l) """

  """
  +
   NAME:
         deriv_growth
  
   PURPOSE:
         For a given redshift (and cosmology passed as kwargs) return the derivative of the linear growth factor
   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy

   CALLING SEQUENCE:
      from comma import *
      Result = deriv_growth(z, **cosmo)
  
   INPUTS:
         z: Redshift  
         cosmo: dictionary of cosmology
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS: 

   OUTPUTS:
          A float 

   MODIFICATION HISTORY (by Alan Duffy):
          
          1/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """

  inv_h = (cosmo['omega_M_0']*(1.+z)**3.+cosmo['omega_lambda_0'])**(-0.5)
  fz = (1.+z) * inv_h**3.

  return growthfactor(z,norm=True, **cosmo)*(inv_h**2.)*1.5*cosmo['omega_M_0']*(1.+z)**2. - fz*growthfactor(z,norm=True, **cosmo)/int_growth(z, **cosmo) 

def growthfactor(z, norm=True, **cosmo):
  """ Returns the linear growth factor at z=z, normalised to z=0 """

  """
  +
   NAME:
         growthfactor
  
   PURPOSE:
         For a given redshift (and cosmology passed as kwargs) return the linear growth factor

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy

   CALLING SEQUENCE:
      from comma import *
      Result = growthfactor(z, norm=True, **cosmo)
  
   INPUTS:
         z: Redshift  
         cosmo: dictionary of cosmology
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS (set to True)
         norm: Normalise growth factor at redshift 'z' by the growth factor at z=0

   OUTPUTS:
          A float 

   MODIFICATION HISTORY (by Alan Duffy):
          
          1/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """


  H = np.sqrt(cosmo['omega_M_0']*(1.+z)**3.+cosmo['omega_lambda_0'])
  growthval = H * int_growth(z, **cosmo) 
  if norm == True:
    growthval /= int_growth(0., **cosmo) 
  return growthval

def cduffy(z0, M0, vir='200crit', relaxed=True):
  """ Give a halo mass and redshift calculate the NFW concentration based on Duffy 08 Table 1 for relaxed / vir def """
  if vir == '200crit':
    if relaxed == True:
      params = [6.71, -0.091, -0.44]
    else:
      params = [5.71, -0.084, -0.47]
  elif vir == 'tophat':
    if relaxed == True:
      params = [9.23, -0.090, -0.69]
    else:
      params = [7.85, -0.081, -0.71]
  elif vir == '200mean':
    if relaxed == True:
      params = [11.93, -0.090, -0.99]
    else:
      params = [10.14, -0.081, -1.01]
  else:
    print "Didn't recognise the halo boundary definition provided ", vir

  return params[0] * ((M0/(2e12/0.72))**params[1]) * ((1.+z0)**params[2])

def minimize_c(c, z0=0., M0=1e10, a_tilde=1., b_tilde=-1., Ascaling = 900., omega_M_0=0.25, omega_lambda_0=0.75):
  """ Trial function to solve 2 equations 17 and 18 from Correa et al 2014b for 1 unknown, the concentration """

  """
  +
   NAME:
         minimize_c
  
   PURPOSE:
         Function to minimize to solution for the concentration of a halo mass M0 at z0, solving when f1 = f2 in Correa et al 2014b (eqns 17 & 18)

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import scipy

   CALLING SEQUENCE:
      from comma import *
      Result = scipy.optimize.brentq(minimize_c, 2., 1000., args=(z0,M0,a_tilde,b_tilde,cosmo['A_scaling'],cosmo['omega_M_0'],cosmo['omega_lambda_0'])) 
  
   INPUTS:
         z0: Redshift of original halo
         M0: Mass of original halo at redshift z0
         A_scaling: A scaling between cosmologies, set into cosmo dict using getAscaling function
         a_tilde: power law growth rate factor
         b_tilde: exponential growth rate factor
         cosmo: dictionary of cosmology, in particular two properties only Total Matter density as cosmo['omega_M_0'] and DE density as cosmo['omega_lambda_0']

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS (set to True):
         norm: Normalise growth factor at redshift 'z' by the growth factor at z=0

   OUTPUTS:
          A float 

   MODIFICATION HISTORY (by Alan Duffy):
          
          1/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """
  #############################################  
  ## Fn 1 (LHS of Eqn 19)
  #############################################  
  Y1 = np.log(2.) - 0.5
  Yc = np.log(1.+c) - c/(1.+c)
  f1 = Y1/Yc

  #############################################  
  ## Fn 2 (RHS of Eqn 19)
  #############################################  

  ## Eqn 14 - Define the mean inner density
  rho_2 = 200.*(c**3.)*Y1/Yc
  ## Eqn 17 rearranged to solve for Formation Redshift (essentially when universe had rho_2 density)
  zf = ( ((1.+z0)**3. + omega_lambda_0/omega_M_0) * (rho_2/Ascaling) - omega_lambda_0/omega_M_0)**(1./3.) - 1.
  ## RHS of Eqn 19
  f2 = ((1.+zf-z0)**a_tilde) * np.exp((zf-z0)*b_tilde)

  ## LHS - RHS should be zero for the correct concentration!  
  return f1-f2

def calc_ab(z0, M0, **cosmo):
  """ Calculate parameters a_tilde and b_tilde from Eqns 9 and 10 of Correa et al 2014b """

  """
  +
   NAME:
         calc_ab
  
   PURPOSE:
         Function to calculate the power-law and exponential index growth rates from Correa et al 2014b (eqns 9 and 10)

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import scipy

   CALLING SEQUENCE:
      from comma import *
      a_tilde, b_tilde = calc_ab(z0, M0, **cosmo)
  
   INPUTS:
         z0: Redshift of original halo
         M0: Mass of original halo at redshift z0
         cosmo: dictionary of cosmology

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS:

   OUTPUTS:
          Two floats
          a_tilde: power law growth rate factor
          b_tilde: exponential growth rate factor
 

   MODIFICATION HISTORY (by Alan Duffy):
          
          2/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """
  
  # When z0 = 0, the a_tilde becomes alpha and b_tilde becomes beta
  zf = -0.0064*(np.log10(M0))**2. + 0.02373*(np.log10(M0)) + 1.8837

  q = 10.**(0.6167)*zf**(-0.9476)

  R0_Mass = cp.perturbation.mass_to_radius(M0, **cosmo) 
  Rq_Mass = cp.perturbation.mass_to_radius(M0/q, **cosmo) 

  sig0, err_sig0 = cp.perturbation.sigma_r(R0_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation
  sigq, err_sigq = cp.perturbation.sigma_r(Rq_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation

  f = (sigq**2. - sig0**2.)**(-0.5)  
  ## Eqn 9 a_tilde is power law growth rate from Correa et al 2014b
  a_tilde = (np.sqrt(2./np.pi)*1.686*deriv_growth(z0, **cosmo)/ growthfactor(z0, norm=True, **cosmo)**2. + 1.)*f
  ## Eqn 10 b_tilde is exponential growth rate from Correa et al 2014b
  b_tilde = -f

  return a_tilde, b_tilde

def acc_rate(z, z0, M0, a_tilde=None, b_tilde=None, **cosmo):
  """ Compute Mass Accretion Rate at redshift 'z' given halo of mass M0 at redshift z0, with z0<z always """


  """
  +
   NAME:
         acc_rate
  
   PURPOSE:
         Function to calculate accretion rate using Correa at al 2014 methodology

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import scipy

   CALLING SEQUENCE:
      from comma import *
      dMdt, Mz, a_tilde, b_tilde = calc_ab(z, z0, M0 [, a_tilde=None, b_tilde=None], **cosmo)
  
   INPUTS:
         z: Redshift of halo at which acc_rate and mass is desired
         z0: Redshift of original halo (note that z0 < z always)
         M0: Mass of original halo at redshift z0 (Msol)
         cosmo: dictionary of cosmology

   OPTIONAL INPUTS (set to None):
         a_tilde: power law growth rate factor, if not provided will be calculated
         b_tilde: exponential growth rate factor, if not provided will be calculated
  
   KEYWORD PARAMETERS:

   OUTPUTS:
         4 floats
         dMdt: accretion rate of halo at redshift z (Msol Gyr^-1)
         Mz: mass of halo at redshift z
         a_tilde: power law growth rate factor
         b_tilde: exponential growth rate factor
 

   MODIFICATION HISTORY (by Alan Duffy):
          
          2/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """
  ## Uses parameters a_tilde and b_tilde following eqns 9 and 10 from Correa et al 2014b 

  if a_tilde == None or b_tilde == None:
    a_tilde, b_tilde = calc_ab(z0, M0, **cosmo)

  ## Accretion rate at z, Msol Gyr^-1
  Mz = np.log10(M0) + np.log10( (1.+z-z0)**a_tilde * np.exp(b_tilde *(z-z0)) )
  dMdt = 9. + np.log10(71.59*(cosmo['h']/0.7)*(-a_tilde - b_tilde*(1.+z-z0))*(10.**(Mz-12.))*np.sqrt(cosmo['omega_M_0']*(1.+z-z0)**3.+cosmo['omega_lambda_0']))

  return 10.**dMdt, 10.**Mz, a_tilde, b_tilde

def MAH(z, z0, M0, **cosmo):
  """ Compute Mass Accretion History at redshift 'z' given halo of mass M0 at redshift z0, with z0<z always """

  """
  +
   NAME:
         MAH
  
   PURPOSE:
         Function to compute mass accretion history using acc_rate

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import scipy

   CALLING SEQUENCE:
      from comma import *
      z_array, dMdt_array, Mz_array = MAH(z, z0, M0, **cosmo)
  
   INPUTS:
         z: Redshift of halo at which acc_rate and mass is desired
         z0: Redshift of original halo (note that z0 < z always)
         M0: Mass of original halo at redshift z0
         cosmo: dictionary of cosmology

   OPTIONAL INPUTS:

   KEYWORD PARAMETERS:

   OUTPUTS:
          Array showing full dMdt and Mz growth rate in deltaz steps
          z_array: redshit steps between z0 and z in steps of deltaz
          dMdt_array: accretion rate of halo from z0 to redshift z (Msol Gyr^-1)
          Mz_array: mass of halo growth from z0 to redshift z 

   MODIFICATION HISTORY (by Alan Duffy):
          
          2/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """

  ## Create a full array
  dMdt_array = np.empty(np.size(z))
  Mz_array = np.empty(np.size(z))

  for ival, zval in enumerate(z):
    if ival == 0:
      ## Uses parameters a_tilde and b_tilde following eqns 9 and 10 from Correa et al 2014b 
      dMdt, Mz, a_tilde, b_tilde = acc_rate(zval, z0, M0, **cosmo)
    else:
      dMdt, Mz, a_tilde, b_tilde = acc_rate(zval, z0, M0, a_tilde=a_tilde, b_tilde=b_tilde, **cosmo)

    dMdt_array[ival] = dMdt
    Mz_array[ival] = Mz

  return z, dMdt_array, Mz_array

def COMLOOP(z_array, Mz_array, a_tilde=None, b_tilde=None, **cosmo):
  """ Loop com call over a range of redshifts/masses """

  c_array = np.empty(np.size(z_array))
  sig0_array = np.empty(np.size(z_array))
  nu_array = np.empty(np.size(z_array))
  ival = 0
  for zval, mzval in itertools.izip(z_array, Mz_array):
    c, sig0, nu = COM(zval, mzval, a_tilde=None, b_tilde=None, **cosmo)
    c_array[ival] = c
    sig0_array[ival] = sig0
    nu_array[ival] = nu
    ival += 1 

  return c_array, sig0_array, nu_array

def COM(z0, M0, a_tilde=None, b_tilde=None, **cosmo):
  """ Given a halo mass and redshift calculate the concentration based on equation 17 and 18 from Correa et al 2014b """

  if a_tilde == None or b_tilde == None:
    a_tilde, b_tilde = calc_ab(z0, M0, **cosmo)
  
  ## Use delta to convert best fit constant of proportionality of rho_crit - rho_2 from Correa et al 2014a to this cosmology
  c = brentq(minimize_c, 2.,1000., args=(z0,M0,a_tilde,b_tilde,cosmo['A_scaling'],cosmo['omega_M_0'],cosmo['omega_lambda_0']))

  R0_Mass = cp.perturbation.mass_to_radius(M0, **cosmo) 

  sig0, err_sig0 = cp.perturbation.sigma_r(R0_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation
  nu = 1.686 / (sig0*growthfactor(z0, norm=True, **cosmo))

  return c, sig0, nu

def creategrid(cosmology, filename=None, zgrid = None, zstart=0., zend=50., deltaz=0.1, mgrid = None, mstart=1., mend=16., deltam=0.1, logm = True, com=True, mah=False):
  """ Call "run" over a grid of redshifts and masses and save to a pickle for later interrogation """

  """
  +
   NAME:
         creategrid
  
   PURPOSE:
         This function calls 'run' across a user defined redshift and mass range, this output is saved as a pickle
         file for future fast interpolation and plotting
  
   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import pickle
          import scipy

   CALLING SEQUENCE:
      from commah import *
      Result = creategrid(cosmology [, filename=False, zstart=0., zend=100., deltaz=0.1, logz=False, mstart=1., mend=1e14, deltam=0.1, logm=True])
  
   INPUTS:
         cosmology: Either a name for a cosmology, default WMAP7 (aka DRAGONS), such as DRAGONS, WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck
                or a dictionary like: 
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
         filename:  If passed then output to " filename+'_'+cosmology+'_COM'/'_MAH'+'.npz' "
         zgrid:     User provided redshift spacings (overwrites other redshift choices) if logarithmic set logz=True
         zstart:    Lowest redshift to consider (default is 0.)
         zend:      Highest redshift to consider (default is 100.)
         deltaz:    How fine the grid should be in redshift (default is 0.01)

         mgrid:     User provided halo mass spacings (overwrites other mass choices) if logarithmic set logm=True
         mstart:    Lowest halo mass to consider (default is 1. Msol)
         mend:      Highest halo mass to consider (default is 1e14 Msol)
         deltam:    How fine the grid should be in mass (default is 0.1)

   KEYWORD PARAMETERS:
         logm:      Whether mass is in logarithmic units (default True, i.e. logarithmic spacing)

         com:    Calculate the NFW properties for haloes in this cosmology (i.e. the c-M relation as a function of z). Set for default True
         mah:    Creat the mass accretion histories for haloes of a given mass, M0, to redshift, z0. Set for default False (it's more limited than com)

   OUTPUTS:
          A numpy pickle with com or mah output to disk if filename passed, or (future version HDF5 output)
          that you can access arrays from by using: 
          filein = np.load(filename)
          filein['c']
          Ultimately will be used to create look up file that c, M, Mdot, z can be requested on per halo basis
  
   RESTRICTIONS:

          
   PROCEDURE:
  
  
   EXAMPLE:
      
   MODIFICATION HISTORY (by Alan Duffy):
          8/01/15 Function created
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """

  ## Establish redshift range
  if zgrid == None:
    zgrid = np.arange(zstart,zend,deltaz)

  ## Establish mass range
  if mgrid == None:
    mgrid = np.arange(mstart,mend,deltam)
  if logm:
    mgrid = 10.**(mgrid)

  if com:
    dataset = np.zeros( (len(mgrid),len(zgrid)), dtype=[('dMdt',float),('c',float),('sig0',float),('nu',float)] )
  elif mah:
    dataset = np.zeros( (len(mgrid),len(zgrid)), dtype=[('dMdt',float)] )

  if com:
    for mind, mloop in enumerate(mgrid):
      print "Solve for Mhalo = ",mloop," Msol, which is mass step ",mind," of ",len(mgrid)
      z_array, dMdt_array, Mz_array, c, sig0, nu = run(cosmology, z=zgrid, z0=min(zgrid), M0=mloop, com=com, mah=mah)  
      dataset[mind,:]['dMdt'] = dMdt_array
      dataset[mind,:]['c'] = c
      dataset[mind,:]['sig0'] = sig0
      dataset[mind,:]['nu'] = nu
  elif mah:
    for mind, mloop in enumerate(mgrid):
      print "Solve for Mhalo = ",mloop," Msol, which is mass step ",mind," of ",len(mgrid)
      z_array, dMdt_array, Mz_array = run(cosmology, z=zgrid, z0=min(zgrid), M0=mloop, com=com, mah=mah)  
      dataset[mind,:]['dMdt'] = dMdt_array

  if filename:
    output = filename+'_'
  else:
    output = ''

  if com:
    output += cosmology+'_COM'+'.pkl'
  elif mah:
    output += cosmology+'_MAH'+'.pkl'

  print "Output to ",output
#  np.savez(output, dataset, M=mgrid, z=zgrid)
  with open(output,'wb') as fout:
    pkl.dump(dict(dataset=dataset, Mhalo=mgrid, Redshift=zgrid), fout)
  
  return "Done"

def run(cosmology, z0=0., M0=1e12, z=[0.,1.,2.,3.,4.,5.], com=True, mah=False):
  """ run commah on a given halo mass 'M0' at a redshift 'z0' at higher redshifts 'z' """

  """
  +
   NAME:
         run
  
   PURPOSE:
         This function determines whether a user wants the mass accretion histories / NFW profiles,
         loads the correct cosmology, and creates lookup tables if required. It will then ensure 
         outputs are in the correct directory
  
   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy

   CALLING SEQUENCE:
      from commah import *
      Result = run(cosmology [z0=0., M0=1e12, mah=True, com=True] )
  
   INPUTS:
         cosmology: Either a name for a cosmology, default WMAP7 (aka DRAGONS), such as DRAGONS, WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck
                or a dictionary like: 
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
         M0:     The mass of the halo that the user wants the accretion history for (can be array)
         z0:     The redshift of the halo when it has mass M0 (can be array)
         z:      Outputs at redshift 'z' (can be array)

   KEYWORD PARAMETERS (set to True)
         com:    Calculate the NFW properties for haloes in this cosmology (i.e. the c-M relation as a function of z) com requires mah
         mah:    Creat the mass accretion histories for haloes of a given mass, M0, to redshift, z0

   OUTPUTS:
          A numpy pickle with com or mah output to disk if filename passed, or (future version HDF5 output)
          that you can access arrays from by using: 
          filein = np.load(filename)
          filein['c']
          Ultimately will be used to create look up file that c, M, Mdot, z can be requested on per halo basis
  
   RESTRICTIONS:

          
   PROCEDURE:
  
  
   EXAMPLE:
      
   MODIFICATION HISTORY (by Alan Duffy):
          23/12/14 Added pickle interogation to the file
          28/10/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """


#########################################################
##
## Check user choices...
##
#########################################################

  if com == False:
    if mah == False:
      print "User has to choose com=True and / or mah=True "

  ## Get the cosmological parameters for the given cosmology
  cosmo = getcosmo(cosmology)

  if mah:
    z_array, dMdt, Mz = MAH(z, z0, M0, **cosmo)
    if np.allclose(z_array,z) != True:
      print len(z_array), len(z)  
    if com:
      ## Loop over the COM routine in redshift (and hence diminishing mass at that z) for the initial mass M0 at z0
      c, sig0, nu = COMLOOP(z_array, Mz, a_tilde=None, b_tilde=None, **cosmo)

      return z_array, dMdt, Mz, c, sig0, nu
    else:
      return z_array, dMdt, Mz
  elif com:
    z_array, dMdt_array, Mz_array = MAH(z, z0, M0, **cosmo)
    if np.allclose(z_array,z) != True:
      print len(z_array), len(z)  
    c, sig0, nu = COMLOOP(z_array, Mz_array, a_tilde=None, b_tilde=None, **cosmo)
  
    return z_array, dMdt_array, Mz_array, c, sig0, nu

def loadval(cosmology=None, filename=None, z=0., M=1e12, val='c'):
  """ Shortcut to interrogate commah datasets for user requested masses and redshifts """


  """
  +
   NAME:
         loadval
  
   PURPOSE:
         This function interpolates over commah pickle and returns concentration and / or accretion rate 
  
   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy

   CALLING SEQUENCE:
          from commah import *
          c = loadval(filename='Full_WMAP5_COM.pkl', z=[0], M=[1e10,1e11], val='c')
          or
          c = loadval(cosmology='WMAP5', z=[0], M=[1e10,1e11], val='c')
  
   INPUTS:
          cosmology: A named cosmology, default WMAP7 (aka DRAGONS), such as DRAGONS, WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck
                     that has a pickle file Full_WMAP7_COM.pkl for example
          filename: Pass an existing pickle file by hand
   OPTIONAL INPUTS:
          z: Array of redshifts for desired outputs
          M: Array of halo mass for desired outputs
          val:  Provide a str value to interrogate pickle file i.e. 'dMdt' for acc rate, 'c' for concentration, 
                'sig0' for sigma fluctuation, 'nu' for dimensionless rareness of fluctuation
   KEYWORD PARAMETERS: 

   OUTPUTS:
          An array of value specified by 'val' for the given redshift and mass range
  
   RESTRICTIONS:
          
   PROCEDURE:
  
   EXAMPLE:
      
   MODIFICATION HISTORY (by Alan Duffy):
          12/01/15 Created example plots
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """


  ## Load the file
  if filename and cosmology:
    print "Which files do you want to open? Currently requesting a named file and the provided cosmology files"
    return -1
  elif cosmology:    
    filename = 'Full_'+cosmology+'_COM.pkl'
  else:
    filename = filename

  with open(filename, 'rb') as fin:
      data = pkl.load(fin)
      dataset = data.get("dataset",[])
      Mhalo = data.get("Mhalo",[])
      Redshift = data.get("Redshift",[])

  val_list = ('c','dMdt','sig0','nu') 
  if val in valist:
    ## Interpolate the output from commah
    interp = RectBivariateSpline(Mhalo, Redshift, dataset[val])
    return interp(M,z)
  else:
    print "You requested val= ",val, " the choices are ",val_list
    return -1

def runexamples(cosmology='WMAP5'):
  """ Example interface commands """

  ## Return the WMAP5 cosmology concentration predicted for z=0 range of masses
  M = [1e8, 1e9, 1e10]
  z = 0.
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'c')

  ## Return the WMAP5 cosmology concentration predicted for MW mass (2e12 Msol) across redshift
  M = 2e12
  z = [0.,0.5,1.,1.5,2.,2.5]
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'c')

  ## Return the WMAP5 cosmology concentration and rarity of high-z cluster
  M = 2e14
  z = 6.
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'c')
  print "Fluctuation sigma of haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'sig0')
  print "Rarity for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'nu')    

  ## Return the WMAP5 cosmology accretion rate prediction for haloes at range of redshift and mass
  M = [1e8, 1e9, 1e10]
  z = [0.,0.5,1.,1.5,2.,2.5]
  print "Concentrations for haloes of mass ", M, " at z=",z
  print loadval(cosmology = cosmology, z=z, M=M, val = 'dMdt')

  return "Done"

def plotexamples(filename='Full_WMAP5_COM.pkl', plotout=None):
  """ Example plot commands """

  import matplotlib.transforms as mtransforms
  import matplotlib.pyplot as plt
  import matplotlib.patches as patches
  import matplotlib.cm as cm
  
  """ Example ways to interrogate the dataset and plot the commah output """

  ## Check that the file type is indeed a pickle
  if filename.rfind('.pkl') >= 0:
    print filename + " is a pickle file"
    ## Load the file

    with open(filename, 'rb') as fin:
      data = pkl.load(fin)
      dataset = data.get("dataset",[])
      Mhalo = data.get("Mhalo",[])
      Redshift = data.get("Redshift",[])

  if filename.rfind('_COM') >= 0:
    print filename + " has Concentration - Mass relation "
    com = True
    mah = False

  if filename.rfind('_MAH') >= 0:
    print filename + " has Mass Accretion and Halo Concentration history "
    com = False
    mah = True

  ## Plot specifics
  plt.rc('text', usetex=True)

  ## Plot the Concentration Mass Relation
  if com:
  ## Plot the c-M relation as a function of redshift
    xval = 'M'
    xarray = 10.**(np.arange(1.,15.,1.))
    yval = 'c'

    ## Interpolate the output from commah
    interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])
 
    ## Specify the redshift range
    zarray = np.arange(0.,5.,0.5)

    xtitle = r"Halo Mass h$^{-1}$ M$_{sol}$"
    ytitle = r"Concentration"
    linelabel = "z="

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(xtitle)            
    ax.set_ylabel(ytitle)   
    colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

    for zind, zval in enumerate(zarray):
      ## Interpolate to the desired output values
      yarray = interp(xarray,zarray[zind])

      ## Plot each line in turn with different colour   
      ax.plot(xarray, yarray, label=linelabel+str(zval), color=colors[zind],)
      ## Overplot the D08 predictions in black
      ax.plot(xarray, cduffy(zval, xarray), color="black")

    ax.set_xscale('log')
    ax.set_yscale('log')

    leg = ax.legend(loc=1)
    leg.get_frame().set_alpha(0) # this will make the box totally transparent
    leg.get_frame().set_edgecolor('white')
    for label in leg.get_texts():                
      label.set_fontsize('small') # the font size   
    for label in leg.get_lines():
      label.set_linewidth(4)  # the legend line width
  
    if plotout:
      fig.tight_layout(pad=0.2)
      print "Plotting to ",plotout+"_CM_relation.png"
      fig.savefig(plotout+"_CM_relation.png", dpi=fig.dpi*5) 
    else:
      plt.show()

## Plot the dMdt-M relation as a function of redshift
  xval = 'M'
  xarray = 10.**(np.arange(1.,15.,1.))
  yval = 'dMdt'

  ## Interpolate the output from commah
  interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])

  ## Specify the redshift range
  zarray = np.arange(0.,5.,0.5)

  xtitle = r"Halo Mass h$^{-1}$ M$_{sol}$"
  ytitle = r"Accretion Rate h$^{-1}$ M$_{sol}$ Gyr$^{-1}$"
  linelabel = "z="

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(xtitle)            
  ax.set_ylabel(ytitle)   
  colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

  for zind, zval in enumerate(zarray):
    ## Interpolate to the desired output values
    yarray = interp(xarray,zarray[zind])

    ## Plot each line in turn with different colour   
    ax.plot(xarray, yarray, label=linelabel+str(zval), color=colors[zind],)

  ax.set_xscale('log')
  ax.set_yscale('log')

  leg = ax.legend(loc=1)
  leg.get_frame().set_alpha(0) # this will make the box totally transparent
  leg.get_frame().set_edgecolor('white')
  for label in leg.get_texts():                
    label.set_fontsize('small') # the font size   
  for label in leg.get_lines():
    label.set_linewidth(4)  # the legend line width

  if plotout:
    fig.tight_layout(pad=0.2)
    print "Plotting to ",plotout+"_MAH_M_relation.png"
    fig.savefig(plotout+"_MAH_M_relation.png", dpi=fig.dpi*5) 
  else:
    plt.show()


## Plot the (dM/M)dt-M relation as a function of redshift
  xval = 'M'
  xarray = 10.**(np.arange(1.,15.,1.))
  yval = 'dMdt'

  ## Interpolate the output from commah
  interp = RectBivariateSpline(Mhalo, Redshift, dataset[yval])

  ## Specify the redshift range
  zarray = np.arange(0.,5.,0.5)

  xtitle = r"Halo Mass h$^{-1}$ M$_{sol}$"
  ytitle = r"Specific Accretion Rate Gyr$^{-1}$"
  linelabel = "z="

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_xlabel(xtitle)            
  ax.set_ylabel(ytitle)   
  colors = cm.rainbow(np.linspace(0, 1, len(zarray)))

  for zind, zval in enumerate(zarray):
    ## Interpolate to the desired output values
    yarray = interp(xarray,zarray[zind])/xarray

    ## Plot each line in turn with different colour   
    ax.plot(xarray, yarray, label=linelabel+str(zval), color=colors[zind],)

  ax.set_xscale('log')
  ax.set_yscale('log')

  leg = ax.legend(loc=1)
  leg.get_frame().set_alpha(0) # this will make the box totally transparent
  leg.get_frame().set_edgecolor('white')
  for label in leg.get_texts():                
    label.set_fontsize('small') # the font size   
  for label in leg.get_lines():
    label.set_linewidth(4)  # the legend line width

  if plotout:
    fig.tight_layout(pad=0.2)
    print "Plotting to ",plotout+"_MAH_M_relation.png"
    fig.savefig(plotout+"_MAH_M_relation.png", dpi=fig.dpi*5) 
  else:
    plt.show()

  return "Done"