##!/usr/bin/env ipython
# -*- coding: utf-8 -*-

"""Routine for creating Mass Accretion Histories and NFW profiles."""

__author__ = 'Camila Correa and Alan Duffy'
__email__ = 'mail@alanrduffy.com'
__version__ = '0.1.0'

from scipy.integrate import ode, quad
from scipy.optimize import brent

import numpy as np
import cosmolopy as cp

def cosmotodict(cosmo=None):
  """ Convert astropy cosmology to cosmolopy cosmology dict """

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
      from comma import *
      Result = comma.delta_sigma(**cosmo)
  
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
         This function takes in a recognised string cosmology (such as WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck) 
         or a cosmo dict that requires the sigma_8, n_s 
   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy

   CALLING SEQUENCE:
      from comma import *
      Result = comma.getAscaling(cosmology, [newcosmo=False] )
  
   INPUTS:
         cosmology: Either a name for a cosmology, default WMAP7, such as WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck
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
    defaultcosmologies = {'wmap1' : 787.01, 'wmap3' : 850.37, 'wmap5' : 903.75, 'wmap9' : 820.37, 'planck' : 798.82} # Values from Correa 14a

    if cosmology.lower() in defaultcosmologies.keys():
      A_scaling = defaultcosmologies[cosmology.lower()]    
    else:
      print "Error, don't recognise your cosmology for A_scaling ", cosmology

  else:
    # Scale from default WMAP5 cosmology using Correa et al 14b eqn C1 
    Mpivot = 1e12
    A_scaling = defaultcosmologies['wmap5'] * delta_sigma(Mpivot, **cosmology)

  return A_scaling

def int_growth(z, **cosmo):
  """ Returns the integral of the linear growth factor from z=200 to z=z """
  zmax = 200.
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
      Result = comma.deriv_growth(z, **cosmo)
  
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
      Result = comma.growthfactor(z, norm=True, **cosmo)
  
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

def minimize_c(c, z0=0., M0=1e10, a_tilde=1., b_tilde=-1., bestfit_param = 900., omega_M_0=0.25, omega_lambda_0=0.75):
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

      Result = scipy.optimize.brent(minimize_c, brack=(0.1,100.), args=(z0,M0,a_tilde,b_tilde,bestfit_param,cosmo['omega_M_0'],cosmo['omega_lambda_0'])) 
  
   INPUTS:
         z0: Redshift of original halo
         M0: Mass of original halo at redshift z0
         bestfit_param: A scaling between cosmologies, use getAscaling function
         a_tilde: power law growth rate factor
         b_tilde: exponential growth rate factor
         cosmo: dictionary of cosmology, in particular two properties only Total Matter density as cosmo['omega_M_0'] and DE density as cosmo['omega_lambda_0']

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS (set to True)
         norm: Normalise growth factor at redshift 'z' by the growth factor at z=0

   OUTPUTS:
          A float 

   MODIFICATION HISTORY (by Alan Duffy):
          
          1/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """
  
  ## Fn 1:
  Y1 = np.log(2.) - 0.5
  Yc = np.log(1.+c) - c/(1.+c)
  f1 = np.log(Y1/Yc)

  ## Fn 2:
  rho_2 = 200.*(c**3.)*Y1/Yc
  zf = ( (rho_2/bestfit_param)/omega_M_0 - omega_lambda_0/omega_M_0 )**(1./3.) - 1.
  f2 = a_tilde * np.log(1.+zf-z0) + b_tilde*(zf-z0)

  return np.abs(f1-f2)

def calc_ab(z0, M0, **cosmo):
  """ Calculate parameters a_tilde and b_tilde from Eqns 9 and 10 of Correa et al 2014b """

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

      Result = scipy.optimize.brent(minimize_c, brack=(0.1,100.), args=(z0,M0,a_tilde,b_tilde,bestfit_param,cosmo['omega_M_0'],cosmo['omega_lambda_0'])) 
  
   INPUTS:
         z0: Redshift of original halo
         M0: Mass of original halo at redshift z0
         bestfit_param: A scaling between cosmologies, use getAscaling function
         a_tilde: power law growth rate factor
         b_tilde: exponential growth rate factor
         cosmo: dictionary of cosmology, in particular two properties only Total Matter density as cosmo['omega_M_0'] and DE density as cosmo['omega_lambda_0']

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS (set to True)
         norm: Normalise growth factor at redshift 'z' by the growth factor at z=0

   OUTPUTS:
          A float 

   MODIFICATION HISTORY (by Alan Duffy):
          
          1/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
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
  ## Uses parameters a_tilde and b_tilde following eqns 9 and 10 from Correa et al 2014b 
  assert(z0 < z);

  if a_tilde == None or b_tilde == None:
    a_tilde, b_tilde = calc_ab(z0, M0, **cosmo)

  ## Accretion rate at z
  Mz = np.log10(M0) + np.log10( (1.+z-z0)**a_tilde * np.exp(b_tilde *(z-z0)) )
  dMdt = np.log10(71.59*(cosmo['h']/0.7)*(-a_tilde - b_tilde*(1.+z-z0))*(10.**(Mz-12.))*np.sqrt(cosmo['omega_M_0']*(1.+z-z0)**3.+cosmo['omega_lambda_0']))

  return 10.**dMdt, 10.**Mz, a_tilde, b_tilde

def MAH(z, z0, M0, **cosmo):
  """ Compute Mass Accretion History at redshift 'z' given halo of mass M0 at redshift z0, with z0<z always """
  ## Uses parameters a_tilde and b_tilde following eqns 9 and 10 from Correa et al 2014b 
  dMdt, Mz = acc_rate(z, z0, M0, **cosmo)

  return Mz

def COM(z0, M0, a_tilde=None, b_tilde=None, **cosmo):
  """ Give a halo mass and redshift calculate the concentration based on equation 17 and 18 from Correa et al 2014b """

  if a_tilde == None or b_tilde == None:
    a_tilde, b_tilde = calc_ab(z0, M0, **cosmo)

  ## Minimize the scalar
  c0 = 3.

  ## Use detal to convert best fit constant of proportionality of rho_crit - rho_2 from Correa et al 2014a to this cosmology
  bestfit_param = 900. * delta_sigma(M0,**cosmo) / (cosmo['omega_M_0']*(1.+z0)**3. + cosmo['omega_lambda_0'])
#  print "bestfit_param is ",bestfit_param

  c = brent(minimize_c, brack=(0.1,100.), args=(z0,M0,a_tilde,b_tilde,bestfit_param,cosmo['omega_M_0'],cosmo['omega_lambda_0'])) 

  R0_Mass = cp.perturbation.mass_to_radius(M0, **cosmo) 
  sig0, err_sig0 = cp.perturbation.sigma_r(R0_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation
  nu = 1.686 / (sig0*growthfactor(z0, norm=True, **cosmo))

  return c, sig0, nu

def run(cosmology, com=True, mah=True):
  """ Determine whether user wants to work out Mass Accretion Histories (MAH) and/or NFW profiles """

  """
  +
   NAME:
         commah.run
  
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
      from comma import *
      Result = comma.run(cosmology [, mah=True, com=True] )
  
   INPUTS:
         cosmology: Either a name for a cosmology, default WMAP7, such as WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck
                or a dictionary like: 
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
  
   KEYWORD PARAMETERS (set to True)
         com:     Calculate the NFW properties for haloes in this cosmology (i.e. the c-M relation as a function of z)
         mah:     Creat the mass accretion histories for haloes of a given mass, M0, to redshift, z0
  
   OUTPUTS:
          A numpy pickle with com or mah output to disk, or dat file if requested, or (future version HDF5 output)
          Ultimately will be used to create look up file that c, M, Mdot, z can be requested on per halo basis
  
   RESTRICTIONS:

          
   PROCEDURE:
  
  
   EXAMPLE:
      

   MODIFICATION HISTORY (by Alan Duffy):
          
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


  defaultcosmologies = {'wmap1' : cp.cparam.DRAGONS(), 'wmap5' : cp.cparam.WMAP5_ML()}#, WMAP5, WMAP7, WMAP9, Planck}
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

  z0 = 0. ## doesn't matter
  print "z0, log10(M0), z, log10(Mz), log10(dMdt), c"

  ## Scan for Mvir(z=0) = 1e10, 1e11, 1e12, 1e13, 1e14 their formation history from 0 -> 5
  for M in np.arange(10.,14.,1.):
    M0 = 10.**M

    #zf = -0.0064*(np.log10(M0))**2.+0.02373*np.log10(M0)+1.8837
    #q = (10.**(0.6167))*zf**(-0.9476)

    #R0_Mass = cp.perturbation.mass_to_radius(M0, **cosmo) 
    #Rq_Mass = cp.perturbation.mass_to_radius(M0/q, **cosmo) 

    #sig0, err_sig0 = cp.perturbation.sigma_r(R0_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation
    #sigq, err_sigq = cp.perturbation.sigma_r(Rq_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation

    #f = (sigq**2. - sig0**2.)**(-0.5)

    for z in np.arange(0.001,5.,1.):
      dMdt, Mz, a_tilde, b_tilde = acc_rate(z, z0, M0, **cosmo)
      c, sig0, nu = COM(z, Mz, a_tilde=a_tilde, b_tilde=b_tilde, **cosmo)

      print z0, np.log10(M0), z, np.log10(Mz), np.log10(dMdt), c

    #print "M, Mh, zf, sig0**2, f"
    #print M, M0, zf, sig0**2., f
    
    ## Assumes omega_Lamda_0 = 1-omega_M_0 TODO fix to include extra parameter
   # Dz0 = cp.perturbation.fgrowth(z, cosmo['omega_M_0'])
   # print Dz0, growthfactor(z, **cosmo)

#  for z in np.arange(0.,10.,0.1):
#    print cp.perturbation.fgrowth(z, cosmo['omega_M_0']) / growthfactor(z, **cosmo)
#    print z, deriv_growth(z,**cosmo)

  ## Acretion rate at z for object of mass M when at zi



#  M = 1e14
#  zi = 0.
#  z = 2.

#  print M, zi, z, acc_rate(M,zi,z, **cosmo)
## Romberg             return romberg(self.CosmoTint, self.EPSCOSMO, dExp**1.5, tol=0.01*self.EPSCOSMO)



