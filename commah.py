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

def checkinput(zi,Mi,z=None, verbose=None):
  """ Ensure user input be it scalar, array or numpy array don't break the code """

  ## How many halo redshifts provided?
  if hasattr(zi, "__len__"):
    lenz = 0
    for test in zi:
      lenz += 1
    zi = np.array(zi)
  else:
    lenz = 1
    zi = np.array([zi]) ## Ensure itertools can work
  
  ## How many halo masses provided?
  if hasattr(Mi, "__len__"):
    lenm = 0
    for test in Mi:
      lenm += 1  
    Mi = np.array(Mi)
  else:
    lenm = 1
    Mi = np.array([Mi]) ## Ensure itertools can work

  ## How many redshifts to output?
  if z:
    if hasattr(z, "__len__"):
      lenzout = 0
      for test in z:
        lenzout += 1  
      z = np.array(z)        
    else:
      lenzout = 1
      z = np.array([z]) ## Ensure itertools can work
  else:
  ## Just output at the input redshifts, in which case Mz = Mi for example
    lenzout = lenz
    z = zi

  ## Check the input sizes for zi and Mi make sense, if not then exit unless 
  ## one axis is length one, in which case replicate values to the size of the other
  if (lenz > 1) & (lenm > 1):
    if lenz != lenm:
      print "Error ambiguous request"
      print "Need individual redshifts for all haloes provided or all haloes at same redshift "
      return -1
  elif (lenz == 1) & (lenm > 1):
    if verbose:
      print "Assume zi is the same for all Mi halo masses provided"
    ## Replicate redshift for all halo masses
    tmpz = zi
    zi = np.empty(lenm)
    if hasattr(tmpz, "__len__") == False:
      zi.fill(tmpz)
    else:
      zi.fill(tmpz[0])
    lenz=lenm
  elif (lenm == 1) & (lenz > 1):
    if verbose:
      print "Assume Mi halo masses are the same for all zi provided"
    ## Replicate redshift for all halo masses
    tmpm = Mi
    Mi = np.empty(lenz)
    if hasattr(tmpm, "__len__") == False:
      Mi.fill(tmpm)
    else:
      Mi.fill(tmpm[0])      
    lenm=lenz
  else:
    if verbose:
      print "A single Mi and zi provided"
    if hasattr(zi, "__len__") == False:
      zi = np.array(zi)
    if hasattr(Mi, "__len__") == False:
      Mi = np.array(Mi)

  return zi, Mi, z, lenz, lenm, lenzout

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
          import scipy

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
          import scipy

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
    'wmap5' : 903.75, 'wmap7' : 850., 'wmap9' : 820.37, 'planck' : 798.82} # Values from Correa 14a

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
          import scipy

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
          import scipy

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

def minimize_c(c, z=0., a_tilde=1., b_tilde=-1., Ascaling = 900., omega_M_0=0.25, omega_lambda_0=0.75, halodef=200.):
  """ Trial function to solve 2 equations 17 and 18 from Correa et al 2014b for 1 unknown, the concentration """

  """
  +
   NAME:
         minimize_c
  
   PURPOSE:
         Function to minimize to solution for the concentration of a halo mass M at z, solving when f1 = f2 in Correa et al 2014b (eqns 17 & 18)

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import scipy

   CALLING SEQUENCE:
      from comma import *
      Result = scipy.optimize.brentq(minimize_c, 2., 1000., args=(z,a_tilde,b_tilde,cosmo['A_scaling'],cosmo['omega_M_0'],cosmo['omega_lambda_0'])) 
  
   INPUTS:
         z: Redshift of halo in question
         M: Mass of halo at redshift z
         A_scaling: A scaling between cosmologies, set into cosmo dict using getAscaling function
         a_tilde: power law growth rate factor (set by M,z)
         b_tilde: exponential growth rate factor (set by M,z)
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
  rho_2 = halodef*(c**3.)*Y1/Yc
  ## Eqn 17 rearranged to solve for Formation Redshift (essentially when universe had rho_2 density)
  zf = ( ((1.+z)**3. + omega_lambda_0/omega_M_0) * (rho_2/Ascaling) - omega_lambda_0/omega_M_0)**(1./3.) - 1.
  ## RHS of Eqn 19
  f2 = ((1.+zf-z)**a_tilde) * np.exp((zf-z)*b_tilde)

  ## LHS - RHS should be zero for the correct concentration!  
  return f1-f2

def formationz(c, z, Ascaling = 900., omega_M_0=0.25, omega_lambda_0=0.75, halodef=200.):
  """ Use equations 17 from Correa et al 2015a to return zf for concentration 'c' at redshift 'z0' """

  Y1 = np.log(2.) - 0.5
  Yc = np.log(1.+c) - c/(1.+c)
  rho_2 = halodef*(c**3.)*Y1/Yc

  return ( ((1.+z)**3. + omega_lambda_0/omega_M_0) * (rho_2/Ascaling) - omega_lambda_0/omega_M_0)**(1./3.) - 1.


def calc_ab(zi, Mi, **cosmo):
  """ Calculate parameters a_tilde and b_tilde from Eqns 9 and 10 of Correa et al 2015a """

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
      a_tilde, b_tilde = calc_ab(zi, Mi, **cosmo)
  
   INPUTS:
         zi: Redshift of halo
         Mi: Mass of halo at redshift zi
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
  
  # When zi = 0, the a_tilde becomes alpha and b_tilde becomes beta
  zf = -0.0064*(np.log10(Mi))**2. + 0.02373*(np.log10(Mi)) + 1.8837

  q = 10.**(0.6167)*zf**(-0.9476)

  R_Mass = cp.perturbation.mass_to_radius(Mi, **cosmo) 
  Rq_Mass = cp.perturbation.mass_to_radius(Mi/q, **cosmo) 

  sig, err_sig = cp.perturbation.sigma_r(R_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation
  sigq, err_sigq = cp.perturbation.sigma_r(Rq_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation

  f = (sigq**2. - sig**2.)**(-0.5)  
  
  ## Eqn 9 a_tilde is power law growth rate from Correa et al 2014b
  a_tilde = (np.sqrt(2./np.pi)*1.686*deriv_growth(zi, **cosmo)/ growthfactor(zi, norm=True, **cosmo)**2. + 1.)*f
  ## Eqn 10 b_tilde is exponential growth rate from Correa et al 2014b
  b_tilde = -f

  return a_tilde, b_tilde

def acc_rate(z, zi, Mi, **cosmo):
  """ Compute Mass Accretion Rate at redshift 'z' given halo of mass Mi at redshift zi, with zi<z always """


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
      dMdt, Mz, a_tilde, b_tilde = arr_rate(z, zi, Mi, **cosmo)
  
   INPUTS:
         z: Redshift of halo at which acc_rate and mass is desired
         zi: Redshift of original halo (note that zi < z always)
         Mi: Mass of original halo at redshift zi (Msol)
         cosmo: dictionary of cosmology

   OPTIONAL INPUTS:

   KEYWORD PARAMETERS:

   OUTPUTS:
         4 floats
         dMdt: accretion rate of halo at redshift z (Msol yr^-1)
         Mz: mass of halo at redshift z
         a_tilde: power law growth rate factor
         b_tilde: exponential growth rate factor
 

   MODIFICATION HISTORY (by Alan Duffy):
          
          2/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """
  ## Uses parameters a_tilde and b_tilde following eqns 9 and 10 from Correa et al 2014b 
  a_tilde, b_tilde = calc_ab(zi, Mi, **cosmo)

  ## Accretion rate at z, Msol yr^-1
  Mz = np.log10(Mi) + np.log10( (1.+z-zi)**a_tilde * np.exp(b_tilde *(z-zi)) )
  dMdt = np.log10(71.59*(cosmo['h']/0.7)*(-a_tilde - b_tilde*(1.+z-zi))*
    (10.**(Mz-12.))*np.sqrt(cosmo['omega_M_0']*(1.+z-zi)**3.+cosmo['omega_lambda_0']))

  return 10.**dMdt, 10.**Mz

def MAH(z, zi, Mi, **cosmo):
  """ Compute Mass Accretion History at redshift 'z' given halo of mass Mi at redshift zi, with zi<z always """

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
      dMdt_array, Mz_array = MAH(z, zi, Mi, **cosmo)
  
   INPUTS:
         z: Redshift of halo at which acc_rate and mass is desired
         zi: Redshift of original halo (note that zi < z always)
         Mi: Mass of original halo at redshift zi
         cosmo: dictionary of cosmology

   OPTIONAL INPUTS:

   KEYWORD PARAMETERS:

   OUTPUTS:
          dMdt_array: accretion rate of halo with mass Mi at redshift zi at redshifts z (Msol yr^-1)
          Mz_array: Mass of halo at redshift z when it was a halo with mass Mi at redshift zi (Msol)

   MODIFICATION HISTORY (by Alan Duffy):
          
          2/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """

  ## Create a full array
  dMdt_array = np.empty(np.size(z))
  Mz_array = np.empty(np.size(z))

  for i_ind, zval in enumerate(z):
    ## Uses parameters a_tilde and b_tilde following eqns 9 and 10 from Correa et al 2014b 
    dMdt, Mz = acc_rate(zval, zi, Mi, **cosmo)

    dMdt_array[i_ind] = dMdt
    Mz_array[i_ind] = Mz

  return dMdt_array, Mz_array

def COM(z, M, **cosmo):
  """ Given a halo mass and redshift calculate the concentration based on equation 17 and 18 from Correa et al 2014b """

  ## Create array
  c_array = np.empty(np.size(z))
  sig_array = np.empty(np.size(z))
  nu_array = np.empty(np.size(z))
  zf_array = np.empty(np.size(z))

  i_ind = 0
  for zval, Mval in itertools.izip(z, M):

    a_tilde, b_tilde = calc_ab(zval, Mval, **cosmo)
    
    ## Use delta to convert best fit constant of proportionality of rho_crit - rho_2 from Correa et al 2014a to this cosmology
    c = brentq(minimize_c, 2.,1000., args=(zval,a_tilde,b_tilde,cosmo['A_scaling'],cosmo['omega_M_0'],cosmo['omega_lambda_0']))

    if np.isclose(c,0.):
      print "Error solving for concentration at this redshift with (probably) too small a mass "
      c = -1.
      sig = -1.
      nu = -1.
      zf = -1.
    else:
      zf = formationz(c, zval, Ascaling = cosmo['A_scaling'], omega_M_0=cosmo['omega_M_0'], omega_lambda_0=cosmo['omega_lambda_0'])

      R_Mass = cp.perturbation.mass_to_radius(Mval, **cosmo) 

      sig, err_sig = cp.perturbation.sigma_r(R_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation
      nu = 1.686 / (sig*growthfactor(zval, norm=True, **cosmo))

    c_array[i_ind] = c
    sig_array[i_ind] = sig
    nu_array[i_ind] = nu
    zf_array[i_ind] = zf

    i_ind +=1

  return c_array, sig_array, nu_array, zf_array

def run(cosmology, zi=0., Mi=1e12, z=None, val=None, com=True, mah=True, verbose=None):
  """ run commah on a given halo mass 'Mi' at a redshift 'zi' solving for higher redshifts 'z' """

  """
  +
   NAME:
         run
  
   PURPOSE:
         Take user requested cosmology along with halo mass at a given redshift as well as desired output
         redshifts if provided to output concentration / sig / nu / formation redshift
  
   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import scipy

   CALLING SEQUENCE:
      from commah import *
      Result = run(cosmology [ zi=0., Mi=1e12, z=[0.,1.,2.] ] )
  
   INPUTS:
         cosmology: Either a name for a cosmology, default WMAP7 (aka DRAGONS), such as DRAGONS, WMAP1, WMAP3, WMAP5, WMAP7, WMAP9, Planck
                or a dictionary like: 
                {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,
                'omega_n_0': 0.0, 'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}

   OPTIONAL INPUTS:
         Mi:     The mass of the halo that the user wants the accretion history for (can be array)
         zi:     The redshift of the halo when it has mass Mi (can be array) Note that zi \le z always
         z:      Outputs at redshift 'z' (can be array) if left as None then z=zi
         val:    If specified only output the value requested (concentration (c), mass variance (sig), 
                 fluctuation height (nu), formation redshift (zf), accretion rate (dM/dt) [Msol / yr] 
                 and halo mass at redshift 'z' (Mz))
         verbose: Whether to output comments, set to None by default

   KEYWORD PARAMETERS (set to True):
         com:   If set to True return concentration (c), mass variance (sig), fluctuation height (nu) and formation redshift (zf)
         mah:   If set to True, return accretion rate (dM/dt) [Msol / yr] and halo mass at redshift 'z' (Mz) [Msol] for object mass Mi at redshift zi
   OUTPUTS:
          A numpy pickle with com or mah output to disk if filename passed, or (future version HDF5 output)
          that you can access arrays from by using: 
          filein = np.load(filename)
          filein['c']
          Ultimately will be used to create look up file that c, M, Mdot, z can be requested on per halo basis
  
   RESTRICTIONS:

          
   PROCEDURE:
  
  
   EXAMPLE:
        import commah
      
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

  ## Convert arrays / lists to np.array and inflate redshift / mass axis to match each other
  zi, Mi, z, lenz, lenm, lenzout = checkinput(zi,Mi,z=z,verbose=verbose)
  ## At this point we will have lenm objects to iterate over, finding lenzout solutions

  ## Get the cosmological parameters for the given cosmology
  cosmo = getcosmo(cosmology)

  if val:
    valname = val
  else:
    if mah and com:
      valname = ('z', 'dMdt', 'Mz', 'c', 'sig', 'nu', 'zf')

  if mah and com:
    if verbose:
      print "Output requested is zi, Mi, z, dMdt, Mz, c, sig, nu, zf"
    dataset = np.zeros( (lenm,lenzout), dtype=[('zi',float),('Mi',float),('z',float),('dMdt',float),('Mz',float),('c',float),('sig',float),('nu',float),('zf',float)] )
  elif mah:
    if verbose:
      print "Output requested is zi, Mi, z, dMdt, Mz"
    dataset = np.zeros( (lenm,lenzout), dtype=[('zi',float),('Mi',float),('z',float),('dMdt',float),('Mz',float)] )
  else:
    if verbose:    
      print "Output requested is zi, Mi, c, sig, nu, zf"
    dataset = np.zeros( (lenm,lenzout), dtype=[('zi',float),('Mi',float),('z',float),('c',float),('sig',float),('nu',float),('zf',float)] )

  i_ind = 0
  for zval, Mval in itertools.izip(zi, Mi):    
    ## Now run the script for each zval and Mval combination
    if verbose:
      print "Output Halo of Mass Mi=",Mval," at zi=",zval    
    if mah and com:
      ## For a given halo mass Mi at redshift zi need to know the output redshifts 'z'
      ## Check that all requested redshifts are greater than the input redshift
      ztemp = np.array(z[z >= zval])
      if ztemp.size > 0:  
        dMdt, Mz = MAH(ztemp, zval, Mval, **cosmo) 
        ## Return accretion rates and halo mass progenitors at redshifts 'z' for object of mass Mi at zi
        c, sig, nu, zf = COM(ztemp, Mz, **cosmo)
        for j_ind, j_val in enumerate(ztemp):
          dataset[i_ind,j_ind] = zval, Mval, ztemp[j_ind], dMdt[j_ind], Mz[j_ind], c[j_ind], sig[j_ind], nu[j_ind], zf[j_ind]
    elif mah:
      ## Check that all requested redshifts are greater than the input redshift
      ztemp = np.array(z[z >= zval])
      if ztemp.size > 0:  
        for j_ind, j_val in enumerate(ztemp):
          dataset[i_ind,j_ind] = zval, Mval, ztemp[j_ind], dMdt[j_ind], Mz[j_ind]
    else:
      ## For any halo mass Mi at redshift zi solve for c, sig, nu and zf
      dataset[i_ind,:] = zval, Mval, zval, COM(zval, Mval, **cosmo) 
    i_ind += 1

  return dataset
