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

# This izip routine is from itertools, saves user having to download extra package for this short routine
def izip(*iterables):
    # izip('ABCD', 'xy') --> Ax By
    iterators = map(iter, iterables)
    while iterators:
        yield tuple(map(next, iterators))

def checkinput(zi,Mi,z=False, verbose=None):
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

  ## Complex test for size / type of incoming array, just in case numpy / list given
  try:
    z.size
  except Exception, e:
    if not z:
    ## Didn't pass anything, set zi = z
     lenzout = 1   
    else:
    ## Passed something, not a numpy array, probably list
      if hasattr(z, "__len__"):
        lenzout = 0
        for test in z:
          lenzout += 1  
        z = np.array(z)        
  else:
    ## Passed an array, what size is it?
    if hasattr(z, "__len__"):
      lenzout = 0
      for test in z:
        lenzout += 1  
      z = np.array(z)      
    else:
    ## Just one entry, force as array
      lenzout = 1
      z = np.array([z]) ## Ensure itertools can work

  return zi, Mi, z, lenz, lenm, lenzout

def getcosmo(cosmology):
  """ Find the cosmological parameters for user provided named cosmology using cosmology.py list """

  defaultcosmologies = {'dragons' : cg.DRAGONS(), 'wmap1' : cg.WMAP1_Mill(), 
  'wmap3' : cg.WMAP3_ML(), 'wmap5' : cg.WMAP5_mean(), 'wmap7' : cg.WMAP7_ML(), 
  'wmap9' : cg.WMAP9_ML(), 'planck' : cg.Planck_2013()}
  if isinstance(cosmology,dict):
    ## User providing their own variables
    cosmo = cosmology
    if ('a_scaling').lower() not in cosmology.keys():
      A_scaling = getAscaling(cosmology, newcosmo=True)
      cosmo.update({'A_scaling':A_scaling})

      ## Add extra variables by hand that cosmolopy requires but that aren't used (set to zero)
      for paramnames in cg.WMAP5_mean().keys():
        if paramnames not in cosmology.keys():
          cosmo.update({paramnames:0.})
  elif cosmology.lower() in defaultcosmologies.keys():
    ## Load by name of cosmology instead
    cosmo = defaultcosmologies[cosmology.lower()]
    A_scaling = getAscaling(cosmology.lower())
    cosmo.update({'A_scaling':A_scaling})
  else:
    print "You haven't passed a dict of cosmological parameters OR a recognised cosmology, you gave ",cosmology
  cosmo = cp.distance.set_omega_k_0(cosmo) ## No idea why this has to be done by hand but should be O_k = 0

  ## Use the cosmology as **cosmo passed to cosmolopy routines
  return cosmo

def getcosmoheader(cosmo):
  """ Output the cosmology to a string for writing to file """

  cosmoheader = "# Cosmology (flat) Om:"+"{0:.3f}".format( cosmo['omega_M_0'] )+ \
    ", Ol:"+"{0:.3f}".format( cosmo['omega_lambda_0'] )+", h:"+"{0:.2f}".format( cosmo['h'] )+ \
    ", sigma8:"+"{0:.3f}".format( cosmo['sigma_8'] )+", ns:"+"{0:.2f}".format( cosmo['n'] )

  return cosmoheader

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
  
def delta_sigma(**cosmo):
  """ Perturb best-fit constant of proportionality Ascaling for rho_crit - rho_2 relation for unknown cosmology (Correa et al 2015c) """

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

  defaultcosmologies = {'dragons' : 887., 'wmap1' : 853., 'wmap3' : 850., 
    'wmap5' : 887., 'wmap7' : 887., 'wmap9' : 950., 'planck' : 880.} # Values from Correa 15c

  if newcosmo:
    # Scale from default WMAP5 cosmology using Correa et al 14b eqn C1 
    A_scaling = defaultcosmologies['wmap5'] * delta_sigma(**cosmology)
  else:
    if cosmology.lower() in defaultcosmologies.keys():
      A_scaling = defaultcosmologies[cosmology.lower()]    
    else:
      print "Error, don't recognise your cosmology for A_scaling ", cosmology

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

def minimize_c(c, z=0., a_tilde=1., b_tilde=-1., Ascaling = 900., omega_M_0=0.25, omega_lambda_0=0.75):
  """ Trial function to solve 2 equations 18 from Correa et al 2015c for 1 unknown, i.e. the concentration """

  """
  +
   NAME:
         minimize_c
  
   PURPOSE:
         Function to minimize to solution for the concentration of a halo mass M at z, solving when f1 = f2 in Correa et al 2015c (eqns 19 and 20)

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
         a_tilde: power law growth rate factor (set by M,z)
         b_tilde: exponential growth rate factor (set by M,z)
         A_scaling: A scaling between cosmologies, set into cosmo['A_scaling']
         omega_M_0: Total Matter density as cosmo['omega_M_0'] 
         omega_lambda_0: Total DE density as cosmo['omega_lambda_0']

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
  ## Fn 1 (LHS of Eqn 18)
  #############################################  
  Y1 = np.log(2.) - 0.5
  Yc = np.log(1.+c) - c/(1.+c)
  f1 = Y1/Yc

  #############################################  
  ## Fn 2 (RHS of Eqn 18)
  #############################################  

  ## Eqn 14 - Define the mean inner density
  rho_2 = 200.*(c**3.)*Y1/Yc
  ## Eqn 17 rearranged to solve for Formation Redshift (essentially when universe had rho_2 density)
  zf = ( ((1.+z)**3. + omega_lambda_0/omega_M_0) * (rho_2/Ascaling) - omega_lambda_0/omega_M_0)**(1./3.) - 1.
  ## RHS of Eqn 19
  f2 = ((1.+zf-z)**a_tilde) * np.exp((zf-z)*b_tilde)

  ## LHS - RHS should be zero for the correct concentration!  
  return f1-f2

def formationz(c, z, Ascaling = 900., omega_M_0=0.25, omega_lambda_0=0.75):
  """ Rearrange eqn 18 from Correa et al 2015c to return zf for concentration 'c' at redshift 'z0' """

  Y1 = np.log(2.) - 0.5
  Yc = np.log(1.+c) - c/(1.+c)
  rho_2 = 200.*(c**3.)*Y1/Yc

  return ( ((1.+z)**3. + omega_lambda_0/omega_M_0) * (rho_2/Ascaling) - omega_lambda_0/omega_M_0)**(1./3.) - 1.

def calc_ab(zi, Mi, **cosmo):
  """ Calculate growth rate indices a_tilde and b_tilde  """

  """
  +
   NAME:
         calc_ab
  
   PURPOSE:
         Function to calculate the power-law and exponential index growth rates from Correa et al 2015a (eqns 19-23)
         generalised to arbitrary redshift in Correa et al 2015c (eqn 9 and 10)

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
         Mi: Mass of halo at redshift zi [Msol]
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

  ## Eqn 23 of Correa et al 2015a (analytically solve from Eqn 16 and 17) 
  ## Arbitray formation redshift, z_-2 in COM is more physically motivated
  zf = -0.0064*(np.log10(Mi))**2. + 0.0237*(np.log10(Mi)) + 1.8837

  ## Eqn 22 of Correa et al 2015a
  q = 4.137*zf**(-0.9476)

  ## Radius of a mass Mi
  R_Mass = cp.perturbation.mass_to_radius(Mi, **cosmo) # [Mpc] 
  ## Radius of a mass Mi/q
  Rq_Mass = cp.perturbation.mass_to_radius(Mi/q, **cosmo) # [Mpc]

  ## Mass variance 'sigma'
  sig, err_sig = cp.perturbation.sigma_r(R_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation
  sigq, err_sigq = cp.perturbation.sigma_r(Rq_Mass, 0., **cosmo) ## evalulate at z=0 to a good approximation

  f = (sigq**2. - sig**2.)**(-0.5)  
  
  ## Eqn 9 and 10 from Correa et al 2015c (generalised to zi from Correa et al 2015a's z=0 special case)
  # a_tilde is power law growth rate 
  a_tilde = (np.sqrt(2./np.pi)*1.686*deriv_growth(zi, **cosmo)/ growthfactor(zi, norm=True, **cosmo)**2. + 1.)*f
  ## b_tilde is exponential growth rate 
  b_tilde = -f

  return a_tilde, b_tilde

def acc_rate(z, zi, Mi, **cosmo):
  """ Compute Mass Accretion Rate at redshift 'z' given halo of mass Mi at redshift zi, with zi<z always """


  """
  +
   NAME:
         acc_rate
  
   PURPOSE:
         Function to calculate accretion rate using Correa at al 2015a (arbitrary redshift generalisation in 2015c)

   CATEGORY:
         Function
  
   REQUIREMENTS:
          import numpy
          import cosmolopy
          import scipy

   CALLING SEQUENCE:
      from comma import *
      dMdt, Mz = arr_rate(z, zi, Mi, **cosmo)
  
   INPUTS:
         z: Redshift of halo at which acc_rate and mass is desired
         zi: Redshift of original halo (note that zi < z always)
         Mi: Mass of original halo at redshift zi [Msol]
         cosmo: dictionary of cosmology

   OPTIONAL INPUTS:

   KEYWORD PARAMETERS:

   OUTPUTS:
         4 floats
         dMdt: accretion rate of halo at redshift z (Msol yr^-1)
         Mz: mass of halo at redshift z

   MODIFICATION HISTORY (by Alan Duffy):
          
          2/12/14 IDL version by Camila Correa translated to Python by Alan Duffy
          Any issues please contact Alan Duffy on mail@alanrduffy.com or (preferred) twitter @astroduff
  """
  ## Uses parameters a_tilde and b_tilde following eqns 9 and 10 from Correa et al 2015c
  a_tilde, b_tilde = calc_ab(zi, Mi, **cosmo)

  ## Halo mass at z, Msol Eqn 8 in Correa et al. 2015c
  Mz = Mi * (1.+z-zi)**a_tilde * np.exp(b_tilde *(z-zi))

  ## Accretion rate at z, Msol yr^-1 in Eqn 11 from Correa et al. 2015c
  dMdt = 71.588*(cosmo['h']/0.7)*(-a_tilde - b_tilde*(1.+z-zi))*(Mz/1e12)*np.sqrt(cosmo['omega_M_0']*(1.+z-zi)**3.+cosmo['omega_lambda_0'])

  return dMdt, Mz

def MAH(z, zi, Mi, **cosmo):
  """ Compute Mass Accretion History at redshift 'z' given halo of mass Mi [Msol] at redshift zi, with zi<z always """

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
         Mi: Mass of original halo at redshift zi [Msol]
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
  """ Given a halo mass [Msol] and redshift calculate the concentration based on equation 19 and 20 from Correa et al 2015c """

  ## Check that z and M are arrays
  if not hasattr(z, "__len__"):
  	z = np.array([z])
  if not hasattr(M, "__len__"):
  	M = np.array([M])

  ## Create array
  c_array = np.empty(np.size(z))
  sig_array = np.empty(np.size(z))
  nu_array = np.empty(np.size(z))
  zf_array = np.empty(np.size(z))

  i_ind = 0
  for zval, Mval in izip(z, M):
    ## Evaluate the indices at each redshift and mass combination that you want a concentration for, different to MAH which uses
    ## one a_tilde and b_tilde at the starting redshift
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
      ## Calculate formation redshift for this concentration, redshift at which the scale radius = virial radius: z_-2 
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

def run(cosmology, zi=0., Mi=1e12, z=False, com=True, mah=True, verbose=None, filename=None):
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
                {'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,'omega_b_0': 0.0458,'omega_lambda_0': 0.725,'sigma_8': 0.816}

   OPTIONAL INPUTS:
         Mi:     The mass of the halo that the user wants the accretion history for (can be array) [Msol]
         zi:     The redshift of the halo when it has mass Mi (can be array) Note that zi \le z always
         z:      Outputs at redshift 'z' (can be array) if left as None then z=zi
         verbose: Whether to output comments, set to None by default

   KEYWORD PARAMETERS (set to True):
         com:   If set to True return concentration (c), mass variance (sig), fluctuation height (nu) and formation redshift (zf)
         mah:   If set to True, return accretion rate (dM/dt) [Msol / yr] and halo mass at redshift 'z' (Mz) [Msol] for object mass Mi at redshift zi
   OUTPUTS:
          A data structure that you can access the desired values with keywords: concentration (c), mass variance (sig), 
                 fluctuation height (nu), formation redshift (zf), accretion rate (dM/dt) [Msol / yr] and halo mass at redshift 'z' (Mz))

   RESTRICTIONS:
          Only correct for halo definition of 200xCrit (due to EPS formalism)
          
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
      return -1

  ## Convert arrays / lists to np.array and inflate redshift / mass axis to match each other
  zi, Mi, z, lenz, lenm, lenzout = checkinput(zi,Mi,z=z,verbose=verbose)
  ## At this point we will have lenm objects to iterate over, finding lenzout solutions

  ## Get the cosmological parameters for the given cosmology
  cosmo = getcosmo(cosmology)

#########################################################
##
## Create output dataset data structure
##
#########################################################

  if filename:
    print "Output to file ",filename
    fout = open(filename, 'wb')

  if mah and com:
    if verbose:
      print "Output requested is zi, Mi, z, dMdt, Mz, c, sig, nu, zf"
    if filename:
      if verbose:
        print "Output to file is z, c, Mz, sig, nu, dMdt"
      fout.write(getcosmoheader(cosmo)+'\n')
      fout.write("# Initial z - Initial Halo  - Output z - concentration -   Final Halo    -   Mass    - Peak    - Accretion"+'\n')
      fout.write("#           -   mass        -          -               -    mass         - Variance  - Height  -   rate"+'\n')
      fout.write("#           -   (M200)      -          -     (c200)    -   (M200)        - (sigma)   - (nu)    -  (dM200/dt)"+'\n')
      fout.write("#           -   [Msol]      -          -               -   [Msol]        -           -         -  [Msol/yr]"+'\n')
    dataset = np.zeros( (lenm,lenzout), dtype=[('zi',float),('Mi',float),('z',float),('dMdt',float),('Mz',float),('c',float),('sig',float),('nu',float),('zf',float)] )
  elif mah:
    if verbose:
      print "Output requested is zi, Mi, z, dMdt, Mz"
    if filename:
      if verbose:
        print "Output to file is z, Mz, dMdt"
      fout.write(getcosmoheader(cosmo)+'\n')
      fout.write("# Initial z - Initial Halo  - Output z -   Final Halo - Accretion"+'\n')
      fout.write("#           -   mass        -          -   mass       -   rate"+'\n')
      fout.write("#           -   (M200)      -          -   (M200)     -  (dM200/dt)"+'\n')
      fout.write("#           -   [Msol]      -          -   [Msol]     -  [Msol/yr]"+'\n')
    dataset = np.zeros( (lenm,lenzout), dtype=[('zi',float),('Mi',float),('z',float),('dMdt',float),('Mz',float)] )
  else:
    if verbose:    
      print "Output requested is zi, Mi, c, sig, nu, zf"
    if filename:
      if verbose:
        print "Output to file is z, c, Mz, sig, nu"
      fout.write(getcosmoheader(cosmo)+'\n')
      fout.write("# z - concentration -   Halo    -   Mass    - Peak"+'\n')
      fout.write("#   -               -   mass    - Variance  - Height"+'\n')
      fout.write("#   -     (c200)    -   (M200)  - (sigma)   - (nu)"+'\n')
      fout.write("#   -               -   [Msol]  -           -     "+'\n')
    dataset = np.zeros( (lenm,lenzout), dtype=[('zi',float),('Mi',float),('z',float),('c',float),('sig',float),('nu',float),('zf',float)] )

  i_ind = 0
  for zval, Mval in izip(zi, Mi):    
    ## Now run the script for each zval and Mval combination
    if verbose:
      print "Output Halo of Mass Mi=",Mval," at zi=",zval    
    ## For a given halo mass Mi at redshift zi need to know the output redshifts 'z'
    ## Check that all requested redshifts are greater than the input redshift, except
    ## if z is None, in which case only solve z at zi, i.e. remove a loop
    try:
       z.size
    except Exception, e:
      if not z:
      ## Didn't pass anything, set as redshift
        ztemp = np.array([zval])
      else:
      ## Passed list perhaps?
        ztemp = np.array(z[z >= zval])
    else:
      ztemp = z
 
    ## Loop over the output redshifts
    if ztemp.size > 0:  
      ## Return accretion rates and halo mass progenitors at redshifts 'z' for object of mass Mi at zi
      dMdt, Mz = MAH(ztemp, zval, Mval, **cosmo) 
      if com:
        ## More computational intensive to return concentrations
        c, sig, nu, zf = COM(ztemp, Mz, **cosmo)
      if mah and com:
        ## Save all arrays
        for j_ind, j_val in enumerate(ztemp):
          dataset[i_ind,j_ind] = zval, Mval, ztemp[j_ind], dMdt[j_ind], Mz[j_ind], c[j_ind], sig[j_ind], nu[j_ind], zf[j_ind]
          fout.write("{}, {}, {}, {}, {}, {}, {}, {} \n".format(zval, Mval, ztemp[j_ind], c[j_ind], Mz[j_ind], sig[j_ind], nu[j_ind], dMdt[j_ind]))
      elif mah:
        ## Save only MAH arrays
        for j_ind, j_val in enumerate(ztemp):
          dataset[i_ind,j_ind] = zval, Mval, ztemp[j_ind], dMdt[j_ind], Mz[j_ind]
          fout.write("{}, {}, {}, {}, {} \n".format(zval, Mval, ztemp[j_ind], Mz[j_ind], dMdt[j_ind]))
      else:
        c, sig, nu, zf = COM(zval, Mval, **cosmo) 
      ## For any halo mass Mi at redshift zi solve for c, sig, nu and zf
        dataset[i_ind,:] = zval, Mval, zval, c, sig, nu, zf
        fout.write("{}, {}, {}, {}, {} \n".format(zval, c, Mval, sig, nu))

      i_ind += 1

  if filename:
    fout.close()
    return dataset        
  else:
    return dataset
