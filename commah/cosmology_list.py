""" Some pre-defined sets of cosmological parameters (e.g. from WMAP)
    copied and expanded from the cosmolopy list that's no longer updated.
"""


def add_extras(cosmo):
    """Sets neutrino number N_nu = 0, neutrino density
       omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.
    """
    extras = {'omega_n_0': 0.0,
              'N_nu': 0,
              'Y_He': 0.24}

    cosmo.update(extras)
    return cosmo


def DRAGONS(flat=False, extras=True):
    """DRAGONS cosmology assumes WMAP7 + BAO + H_0 mean from
    Komatsu et al. (2011) ApJS 192 18K (arxiv:1001.4538v1)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.2292
    omega_b_0 = 0.0458
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.725,
             'h': 0.702,
             'n': 0.963,
             'sigma_8': 0.816,
             'tau': 0.088,
             'z_reion': 10.6,
             't_0': 13.76,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP1_2dF_mean(flat=False, extras=True):
    """WMAP1 with 2dF and ACBAR results from
    Spergel et al. (2003) ApJS 148 175S (arXiv:astro-ph/0302209)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.206
    omega_b_0 = 0.044
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.75,
             'h': 0.73,
             'n': 0.97,
             'sigma_8': 0.9,
             'tau': 0.148,
             'z_reion': 17.,
             't_0': 13.7,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP1_Mill(flat=False, extras=True):
    """WMAP1 Millennium cosmology

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.206
    omega_b_0 = 0.044
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.75,
             'h': 0.73,
             'n': 1.0,
             'sigma_8': 0.9,
             'tau': 0.148,
             'z_reion': 17.,
             't_0': 13.7,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP3_mean(flat=False, extras=True):
    """WMAP3 Maximum Liklihood from Spergel et al. (2007) ApJS 170 377-408
    (arXiv:astro-ph/0603449)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.196
    omega_b_0 = 0.041
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.763,
             'h': 0.73,
             'n': 0.954,
             'sigma_8': 0.756,
             'tau': 0.091,
             'z_reion': 11.3,
             't_0': 13.73,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP3_ML(flat=False, extras=True):
    """WMAP3 Maximum Liklihood from Spergel et al. (2007) ApJS 170 377-408
    (arXiv:astro-ph/0603449)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.1959
    omega_b_0 = 0.0411
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.763,
             'h': 0.732,
             'n': 0.954,
             'sigma_8': 0.756,
             'tau': 0.091,
             'z_reion': 11.3,
             't_0': 13.73,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP5_BAO_SN_mean(flat=False, extras=True):
    """WMAP5 + BAO + SN parameters from Komatsu et al. (2009ApJS..180..330K).

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

    Notes
    -----

    From the abstract of the paper:

      The six parameters and the corresponding 68% uncertainties,
      derived from the WMAP data combined with the distance
      measurements from the Type Ia supernovae (SN) and the Baryon
      Acoustic Oscillations (BAO) in the distribution of galaxies,
      are:

      Omega_B h^2 = 0.02267+0.00058-0.00059,
      Omega_c h^2 = 0.1131 +/- 0.0034,
      Omega_Lambda = 0.726 +/- 0.015,
      n_s = 0.960 +/- 0.013,
      tau = 0.084 +/- 0.016, and
      Delata^2 R = (2.445 +/- 0.096) * 10^-9 at k = 0.002 Mpc^-1.

      From these, we derive

      sigma_8 = 0.812 +/- 0.026,
      H0 = 70.5 +/- 1.3 km s^-11 Mpc^-1,
      Omega_b = 0.0456 +/- 0.0015,
      Omega_c = 0.228 +/- 0.013,
      Omega_m h^2 = 0.1358 + 0.0037 - 0.0036,
      zreion = 10.9 +/- 1.4, and
      t0 = 13.72 +/- 0.12 Gyr

      """
    omega_c_0 = 0.228
    omega_b_0 = 0.0456
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.726,
             'h': 0.706,
             'n': 0.960,
             'sigma_8': 0.812,
             'tau': 0.084,
             'z_reion': 10.9,
             't_0': 13.72
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP5_ML(flat=False, extras=True):
    """WMAP5 parameters (using WMAP data alone) from Komatsu et
    al. (2009ApJS..180..330K).

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

    Notes
    -----

    Values taken from "WMAP 5 Year ML" column of Table 1 of the paper.

      """
    omega_c_0 = 0.206
    omega_b_0 = 0.043
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.751,
             'h': 0.724,
             'n': 0.961,
             'sigma_8': 0.787,
             'tau': 0.089,
             'z_reion': 11.2,
             't_0': 13.69
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP5_mean(flat=False, extras=True):
    """WMAP5 parameters (using WMAP data alone) from Komatsu et
    al. (2009ApJS..180..330K).

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

    Notes
    -----

    Values taken from "WMAP 5 Year Mean" of Table 1 of the paper.

    """
    omega_c_0 = 0.214
    omega_b_0 = 0.044
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.742,
             'h': 0.719,
             'n': 0.963,
             'sigma_8': 0.796,
             'tau': 0.087,
             'z_reion': 11.0,
             't_0': 13.69
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP7_ML(flat=False, extras=True):
    """WMAP7 ML parameters from Komatsu et al. (2011) ApJS 192 18K
    (arxiv:1001.4538v1)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.217
    omega_b_0 = 0.0445
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.738,
             'h': 0.714,
             'n': 0.969,
             'sigma_8': 0.803,
             'tau': 0.086,
             'z_reion': 10.3,
             't_0': 13.71,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP7_BAO_H0_mean(flat=False, extras=True):
    """WMAP7 + BAO + H_0 parameters from Komatsu et al. (2011) ApJS 192 18K
    (arxiv:1001.4538v1)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.227  # 0.228
    omega_b_0 = 0.0456  # 0.0456
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.728,  # 0.726,
             'h': 0.704,  # 0.706,
             'n': 0.963,  # 0.960,
             'sigma_8': 0.809,  # 0.812,
             'tau': 0.087,  # 0.084,
             'z_reion': 10.4,  # 10.9,
             't_0': 13.75,  # 13.72
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def WMAP9_ML(flat=False, extras=True):
    """WMAP Maximum Likelihood from Hinshaw et al. (2013) ApJS 208 19
    (arxiv:1212.5226v3)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.235
    omega_b_0 = 0.0465
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.7185,
             'h': 0.693,
             'n': 0.971,
             'sigma_8': 0.820,
             'tau': 0.0851,
             'z_reion': 10.36,
             't_0': 13.76,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def Planck_2013(flat=False, extras=True):
    """Planck 2013 XVI: Scalar Perturbations only (Maximum Likelihood)
    from Ade et al. (2013) A&A 571 16 (arxiv:1303.5076)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_c_0 = 0.267
    omega_b_0 = 0.05
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': omega_b_0 + omega_c_0,
             'omega_lambda_0': 0.683,
             'h': 0.671,
             'n': 0.9624,
             'sigma_8': 0.82344,
             'tau': 0.0925,
             'z_reion': 11.35,
             't_0': 13.82,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo


def Planck_2015(flat=False, extras=True):
    """Planck 2015 XII: Cosmological parameters Table 4
    column Planck TT, TE, EE + lowP + lensing + ext
    from Ade et al. (2015) A&A in press (arxiv:1502.01589v1)

    Parameters
    ----------

    flat: boolean

      If True, sets omega_lambda_0 = 1 - omega_M_0 to ensure omega_k_0
      = 0 exactly. Also sets omega_k_0 = 0 explicitly.

    extras: boolean

      If True, sets neutrino number N_nu = 0, neutrino density
      omega_n_0 = 0.0, Helium mass fraction Y_He = 0.24.

      """
    omega_b_0 = 0.02230/(0.6774**2.)
    cosmo = {'omega_b_0': omega_b_0,
             'omega_M_0': 0.3089,
             'omega_lambda_0': 0.6911,
             'h': 0.6774,
             'n': 0.9667,
             'sigma_8': 0.8159,
             'tau': 0.066,
             'z_reion': 8.8,
             't_0': 13.799,
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        add_extras(cosmo)
    return cosmo
