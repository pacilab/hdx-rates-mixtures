#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last modified on Jun 5th, 2025

Author: Antonio Grimaldi

"""

#%%

import numpy as np

#%%

def get_effective_acidity(T, x, pH_read):
    """
    Parameters
    ----------
    T : tempertature (K)
    x : D isotopic abundance in solution
    pH_read : measured pH (glass electrode) for the H2O-D2O mixture
    
    Returns
    -------
    a dictionary of keywords
        pKW : dissociation constant of pure H2O at temperature T
        pKx : effective dissociation constant of the mixture (x,T)
        pL :  effective acidity  pL = - log10[L+], L={H,D}
    
    The dissociation constant of pure H2O is computed according to
    Sweeton, Fo H., R. E. Mesmer, and C. F. Baes. "Acidity measurements at elevated temperatures. VII. Dissociation of water." Journal of Solution Chemistry 3 (1974): 191-214.    
    using the relation
    
    pKW_H2O = 670.86 - 34691.6 / T - 105.15 ln(T) + 0.10757 * T + 2358000 / T^2
           := pKW(0)
    
    The effective dissociation constant is computed according to
    Gold, V., and B. M. Lowe. "Measurement of solvent isotope effects with the glass electrode. Part I. The ionic product of D 2 O and D 2 O–H 2 O mixtures." Journal of the Chemical Society A: Inorganic, Physical, Theoretical (1967): 936-943.
    using the relation
    
    pKW(x) = pKW(0) + delta_pKW(x),    where
    
    delta_pKW(x) = 0.7282 * x + 0.05122 * x**2 + 0.0826 * x**3
    
    The effective acidity       pL = - log10[L+] (L={H,D})
    is approximated assuming    pL(x) = pH_read * pKW(x) /pKW(0)
                                      = pH_read (1 + delta_pKW(x) / pKW(0))
    
    """
    # pure water dissociation constant (pKW)
    pKW = 670.86 - 34691.6 / T - 105.15 * np.log(T) + 0.10757 * T + 2358000 / (T**2)
    
    # mixture effective dissociation constant (pKx)
    delta_pKW = 0.7282 * x + 0.05122 * (x**2) + 0.0826 * (x**3)
    pKx = pKW + delta_pKW
    
    # mixture effective acidity (pL)
    pL = pH_read * pKx / pKW 
    
    return {'pKW' : pKW,
            'pKx' : pKx,
            'pL' : pL}

def get_ions_isotopic_abundance(x, phi_L = 0.69, phi_OL = 0.47):
    """
    Parameters
    ----------
    x : D atom fraction in solution

    Returns
    -------
    dictionary of keywords
        frac_L : isotopic abundance of D in L+ ions
        frac_OL : isotopic abundance of D in OL- ions
    calculated according to definitions given in
    Gold, V. "Protolytic processes in H2O-D2O mixtures." Advances in Physical Organic Chemistry. Vol. 7. Academic Press, 1969. 259-331.

    """
    
    return {'frac_L' : x * phi_L / (1 - x + x * phi_L),
            'frac_OL' : x * phi_OL / (1 - x + x * phi_OL)}

#%%

def get_pH_read_from_pL(T, x, pL):
    """
    Parameters
    ----------
    T : tempertature (K)
    x : D isotopic abundance in solution
    pL : effective acidity  pL = - log10[L+], L={H,D}

    Returns
    -------
    a dictionary of keywords
        pKW : dissociation constant of pure H2O at temperature T
        pKx : effective dissociation constant of the mixture (x,T)
        pH_read : measured pH (glass electrode) for the H2O-D2O mixture
    """

    pKW =  670.86 - 34691.6 / T - 105.15 * np.log(T) + 0.10757 * T + 2358000 / (T**2)
    delta_pKW = 0.7282 * x + 0.05122 * (x**2) + 0.0826 * (x**3)
    pKx = pKW + delta_pKW

    pH_read = pL * pKW / pKx

    return pH_read
