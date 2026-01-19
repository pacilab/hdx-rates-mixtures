#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last modified on Jan 12th, 2026

Author: Antonio Grimaldi

"""

#%%

# Import required packages, modules, functions

import os
import numpy as np
import pandas as pd

from utils_pL import get_effective_acidity, get_ions_isotopic_abundance

#%%

file_dir = os.path.dirname(os.path.abspath(__file__))

file_kref = os.path.join(file_dir, 'kint_ref_values.csv')
file_table = os.path.join(file_dir, 'kint_Bai_factors.csv')

#%%

def read_seq(input_seq):
    """
    Parameters
    ----------
    input_seq : input sequence (str or path-like object)

    Returns
    -------
    seq : protein sequence (str)

    """
    try:
        with open(input_seq, 'r') as file:
            seq = file.readline().strip()
    except:
        seq = str(input_seq)
    
    return seq

#%%

# temperature dependence of intrinsic rates
    # A : acid-catalyzed
    # B : base-catalyzed
    # W : water-catalyzed
    
R = 1.987 # cal/(mol K)

Arrhenius = lambda T, T_ref, Ea : np.exp(-Ea * (1/T - 1/T_ref) / R)    

def get_temperature_factors(T, T_ref=293):
    
    Ea_A = 14000 # cal/mol
    Ea_B = 17000 # cal/mol
    Ea_W = 19000 # cal/mol
    
    fTA = Arrhenius(T, T_ref, Ea_A)
    fTB = Arrhenius(T, T_ref, Ea_B)
    fTW = Arrhenius(T, T_ref, Ea_W)
    
    return fTA, fTB, fTW

#%%

def get_kref(T, ref='PDLA'):
    """
    Parameters
    ----------
    T : temperature (in K)
    ref : reference values to be used ('PDLA' or '3Ala'). Default to 'PDLA'.  

    Returns
    -------
    kref : dataframe containing reference rates (/s) for (A) acid-, (B) base-,
        and (W) water-catalyzed exchange reactions at temperature T [1/s]

    """
    constants = pd.read_csv(file_kref, index_col=0)
    
    kref = pd.DataFrame(constants.loc[['HH', 'DH', 'HD', 'DD']])

    kref = 10 ** kref # rates
    kref /= 60 # convert from 1/min to 1/s
    
    # Temperature dependence
    
    kref *= get_temperature_factors(T, T_ref=293)
    kref = kref.round(6)
    
    # Adjust to the selected reference
    
    if ref == '3Ala':
        kref = kref * constants.loc['3Ala']
    elif ref == 'PDLA':
        pass
    
    return kref

#%%

get_pKc = lambda T, T_ref, Ea, cc : - np.log10(
    (10 ** cc) * Arrhenius(T, T_ref, Ea)
    ).round(6)

def get_pKc_Asp(T, T_ref=278):
    
    Ea_HD = 1000 # cal/mol
    Ea_DH =  960 # cal/mol
    
    c_HD = -4.48
    c_DH = -3.87
    
    return {'D' : get_pKc(T, T_ref, Ea_HD, c_HD),   # 'D' means in D2O
            'H' : get_pKc(T, T_ref, Ea_DH, c_DH)}   # 'H' means in H2O

def get_pKc_Glu(T, T_ref=278):
    
    Ea_HD = 1083 # cal/mol
    Ea_DH = 1083 # cal/mol
    
    c_HD = -4.93
    c_DH = -4.33
    
    return {'D' : get_pKc(T, T_ref, Ea_HD, c_HD),   # 'D' means in D2O
            'H' : get_pKc(T, T_ref, Ea_DH, c_DH)}   # 'H' means in H2O

def get_pKc_His(T, T_ref=278):
    
    Ea_HD = 7500 # cal/mol
    Ea_DH = 7500 # cal/mol
    
    c_HD = -7.42
    c_DH = -7
    
    return {'D' : get_pKc(T, T_ref, Ea_HD, c_HD),   # 'D' means in D2O
            'H' : get_pKc(T, T_ref, Ea_DH, c_DH)}   # 'H' means in H2O

def get_AL_CT(pKc_Glu,pL):#(T, pL):
    
    #pKc_Glu = get_pKc_Glu(T)
    
    log_den_CT = lambda pKc_Glu, pL : 10**(-pKc_Glu) + 10**(-pL)
    log_arg_CT = lambda pKc_Glu, pL : (10**(0.05-pL) + 10**(0.96-pKc_Glu))/log_den_CT(pKc_Glu, pL)
    
    return {'D' : np.log10(log_arg_CT(pKc_Glu['D'], pL)), # 'D' : in D2O
            'H' : np.log10(log_arg_CT(pKc_Glu['H'], pL))} # 'H' : in H2O

def get_Asp_line(table, pL, pKc_Asp):
    
    log_den = lambda pKc, pL : 10**(-pKc)+10**(-pL)
    log_num = lambda pKc, pL, table : 10**(table.loc['D+']-pL) + 10**(table.loc['D0']-pKc)
    
    log_arg = lambda pKc, pL, table : log_num(pKc,pL,table)/log_den(pKc,pL)
    
    return np.log10(log_arg(pKc_Asp, pL, table))

def get_Glu_line(table, pL, pKc_Glu):
    
    log_den = lambda pKc, pL : 10**(-pKc)+10**(-pL)
    log_num = lambda pKc, pL, table : 10**(table.loc['E+']-pL) + 10**(table.loc['E0']-pKc)
    
    log_arg = lambda pKc, pL, table : log_num(pKc,pL,table)/log_den(pKc,pL)
    
    return np.log10(log_arg(pKc_Glu, pL, table))

def get_His_line(table, pL, pKc_His):
    
    log_den = lambda pKc, pL : 10**(-pKc)+10**(-pL)
    log_num = lambda pKc, pL, table : 10**(table.loc['H+']-pL) + 10**(table.loc['H0']-pKc)
    
    log_arg = lambda pKc, pL, table : log_num(pKc,pL,table)/log_den(pKc,pL)
    
    return np.log10(log_arg(pKc_His, pL, table))

# to do: -NHMe

#%%

def make_Bai_tables(T, pL):
    
    isotope_pairs = ['HH', 'DH', 'HD', 'DD']

    tables = {}
    
    pKc_Asp = get_pKc_Asp(T)
    pKc_Glu = get_pKc_Glu(T)
    pKc_His = get_pKc_His(T)
    
    AL_CT = get_AL_CT(pKc_Glu, pL)
    
    for pair in isotope_pairs:
        
        tables[pair] = pd.read_csv(file_table, index_col=0)#.fillna(0)
        
        if pair.endswith('H'):
            key = 'H'
        elif pair.endswith('D'):
            key = 'D'
        
        # CT
        tables[pair].loc['CT', 'AL'] = AL_CT[key].round(2)
        # Asp
        tables[pair].loc['D'] = get_Asp_line(table=tables[pair],
                                             pL=pL,
                                             pKc_Asp=pKc_Asp[key]).round(2)
        # Glu
        tables[pair].loc['E'] = get_Glu_line(table=tables[pair],
                                             pL=pL,
                                             pKc_Glu=pKc_Glu[key]).round(2)
        # His
        tables[pair].loc['H'] = get_His_line(table=tables[pair],
                                             pL=pL,
                                             pKc_His=pKc_His[key]).round(2)
        
        
    return tables
        
#%%

# get primary sequence coefficients from Bai table 

def acid(resn, seq, table):
    
    res = seq[resn-1]
    
    if resn == 1:
        # multiply by 0 for first residue
        return 0   
    
    elif res == ('P' or 'Pc'):
        # multiply by 0 for first residue
        return 0
    
    else:
        resL = res
        resR = seq[resn-2]
        
        AL = table.loc[resL, 'AL']
        AR = table.loc[resR, 'AR']
        
        A = AL + AR
        
        if resn == 2:
            # N-terminal
            A += table.loc['NT', 'AR']
        elif resn == len(seq):
            # C-terminal
            A += table.loc['CT', 'AL']
            
        return 10 ** A

def base(resn, seq, table):
    
    res = seq[resn-1]
    
    if resn == 1:
        # multiply by 0 for first residue
        return 0   
    
    elif res == ('P' or 'Pc'):
        # multiply by 0 for first residue
        return 0
    
    else:
        resL = res
        resR = seq[resn-2]
        
        BL = table.loc[resL, 'BL']
        BR = table.loc[resR, 'BR']
        
        B = BL + BR
        
        if resn == 2:
            # N-terminal
            B += table.loc['NT', 'BR']
        elif resn == len(seq):
            # C-terminal
            B += table.loc['CT', 'BL']
            
        return 10 ** B

def calculate_res_Bai_factors(resn, seq, table):
    
    A = acid(resn, seq, table)
    B = base(resn, seq, table)
    
    return {'A' : A,
            'B' : B}
        
#%%

# calculate intrinsic rate of a residue (base-catalyzed only, pH > 4)

def calculate_kint(seq, T, x, pH_read, ref='PDLA', shift=0):
    """
    Parameters
    ----------
    seq : sequence (str)
    T : absolute temperature (K)
    x : isotopic abundance of D in solvent
    pH_read : measured pH (glass electrode) for the H2O-D2O mixture
    
    Returns
    -------
    kint : dataframe with columns
        index : residue index
        res : residue symbol (1-letter code)
        kforw : estimated intrinsic rate for forward exchange reaction
        kback : estimated intrisnic rate for back exchange reaction
        
    Unless otherwise specified, intrinsic rates are given in 1/s
    
    """
    
    isotope_pairs = ['HH', 'DH', 'HD', 'DD']
    
    # reference intrinsic rates, corrected for temperature
    kref = get_kref(T, ref=ref)
    
    # get mixture effective ionization constant and acidity
    mixture_params = get_effective_acidity(T, x, pH_read)
    pKW_mix = mixture_params['pKx']
    pL = mixture_params['pL']
    concentration_L = 10 ** (-pL)
    concentration_OL = 10 ** (pL - pKW_mix)

    # generate tables of Bai coefficients for all isotope combinations
    Bai_tables = make_Bai_tables(T, pL)
    
    # isotopic abundance of D+ in L+
    # [D+] = [L+] * fD, [H+] = [L+] * (1-fD)
    fD = get_ions_isotopic_abundance(x)['frac_L']
    # isotopic abundance of OD- in OL-
    # [OD-] = [OL-] * fOD, [OH-] = [OL-] * (1-fOD)
    fOD = get_ions_isotopic_abundance(x)['frac_OL']
    
    kint = pd.DataFrame(columns=['res', 'kforw', 'kback'])
            
    for resn in range(1, len(seq)+1):
        
        kforw, kback = 0, 0
        
        for pair in isotope_pairs:
            
            table = Bai_tables[pair]
            B = calculate_res_Bai_factors(resn, seq, table)['B']

            if pair == 'HH':
                # base
                kforw += x * kref.loc[pair, 'kB'] * B * concentration_OL * (1-fOD) # removal of H by OH-
            if pair == 'HD':
                # acid
                A = calculate_res_Bai_factors(resn, seq, table)['A']
                kforw += kref.loc[pair, 'kA'] * A * concentration_L * fD # protonation by D+
                # base
                kforw += x * kref.loc[pair, 'kB'] * B * concentration_OL * fOD # removal of H by OD-
                # water
                kforw += x * kref.loc[pair, 'kW'] * B
            if pair == 'DH':
                # acid
                A = calculate_res_Bai_factors(resn, seq, table)['A']
                kback += kref.loc[pair, 'kA'] * A * concentration_L * (1-fD) # protonation by H+
                # base
                kback += (1-x) * kref.loc[pair, 'kB'] * B * concentration_OL * (1-fOD) # removal by OH-
                # water
                kback += (1-x) * kref.loc[pair, 'kW'] * B
            if pair == 'DD':
                # base
                kback += (1-x) * kref.loc[pair, 'kB'] * B * concentration_OL * fOD # removal by OD-
        
        kint.loc[resn+shift] = (seq[resn-1], kforw, kback)

    return kint

#%%

def convert_rate_units(kint, units):
    if units == "h":
        factor = 3600
    elif units == "m":
        factor = 60
    elif units == "s":
        factor = 1
    else:
        print("Invalid time unit chosen. Rates given as by default - 1/s.")
        factor = 1
        
    kint[['kforw', 'kback']] *= factor
    
    return kint

#%%
