#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last modified on Jan 12th, 2026

Author: Antonio Grimaldi

"""

import os
import argparse
from utils_kint import calculate_kint,\
    read_seq, convert_rate_units

#%%
    
def float_01(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} not a floating-point literal")
    
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x

def run():
    # Set up argument parsing
    parser = argparse.ArgumentParser()
    # Required arguments
    parser.add_argument("--seq",
                        required=True,
                        help="Sequence (str or path-like object)"
                        )
    parser.add_argument("--pH",
                        required=True, type=float,
                        help='pH read (glass electrode) for the mixture'
                        )
    parser.add_argument("--temp",
                        required=True, type=float,
                        help='Absolute temperature (K)'
                        )
    parser.add_argument("--deut",
                        required=False, default=1.0, type=float_01,
                        help='solvent deuteration level; must be in [0,1]'
                        )
    parser.add_argument("--ref",
                        required=False, type=str, default="PDLA",
                        help='Reference for computation of rates; choose "PDLA" or "3Ala"')
    parser.add_argument("--time",
                        required=False, type=str, default="s",
                        help='Rate (inverse) units; choose "h" or "m" or "s"')
    parser.add_argument("--out",
                        required=False,  type=str,
                        help='Output file name (path-like object)'
                        )
    parser.add_argument("--shift",
                        required=False, type=int, default=0,
                        help='Shift in starting residue index'
                        )
                        
    #
    args = parser.parse_args()
    
    # Read sequence
    
    seq = read_seq(args.seq)
    
    # Set parameters
    
    pH_read = args.pH
    T = args.temp
    x = args.deut
    
    # Output file
    data_dir = os.getcwd()
    if args.out:
        out_file = os.path.join(data_dir, args.out)
    else:
        filename = "rates_%3iD_%3iK_pH%4i.csv" % (int(x*100), int(T), int(pH_read*100))
        out_file = os.path.join(data_dir, filename)
    
    kint = calculate_kint(seq=seq, T=T, x=x, pH_read=pH_read,
                          ref=args.ref, shift=args.shift)
    
    kint = convert_rate_units(kint=kint, units=args.time)
    
    kint[['kforw', 'kback']].round(6)
    
    kint.to_csv(out_file)
    
if __name__ == "__main__":
    run()
    
    
    
    