"""
Extract best ORF from each orf_group in the output from get_phasing.py and 
computes bio-type and orf-type specific FDRs
By default, all pseudogene biotypes are removed (this is can be avoided by the --keep-pseudo flag)

usage:
python get_results.py --phasing <ribofy_phasing.txt-file> --output <output-file>

"""

import sys
import os
import argparse
import pandas as pd
import numpy as np


def p_adjust_bh (p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    found here: https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
    """
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def get_results (phasing, output, keep_pseudo=False, p_method = "wilcox"):

    print ("### get_results ###")

    pcol = "p_" + p_method


    pd_phasing = pd.read_csv (phasing, sep="\t")


    # remove None
    pd_phasing = pd_phasing.dropna(subset=[pcol])

    if not keep_pseudo:
        pd_phasing = pd_phasing[~pd_phasing.bio_type.str.contains("pseudo")]
        
        
    # group by orf_groups - only keep the most significant orf for each orf_group
    #pd_collapse = pd_phasing.sort_values([pcol], ascending=True) \
    pd_collapse = pd_phasing.sort_values(["total_counts"], ascending=False) \
        .groupby(['orf_group']) \
        .head(1) 

    # calculate group-specific FDR
    pd_collapse['fdr_orftype'] = pd_collapse.groupby('orf_type')[pcol].transform (lambda x: p_adjust_bh(x.astype("float"))) #.reset_index(name='fdr')
    pd_collapse['fdr_biotype'] = pd_collapse.groupby('bio_type')[pcol].transform (lambda x: p_adjust_bh(x.astype("float"))) #.reset_index(name='fdr')
    
    pd_collapse['fdr'] = pd_collapse[pcol].transform (lambda x: p_adjust_bh(x.astype("float"))) #.reset_index(name='fdr')

    print (pd_collapse)
    
    pd_collapse.to_csv (output, sep="\t")

    print ("### Done ###")


def ribofy_results ():

    parser = argparse.ArgumentParser(description='collapse ORFs and calculate FDR')
    parser._action_groups.pop()

    # required
    required = parser.add_argument_group('required arguments')
    required.add_argument("--phasing", dest='phasing', required=True, help="phasing - generated by get_phasing.py")
    required.add_argument("--output", dest='output', default="ribofy.results.txt", help="output; default=ribofy.results.txt")

    # optional    
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--p_method", dest='p_method', default="wilcox", choices = ["wilcox", "binom"], help="statistics used for enrichment of phased reads; either wilcox or binom")
    optional.add_argument("--keep_pseudo", dest='keep_pseudo', action='store_true', default=False, help="Keep pseudogenes in analysis?")

    args = parser.parse_args()

    get_results (args.phasing, args.output, args.keep_pseudo, args.p_method)


if __name__ == "__main__":

    ribofy_results ()
