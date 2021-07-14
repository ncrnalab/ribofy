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
    adapted from here: https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python to allow NaNs
    """
    p = np.asfarray(p)
    
    nna = ~np.isnan (p)
    q = np.empty ((len(p)))
    q[:] = np.nan
    pnna = p[nna]

    by_descend = pnna.argsort()[::-1]
    by_orig = by_descend.argsort()

    n = len(pnna) #[~np.isnan (p)])
    i = np.arange(len(pnna), 0, -1)
    q[nna] = np.minimum(1, np.fmin.accumulate((float (n)/i) * pnna[by_descend]))[by_orig]
    return q
        



def get_filtered_padj (s, pcol="p_glm", name="filtered_padj"):
    """Adapted from DESeq2; filtering by expression, the BH padjustment if performed solely on ORFs exceeding the expression threshold. 
    Then, the threshold that maximized number of rejections (i.e. significant ORFs) are used and the rest obtain NaN fdr values. The maximization is
    not based on lowess regression, but simply the cutoff with max rejections (lowess implementation TODO)
    """
    
    filter=np.array (s['n'])
    p=np.array(s[pcol])
    nrows = s.shape[0]

    if nrows < 50:        
        s[name] = p_adjust_bh(p)   
        return (s)

    lq = np.mean(filter == 0)
    uq = .95 if lq < .95 else 1

    r = np.array (np.linspace(start=lq, stop=uq, num=50))

    cutoffs = np.quantile (filter, r)

    result = np.empty((nrows,len(cutoffs)))
    result[:] = np.nan

    for i in range (len (cutoffs)):
        
        use = filter >= cutoffs[i]    
        
        if (np.any(use)):
            
            use_p = p[use]        
            result[use, i] = p_adjust_bh(use_p)        


    numRej = np.sum (result < 0.05, axis=0)
    j = np.argmax(numRej)

    s[name] = result[:,j]
    return (s)



def get_results (phasing, output, keep_pseudo=False, p_methods = ["glm"]):
    """Main function - loads the output from phasing analysis, removes pseudogenes (except if --keep_pseudo is set) and 
    performed multiple testing corrections 
    phasing: path/to/phasing file
    output: path/to/output file  
    keep_pseudo: flag to allow pseudogenes in the final output
    p_methods: the p methods to use in multiple testing correction, should be a combination if glm, binom, wilcox and taper, 
    and the chosen statistics should be been performed in the phasing analysis
    """

    print ("### get_results ###")

    pd_phasing = pd.read_csv (phasing, sep="\t")

    # remove pseudo genes    
    if not keep_pseudo:
        pd_phasing = pd_phasing[~pd_phasing.bio_type.str.contains("pseudo")]
        
    # group by orf_groups - only keep the most significant orf for each orf_group
    #pd_collapse = pd_phasing.sort_values([pcol], ascending=True) \
    pd_collapse = pd_phasing.sort_values(["total_counts"], ascending=False) \
        .groupby(['orf_group']) \
        .head(1) 

    # calculate group-specific FDR
    for p in p_methods:

        pcol = "p_" + p            
        #pd_collapse = pd_collapse.dropna(subset=[pcol])

        # combine types
        pd_collapse['type'] = pd_collapse[['bio_type', 'orf_type']].agg('_'.join, axis=1)

        # non-filtered padj, i.e. all detected ORFs with p-values are included
        pd_collapse['fdr_' + p] = pd_collapse[pcol].transform (lambda x: p_adjust_bh(x.astype("float"))) 
        pd_collapse['fdr_type_' + p] = pd_collapse.groupby('type')[pcol].transform (lambda x: p_adjust_bh(x.astype("float")))
        pd_collapse['fdr_orftype_' + p] = pd_collapse.groupby('orf_type')[pcol].transform (lambda x: p_adjust_bh(x.astype("float")))
        pd_collapse['fdr_biotype_' + p] = pd_collapse.groupby('bio_type')[pcol].transform (lambda x: p_adjust_bh(x.astype("float")))        
                
        # filtered padj
        pd_collapse = get_filtered_padj (pd_collapse, pcol=pcol, name="filtered_fdr_" + p)    
        pd_collapse = pd_collapse.groupby('type').apply (get_filtered_padj, pcol=pcol, name="filtered_fdr_type_" + p)
        pd_collapse = pd_collapse.groupby('bio_type').apply (get_filtered_padj, pcol=pcol, name="filtered_fdr_biotype_" + p)
        pd_collapse = pd_collapse.groupby('orf_type').apply (get_filtered_padj, pcol=pcol, name="filtered_fdr_orftype_" + p)
        

    print (pd_collapse)
    
    pd_collapse.to_csv (output, sep="\t")

    print ("### Done ###")


def ribofy_results ():

    parser = argparse.ArgumentParser(description='collapse ORFs and calculate FDR')
    parser.add_argument('results', nargs='?', help='') # dummy argument

    # required
    parser.add_argument("--phasing", dest='phasing', required=True, help="phasing - generated by get_phasing.py")
    parser.add_argument("--output", dest='output', default="ribofy.results.txt", help="output; default=ribofy.results.txt")

    # optional    
    parser.add_argument("--p_methods", dest='p_methods', nargs='*', default=["glm"], help="statistics used for enrichment of phased reads; either wilcox or binom")
    parser.add_argument("--keep_pseudo", dest='keep_pseudo', action='store_true', default=False, help="Keep pseudogenes in analysis?")

    args = parser.parse_args()

    get_results (args.phasing, args.output, args.keep_pseudo, args.p_methods)


if __name__ == "__main__":

    ribofy_results ()
