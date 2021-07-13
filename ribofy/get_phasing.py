"""
For each ORF, get_phasing determines P-site positions for all mapped reads and computes statistics for on-frame enrichment

usage:
python get_phasing.py --bam <bam-file> --orfs <ribofy orfs-file> --offsets <ribofy offsets-file> --output <output-file>

"""

import argparse
import pysam
import pandas as pd
import numpy as np
import warnings
from tqdm import tqdm
from scipy.stats import wilcoxon, binomtest, f
import statsmodels.formula.api as smf 
import statsmodels.api as sm
import statsmodels.tools.sm_exceptions as sme
from .bam_utils import get_tid_info

libmtspec = True
try:
    from mtspec import mtspec
except ModuleNotFoundError:
    libmtspec = False
    

"""
converts from 1D array to 2D matrix
"""
def get_phasing_matrix (psites):
    
    mat = np.reshape (psites, (int (len(psites)/3), 3))
    
    return (mat)


"""
counts number of p-sites in each frame
"""
def get_counts (psites):

    mat = get_phasing_matrix (psites)

    return ({
        'total' : np.sum(mat),
        'frame0' : np.sum(mat[:,0]),
        'frame1' : np.sum(mat[:,1]),
        'frame2' : np.sum(mat[:,2])
    })


def get_taper (psites, time_bandwidth = 3, ntapers = "default", nfft = "default"):
    """Performs multitaper analysis (as in ribotaper) with Ftest statistics for 1/3 frequency
    psites: 1D array (ORF length) with P-site counts ncodons
    returns: p-value
    """

    if sum (psites) == 0:
        return (np.nan)

    if nfft == "default":
        nfft = int(2 * 2**np.ceil(np.log2(len(psites))))

    if ntapers == "default":
        ntapers = int(2*time_bandwidth) - 1

    # Calculate the spectral estimation.
    spec, freq, jackknife, fstatistics, _ = mtspec(data=np.array(psites), delta = 1, time_bandwidth = time_bandwidth, number_of_tapers=ntapers, nfft=nfft, statistics=True, rshape=0)

    m = int(np.round (nfft/3))
    sf = f.sf (fstatistics[m],dfn=2,dfd=(2*ntapers)-2)
    return (sf)


def get_wilcox (mat):
    """Paired wilcoxon-test for frame0 > mean (frame1, frame2)
    mat: 2D matrix with shape (3, ncodons)
    returns: p-value
    """

    frame0 = mat[:,0]
    frame12 = np.mean (mat[:,1:3], axis=1)

    wilcox_stat, wilcox_p = wilcoxon(frame0, frame12, alternative="greater", paired=True) if not np.all (frame0-frame12==0) else (np.nan, np.nan)
    return (wilcox_p)



def get_binom (mat):
    """Perform binomial-test for n(frame0 > frame1 and frame0 > frame2). Adding random noise to reduce draw-bias, otherwise on-frame is max on draw
    mat: 2D matrix with shape (3, ncodons)
    returns: p-value"""


    mat = mat + np.random.uniform(low=0.0, high=0.99, size=mat.shape)

    index_max = np.argmax (mat, axis=1)
    binom_p = binomtest (k=np.sum (index_max == 0), n=len(index_max), p=1/3, alternative="greater").pvalue if len (index_max) > 0 else np.nan
    return (binom_p)


def get_glm (mat):
    """Fits a negative binomial GLM to the p-sites with a two-class frame feature (on or off-frame) and extracts the parameter for the frame coefficient. 
    mat: 2D matrix with shape (3, ncodons)
    returns: p-value"""

    df_glm = pd.DataFrame ({
        'counts' : mat.reshape (-1),
        'frame' : ['onframe', 'offframe', 'offframe'] * mat.shape[0]
    })



    try:
        model = smf.glm(formula = "counts ~ frame", data=df_glm, family=sm.families.NegativeBinomial()).fit()
        #r = np.array([1, -.5, -.5])
                
        #r = np.array([1, -1, -1])
        #glm_ttest = model.t_test(r)
        
        glm_p = model.pvalues[1] # glm_ttest.pvalue

        if glm_p < 1e-8:            
            print (model.params, glm_p)

        # converting to one-tailed
        if model.params[1] > 0: #== max (model.params):
            glm_p_onetailed = glm_p/2
        else:
            glm_p_onetailed = 1-glm_p/2

        if glm_p < 1e-8:            
            print (glm_p_onetailed)

        return (glm_p_onetailed)


    except sme.PerfectSeparationError:
        return (np.nan)


    
def get_phasing_stats (input, p_methods = ['wilcox', 'binom', 'glm']):
    """calculates phasing statistics from input (either 1D array with p-site counts or 2D matrix with frame-stratified counts (one column per frame))
    input: 1D or 2D object with p-site distribtion
    returns: dict with p-values for the specified tests
    """
    
    if len (input.shape) == 1: #1D input
        psites = input
        mat = get_phasing_matrix (psites)
    else:  #2D input
        mat = input
        psites = mat.reshape(-1)

    if isinstance (p_methods, str):
        p_methods = [p_methods]

    # remove codons without any signal
    mat = mat[~np.all(mat == 0, axis=1)]

    output_stats = {'n' : mat.shape[0]}

    for m in p_methods:
        output_stats[m] = np.nan
    
    if mat.shape[0] == 0:
        return (output_stats)


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if 'glm' in p_methods:
            output_stats['glm'] = get_glm (mat)
            
        if 'wilcox' in p_methods:
            output_stats['wilcox'] = get_wilcox (mat) 

        if 'binom' in p_methods:
            output_stats['binom'] = get_binom (mat)

        if 'taper' in p_methods:
            output_stats['taper'] = get_taper (psites)

    return (output_stats)




def get_psites (orf_tid, start, stop, bamfiles, pd_offsets, tid2ref_dict = None, bam_dict = None):
    """Iterates through bam-files and orfs to collect all p-sites"""

    if tid2ref_dict == None or bam_dict == None:
        
        tid2ref_dict, bam_dict = {},{}
        for bamfile in bamfiles:

            bam_dict[bamfile] = pysam.Samfile (bamfile)

            _, tid2ref = get_tid_info (bamfile)
            tid2ref_dict[bamfile] = tid2ref


    psites = np.zeros ((stop-start))

    # fetch read from all bams mapping on tid
    for bamfile in bamfiles:

        dtid2ref = tid2ref_dict[bamfile]

        if not orf_tid in dtid2ref:
            continue

        tid = dtid2ref[orf_tid]

        bam_offsets = pd_offsets[pd_offsets.bam == bamfile]
        doffsets = bam_offsets[["read_length", "offset_key"]].set_index('read_length').to_dict ()['offset_key']

        bam = bam_dict[bamfile]
        

        if not orf_tid in dtid2ref:
            print (f"ERROR: {orf_tid} not found in bam... skipping")
            continue
        
        # array of p-site position counts
        # if len(cds)%3 != 0:
        #     print (f"ERROR: CDS-length invalid {tid}:{start}-{stop}")
        #     continue

        for read in bam.fetch (tid, start, stop):

            if read.is_reverse:
                continue
                    
            read_length = read.infer_read_length () 
            
            if not read_length in doffsets:
                continue

            length_offset = int(doffsets[read_length])
            
            offset_pos = read.pos + length_offset
                    
            if offset_pos >= start and offset_pos < stop: 
                try:
                    psites[offset_pos-start] += 1
                except IndexError:
                    print ("error", offset_pos, start, len(psites))


    return (psites)


"""
get_phasing: Main function. Iterates through all ORF

"""   
def get_phasing (bamfiles, orfs, offsets, output, percentile=0.9, alpha = 0.01, p_methods = ["glm"], shuffle = False, multiplier = [1]):

    
    print ("### get_phasing ###")

    if isinstance (p_methods, str):
        p_methods = [p_methods]
    
    
    if not libmtspec and 'taper' in p_methods:
        print("module 'mtspec' is required for taper statistics")
        p_methods.remove ("taper")

        if len (p_methods) == 0:
            print ("no methods left... exiting")
            return ()

    

    bam_dict = {}
    tid2ref_dict = {}
        
    print ("loading bam...")
    
    for bamfile in bamfiles:

        bam_dict[bamfile] = pysam.Samfile (bamfile)

        _, tid2ref = get_tid_info (bamfile)
        tid2ref_dict[bamfile] = tid2ref


    
    print ("loading offsets...")

    pd_offsets = pd.read_csv (offsets, sep="\t")

    # only use read_length showing significant phasing
    
    pd_offsets = pd_offsets[pd_offsets.offset_pct >= percentile]
    
    for s in p_methods:
        pd_offsets = pd_offsets[pd_offsets["p_" + s] <= alpha]
    
    header = {}

    # output file
    fout = open (output, "w")

    # checking orfs - line by line
    print ("checking orfs...")

    num_lines = sum(1 for line in open(orfs,'r'))

    with open(orfs, 'r') as f:

        for i, line in enumerate (tqdm(f, total=num_lines)):

            columns = line.strip ("\n").split ("\t")

            if i == 0:
                for icol, col in enumerate (columns):
                    header[col] = icol

                # printing output header

                header_col = columns + ["total_counts", "frame0", "frame1", "frame2"]                                
                
                p_headers = ["p_"+s for s in p_methods]

                for f in multiplier:                    
                    header_col += [s + "_" + str(f) if f != 1 else s for s in p_headers + ['n']]                    

                print ("\t".join (header_col), file=fout)
                continue

            orf_id = columns[header['orf_id']] 
            orf_tid = columns[header['tid']] 
            start, stop = int(columns[header['start']]), int(columns[header['stop']])+3

            psites = get_psites (orf_tid, start, stop, bamfiles, pd_offsets, tid2ref_dict=tid2ref_dict, bam_dict=bam_dict)

            if shuffle:
                np.random.shuffle (psites)


            counts = get_counts (psites)


            output = columns + [str (counts['total']), str(counts['frame0']), str(counts['frame1']), str(counts['frame2'])]

            # the multiplier is only relevent to show signal-length vs statistics
            for f in multiplier:
                    
                if f < 1:
                    new_l = int(len (psites)*f)
                    new_l -= new_l%3
                    sub_psites = psites[0:new_l]
                else:
                    sub_psites = psites * int(f)

                output_stats = get_phasing_stats (sub_psites, p_methods)
                
                output += [str(output_stats[p]) for p in p_methods + ['n']]


            print ("\t".join (output), file=fout)

            # testing
            if shuffle and output_stats['glm'] < 1e-8:
                pd_out = pd.DataFrame (get_phasing_matrix (psites), columns=["frame0", "frame1", "frame2"])
                pd_out.to_csv ("test.txt", sep="\t")
                print (output)


    print ("### Done ###")


def ribofy_phasing ():

    parser = argparse.ArgumentParser(description='get phasing')
    parser.add_argument("--bam", dest='bam', required=True, nargs="+", help="Bam file - sorted and indexed")
    parser.add_argument("--orfs",   dest='orfs', required=True, help="orfs - generated by get_ORFs.py")
    parser.add_argument("--offsets", dest='offsets', required=True, help="offsets - generated by get_offsets.py")
    parser.add_argument("--output", dest='output', default = "ribofy_phasing.txt", help="output")

    parser.add_argument("--percentile", dest='percentile', default = 0.9, help="Percentile of consistent offset-determinants")
    parser.add_argument("--alpha", dest='alpha', default = 0.01, help="cutoff p-value for phase-detection at percentile")
    parser.add_argument("--p_methods", dest='p_methods', nargs="*", default = ["binom", "wilcox", "glm"], help="Statistics: possibilities: binom, wilcox, glm, and taper")
    parser.add_argument("--multiplier", dest='multiplier', nargs="*", default = [1], type=float, help="multipler")


    args = parser.parse_args()

    get_phasing (args.bam, args.orfs, args.offsets, args.output,
                 percentile=args.percentile,
                 alpha=args.alpha,
                 p_methods=args.p_methods,
                 multiplier=args.multiplier)



if __name__ == "__main__":
        
    ribofy_phasing ()