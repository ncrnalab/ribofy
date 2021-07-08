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
    
def get_phasing_matrix (cds):
    
    mat = np.array ([[s for i, s in enumerate (cds) if i%3 == 0], # on-frame
                     [s for i, s in enumerate (cds) if i%3 == 1],
                     [s for i, s in enumerate (cds) if i%3 == 2]]).T

    
    return (mat)


def get_phasing_stats (mat):

    # remove codons without any signal
    mat = mat[~np.all(mat == 0, axis=1)]

    if mat.shape[0] < 10:
        return (np.nan, np.nan, np.nan)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # nb glm
        
        df_glm = pd.DataFrame ({
            'counts' : mat.reshape (-1),
            'frame' : ['f1', 'f2', 'f3'] * mat.shape[0]
        })

        try:
            model = smf.glm(formula = "counts ~ frame", data=df_glm, family=sm.families.NegativeBinomial()).fit()
            r = np.array([1,-.5, -.5])
            glm_ttest = model.t_test(r)
            glm_p = glm_ttest.pvalue
            # converting to one-tailed
            if model.params[0] >= max (model.params):
                glm_p = glm_p/2
            else:
                glm_p = 1-glm_p/2


        except sme.PerfectSeparationError:
            glm_p = np.nan


        # wilcoxon-test for frame0 > mean (frame1, frame2)
        wilcox_stat, wilcox_p = wilcoxon(mat[:,0], np.mean (mat[:,1:3], axis=1), alternative="greater") if mat.shape[0] >= 10 else (np.nan, np.nan)

        # binomial-test for n(frame0 > frame1 and frame0 > frame2)
        # add random noise to reduce draw-bias
        mat = mat + np.random.uniform(low=0.0, high=0.99, size=mat.shape)

        index_max = np.argmax (mat, axis=1)
        binom_p = binomtest (k=np.sum (index_max == 0), n=len(index_max), p=1/3, alternative="greater").pvalue if len (index_max) > 0 else np.nan

    return (wilcox_p, binom_p, glm_p)

def get_spectrum (cds, time_bandwidth = 3, ntapers = "default", nfft = "default"):
    
    if sum (cds) == 0:
        return (np.nan)

    if nfft == "default":
        nfft = int(2 * 2**np.ceil(np.log2(len(cds))))

    if ntapers == "default":
        ntapers = int(2*time_bandwidth) - 1

    # Calculate the spectral estimation.
    spec, freq, jackknife, fstatistics, _ = mtspec(data=np.array(cds), delta = 1, time_bandwidth = time_bandwidth, number_of_tapers=ntapers, nfft=nfft, statistics=True, rshape=0)

    m = int(np.round (nfft/3))
    sf = f.sf (fstatistics[m],dfn=2,dfd=(2*ntapers)-2)
    return (sf)


def get_psites (orf_tid, start, stop, bamfiles, pd_offsets, tid2ref_dict = None, bam_dict = None):

    if tid2ref_dict == None or bam_dict == None:
        
        tid2ref_dict, bam_dict = {},{}
        for bamfile in bamfiles:

            bam_dict[bamfile] = pysam.Samfile (bamfile)

            _, tid2ref = get_tid_info (bamfile)
            tid2ref_dict[bamfile] = tid2ref


    cds = [0] * (stop-start)

    # fetch read from all bams mapping on tid
    for bamfile in bamfiles:

        dtid2ref = tid2ref_dict[bamfile]

        if not orf_tid in dtid2ref:
            continue

        tid = dtid2ref[orf_tid]

        bam_offsets = pd_offsets[pd_offsets.bam == bamfile]
        doffsets = bam_offsets[["read_length", "offset_table_key"]].set_index('read_length').to_dict ()['offset_table_key']

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

            length_offset = doffsets[read_length]
            
            offset_pos = read.pos + length_offset
                    
            if offset_pos >= start and offset_pos < stop: 
                cds[offset_pos-start] += 1


    return (cds)

def get_phasing (bamfiles, orfs, offsets, output, percentile=0.9, alpha = 0.01, devel = False):

    
    print ("### get_phasing ###")
    
    if not libmtspec:
        print("module 'mtspec' is required for devel mode - setting devel = False")
        devel = False

    if devel:
        multiplier = [.5, 1, 2]
    else:
        multiplier = [1]


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
    
    pd_offsets = pd_offsets[(pd_offsets.offset_table_pct >= percentile) & 
                            (pd_offsets.p_wilcox_percentile <= alpha) & 
                            (pd_offsets.p_binom_percentile <= alpha)]


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

                stats = ["p_wilcox", "p_binom", "p_glm"]
                
                if devel:
                    stats += ["p_taper"]

                for f in multiplier:                    
                    header_col += [s + "_" + str(f) if f != 1 else s for s in stats]
                    

                print ("\t".join (header_col), file=fout)
                continue

            orf_tid = columns[header['tid']] 
            start, stop = int(columns[header['start']]), int(columns[header['stop']])+3

            cds = get_psites (orf_tid, start, stop, bamfiles, pd_offsets, tid2ref_dict=tid2ref_dict, bam_dict=bam_dict)

            mat = get_phasing_matrix (cds)

            output = columns + [str(sum(sum(mat))), str(sum(mat[:,0])), str(sum(mat[:,1])), str(sum(mat[:,2]))]

            # the multiplier is only relevent in devel-mode to show signal-length vs statistics
            for f in multiplier:
                    
                if f < 1:
                    new_l = int(len (cds)*f)
                    new_l -= new_l%3
                    sub_cds = cds[0:new_l]
                else:
                    sub_cds = cds * int(f)

                mat = get_phasing_matrix (sub_cds)
                wilcox_p, binom_p, glm_p = get_phasing_stats (mat)
                
                output += [str(wilcox_p), str(binom_p), str(glm_p)]

                if devel:
                    output += [str(get_spectrum (sub_cds))]
            

            print ("\t".join (output), file=fout)


    print ("### Done ###")


def ribofy_phasing ():

    parser = argparse.ArgumentParser(description='get phasing')
    parser.add_argument("--bam", dest='bam', required=True, nargs="+", help="Bam file - sorted and indexed")
    parser.add_argument("--orfs",   dest='orfs', required=True, help="orfs - generated by get_ORFs.py")
    parser.add_argument("--offsets", dest='offsets', required=True, help="offsets - generated by get_offsets.py")
    parser.add_argument("--output", dest='output', default = "ribofy_phasing.txt", help="output")

    parser.add_argument("--percentile", dest='percentile', default = 0.9, help="Percentile of consistent offset-determinants")
    parser.add_argument("--alpha", dest='alpha', default = 0.01, help="cutoff p-value for phase-detection at percentile")



    args = parser.parse_args()

    get_phasing (args.bam, args.orfs, args.offsets, args.output,
                 percentile=args.percentile,
                 alpha=args.alpha)



if __name__ == "__main__":
        
    ribofy_phasing ()