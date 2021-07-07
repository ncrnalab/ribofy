"""
get_offset determines the optimal p-site offset for each read-length on the top 10 most abundant ORFs in the bam-file


usage:
python get_offset.py --bam <bam-file> --orfs <ribofy orfs-file> --output <output-file>

By default, get_offset analyses reads between 25 and 35 nt, 
but this is customizable with the --min_read_length and --max_read_length options

And change number of ORFs used in offset calculation by the --norfs option

"""

import pysam
import pandas as pd
import numpy as np
from collections import Counter
from .argparse2 import argparse2
from .get_phasing import get_phasing_stats, get_phasing_matrix
from .bam_utils import get_tid_info


def agg_percentile(n, postfix = ""):
    def percentile_(x):
        return x.quantile(n)
    percentile_.__name__ = ('percentile_' + postfix).strip ("_")
    return percentile_

def agg_table(ret = "key"): #key/value/pct
    def table_(x):
        ctr = Counter(x).most_common()
        if ret == "key":            
            return ([k for (k,v) in ctr][0])
        elif ret == "value":
            return ([v for (k,v) in ctr][0])
        elif ret == "pct":
            return ([v/len(x) for (k,v) in ctr][0])


    table_.__name__ = f"table_{ret}"
    return table_


def get_offset (bamfiles, orfs, output, norfs=20, min_read_length=25, max_read_length=35, percentile = .9, devel=False):

    print ("### get_offset ###")
    print ("loading orfs...")

    pd_orfs = pd.read_csv (orfs, sep="\t")

    pd_annot = pd_orfs[(pd_orfs.orf_type == "annotated") & (pd_orfs.orf_length >= 500)] \
        .groupby ("orf_group") \
        .head (1) 

    pd_combined = pd.DataFrame ()
    pd_devel = pd.DataFrame ()

    for bamfile in bamfiles:

        # get transcripts with most counts
        print ("get transcript counts from bam...")

        # load bam
        bam = pysam.Samfile (bamfile)
        
        dtid2count, dtid2ref = get_tid_info (bamfile)
        
        # add read counts to dataframe
        pd_annot['total_reads'] = pd_annot['tid'].transform (lambda x: dtid2count[x] if x in dtid2count else 0)
        
        pd_annot = pd_annot.sort_values ('total_reads', ascending=False)

        print ("determining p-site offsets using the following transcripts..:")
        print (pd_annot[['tid', 'total_reads']].head (norfs))


        # initialize count_offsets
        length_range = range (min_read_length, max_read_length+1) 
        
        count_offsets = {}

        for x in length_range:
            count_offsets[x] = [0,0,0]

        off_conv = {0:0, 1:2, 2:1}

        offset_stats = []

        for i, row in pd_annot.head (norfs).iterrows ():

            tid, start, end = dtid2ref[row['tid']], int(row['start']), int(row['stop'])+3

            dcds = {}
            for lr in length_range:
                dcds[lr] = [0] * (end-start)
                    
            for read in bam.fetch (tid, start, end):

                init_offset_pos = read.pos + 12

                read_length = read.infer_read_length () 

                if read_length < min_read_length or read_length > max_read_length:
                    continue

                if init_offset_pos >= start and init_offset_pos < end: 

                    init_rel_pos = init_offset_pos - start            
                    offset = off_conv[init_rel_pos % 3]
                    
                    rel_pos = (12 + offset + read.pos) - start

                    if rel_pos % 3 != 0:
                        print ("something wrong with offset")
        
                    count_offsets[read_length][offset] += 1
                    
                    if init_rel_pos >= 0 and init_rel_pos < len (dcds[read_length]): 
                        dcds[read_length][init_rel_pos] += 1

            
            # check phasing for each read length individually
            for lr in length_range:
                
                mat = get_phasing_matrix (dcds[lr])

                frame_sort = np.argsort(mat.sum(axis=0))[::-1]

                # set the frame with most count as onframe (index 0)
                mat = mat[:,frame_sort]

                # test whether onframe counts are significantly enriched across CDS
                wilcox_p, binom_p, glm_p = get_phasing_stats (mat)

                offset_stats.append ({
                    'read_length' : lr,
                    'tid' : tid,
                    'offset' : 12 + off_conv[frame_sort[0]], # best offset
                    'p_wilcox' : wilcox_p,  # stats for best offset
                    'p_binom' : binom_p,   # stats for best offset
                    'p_glm' : glm_p,   # stats for best offset
                    'onframe' : np.sum (mat[:,0]),  # sum of reads with onframe psites
                    'offframe1' : np.sum (mat[:,1]),
                    'offframe2' : np.sum (mat[:,2])                
                })


        pd_stats = pd.DataFrame (offset_stats)

        if devel:
            pd_devel = pd.concat ([pd_devel,pd_stats])
            
        # for each read-length aggregate the results for the analysed transcripts
        pd_stats = pd_stats \
            .dropna().groupby ('read_length') \
            .agg ({'p_wilcox' : [agg_percentile (percentile)],
                'p_binom' : [agg_percentile (percentile)],
                'offset' : [agg_table("key"), agg_table("value"), agg_table("pct")], # [percentile (.05), percentile (.95)],
                'onframe' : [np.sum],
                'offframe1' : [np.sum],
                'offframe2' : [np.sum]})

        pd_stats.columns = ['_'.join(col).strip("_") for col in pd_stats.columns.values]
        
        pd_stats['bam'] = bamfile

        pd_combined = pd.concat ([pd_combined, pd_stats])


    print ("extracted offsets:")
    print (pd_combined[['offset_table_key', 'offset_table_pct']])

    # save to output
    pd_combined.to_csv (output, sep="\t")

    if devel:
        pd_devel.to_csv (output + "_devel.txt", sep="\t")

    print ("### Done ###")
    
   

def ribofy_offset ():

    parser = argparse2 (
        description="",
        usage="",
        help=""
    )

    parser.add_argument('detect', nargs='?', help='') # dummy argument
    parser._action_groups.pop()

    parser.add_argument("--bam",   dest='bam', nargs="+", required=True, help="bam file - sorted and indexed")
    parser.add_argument("--orfs",  dest='orfs', required=True, help="orfs - generated by get_ORFs.py")
    parser.add_argument("--output", dest='output', default = "ribofy_offsets.txt", help="output")
    
    #optional
    parser.add_argument('--norfs', dest='norfs', default = 20, type = int, help="number of distinct orfs to build offsets")
    parser.add_argument('--min_read_length', dest='min_read_length', default = 25, type = int, help="minimum read length used in analysis")
    parser.add_argument('--max_read_length', dest='max_read_length', default = 35, type = int, help="maximum read length used in analysis")
    parser.add_argument("--percentile", dest='percentile', default = 0.9, help="Percentile of consistent offset-determinants")
    
    args = parser.parse_args()

    get_offset (args.bam, args.orfs, args.output, 
                norfs=args.norfs, 
                min_read_length = args.min_read_length,
                max_read_length = args.max_read_length,
                percentile=args.percentile)


if __name__ == "__main__":
        
    ribofy_offset ()