"""
get_offset determines the optimal p-site offset for each read-length on the top 10 most abundant ORFs in the bam-file


usage:
python get_offset.py --bam <bam-file> --orfs <ribofy orfs-file> --output <output-file>

By default, get_offset analyses reads between 25 and 35 nt, 
but this is customizable with the --min_read_length and --max_read_length options

And change number of ORFs used in offset calculation by the --norfs option

"""


import sys
import pysam
import pandas as pd
import numpy as np
from .argparse2 import argparse2
from .get_phasing import get_phasing_stats, get_phasing_matrix
from .bam_utils import get_tid_info


def percentile(n):
    def percentile_(x):
        return x.quantile(n)
    percentile_.__name__ = 'percentile_{:2.0f}'.format(n*100)
    return percentile_

def get_offset (bamfile, orfs, output, norfs=10, min_read_length=25, max_read_length=35):

    print ("### get_offset ###")

    print ("loading orfs...")

    pd_orfs = pd.read_csv (orfs, sep="\t")
    pd_orfs = pd_orfs[(pd_orfs.orf_type == "annotated")] \
        .groupby ("orf_group") \
        .head (1) 

    # get transcripts with most counts
    print ("get transcript counts from bam...")

    # load bam
    bam = pysam.Samfile (bamfile)
    
    dtid2count, dtid2ref = get_tid_info (bamfile)

    
    # add read counts to dataframe
    pd_orfs['total_reads'] = pd_orfs['tid'].transform (lambda x: dtid2count[x])
    pd_orfs = pd_orfs.sort_values ('total_reads', ascending=False)

    print ("using the following transcripts..:")
    print (pd_orfs[['tid', 'total_reads']].head (norfs))


    # initialize count_offsets
    length_range = range (min_read_length, max_read_length+1) 
    
    count_offsets = {}

    for x in length_range:
        count_offsets[x] = [0,0,0]

    off_conv = {0:0, 1:2, 2:1}

    offset_stats = []

    for i, row in pd_orfs.head (norfs).iterrows ():

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

        #print (dcds)
        for lr in length_range:
            
            mat = get_phasing_matrix (dcds[lr])

            frame_sort = np.argsort(mat.sum(axis=0))[::-1]

            mat = mat[:,frame_sort]

            wilcox_p, binom_p = get_phasing_stats (mat)

            offset_stats.append ({
                'read_length' : lr,
                'tid' : tid,
                'p_wilcox' : wilcox_p,
                'p_binom' : binom_p,
                'offset' : 12 + off_conv[frame_sort[0]],
                'onframe' : np.sum (mat[:,0]),
                'offframe1' : np.sum (mat[:,1]),
                'offframe2' : np.sum (mat[:,2])                
            })

    pd_stats = pd.DataFrame (offset_stats)
    
    pd_agg = pd_stats.dropna().groupby ('read_length').agg ({'p_wilcox' : [percentile (.1), percentile (.9)],
                                          'p_binom' : [percentile (.1), percentile (.9)],
                                          'offset' : [percentile (.1), percentile (.9)],
                                          'onframe' : [np.sum],
                                          'offframe1' : [np.sum],
                                          'offframe2' : [np.sum]})

    pd_agg['offset'] = pd_agg['offset'].astype ("int")
    pd_agg['mean_offset'] = pd_agg['offset'].mean (axis=1)

    pd_agg.columns = ['_'.join(col).strip("_") for col in pd_agg.columns.values]
    
    pd_agg['bam'] = bamfile

    print (pd_agg)

    pd_agg.to_csv ("pd_stats.txt", sep="\t")
    
    
    # # compile dataframe
    # pd_offsets = pd.DataFrame ({        
    #     'read_length' : list (length_range),
    #     'offset' : [str(12+count_offsets[x].index (max(count_offsets[x]))) for x in length_range],
    #     'reads12' : [count_offsets[x][0] for x in length_range],
    #     'reads13' : [count_offsets[x][1] for x in length_range],
    #     'reads14' : [count_offsets[x][2] for x in length_range],
    #     'bam' : bamfile
    # })


    pd_agg.to_csv (output, sep="\t")

    print (pd_agg.columns)

    print ("extracted offsets:")
    print (pd_agg[['mean_offset']])
    print ("done")


def ribofy_offset ():

    parser = argparse2 (
        description="",
        usage="",
        help=""
    )

    parser.add_argument('detect', nargs='?', help='') # dummy argument
    parser._action_groups.pop()

    parser.add_argument("--bam",   dest='bam', required=True, help="bam file - sorted and indexed")
    parser.add_argument("--orfs",  dest='orfs', required=True, help="orfs - generated by get_ORFs.py")
    parser.add_argument("--output", dest='output', default = "ribofy_offsets.txt", help="output")
    
    #optional
    parser.add_argument('--norfs', dest='norfs', default = 10, type = int, help="number of distinct orfs to build offsets")
    parser.add_argument('--min_read_length', dest='min_read_length', default = 25, type = int, help="minimum read length used in analysis")
    parser.add_argument('--max_read_length', dest='max_read_length', default = 35, type = int, help="maximum read length used in analysis")
    
    args = parser.parse_args()

    get_offset (args.bam, args.orfs, args.output, 
                norfs=args.norfs, 
                min_read_length = args.min_read_length,
                max_read_length = args.max_read_length)


if __name__ == "__main__":
        
    ribofy_offset ()