"""
get_ORFs retrieves all putative ORFs in transcriptome. This only has to be done once per gtf/genome


usage:
python get_ORFs.py --gtf <gtf-file> --fa <genome fasta file> --output <output-file>

By default, any ORFs less than 90 nts (30 amino acid) are discarded, but this is set
by the --min_aa_length flag

"""

from ribofy.stats import get_2D_matrix
import sys
import pysam
import re
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from collections import defaultdict

from . import __version__
from .argparse2 import argparse2
from .utils import rev_comp, translate
from .gtf2 import *


libpybigwig = True
try:
    import pyBigWig
except ModuleNotFoundError:
    libpybigwig = False
    

def pickle_save (obj, file):
    """Save as gzipped pickle"""

    try:
        import pickle, gzip
        with gzip.open(file, 'wb') as handle:
            pickle.dump(obj, handle)
        
    except ModuleNotFoundError:
        print("modules 'pickle' and 'gzip' are required for saving, skipped...")


def get_orfs_from_seq (seq, start_codon, stop_codon, min_aa_length):
    """Get ORFs from sequence"""

    seq = seq.upper()

    # searching for start (ATG) and stop-codons            
    start = [m.start() for m in re.finditer(start_codon, seq)]
    stop = [m.start() for m in re.finditer(stop_codon, seq)]

    start_frame = [s % 3 for s in start]
    stop_frame = [s % 3 for s in stop]


    # getting all longest possible ORFS in each frame

    dorf = []

    for frame in [0,1,2]:

        ipos = 0

        while True:

            fstart = [s for i, s in enumerate(start) if start_frame[i] == frame and s >= ipos]
            
            if len (fstart) == 0:
                break

            fstop = [s for i, s in enumerate(stop) if stop_frame[i] == frame and s > fstart[0]]

            if len (fstop) == 0:
                break

            orf_length = abs (fstop[0] - fstart[0])

            ipos = fstop[0]

            if orf_length < min_aa_length*3:
                continue

            dorf.append ({'start':fstart[0], 'stop': fstop[0], 'frame':frame, 'orf_length':orf_length})


    # sorting - not required, but then low orf_id corresponds to longer ORFs
    dorf = sorted(dorf, key = lambda i: i['orf_length'], reverse=True)

    return (dorf)


def get_orf_groups (edge_dict, graph_output):
    """get ORF groups based on chromosome-stratified (to save memory) network (edge_dict)
    
    """

    orf_groups = {}
    group_count = 0    
    group2edge = {}

    for chrom in tqdm(edge_dict):
    
        elist = [(e1, e2) for (e1, e2) in edge_dict[chrom]]
        
        # build graph
        g = nx.Graph()
        g.add_edges_from(elist)
        
        # extract graph groups        
        for ig, cc in enumerate (list(nx.connected_components(g))):

            id_text = f"group_{(group_count+1):05d}"

            for group_id in cc:                
                orf_groups[group_id] = id_text

                if graph_output != "":
                    group2edge[id_text] = (chrom, ig)

            group_count += 1

    return (orf_groups, group_count, group2edge)



def get_orfs (gtf, fa, output, start_codon = "ATG", stop_codon = "TGA|TAA|TAG", min_aa_length=30,  output_fa = False, error_output="", graph_output="", cons_bw = "", ):
    """
    Main function: Extract putative ORFs from GTF, establish ORF-network and assign type to each ORF

    Parameters
    ----------
    gtf: str
        Path/to/gtf-file (only tested with gencode)
    fa: str
        Path/to/fa; genome fasta file (compatible with GTF, and indexed: samtools faidx path/to/fa)
    output: str
        path/to/output (to be used in ribofy detect)
    start_codon: str, default= "ATG"
        start codon sequence used in finding ORFs - for multiple codons, seperate by '|', e.g. 'ATG|CTG' to allow for non-canonical start codons
    stop_codon: str, default= "TGA|TAA|TAG"
        stop codon sequence used in finding ORFs
    min_aa_length: int, default=30
        Minimun length of ORF (in amino acids)
    output_fa: str, optional
        Output the nucleotide and amino-acid sequences of ORFs (<output>.nt.fa and <output>.aa.fa, respectively)

    graph_output: str, optional
        Save graph network to file
    cons_bw: str, optional
        Bigwig file with conservation scores... greatly increases duration
        
    """


    print ("### get_orfs ###")

    start_codon = start_codon.upper().replace ("U", "T")
    stop_codon = stop_codon.upper().replace ("U", "T")

    print ("reading gtf...")
    gtf = gtf2 (gtf)
    fa = pysam.FastaFile (fa)

    if cons_bw != "" and not libpybigwig:
        print ("pyBigwig must be install to extract conservation scores from bw-files")
        cons_bw = ""
    
    if cons_bw != "":
        bw = pyBigWig.open(cons_bw)
    

    ferror = open (error_output, "w") if error_output != "" else None    
    fseq_aa  = open (output + ".aa.fa", "w") if output_fa else None
    fseq_nt  = open (output + ".nt.fa", "w") if output_fa else None

    start_codon = start_codon.upper()

    lorfs = []
    
    cons_dict, edge_dict = {}, {}

    print ("finding ORFs in all transcripts...")

    for tid in tqdm(gtf.get_all_tids ()):

        seq = "" # transcript sequence
        dt2g = {} # dict, converting from transcript position to genomic position
        cons_arr = [] # conservation scores (only in use if cons_bw is set)

        gannot_start = gtf.get_startcodon (tid)-1
        gannot_stop = gtf.get_stopcodon (tid)-1
        gid, symbol, biotype = gtf.get_gene_id (tid), gtf.get_name (tid), gtf.get_type (tid)        
        chrom, strand = gtf.get_chrom (tid), gtf.get_strand (tid)

        annot_start, annot_stop = 0, 0  # relative, annotated start/stop codons

        # iterate through exons
        for (chrom, start, end, strand) in gtf.get_exon_coords (tid):

            for i, g in enumerate (range (start-1, end)):
                dt2g[i+len(seq)] = g
                
            seq += fa.fetch (chrom, start-1, end)
            
            annot_start += max (0, min (end, gannot_start) - (start-1))
            annot_stop +=  max (0, min (end, gannot_stop) - (start-1))

            if cons_bw != "":
                cons_arr.extend (bw.values (chrom, start-1, end))

        # reverse complement '-'-strand transcripts
        if strand == "-":
            seq = rev_comp (seq)
            annot_start, annot_stop = len(seq) - annot_start - 1, len(seq) - annot_stop - 1
            cons_arr = cons_arr[::-1]

        #if cons_bw != "":
        #    cons_dict[tid] = cons_arr

        # no proper annotation
        if gannot_start <= 0 or gannot_stop <= 0:
            annot_start, annot_stop = -1,-1

        dorf = get_orfs_from_seq (seq, start_codon, stop_codon, min_aa_length)

        if len (dorf) == 0:
            continue

        for i, orf in enumerate (dorf):
            
            orf['tid']    = tid
            orf['gid']    = gid
            orf['symbol'] = symbol                        
            orf['chrom']  = chrom
            orf['strand'] = strand            
            orf['gannot_start'] = gannot_start
            orf['gannot_stop']  = gannot_stop
            orf['annot_start']  = annot_start 
            orf['annot_stop']   = annot_stop                     
            orf['bio_type']     = biotype

            orf['orf_id'] = f"{tid}_orf_{(i+1):05d}"

            orf_seq = seq[orf['start']:orf['stop']]

           

            if output_fa:
                #fa_header = 
                #print (f">{orf['orf_id']}|{orf['gid']}|{orf['symbol']}|{orf['bio_type']}|{row.orf_type}")

                print (f">{orf['orf_id']}\n{translate(orf_seq)}", file=fseq_aa)
                print (f">{orf['orf_id']}\n{orf_seq}", file=fseq_nt)


            if cons_bw:
                cons_mat = get_2D_matrix (cons_arr[orf['start']:orf['stop']])
                orf['cons_f0'] = np.mean (cons_mat[:,0])
                orf['cons_f1'] = np.mean (cons_mat[:,1])
                orf['cons_f2'] = np.mean (cons_mat[:,2])

         
            orf['tid_length'] = len(seq)
            
            #using annotated stop
            if orf['stop'] == annot_stop and gannot_start > 0:
                orf['orf_type'] = "annotated"

                if annot_start < orf['start'] and abs(annot_stop-annot_start)%3 != 0:
                                        
                    if error_output != "":
                        ## ERROR: invalid ORF annotation
                        error_data = [str(orf[c]) for c in orf]
                        error_data += [seq]
                        error_data += [seq[orf['start']:orf['stop']]]
                        error_data += [translate (seq[orf['start']:orf['stop']])]

                        print ("\t".join (error_data), file=ferror)

                    # overwrite annot
                    annot_start = orf['start']                     

                else:
                    orf['start'] = annot_start
                    orf['orf_length'] = orf['stop'] - orf['start']



            elif orf['stop'] < annot_stop and orf['start'] > annot_start:
                continue

            elif orf['stop'] > annot_stop and annot_stop != annot_start:
                orf['orf_type'] = "dORF"

            elif orf['stop'] < annot_stop and orf['start'] < annot_start and annot_stop != annot_start:
                orf['orf_type'] = "uORF"
                
            else:
                orf['orf_type'] = "novel"

            orf['gstart'] = dt2g[orf['start']] if strand == "+" else dt2g[len(seq)-orf['start']-1]
            orf['gstop'] = dt2g[orf['stop']] if strand == "+" else dt2g[len(seq)-orf['stop']-1]

            range1 = range (orf['start'], orf['stop'], 3)
            range2 = range (orf['start']+3, orf['stop']+3, 3)

            for pos1, pos2 in zip (range1, range2):

                p1 = dt2g[pos1] if strand == "+" else dt2g[len(seq)-pos1-1]
                p2 = dt2g[pos2] if strand == "+" else dt2g[len(seq)-pos2-1]
      
                e1 = f"{chrom}:{p1}{strand}"
                e2 = f"{chrom}:{p2}{strand}"

                if not orf['chrom'] in edge_dict:
                    edge_dict[orf['chrom']] = {} 

                edge_dict[orf['chrom']][(e1, e2)] = 1

            lorfs.append (orf)



    if error_output != "":
        ferror.close()

    if output_fa:
        fseq_aa.close()
        fseq_nt.close()

    
    print ("infering ORF groups...")

    orf_groups, group_count, group2edge,  = get_orf_groups (edge_dict, graph_output)
    
    print (f"found {group_count} ORF-groups in {len (lorfs)} total ORFs")

    if graph_output != "":

        print (f"saving network edges")
        pickle_save ({'edges' : edge_dict, 'group2edge' : group2edge}, graph_output)
               

    # if cons_bw != "":

    #     print (f"saving conservation scores")
    #     pickle_save (cons_dict, output + ".cons.gz")
        

    print (f"assigning ORF-group to individual ORF")

    connected = 0
    # Assign group to each orf  
    for orf in lorfs:
            
        groupid_from_start = orf_groups[f"{orf['chrom']}:{orf['gstart']}{orf['strand']}"]
        groupid_from_stop = orf_groups[f"{orf['chrom']}:{orf['gstop']}{orf['strand']}"]      

        if groupid_from_stop != groupid_from_start:
            print (f"ERROR: Unconnected network in {orf['tid']} : {groupid_from_start} vs {groupid_from_stop} {orf['strand']}")
        else:
            connected += 1
        groupid = orf_groups[f"{orf['chrom']}:{orf['gstart']}{orf['strand']}"]

        orf['orf_group'] = groupid

    # set group type: annotated > uORF > dORF > novel
    group_score = {}
    group_conv = {'annotated' : 3, 'uORF' : 2, 'dORF' : 1, 'novel' : 0}

    for orf in lorfs:
        
        score = group_score[orf['orf_group']] if orf['orf_group'] in group_score else 0
        group_score[orf['orf_group']] = max (score, group_conv[orf['orf_type']])
  
    
    print ("outputting...")

    

    columns = ["gid", "symbol", "tid", "start", "stop", "tid_length", "annot_start", "annot_stop", "frame",
               "chrom", "gstart", "gstop", "strand", "orf_length", "orf_type", 
               "bio_type", "orf_id", "orf_group"]

    with open (output, "w") as fout:
        
        print ("\t".join (columns), file=fout)
        for orf in lorfs:
                      
            gscore = group_score[orf['orf_group']]
            oscore = group_conv[orf['orf_type']]
            
            if oscore > gscore:
                print ("ERROR: orf_group scores invalid")

            if oscore >= gscore:
                print ("\t".join ([str(orf[col]) for col in columns]), file=fout)

     
    print ("### Done ###")





def ribofy_orfs ():

    info_text = """
        ribofy orfs: extracting ORFs from GTF
    """

    help_text = f"""    
        ribofy orfs - version {__version__}

        required arguments:
        --gtf <file>            GTF file, GENCODE-style
        --fa <file>             Genome Fasta file (indexed with samtools faidx)
        --output <str>          Output filename, default=orfs.txt   

        optional arguments:
        --start_codon <str>     Specify start_codons for ORF detection. default="ATG"
        --stop_codon <str>      Specify stop_codons for ORF detection. default="TGA|TAA|TAG"
        --min_aa_length <INT>   Minimum peptide length, default=30
        --output_fa             If set, outputs nucleotide and amino-acid fasta files (<output>.nt.fa and 
                                <output>.aa.fa, respectively) for all ORFs found
        
        usage: ribofy orfs --gtf GTF --fa FA [--output OUTPUT]\n"""

    parser = argparse2 (
        description=info_text,
        usage=help_text,
        help=help_text
    )
    parser.add_argument('orfs', nargs='?', help='') # dummy positional argument
    
    # required    
    parser.add_argument("--gtf", dest='gtf', required=True)
    parser.add_argument("--fa", dest='fa', required=True)
    parser.add_argument("--output", dest='output', default = "orfs.txt")

    # optional        
    parser.add_argument("--start_codon", dest='start_codon', type=str, default="ATG")
    parser.add_argument("--stop_codon", dest='stop_codon', type=str, default="TGA|TAA|TAG")
    parser.add_argument("--min_aa_length", dest='min_aa_length', type=int, default=30)
    parser.add_argument("--output_fa", dest='output_fa', action="store_true")
    
    parser.add_argument("--error_output", dest='error_output', type=str, default="")
    parser.add_argument("--graph_output", dest='graph_output', type=str, default="")

    parser.add_argument("--cons_bw", dest='cons_bw', type=str, default="")



    args = parser.parse_args()

    get_orfs (args.gtf, args.fa, args.output, 
              start_codon=args.start_codon, stop_codon=args.stop_codon, 
              min_aa_length=args.min_aa_length, output_fa=args.output_fa, 
              error_output=args.error_output,
              graph_output=args.graph_output,
              cons_bw = args.cons_bw)



if __name__ == "__main__":

    ribofy_orfs ()







        




