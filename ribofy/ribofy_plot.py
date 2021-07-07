import sys
import argparse
import pandas as pd

from . import __version__
from .argparse2 import argparse2
from .get_phasing import get_psites

import matplotlib    
matplotlib.use('Agg')    
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



def ribofy_plot ():

    info_text = """
        ribofy orfs: extracting ORFs from GTF
    """

    help_text = f"""    
        ribofy plot - version {__version__}

        TODO
        usage: ribofy plot ..."""

    parser = argparse2 (
        description=info_text,
        usage=help_text,
        help=help_text
    )
    parser.add_argument('plot', nargs='?', help='') # dummy positional argument
    
    
    parser.add_argument("--id", dest='id', type=str, required=True, help="tid/gid/symbol")
    parser.add_argument("--orfs", dest='orfs', type=str, required=True, help="orfs")
    parser.add_argument("--bams", dest='bams', nargs="+", type=str, required=True, help="bams")
    parser.add_argument("--offsets", dest='offsets', nargs="+", type=str, required=True, help="offsets")


    args = parser.parse_args()

    print ("reading orfs...")

    pd_orfs = pd.read_csv (args.orfs, sep="\t")

    pd_gene = pd_orfs[pd_orfs.tid == args.id]
    if pd_gene.shape[0] == 0:
        pd_gene = pd_orfs[pd_orfs.gid == args.id]
    if pd_gene.shape[0] == 0:
        pd_gene = pd_orfs[pd_orfs.symbol == args.id]

    if pd_gene.shape[0] == 0:
        print (f"{args.id} not found")
        sys.exit ()

    print ("reading offsets...")
    pd_offsets = pd.concat ([pd.read_csv (o, sep="\t") for o in args.offsets])
    

    for tid in list(set(pd_gene.tid.values)):

        pd_tid = pd_gene[pd_gene.tid == tid] 
                
        length = pd_tid.tid_length.values[0]
        symbol = pd_tid.symbol.values[0]
        
        print (f"retrieving p-sites ({tid})...")
        cds = get_psites (tid, 0, length, args.bams, pd_offsets)

        print ("plotting...")

        frame_color = {
            0:"green",
            1:"blue",
            2:"red"       
        }

        fig, ax = plt.subplots(figsize=(20, 10)) #, dpi=600)
        
        for i in range (3):

            x = [p for p in range (len (cds)) if p%3==i]
            y = [cds[p] for p in range (len (cds)) if p%3==i]           
            ax.bar(x, y, color=frame_color[i], label=f"frame{i}")

        max_y = max (cds)
        frame_unit = max_y/100
        min_y = -frame_unit*5
        ax.set_ylim([min_y, max_y])
                
        print ("adding orfs...")

        for idx, row in pd_tid.iterrows():  
            start = row['start']
            length = row['orf_length']
            frame = row['frame']
            frame_y_pos = -(frame+2) * frame_unit
            
            ax.add_patch(Rectangle((start, frame_y_pos), length, frame_unit-1, color=frame_color[frame]))

        ax.legend()

        print ("saving...")

        #plt.savefig(tid + ".png")
        plt.savefig(f"{symbol}_{tid}.pdf", format='pdf')


if __name__ == "__main__":

    ribofy_plot ()







        




