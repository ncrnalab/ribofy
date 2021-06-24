"""
Ribofy: One step wrapper for get_offset, get_phasing and get_results in the ribofy pipeline using default setting.

Before usage, retrieve ORFs in genome using the get_ORFs.py script. This generates the 
ribofy ORFs-file used as input (--orfs).

Each sub-script has its own output: <prefix>_offset.txt, <prefix>_phasing.txt and <prefix>_results.txt, respectively.

The <prefix>_results.txt is the final output with typically the 

usage:
python ribofy.py --bam <bamfile> --orfs <ribofy orfs-file> --prefix <prefix for output-files>

"""


import sys
from . import __version__
from .argparse2 import argparse2
from .get_phasing import get_phasing
from .get_offset import get_offset
from .get_results import get_results


def ribofy_detect ():

    info_text = """
        ribofy detect: extracts phased riboseq reads across pre-defined ORFs
    """

    help_text = f"""    
        ribofy detect - version {__version__}

        required arguments:
        --bam <file>            Bam file - sorted and indexed
        --orfs <file>           orfs - generated by ribofy orfs
        --prefix <str>          output prefix, default=ribofy   
        
        usage: ribofy detect --bam BAM --orfs ORFS [--prefix PREFIX]"""

     
    parser = argparse2 (
        description=info_text,
        usage=help_text,
        help=help_text
    )
    parser.add_argument('detect', nargs='?', help='') # dummy argument
    parser._action_groups.pop()

    # required
    required = parser.add_argument_group('required arguments')
    required.add_argument("--bam", dest='bam', required=True, help="Bam file - sorted and indexed")
    required.add_argument("--orfs",   dest='orfs', required=True, help="orfs - generated by ribofy orfs")
    required.add_argument("--prefix", dest='prefix', default = "ribofy", help="output prefix, default=ribofy")

    args = parser.parse_args()

    offset = f"{args.prefix}_offsets.txt"
    phasing = f"{args.prefix}_phasing.txt"
    result = f"{args.prefix}_results.txt"

    get_offset (args.bam, args.orfs, offset)
    get_phasing (args.bam, args.orfs, offset, phasing)
    get_results (phasing, result)


if __name__ == "__main__":
    ribofy_detect ()
   
