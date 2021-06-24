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
import argparse

from . import __version__
from .ribofy_detect import ribofy_detect
from .ribofy_orfs import ribofy_orfs

def main ():

    text = """
        ribofy <CMD> [arguments] ...\n
        Where <CMD> can be one of:
        \torfs\t\tbuilds a compilation of ORFs from gtf annotation file
        \tdetect\t\tdetects phased riboseq reads across pre-established ORFs (from riboseq orfs)
    """


    if len (sys.argv) < 2:
        print (text)
        sys.exit (0)

    if sys.argv[1] == "orfs":
        ribofy_orfs ()
    
    elif sys.argv[1] == "detect":
        ribofy_detect ()

    elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
        print ("ribofy version", __version__)

    else:
        print (text)
        sys.exit (0)


if __name__ == "__main__":
    main ()
   
