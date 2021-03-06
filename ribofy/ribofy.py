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
from .ribofy_plot import ribofy_plot
from .ribofy_raw import ribofy_raw
from .get_offset import ribofy_offset
from .get_phasing import ribofy_phasing
from .get_results import ribofy_results



def main ():

    text = f"""
        ribofy - version {__version__}

        usage: ribofy <CMD> [arguments] ...

        Where <CMD> can be one of:
            orfs            builds ORFs from gtf annotation
            detect          detects phased riboseq reads across pre-build ORFs
            plot            plots the p-site distibution for specified transcripts\n\n"""

    if len (sys.argv) < 2:
        print (text)
        sys.exit (0)

    if sys.argv[1] == "orfs":
        ribofy_orfs ()
    
    elif sys.argv[1] == "detect":
        ribofy_detect ()

    elif sys.argv[1] == "offset" or sys.argv[1] == "offsets":
        ribofy_offset ()

    elif sys.argv[1] == "phasing":
        ribofy_phasing ()

    elif sys.argv[1] == "results":
        ribofy_results ()

    elif sys.argv[1] == "plot":
        ribofy_plot ()

    elif sys.argv[1] == "raw":
        ribofy_raw ()

    elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
        print ("ribofy version", __version__)

    else:
        print (text)
        sys.exit (0)


if __name__ == "__main__":
    main ()
   
