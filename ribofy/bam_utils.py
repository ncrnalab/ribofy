import pysam


def get_tid_info (bamfile):

    dtid2count = {}
    dtid2ref = {}

    for idx in pysam.idxstats(bamfile).split ("\n"):
        idx_split = idx.strip ("\n").split ("\t")    
        if len (idx_split) >= 3:
            tid = idx_split[0].split("|")[0]

            dtid2count[tid] = int(idx_split[2]) + int(idx_split[3])
            dtid2ref[tid] = idx_split[0]

    return (dtid2count, dtid2ref)
    