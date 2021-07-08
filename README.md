
# ribofy: ORF detection using RiboSeq data

Ribofy is a fast and simple python-based tool for detection of phased p-sites across open-reading-frames (ORFs)

  

## Installation

* pip (soon)
``` 
pip install ribofy
```

* from source
``` 
git clone https://github.com/ncrnalab/ribofy.git
cd ribofy
python setup.py install
```

## Running ribofy
  

First, all ORFs are assembled from an annotation file (preferably [gencode](https://www.gencodegenes.org/) GTF) and the corresponding genome fasta (should not take more than 5-10 minutes). This is only required once per genome/annotation:

```
ribofy orfs --gtf <path/to/gtf> --fa <path/to/fasta>
```
 
The genome fasta-file must be indexed prior to ORF assembly:
```
samtools faidx <path/to/fasta>
```
 
Currently, ribofy is compatible with STAR, kallisto and salmon mapped reads. Recommended mapping commands:

* STAR
```
STAR --genomeDir <path/to/STAR_index> --outSAMtype BAM SortedByCoordinate 
--readFilesIn <path/to/fastq_files> --readFilesCommand zcat --outFileNamePrefix <sample_prefix>. 
```
* salmon  
```
salmon quant -i </path/to/salmon_index -l A </path/to/fastq_files> --gcBias --validateMappings {additional_params} --writeMappings=</path/to/output.bam> -o </path/to/output>
```
* kallisto
```
kallisto quant -i </path/kallisto_index> --bias -o </path/to/output> --single --pseudobam --fr-stranded -l 30 -s 2 </path/to/fastq_files>
```

Note that for kallisto and salmon, genome indexing should be performed with reduced k-mer value to allow mapping of <30nt ribosome-protected fragments.

Before running ribofy, bam-files should be sorted and indexed:
```
samtools sort </path/to/bamfile> > </path/to/sorted/bam_file>
samtools index </path/to/sorted/bam_file>
```
 

Then, run ribofy:

```
ribofy detect --orfs <path/to/orf/assembly> --bams </path/to/bamfiles> --prefix <prefix>
```
 

## Under the hood

1) Ribofy infers the p-site offsets for read-lengths between 25 and 35 (although this can be customized) and outputs the \<prefix>.offset.txt

2) Then, for each ORF, ribofy counts the p-sites and evaluates the statistical enrichment of in-frame p-sites. This outputs the \<prefix>.phasing.txt

3) Finally, Ribofy collects the individual ORFs into ORF-groups (collapsing overlapping and correlating ORFs), preserving only the highest expressed ORF (based on overall coverage), performs ORF-type specific FDR corrections and outputs the final \<prefix>.results.txt

## Citation
*in preparation*

## Contact
Thomas Hansen (tbh@mbg.au.dk)