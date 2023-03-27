## Virus agnostic (but mostly dengue) analysis pipeline

Created by Verity Hill and Chrispin Chaguza, Grubaugh Lab

This pipeline takes raw read data in the form of fastq files, maps them against provided bed files and then provides a series of outputs for further analysis including consensus sequences.

It calls input files as a virus type if it has more than 50% coverage of the reference genome provided.

### Installation instructions

- Clone this repo
- conda env create -f environment.yml
- pip install .

If step 2 fails on a server because of a "bus error", then first run the command "salloc" to request more memory


### Main inputs

As a minimum you need fastQ files to analyse. There are two ways you can provide these:

1. If you are running on the Yale cluster using their symlinks, simply provide the symlink emailed to you by YCRC (the second half of the link) using ``--symlink`` and the pipeline will deal with it for you.
2. Otherwise, must be separately in their own folder named by sample. In each file, must have the forward and reverse fastq files, defined by them containing "R1" and "R2" somewhere in the name. See example input file for more information.


- For the second option, you can use ``--indir`` to indicate where the folders of samples are kept. Default is the same as the output directory (this is for when you already have a input/output folder where you're working)
- NB sample names are found by looping through directories in the input directory. The script ignores temporary and download folders, but will get confused if there are non-sample directories in there other than that.

To get consensus files for dengue, you don't need anything else. If you want to try other viruses you need some bed files to compare the sequencing data to:

- Here, we provide each of the four dengue virus serotypes as package data, and these are used as default. 
- Otherwise, please use the same format and provide the stem of the files in a "refs.txt" text file in the same folder. Use the ``--primer-directory`` option to provide the path to the directory.


### Main outputs

Specify the main output directory using ``--outdir``. Default is the generation of a new folder with today's date.

Within this, there will be:

1. results:
	- DENV.serotype.calls.final.tsv: Contains virus calls per sample which have more than 50% coverage
	- DENV.top.serotype.cals.all.samples.tsv: Contains all top calls per sample, regardless of coverage
	- summary.all.samples.tsv: contains information about all possible options provided per sample
	- alignment - contains trimmed and untrimmed alignments by virus type NB not to be used for phylogenetics because it is only rough for estimating coverage
	- consensus - consensus sequences of the called virus for each sample
	- depth - text files of each position of the genome and their sequencing depth by sample
	- variants - contains a summary file of the number of variants by sample, and then a file for each sample containing additional information about variants

2. temporary_files (if option ``--temp`` is used):

	All files produced in the run of the pipeline that don't make it to results.
	This includes intermediate files produced by software, and bam/fasta/variant files for virus comparisons which did not meet the threshold for the sample to be called as that virus
	
	Mostly for debugging purposes
	
3. downloads (if option ``--download`` is used):

	The same as results, but without the bam files. 
	
	For download and storage in places with less data storage capacity (eg dropbox)
	
4. Within results, there are optional QC plots produced.

	- ct_plot: plots genome coverage against Ct value coloured by call. Requires additional input of csv or tsv containing sample name and ct value. 
	- variant_plot: plots sample id against the number of variants (i.e. mutations compared to the reference genome) in between 20% and 80% of the reads. Useful for detecting co-infections


### All options

``--outdir`` location where files will be stored

``--indir`` directory containing samples. Each sample must be a folder with the forward and reverse runs in. Default is same as output directory.


``--primer-directory`` or ``-pd`` location where bed files etc for references are. Default is the dengue directory provided as package information

``--depth`` minimum depth to call consensus. Default is 20


``--temp`` keep temporary files

``--tempdir`` location of temporary files. Default is folder in output directory called "temporary_files"

``--download`` produce downloads directory (see above)


``--overwrite`` delete old results files and make new ones in the run

``--verbose`` print more stuff to screen

``--help`` print help message to screen and quit


#### Yale specific options

``--symlink`` for use on Yale's HPC. Provide the second half of the link emailed by genomics core. 

``--slurm`` flag for use on Yale's HPC. This will parallelise the analysis. This uses DSQ (dead simple queue) to run the jobs text file. If your HPC uses DSQ and slurm, this should also work for you.

    
    
