## Virus agnostic (but mostly dengue) analysis pipeline

Created by Verity Hill and Chrispin Chaguza, Grubaugh Lab

Still a bit beta - please feel free to have a go and report any issues you get! Please provide a screenshot of the error and all log files produced so that I can see what's gone wrong.

**This pipeline is still in active development including argument names and input files. Use ```git pull``` and reinstall using ```pip install .``` every time you run it**

This pipeline takes raw read data in the form of fastq files, maps them against provided bed files and then provides a series of outputs for further analysis including consensus sequences. 

IMPORTANT: the bed files must correspond to the wet lab protocol that you used and the reference sequence used to generate them otherwise the sequences generated will be incorrect. 

It estimates input files virus types based on the minhash algorith as implemented in [Mash](https://github.com/marbl/Mash). This is confirmed if it has more than 50% coverage of the reference genome provided.

If running on a server, it is highly recommended to run using screen/tmux or similar.
Running with apptainer / singularity is recommended, but local install instructions are also given 

# Installation instructions (apptainer / singularity)

1. Clone this repo by going to your command line interface and typing ```git clone https://github.com/sethnr/pathag_pipeline.git```
2. Navigate to your local copy of the repo by typing ```cd pathag_pipeline```
3. Install packages from the repo by typing ```pip install .```

# Installation instructions (local)

1. Clone this repo by going to your command line interface and typing ```git clone https://github.com/sethnr/pathag_pipeline.git```
2. Navigate to your local copy of the repo by typing ```cd pathag_pipeline```
3. Create the conda environment by typing ```conda env create -f environment.yml``` You will need conda to be installed first, installation instructions found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
4. Activate the environment by typing ```conda activate analysis_env```
5. Install packages from the repo by typing ```pip install .```


# Running the pipeline

```pathag_pipeline -h``` will give you all of the usage instructions (information on these options at the bottom of this readme)

### Example commands

```pathag_pipeline --indir <input-directory> --reference-directory <reference-directory> --singularity --slurm```

### Main inputs

#### 1. Fastq files
As a minimum you need fastQ files to analyse. There are two ways you can provide these:

**Option A** (most people!): The fastq files must be separately in their own folder named by sample. In each file, must have the forward and reverse fastq files, defined by them containing "R1" and "R2" somewhere in the name and with the sample name at the start of the file name. See example input file for more information.

For the second option, you can use ``--indir`` to indicate where the folders of samples are kept. Default is the same as the output directory (this is for when you already have a input/output folder where you're working)

NB sample names are found by looping through directories in the input directory. The script ignores temporary and download folders, but will get confused if there are non-sample directories in there other than that.

**Option B** (for Yale users): If you are running on the Yale cluster using their symlinks, simply provide the symlink emailed to you by YCRC (the second half of the link) using ``--symlink`` and the pipeline will deal with it for you.
NB you must run ```module load ycga-public``` in the command line before using this option

#### 2. Reference sequences and bed files
To get consensus files for dengue if you are using our sequencing protocol (https://www.protocols.io/view/dengueseq-a-pan-serotype-whole-genome-amplicon-seq-kqdg39xxeg25/v2), you don't need anything else. 

If you want to try other viruses, or use your own reference and bed files:

- Look at our directory "DENV_primers_and_refs" for formatting file names etc
- Provide the stem of each file in the "refs.txt" tesxt file in the same folder
- Use the ``--reference-directory`` option to provide the path to the directory. 


#### 3. Config file
You can also provide all of your arguments using a config file. This is specified using the ```--config``` option. An template can be found in the main repository. Any command line arguments you specify will overwrite config file arguments.


### Main outputs

Specify the main output directory using ``--outdir``. Default is the generation of a new folder with today's date.

Within this, there will be separate folders containing:
 - <outdir>/results/mash      -> target genome estimates
 - <outdir>/results/bams      -> alignment bams
 - <outdir>/results/consensus -> consensus sequence
 - <outdir>/results/align     -> aligned consensus
 - <outdir>/results/variants  -> sample SNVs
 - <outdir>/results/summary   -> summary files & QC
   

1. results:
	- virus_calls.tsv: Contains virus calls per sample. I.e. those viruses with which the sample has more than 50% coverage
	- top_virus_all_samples.tsv: Contains the highest coverage virus per sample, regardless of coverage
	- summary_all_samples.tsv: contains information about all possible virus options provided per sample
	- alignment - contains alignments by virus type NB not to be used for phylogenetics because it is only rough for estimating coverage. If a trimmed bed file was provided then these are trimmed down, otherwise they are only untrimmed.
	- consensus - consensus sequences of the called virus for each sample
	- depth - text files of each position of the genome and their sequencing depth by sample. NB positions are relative to reference sequence
	- variants - contains a summary file of the number of variants by sample, and then a file for each sample containing additional information about variants

2.  Within results, there are some QC plots:

	- variant_plot: plots sample id against the number of variants (i.e. mutations compared to the reference genome) in between 20% and 80% of the reads. Useful for detecting co-infections
	- ct_plot: plots genome coverage against Ct value coloured by call.
	
	 To make the Ct plot, you must provide:
	 	- File containing the Ct values using ``--ct-file``
	 	- the name of the column containing Ct values with ``--ct-column``
	 	- the name of the column containing IDs which match the ID names on the fastq files/directories using ``--id-column``

3. Log_files:

	Log files for each step of the pipeline, named by snakemake rule and sample - mostly useful for debugging.


	
	




### All options

``--config`` name of optional config file to specify all following options

``--dry-run`` just run set up and not the full pipeline to check that it all works properly


``--outdir`` location where files will be stored

``--indir`` directory containing samples. Each sample must be a folder with the forward and reverse runs in. Default is same as output directory.


``--reference-directory`` or ``-rd`` location where the bed files for primer trimming and associated reference sequences are. Default is the dengue directory provided as package information

``--depth`` minimum depth to call consensus. Default is 10


``--temp`` keep temporary files

``--tempdir`` location of temporary files. Default is folder in output directory called "temporary_files"

``--download`` produce downloads directory (see above)


``--slurm`` parallelise on a cluster which uses slurm. The outer pipeline will run on the login node, but all actual analysis will be submitted to compute nodes.

``--slurm-cores`` number of slurm cores to use default is 10


``--overwrite`` delete old results files and make new ones in the run

``--verbose`` print more stuff to screen

``--help`` print help message to screen and quit


#### Yale specific options

``--symlink`` for use on Yale's HPC. Provide the second half of the link emailed by genomics core. 

    
    
