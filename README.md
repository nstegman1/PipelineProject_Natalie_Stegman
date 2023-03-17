# PipelineProject_Natalie_Stegman

Hello everyone! This is track 2 of the Pipeline Project.

These are the tools that must be installed before running this pipeline:
  - SRA Tool Kit
  - Bowtie2
  - Spades
  - Linux/Unix
  - Python3

# Before you run:
  - Included in this repo are the 8 truncated fastq files
  - You must download these files to run the test data
  - These reads include the first 70,000 reads of each fastq file

# Command to run:

- BEFORE RUNNING, you must be in the location where the trimmed_fastq file folder is, which is PipelineProject_Natalie_Stegman
- **git clone https://github.com/nstegman1/PipelineProject_Natalie_Stegman.git**
- **cd PipelineProject_Natalie_Stegman**
- **python3 pipeline_project.py**
- To view the log file of the test data, make sure you are in the PipelineProject_Natalie_Stegman/PipelineProject_Natalie_Stegman directory.
            
  
# The files in the produced file include:
  - spades_raw_assembly: This is the assembly produced by spades using the fastq files
  - mapped_strains: These are the .sam files
  - bowtie_indexes: These files are the bowtie index
  - blast: This includes the data needed to blast the longest contig in the assembly against the Betaherpesvirinae family
  - aligned_reads: These are the filttered fastq files with reads that map to the bowtie index
  - aligned_reads_unzipped:These are the aligned reads unzipped
  - PipelineProject.log: This includes the output asked for in the project description
  - bowtie_fasta: This is the fasta file used for the Bowtie index
