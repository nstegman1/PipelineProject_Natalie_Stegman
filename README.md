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
  - These reads include the first 75,000 reads of each fastq file

# Command to run:
**git clone insert link here**
- after you download the truncated dataset and did git clone, you can now run the pipeline:
            **python3 pipeline_project.py**
  
# The files in the produced file include:
  - Trimmed_fastq_files: These are the truncated fastq files so that the pipeline runs under 2 minutes
  - spades_raw_assembly: This is the assembly produced by spades using the fastq files
  - mapped_strains: These are the .sam files
  - bowtie_indexes: These files are the bowtie index
  - blast: This includes the data needed to blast the longest contig in the assembly against the Betaherpesvirinae family
  - aligned_reads: These are the filttered fastq files with reads that map to the bowtie index
  - aligned_reads_unzipped:These are the aligned reads unzipped
  - PipelineProject.log: This includes the output asked for in the project description
  - bowtie_fasta: This is the fasta file used for the Bowtie index
