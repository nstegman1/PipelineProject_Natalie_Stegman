import os
from Bio import Entrez
from Bio import SeqIO
import glob
import subprocess


#This hold the base cwd that we put the project folders in 
foundation = os.getcwd()

#this gets the path the pipeline project folder
base_path = os.getcwd() 
base_path = base_path+'/PipelineProject_Natalie_Stegman'
os.makedirs(base_path)

#this opens the log file that will be written into
log_file = open(base_path+'/PipelineProject.log', 'w+')



#This gets the bowtie index of HCMV 
def get_index(cwd):
    
    #this gets the location path of folder of the bowtie index
    file_location = cwd+'/bowtie_fasta'
    
    #This makes the directory that will hold the index files
    os.makedirs(cwd+'/bowtie_indexes')
    index_file_location = cwd+'/bowtie_indexes'
    
    #this searched for the for the HCMV fasta file. and 
    #writes it into the file location
    handle = Entrez.efetch(db="nucleotide", id= 'NC_006273.2', rettype='fasta')
    records = list(SeqIO.parse(handle, 'fasta'))   
    SeqIO.write(records, file_location, "fasta")
    
    #this calls the command to build the bowtie index based on the HCMV fasta file
    print('bowtie2-build '+file_location+' '+index_file_location+'/HCMV')
    os.system('bowtie2-build '+file_location+' '+index_file_location+'/HCMV')
    
    
def map_strains(cwd):
    
    #this will hold the trimmed fastq reads
    reads = []
    
    #holds the path of the folder of the trimmed fastq files
    path = os.getcwd()+'/trimmed_fastq'
    
    #makes the folders to hold the mapped and aligned strains
    os.makedirs(cwd+'/mapped_strains')
    os.makedirs(cwd+'/aligned_reads')
 
    #for each file in the trimmed fastq file folder
    for filename in os.listdir(path):
        
        #make the entire path of the file name 
        f = os.path.join(path, filename)
        
        # checking if it is a file
        if os.path.isfile(f):
            
            #adds the file name to the reads list
            reads.append(filename)
    
    #sorts the reads list to that the file names are in order
    reads = sorted(reads)
    print(reads)
    
    #for each pair of reads, 
    for i in range(0, len(reads),2):
        
        #makes the name of the sam file
        file_name = cwd+'/mapped_strains/'+reads[i][0:10]+'.sam'
        
        #this is the path to the HCMV indexes
        index_path = cwd+'/bowtie_indexes/HCMV'        
        print(file_name)
        
        #This makes the bowtie command and calls it 
        #command = 'bowtie2 --quiet -x '+index_path+' -1 '+cwd+'/trimmed_fastq/'+reads[i]+' -2 '+cwd+'/trimmed_fastq_files/'+reads[i+1]+' -S '+file_name+' --al-conc-gz '+cwd+'/aligned_reads/'+reads[i][0:10]+'_mapped_%.fq.gz'
        command = 'bowtie2 --quiet -x '+index_path+' -1 '+path+'/'+reads[i]+' -2 '+path+'/'+reads[i+1]+' -S '+file_name+' --al-conc-gz '+cwd+'/aligned_reads/'+reads[i][0:10]+'_mapped_%.fq.gz'
        print(command)
        os.system(command)
        
#this counts the reads before and after the bowtie filtering
def count_reads(cwd):
    
    #list of files before 
    before_path = os.getcwd()+'/trimmed_fastq'
    before_files = []
    before_numbers = []
    
    #This adds the trimmed fastq file names to a list
    for filename in os.listdir(before_path):
        f = os.path.join(before_path, filename)
        before_files.append(filename)
        
    #this sorted the file names so that they are in order
    before_files = sorted(before_files)

    #for each file name before they are filtered     
    for i in before_files:
        
        #get the path of the file
        file = before_path+'/'+i
        
        #open the file and append each read to a list
        with open(file, 'r') as f:
            reads = []
            for line in f:
                if line.startswith("@"):
                    reads.append(line)
        
        #append the number of reads in each file to a list of length reads
        before_numbers.append([i, len(reads)])       
    print(before_numbers)
    
    #list of files after alignment to index
    after_path = cwd+'/aligned_reads'
    
    #first you have to unzip the files
    os.makedirs(cwd+'/aligned_reads_unzipped')
    
    #this iterates through the aligned reads and puts the un
    #zipped filed in a new folder
    for filename in os.listdir(after_path):
        f =  os.path.join(after_path, filename)
        
        #This calls the unzip command
        #print('gunzip -c '+f+' > '+cwd+'/aligned_reads_unzipped/'+filename[0:10]+filename[-8:-6]+'_unzipped.fastq')
        os.system('gunzip -c '+f+' > '+cwd+'/aligned_reads_unzipped/'+filename[0:10]+filename[-8:-6]+'_unzipped.fastq')
    
    #makes an empty list for the filtered file names
    after_files = []
    
    #initialized a list for the read counts after filtering
    after_numbers = []
    
    #gets the path for where the unzipped filed are
    unzipped_path = cwd+'/aligned_reads_unzipped'
    
    #gets list of filenames after filtering
    for filename in os.listdir(unzipped_path):
        f = os.path.join(unzipped_path, filename)
        after_files.append(filename)
    
    #this sorted the file names so they are in order
    after_files = sorted(after_files)

    #gets the number of reads of the filtered file
    for i in after_files:
        file = unzipped_path+'/'+i
        #print(file)
        with open(file, 'r') as f:
            reads = []
            for line in f:
                if line.startswith("@"):
                    reads.append(line)
        
        after_numbers.append([i, len(reads)])
        
    print(after_numbers)
    
    #this writes out the number of reads in the instructed way 
    log_file.write('Donor 1 (2dpi) had '+str(before_numbers[0][1])+' read pairs before Bowtie2 filtering and '+str(after_numbers[0][1])+' read pairs after.'+'\n')
    log_file.write('Donor 1 (6dpi) had '+str(before_numbers[2][1])+' read pairs before Bowtie2 filtering and '+str(after_numbers[2][1])+' read pairs after.'+'\n')
    log_file.write('Donor 3 (2dpi) had '+str(before_numbers[4][1])+' read pairs before Bowtie2 filtering and '+str(after_numbers[4][1])+' read pairs after.'+'\n')
    log_file.write('Donor 3 (6dpi) had '+str(before_numbers[6][1])+' read pairs before Bowtie2 filtering and '+str(after_numbers[6][1])+' read pairs after.'+'\n')
 


def assemble_transcripts(cwd):
    
    #stores the path name of unzipped aligned reads
    path = cwd+'/aligned_reads_unzipped'
    files = []
     
    #add the paths to the file to a list of strings
    for filename in os.listdir(path):
        f = os.path.join(path, filename)
        files.append(f)
    
    #sort this list so that it is in order
    files = sorted(files)
    #print(files)
    
    #FIX: make link
    command_path = '/home/nataliestegman/bioi_tools/SPAdes-3.15.5-Linux/bin/rnaspades.py' 
    output_path = cwd+'/spades_raw_assembly' 

    #this assembles the aligned reads using spades    
    command = 'python3 '+command_path+' --pe-1 1 '+files[0]+' --pe-2 1 '+files[1]+' --pe-1 2 '+files[2]+' --pe-2 2 '+files[3]+' --pe-1 3 '+files[4]+' --pe-2 3 '+files[5]+' --pe-1 4 '+files[6]+' --pe-2 4 '+files[7]+' -o '+output_path
    #command = 'python3 rnaspades.py --pe-1 1 '+files[0]+' --pe-2 1 '+files[1]+' --pe-1 2 '+files[2]+' --pe-2 2 '+files[3]+' --pe-1 3 '+files[4]+' --pe-2 3 '+files[5]+' --pe-1 4 '+files[6]+' --pe-2 4 '+files[7]+' -o '+output_path

    #print(command)
    os.system(command)
    
    #this writes out the command string into the log file
    log_file.write(command+'\n')
    

def num_contigs(cwd):
    
    #this stores the path of the assembly of the 8 reads
    filename = cwd+'/spades_raw_assembly/transcripts.fasta'
    
    #this reads in the fasta file
    records = SeqIO.parse(filename, "fasta")
    base_count = 0
    contig_count = 0
    
    #this counts how many contigs>1000bp that are in the assembly
    #and how many total bp are in the assembly
    for record in records:       
        if len(record.seq) > 1000:
            contig_count+=1
            base_count+=len(record.seq)
            
    #this writes out these numbers nicely
    log_file.write('There are '+str(contig_count)+' contigs > 1000 bp in the assembly.'+ '\n')
    log_file.write('There are '+str(base_count)+' bp in the assembly.'+'\n')
    
    
def blast(cwd):
    
    #this makes the folder to hold the blast data
    os.makedirs(cwd+'/blast')
    
    #path to the spades assembly fasta file
    filename = cwd+'/spades_raw_assembly/transcripts.fasta'
    #this parses the fasta file
    records = SeqIO.parse(filename, "fasta")
    max_len = 0
    max_contig = ''
    
    #this finds the longest contig in transcripts.fasta
    for record in records:
        if len(record.seq) > max_len:
            max_len = len(record.seq)
            max_contig = str(record.seq)
    
    #writes out longest contig into separate fasta file
    longest_contig = open(cwd+'/blast/longest_contig.fasta' ,'w')
    
    #longest_contig.write(max_id)
    #longest_contig.write('\n')
    
    #this writes out the max contig sequence
    longest_contig.write(max_contig)
    longest_contig.close()
    
    #this opens the HCMV fasta file
    HCMV_db = open(cwd + "/blast/HCMV.fasta", "w")
    
    #this gets all the records from the Betaherpesvirinae family
    handle = Entrez.esearch(db = "nucleotide", term = ("Betaherpesvirinae[Organism] OR Betaherpesvirinae[All Fields]"), retmax = 1500) 
    record = Entrez.read(handle)
    
    list_of_ids = record["IdList"]
    handle = Entrez.efetch(db = "nucleotide", id = list_of_ids, rettype = "fasta")
    
    records = list(SeqIO.parse(handle, format = "fasta"))

    for record in records: #for each record, write it to the multifasta file for HCMV db
            HCMV_db.write(">" + str(record.description))
            HCMV_db.write("\n")
            HCMV_db.write(str(record.seq))
            HCMV_db.write("\n")
        
    HCMV_db.close()
    
    
    #DATABASE COMMAND: makes the blast db
    db_command_path = subprocess.check_output('find $(pwd) -name makeblastdb', shell=True)
    print(db_command_path.decode("utf-8"))
    db_path = cwd+'/blast/HCMV.fasta'
    output_db = cwd+'/blast/blastdb/blastdb'
    db_command = 'makeblastdb -in '+db_path+' -out '+output_db+' -dbtype nucl'
    print(db_command)
    os.system(db_command)

    
    #BLAST COMMAND: blasts the longest contig against the made db
    query_path = cwd+'/blast/longest_contig.fasta'
    output_path = cwd+'/blast/blast_outfile'
    blast_command = 'blastn -query '+query_path+' -db '+output_db+" -num_threads 2 -max_hsps 1 -max_target_seqs 10 -out "+output_path+" -outfmt='6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    print(blast_command)
    #subprocess.call(blast_command, shell = True)
    os.system(blast_command)

    #write out to log file
    blast_output = open(cwd+'/blast/blast_outfile','r')
    log_file.write("sacc" + "\t" + "pident" + "\t" + "length" + "\t" +  "qstart" + "\t" + "qend" + "\t" + "sstart" + "\t" + "send" + "\t" + "bitscore" + "\t" + "evalue" + "\t" + "stitle" + "\n")
    log_file.write(blast_output.read())  

get_index(base_path)
map_strains(base_path)
count_reads(base_path)
assemble_transcripts(base_path)
num_contigs(base_path)
blast(base_path)

log_file.close()