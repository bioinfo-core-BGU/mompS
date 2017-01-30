====================================================
Welcome to mompS tool's documentation!
====================================================

**Software author**: Dr. Michal Gordon

**PI**: Prof. Jacob Moran-Gilad

**Citation**: Gordon, M., Yakunin, E., Valinsky, L., Chalifa-Caspi, V., & Moran-Gilad, J. (2017). A Bioinformatics Tool for Ensuring the Backwards Compatibility of Legionella pneumophila Typing in the Genomic Era. Clinical Microbiology and Infection.‏‏

**Date**:   29/1/2017

Why should I need this?
-------------------------

Whole genome sequencing (WGS) has revolutionized the subtyping of L. pneumophila (Lp) but calling the traditional sequence based type from genomic data is hampered by presence of multiple copies of the mompS locus. We propose a novel bioinformatics solution for rectifying that limitation, ensuring the feasibility of WGS for cluster investigation. 

Usage
---------

1. Download the momps.zip folder to your unix
2. Unzip the momps folder
3. In 'config.txt', set the path to blast, bwa, samtools, picard and freebayse
4. cd into the momps folder
5. Run the program as follows::

    perl All_Extract_MLST_Legionella.pl   \
        -f forward_fastq_file.fq   \
        -r reverse_fastq_file.fq   \
        -a assembly_file.fasta   \
        -o output_folder_name/ \
        -p prefix
        
    
Where:

-f        full path to forward fastq file
-r        full path to reverse fastq file
-a        full path to assembly fasta file
-o        output folder for the results
-p        prefix (e.g. sample_1)

Dependencies
-------------


* blast
* bwa
* samtools (version: 1.3.1)
* picard
* freebayse
* perl



