
************************************************************************************************
Data cleaning


i) What is does ?

   The code is use to remove all the non-acgt charaters including annotationas and delimeters from a genome sequence. 

ii) System commands to do the task

   To compile the source code 'data_cleaning.c':    $gcc data_cleaning.c -o data_cleaning -DMAX=0 
  
   (The option -DMAX=0 indicates acgt characters in the sequence are kept unchanged)

   To run the executable code 'data_cleaning':   $./data_cleaning Genome.fna Genome.fna.txt

   'Genome.fna' is input file. It contains raw genome sequence that includes the non-acgt characters, annotationas and delimeters.
    The output will be stored in a file 'Genome.fna.txt' or any user specified filename, which will be the clean sequence containing only the acgt characters. 



****************************************************************************************************
Transforming genomic sequences into their corresponding TRIM vectors


     To compile the source code 'TRIM.c':   $gcc TRIM.c -DMAX=3900000000 -lm -o TRIM

     (The option -DMAX indicates maximum length of genome the program will accomodate in RAM, which is presently 
     the maximum available genome length of the species taken from in NCBI repository under this study. 
     For any genome upto this length, user need not specify length of the experimenting sequence.)
 
     To run the executable file 'TRIM':   $./TRIM 
     
     A typical output
     
     .............................................
     Enter motif/k-mer length: 7
     Enter the number of genomes to convert in TRIM vector: 4
     Input file name: Genome1.fna.txt
     Constructing TRIM vector for genome stored in file: Genome1.fna.txt
     Storing the TRIM vector in: TRIM_vector_Genome1.txt
     Input file name: Genome2.fna.txt
     Constructing TRIM vector for genome stored in file: Genome2.fna.txt
     Storing the TRIM vector in: TRIM_vector_Genome2.txt
     Input file name: Genome3.fna.txt
     Constructing TRIM vector for genome stored in file: Genome3.fna.txt
     Storing the TRIM vector in: TRIM_vector_Genome3.txt
     Input file name: Genome4.fna.txt
     Constructing TRIM vector for genome stored in file: Genome4.fna.txt
     Storing the TRIM vector in: TRIM_vector_Genome4.txt
     
     

     

