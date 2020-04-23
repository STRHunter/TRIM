
************************************************************************************************
Data cleaning


i) What is does ?

   The code is use to remove all the non-acgt charaters including annotationas and delimeters from a genome sequence. 

ii) System commands to do the task

   To compile the source code 'data_cleaning.c':    $gcc data_cleaning.c -o data_cleaning -DMAX=0 
  
   (The option -DMAX=0 indicates acgt characters in the sequence are kept unchanged)

   To run the executable code 'data_cleaning':   $data_cleaning Genome.fna Genome.fna.txt

   'Genome.fna' is input file. It contains raw genome sequence that includes the non-acgt characters, annotationas and delimeters.
    The output will be stored in a file 'Genome.fna.txt' or any user specified filename, which will be the clean sequence containing only the acgt characters. 



****************************************************************************************************
Transforming genomic sequences into their corresponding TRIM vectors


     To compile the source code 'TRIM.c':   $gcc TRIM.c -DMAX=3900000000 -lm -o TRIM

     (The option -DMAX indicates maximum length of genome the program will accomodate in RAM, 
      which is presently the maximum available genome length of the species taken from in NCBI repository under this study. 
      For any genome upto this length, user need not specify length of the experimenting sequence.)
 
     To run the executable file 'TRIM':   $TRIM Genome.fna.txt -k 7 > TRIM_Vector_Genome.txt

     'Genome.fna.txt' is the input file containing the genome sequence with only acgt characters.
      The output of the program is the corresponding TRIM vector which will be stored in the output file 'TRIM_Vector_Genome.txt'.

      The option -k 7 indicates motif length of Tandem Repeat is 7nt.

