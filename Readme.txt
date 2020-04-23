
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
TRIM vector computation 


i) What is does ?

   Transforming genomic sequences into their corresponding TRIM vectors
 
ii)  System commands to do the task

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
     
     

*****************************************************************************************************

Trainig the classifier TRIMEC

i) What is does ?

  This code was used to train the TRIMEC. However, it is not essential to run the taxa identification utility. We have uploaded 
  and provided the instruction about how to use pre-built training in the taxa identification section. So user can skip this 
  steps.
  
ii)  System commands to do the task  
  
To compile the source code 'Bagging_classifier_v1_Building.c': $gcc Bagging_classifier_v1_Building.c -lm -o TRIM_classifier_Bulder

To run the executable file 'TRIM':   $./TRIM_classifier_Bulder

The output of the runtime instance used to create the file 'Bagging_claffifier_k=7.txt':

...........................................
Enter the number of TRIM Vectors to loan in-meory: 218
Enter the motif lengh used to construct TRIM Vectors: 7
Collecting prior-knowledge aout "class-label" of TRIM Vectors
Provide the range of indices of TRIM Vectors for the class Amphibian
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 207 212 4
Provide the range of indices of TRIM Vectors for the class Reptile
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training  213 216 2
Provide the range of indices of TRIM Vectors for the class Bird
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 68 93 20
Provide the range of indices of TRIM Vectors for the class Insect
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 1 3  3
Provide the range of indices of TRIM Vectors for the class Fish
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 4 67 55
Provide the range of indices of TRIM Vectors for the class Mammal
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 94 206 90
Provide the range of indices of TRIM Vectors for the class Non-animal
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 217 218 2
Classifier with motif length %d has been written to file Bagging_claffifier_k=7.txt













*******************************************************************************************************
Taxa identification

i) What is does ?

  It is used to identify taxa of a given genome sequence and do not require any other inforation.
 
ii)  System commands to do the task

 To compile the source code 'Bagging_classifier_v1.c': $gcc Bagging_classifier_v1.c -DAMX=38900000000 lm -o TRIMEC
 
 To run the executable file 'TRIMEC': $./TRIMEC
 
 Supplementary_File_F1.docx contains several runtime instatnces of TRIMEC.
  


