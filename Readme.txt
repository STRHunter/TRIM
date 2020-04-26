
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

Note: Cleaning genome sequence is mandatory for correct working of all other executables involved in following computations. 

****************************************************************************************************
TRIM vector computation 


i) What is does ?

   Transforming genomic sequences into their corresponding TRIM vectors. Here we assume that all the genome sequences 
   have been cleaned by the 'data_cleaning' utility.  
 
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
     
     Note: The utility 'TRIM' autimaically assigns the filename 'TRIM_vector_GenomeN.txt' and store the TRIM vector correspond 
     to the Nth input genome. Input file name can be any user specified filenane that contain a clean genomic sequence.
     
     
     
     We have used the utility to generate 218 TRIM vectos for 296 animal genomes and 2 non-animal genomes at k=7 as 
     
     TRIM_vector_Genome1.txt to TRIM_vector_Genome3.txt : Taxa Insect
     TRIM_vector_Genome4.txt to TRIM_vector_Genome67.txt : Taxa Fish
     TRIM_vector_Genome68.txt to TRIM_vector_Genome93.txt : Taxa Bird
     TRIM_vector_Genome94.txt to TRIM_vector_Genome206.txt : Taxa Mammal
     TRIM_vector_Genome207.txt to TRIM_vector_Genome212.txt : Taxa Amphibian
     TRIM_vector_Genome213.txt to TRIM_vector_Genome216.txt : Taxa Reptile
     TRIM_vector_Genome217.txt to TRIM_vector_Genome218.txt : Taxa Non_animal
                           
                                and
                                
    We have used the utility to generate 91 TRIM vectos for 78 animal genomes and 13 non-animal genomes at k=6 as
    
    TRIM_vector_Genome1.txt to TRIM_vector_Genome59.txt : Taxa Insect
    TRIM_vector_Genome60.txt to TRIM_vector_Genome78.txt : Taxa Fish
    TRIM_vector_Genome79.txt to TRIM_vector_Genome91.txt : Taxa Non-animal

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
Classifier with motif length 7 has been written to file Bagging_claffifier_k=7.txt




The output of the runtime instance used to create the file 'Bagging_claffifier_k=6.txt':

.......................................................
Enter the number of TRIM Vectors to loan in-meory: 91
Enter the motif lengh used to construct TRIM Vectors: 6
Collecting prior-knowledge aout "class-label" of TRIM Vectors
Provide the range of indices of TRIM Vectors for the class Amphibian
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 0 0 0
Provide the range of indices of TRIM Vectors for the class Reptile
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training  0 0 0
Provide the range of indices of TRIM Vectors for the class Bird
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 0 0 0
Provide the range of indices of TRIM Vectors for the class Insect
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 1 59  50
Provide the range of indices of TRIM Vectors for the class Fish
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 60 78 15
Provide the range of indices of TRIM Vectors for the class Mammal
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 0 0 0
Provide the range of indices of TRIM Vectors for the class Non-animal
Enter lower and upper indices along with the number of TRIM Vectors for training: Put zeros to skip training 79 91 10
Classifier with motif length 7 has been written to file Bagging_claffifier_k=6.txt



*******************************************************************************************************
Taxa identification

i) What is does ?

  It is used to identify taxa of a given genome sequence and do not require any other inforation.
 
ii)  System commands to do the task

 To compile the source code 'Bagging_classifier_v1.c': $gcc Bagging_classifier_v1.c -DAMX=38900000000 lm -o TRIMEC
 
 To run the executable file 'TRIMEC': $./TRIMEC
 
 A typical instance of the TRIMEC is as 
 
 .......................................................
 ...........................Welcome to genome classifier utility.....................

Input file name of the genome to identify its taxa: 
...........................System is performing the following tasks................

Constructing TRIM vector for the genome stored in file: D1.txt ... 
at k=6
20276.291245
at k=7
28686.171945
sum 20276.291245 
0.000000 
sum 28686.171945 
0.000000 

Total entropy maximizes at k=7
Locading corresponding classifier....

Instance 0   16384
Instance 1   16384
Instance 2   16384
Instance 3   16384
Instance 4   16384
Instance 5   16384
Instance 6   16384
Instance 7   16384
Instance 8   16384
Instance 9   16384
Instance 10   16384
Instance 11   16384
Instance 12   16384
Instance 13   16384
Instance 14   16384
Instance 15   16384
Instance 16   16384
Instance 17   16384
Instance 18   16384
Instance 19   16384
Instance 20   16384
Instance 21   16384
Instance 22   16384
Instance 23   16384
Instance 24   16384

0.068795 
0.074580 
0.085106 
0.317337 
0.075354 
0.075393 
0.500000 

0.028597 
0.073390 
0.087673 
0.300203 
0.074356 
0.073517 
0.500000 

0.046292 
0.073192 
0.086085 
0.365971 
0.075493 
0.074496 
0.500000 

0.034589 
0.075352 
0.085452 
0.321611 
0.076503 
0.075101 
0.500000 

0.043827 
0.063790 
0.084590 
0.300203 
0.078164 
0.074720 
0.500000 

0.068795 
0.070581 
0.088238 
0.326818 
0.080512 
0.075077 
0.500000 

0.064444 
0.070353 
0.087800 
0.300203 
0.075759 
0.073566 
0.500000 

0.040569 
0.067919 
0.085370 
0.300203 
0.073270 
0.075941 
0.500000 

0.071189 
0.070353 
0.087103 
0.306424 
0.075623 
0.074947 
0.500000 

0.057695 
0.073390 
0.085872 
0.309254 
0.075374 
0.073982 
0.500000 

0.075128 
0.067919 
0.086699 
0.326818 
0.077091 
0.075563 
0.500000 

0.087444 
0.083584 
0.085746 
0.306424 
0.075192 
0.076949 
0.500000 

0.080872 
0.080793 
0.087114 
0.321611 
0.074348 
0.076758 
0.500000 

0.082834 
0.075352 
0.084662 
0.300203 
0.075373 
0.074815 
0.500000 

0.048884 
0.073192 
0.086113 
0.321611 
0.077589 
0.076411 
0.500000 

0.056742 
0.088158 
0.084874 
0.326818 
0.075996 
0.076420 
0.500000 

0.037962 
0.070581 
0.083478 
0.300203 
0.077337 
0.074267 
0.500000 

0.045507 
0.069565 
0.087489 
0.309254 
0.076216 
0.075066 
0.500000 

0.075060 
0.070353 
0.085335 
0.300203 
0.075295 
0.075401 
0.500000 

0.089092 
0.068500 
0.086290 
0.321611 
0.075817 
0.075501 
0.500000 

0.072141 
0.083584 
0.083717 
0.341782 
0.078701 
0.074996 
0.500000 

0.058619 
0.067919 
0.087616 
0.321611 
0.079172 
0.075316 
0.500000 

0.085734 
0.067919 
0.084319 
0.300203 
0.074794 
0.074626 
0.500000 

0.041134 
0.069565 
0.085541 
0.321611 
0.076691 
0.074636 
0.500000 

0.037962 
0.068500 
0.083992 
0.300203 
0.074171 
0.075015 
0.500000 
This is a Amphibian
  


