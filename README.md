#data_cleaning

To compile the source code data_cleaning.c

$gcc data_cleaning.c -o data_cleaning

To run the executable code data_cleaning

$data_cleaning raw_data.fna clean_data.txt

'raw_data.fna' contains genome sequene along with non-acgt characters, annotationas and delimeters.
'clean_data.txt' contains genome sequence, only with acgt characters.
      

# TRIM
To compile the source code TRIM.c

$gcc TRIM.c -DMAX=3900000000 -lm -o TRIM

To run the executable file TRIM

$TRIM clean_data.txt -k 7 > output.txt

'clean_data.txt' contains the genomce sequence 
'output.txt' contains the TRIM vector corresponds to the genome sequence.

