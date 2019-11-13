#data_cleaning
To compile the source code data_cleaning.c

$gcc data_cleaning.c -o data_cleaning

To run the executable code data_cleaning

$data_cleaning unclean.fna input.txt

where unclean.fna contains genome sequene along with non-acgt characters, annotationas and delimeters
      input.txt contains genomes only with acgt characters.

# TRIM
To compile the source code TRIM.c

$gcc TRIM.c -DMAX=3900000000 -lm -o TRIM

To run the executable file TRIM

$TRIM input.txt -k 7 > output.txt

where input.txt contains the genomce sequence 
      output.txt contains the TRIM vector corresponds to the genome sequence.

