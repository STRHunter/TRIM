Experiments consists of four stages. Stages 1 is the preparatory step, stage 2 is implementation of the algorithm for transforming a genomic sequence
into the corresponding TRIM vector. Stage 3 and 4 are independent applications of TRIM vector into phylogeny construction and genome identication respectively.  

************************************************************************************************
Stage 1: Data cleaning


i) What is does ?

   The code is use to remove all the non-acgt charaters including annotationas and delimeters from a genome sequence. 

ii) System commands to do the task

   To compile the source code 'data_cleaning.c':    $gcc data_cleaning.c -o data_cleaning -DMAX=0 
  
   (The option -DMAX=0 indicates acgt characters in the sequence are kept unchanged)

   To run the executable code 'data_cleaning':   $data_cleaning Genome.fna Genome.fna.txt

   'Genome.fna' is input file. It contains raw genome sequence that includes the non-acgt characters, annotationas and delimeters.
    The output will be stored in a file 'Genome.fna.txt' or any user specified filename, which will be the clean sequence containing only the acgt characters. 



****************************************************************************************************
Stage 2: Transforming a genomic sequence into the corresponding TRIM vector


     To compile the source code 'TRIM.c':   $gcc TRIM.c -DMAX=3900000000 -lm -o TRIM

     (The option -DMAX indicates maximum length of genome the program will accomodate in RAM, 
      which is presently the maximum available genome length in NCBI repository. For any genome 
      upto this length, user need not specify length of the experimenting sequence.)
 
     To run the executable file 'TRIM':   $TRIM Genome.fna.txt -k 7 > TRIM_Vector_Genome.txt

     'Genome.fna.txt' is the input file containing the genome sequence with only acgt characters.
      The output of the program is the corresponding TRIM vector which will be stored in the output file 'TRIM_Vector_Genome.txt'.

      The option -k 7 indicates motif length of Tandem Repeat is 7nt.



*******************************************************************************
Stage 3: Phylogeny tree construction of N number of genome sequences based on TRIM vectors

For all N genomes, Stage 1 and Stage 2 will have to be executated sequentially N times, and the output for the sequences are required to be 
stored in files 'TRIM_Vector_Genome1.txt', 'TRIM_Vector_Genome2.txt', ..., 'TRIM_Vector_GenomeN.txt'.


(a) To load the TRIM vectors in Matlab, user can execute the following Matlab script:


   for i=1:1:N
    myfilename = sprintf('TRIM_Vector_Genome%d.txt',i);
    text = fileread(myfilename);
    arr_text = strsplit(text);
    x = char(arr_text); 
    x_n = str2num(x);
    x_n= x_n';
    TRIM_DB(i,:) = x_n; 
   end

'TRIM_DB' is a 2-D array of N rows and 4^k columns to store the TRIM vectors, each row for each TRIM vector. 

(b) Construction of TRIM distance matrix and generating phylogney tree 

The Matlab script for the task is 

 D =0;
 l =1;
 for i = 1 : k-1 
 for j = i+1 : k
  part1 =0; part2 =0;
  p = TRIM_DB(i,1:16384); q = TRIM_DB(j,1:16384); m = (p+q)/2;
  for index =1:16384
   if m(index) ~=0 && p(index) ~=0
     part1 = part1 + p(index) * log2(p(index)/m(index));
    else
     part1 = part1 + 0;
   end
   if m(index)~=0 && q(index) ~=0 
    part2 = part2 + q(index) * log2(q(index)/m(index));
   else
     part2 = part2 +0;
    end
   end
  D(l) = (part1+part2)/2;
  l = l + 1;
 end
 end
% Constructing phylogeny tree from the TRIM distance matrix uising NJ method. 
% Variable 'phy' contains the phylogney tree object
  phy = seqneighjoin(D,'equivar')
% To get the tree in Newick file format 
  getnewckstr(phy)
% To view the tree
  view(phy)
% end of tree construction


Stage 4: Genome identification based on TRIM vectors 
**********************************************************************************************************

(a) To load the TRIM vectors in Matlab as in the previous step, user needs to execute the following Matlab script:
    (This step is not required if user continues with the TRIM vectors already loaded in Matlab) 

   for i=1:1:N
    myfilename = sprintf('TRIM_Vector_Genome%d.txt',i);
    text = fileread(myfilename);
    arr_text = strsplit(text);
    x = char(arr_text); 
    x_n = str2num(x);
    x_n= x_n';
    TRIM_DB(i,:) = x_n; 
   end

'TRIM_DB' is a 2-D array of N rows and 4^k columns to store the TRIM vectors, each row for each TRIM vector. 

(b) Construction of the taxa specific Mean TRIM Vectors:

% Compute Mean TRIM Vectors for six animal taxa and one non-animal taxa by 

for i =1:1:16384
Amphibian(i) = sum(TRIM_DB(1:5,i)); % TRIM_DB(6,:) is left for test data
MTV(1) = Amphibian/5;

Reptile(i) = sum(TRIM_DB(7:9,i)); % TRIM_DB(10,:) is left for test data
MTV(2)= Reptile/3;

Bird(i) = sum(TRIM_DB(11:30,i)); % TRIM_DB(31,:) to TRIM_DB(36,:) is left for test data
MTV(3)= Bird/20;

Insect(i) = sum(TRIM_DB(37:92,i)); % TRIM_DB(93,:) to TRIM_DB(99,:) is left for test data
MTV(4)= Insect/56;

Fish(i) = sum(TRIM_DB(100:173,i)); % TRIM_DB(174,:) to TRIM_DB(183,:) is left for test data
MTV(5)= Fish/74;

Mammal(i) = sum(TRIM_DB(184:271,i)); % TRIM_DB(272,:) to TRIM_DB(296,:) is left for test data
MTV(6)= Mammal/88;

NonAnimal(i) = sum(TRIM_DB(297:301,i)); 
NonAnimal(i) = sum(TRIM_DB(307:311,i)); % TRIM_DB(302,:) to TRIM_DB(306,:) is left for test data
MTV(7)= NonAnimal/10;

end
% end script


%The following script will report taxa of unknown genome using closeness to respective Mean TRIM Vectors

% Input any of the test data from TRIM_DB(6,:), TRIM_DB(10,:), TRIM_DB(31,:) to TRIM_DB(36,:), 
% TRIM_DB(93,:) to TRIM_DB(99,:), TRIM_DB(174,:) to TRIM_DB(183,:), TRIM_DB(272,:) to TRIM_DB(296,:) 
% and TRIM_DB(302,:) to TRIM_DB(306,:) to check efficiency of TRIM to identify unknown genome

prompt='Enter the number of sequences';
N = input(prompt)
for h=1:1:N

	prompt='Enter the sequence number for taxa identification';
        unknown = input(prompt)
	% variable 'unknown' contains the serial number for unknown genome
	% variable 'taxa' contains the name of the taxa from which TRIM distance of the unknown genome is computed

	for k=1:1:7  
	part1 =0; part2 =0;
	p = TRIM_DB(unknonw,:); q = MTV(k); 
	m = (p+q)/2;
	for index =1:16384
	if m(index) ~=0 && p(index) ~=0
		part1 = part1 + p(index) * log2(p(index)/m(index));
	else
		part1 = part1 + 0;
	end
	if m(index)~=0 && q(index) ~=0 
		part2 = part2 + q(index) * log2(q(index)/m(index));
	else
		part2 = part2 +0;
	end
	end
	d = (part1+part2)/2;
	end
	% variable d contains the TRIM distance between unknown geome and the taxa 
	min = d; taxa =k;
	if (d < min)
		min=d; taxa=k;
	end
	end

	if taxa=1
  		print('Amphibian');
	else if taxa=2
  		print('Reptile');
	else if taxa=3
  		print('Bird');
	else if taxa=4
  		print('Insect');
	else if taxa=5
  		print('Fish');
	else if taxa=6
  		print('Mammal');
	else if taxa=7
  		print('NonAnimal');
	end
end


********************************************************************************
User can now create his/her own training and test datasets by dividing the TRIM vectors 
in TRIM_DB(:) to evaluate performance of TRIM in genome identification.  