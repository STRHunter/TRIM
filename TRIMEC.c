#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
unsigned char *dna;
unsigned long dna_len;
unsigned long counter;
double length;
long int **dist;   
double **interval;  
double sum, mean, non_zero_interval_counter,std; 
double entropy; 
double total_entropy_density; double total_entropy; 
unsigned long atomic_counter;
long int *repeats;
int cut_off; 
int k;  
unsigned char motif[7];
double **TRIM_Vector;
FILE *fp;
void TRIM(int, char *); 
void strs(int k);
unsigned long len; 
unsigned int lookup[20]={0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};


int main (int argc, char *argv[]) 
{
    
  unsigned long i,j,k=0; char filename[500];
  double max=0.0, temp=0.0;  double myvariable=0.0;
  ///////////////////////////////////////////////////////////////
  // allocating memory for candidate TRIM vectors for the input genome
  // As the classifier will check the H for different k values, we need to store the TRIM vectors for 
  // the given unknown genome on those k values. Usually the classifier will create TRIM vectors
  // for k =6 and 7 as currently we have training set based on genomes whore H maximizes at either k=6 or 7.  
  TRIM_Vector = (double **)malloc(6*sizeof(double *));
    if (TRIM_Vector == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<=1; i++) { 
   TRIM_Vector[i] = (double *)malloc(pow(4,8)*sizeof(double));
   if (TRIM_Vector[i] == NULL) {
      printf("Error in allocating \n");
      //return 0;
    }
  }
  
  ///////////////////////////////////////////////////////////////
  /// calling TRIM vector construction routine for k=6 and 7 for the input genome
  printf("........................... Welcome to genome classifier utility: TRIMEC .......................\n\n");
  printf("Input file name of the genome to identify its taxa: ");
  scanf("%s", filename);
  // Now cleaning the genome sequence, if already cleanded the sequence will be kept as it it. 
  // declaring necessary variables to read data files
   FILE *fp1,  *fp2;
   
   unsigned char *buffer, *dna;
   float actual=0.0;
   // opening the data files

   fp1 = fopen(filename, "rb");// opening anonated or uncleanded file
   if( fp1 == NULL )  {
      perror ("Error opening file");
      return(-1);
   }
   
   fp2 = fopen("clean_genome.txt", "wb");// opening file to write clean genome
   if( fp2 == NULL )  {
      perror ("Error opening file");
      return(-1);
   }
   
   // collecting length of the data in anonated file
   fseek(fp1, 0, SEEK_END);
   len = ftell(fp1);
   printf("Length of the annotated file is %lu\n",len);    
   // allocating a memory space 'buffer' to hold anotated file content
   buffer = (unsigned char *) malloc(sizeof(unsigned char)*len +1);
   
   // reading the data file by single disk access
   fseek(fp1, 0, SEEK_SET);// going begin of the file
   fread(buffer, len , 1, fp1); // reading file content into buffer 
   fclose(fp1);// closing file
         
   // computing the raw dna size from annotated file 
   j=0;
   while (i <= len-2)
  {
     if (buffer[i]=='a' || buffer[i] =='c' || buffer[i]=='g' || buffer[i]=='t' || buffer[i]=='A' || buffer[i] =='C' || buffer[i]=='G' || buffer[i]=='T') 
        j++;
     i++;
  }// end while 

 printf("Raw genome(acgt only) size from original file %s = %lu\n", filename, j);
 // allocating space for raw dna 
 dna= (unsigned char *) malloc(sizeof(unsigned char) * j +1);
 i=0;j=0;
  while (i <= len -2)
  {
 if (buffer[i]=='a' || buffer[i] =='c' || buffer[i]=='g' || buffer[i]=='t' || buffer[i]=='A' || buffer[i] =='C' || buffer[i]=='G' || buffer[i]=='T') //  collect bases for a dna from annotated file
  { dna[j]= buffer[i]; 
    j++;
  }// end if 
  i++;
  
  }// end while
  

  // changing capital letters to small letters
  i=0;
  while (i < j)
  {
  if (dna[i] == 'A') dna[i] ='a';
  else if (dna[i] == 'C') dna[i] ='c';
  else if (dna[i] == 'G') dna[i] = 'g';
  else if (dna[i] == 'T') dna[i] ='t'; 
  i++;
   
  }// end while

  //printf("Null char index =  %lu\n", j);
  dna[j] = '\0';  
  //printf("The raw dna is as \n %s \n The size is %d\n",dna, (int)strlen(dna));
  
  fwrite(dna,1,(unsigned long)strlen(dna),fp2);
  fclose(fp2);
  len = (unsigned long)strlen(dna); 
  printf("Actual genome length %lu\n", len);
  free(buffer);
  free(dna);
  
  printf("\n........................... System is performing the following tasks ................\n\n");
  printf("Constructing TRIM vector for the genome stored in file: %s ... \n",filename); 
  for (i=6;i<=7; i++)
  {
    
   printf("at k=%lu",i); 
   TRIM(i,"clean_genome.txt"); // motif length, file name
     
  }
  // Now In TRIM_Vector[][] we have 2 different vectors for the genome
  // Based on maximum value of total entropy we will call correcponding classifier 
  max=0.0;      
  for (i=0;i<=1;i++)
  {
     for (j=0;j<pow(4,i+6);j++)
       temp += TRIM_Vector[i][j]; // computing H 
     //printf("sum %lf \n",temp); // Display H 
     for (j=0;j<pow(4,i+6);j++)
       //TRIM_Vector[i][j] = TRIM_Vector[i][j] / temp; // normalizing the vector
     if (temp > max){ max = temp ; k= i+6; }
     //temp=0.0;
     //for (j=0;j<pow(4,i+6);j++)
     //  temp += TRIM_Vector[i][j]; 
     //printf("Normalized sum %lf \n",temp);
     temp =0.0; 
     //printf("%lf \n",temp);
  }
  
 // We need to apply feature selection on the TRIM vectors of the unknown input genome
 // TRIM vector of the input genome at k=6 is trored in TRIM[0][0-4095]
 // TRIM vector of the input genome at k=7 is trored in TRIM[1][0-16383]
 if (k==7)
{
  FILE *fp;
  fp = fopen("first_2500_feature.txt", "r"); 
   if(fp == NULL)  {
      perror ("Error opening file");
      //printf("No such file: %s\n",inputfile);
      return(-1);
     }

 int feature[2500], feature_index[16384];
 for (j=0;j<pow(4,k);j++)
  feature_index[j] =0; 

 for (j=0;j<2500;j++)
 { 
    fscanf(fp,"%lf",&myvariable);
    feature[j] = myvariable;
 }

 // now we have features in the array feature[], 
 // so we will cancel all other features from  TRIM vectors excepts those are in this array
 // and then perfrom normalization of TRIM vectors
 
 for (j=0;j<2500;j++)
  feature_index[feature[j]]=1; 

 for (i=0;i<2;i++) 
 {
   for (j=0;j<pow(4,k);j++)
    {
       if (feature_index[j]==0) 
       TRIM_Vector[i][j] = 0; // In each TRIM vector the jth feature 
    }
 } 


 // normalizing TRIM vectors
 
 for (i=0;i<2;i++) 
 {
  temp =0.0; sum=0.0;
  for (j=0;j<pow(4,k);j++) {
  //printf("%lf ", TRIM_Vectors[i][j]); 
  sum += TRIM_Vector[i][j];
  }
  for (j=0;j<pow(4,k);j++) {
  TRIM_Vector[i][j] = TRIM_Vector[i][j] /sum; 
  temp += TRIM_Vector[i][j];
  }
  //printf("The sum of %dth TRIM vector is %lf ",i,temp);
  //printf("\n"); 
}

}

else if (k==6)
{
  FILE *fp;
  fp = fopen("first_4000_feature.txt", "r"); 
   if(fp == NULL)  {
      perror ("Error opening file");
      //printf("No such file: %s\n",inputfile);
      return(-1);
     }

 int feature[4000], feature_index[4096];
 for (j=0;j<pow(4,k);j++)
  feature_index[j] =0; 

 for (j=0;j<4000;j++)
 { 
    fscanf(fp,"%lf",&myvariable);
    feature[j] = myvariable;
 }

 // now we have features in the array feature[], 
 // so we will cancel all other features from  TRIM vectors excepts those are in this array
 // and then perfrom normalization of TRIM vectors
 
 for (j=0;j<4000;j++)
  feature_index[feature[j]]=1; 

 for (i=0;i<2;i++) 
 {
   for (j=0;j<pow(4,k);j++)
    {
       if (feature_index[j]==0) 
       TRIM_Vector[i][j] = 0; // In each TRIM vector the jth feature 
    }
 } 


 // normalizing TRIM vectors
 
 for (i=0;i<2;i++) 
 {
  temp =0.0; sum=0.0;
  for (j=0;j<pow(4,k);j++) {
  //printf("%lf ", TRIM_Vectors[i][j]); 
  sum += TRIM_Vector[i][j];
  }
  for (j=0;j<pow(4,k);j++) {
  TRIM_Vector[i][j] = TRIM_Vector[i][j] /sum; 
  temp += TRIM_Vector[i][j];
  }
  //printf("The sum of %dth TRIM vector is %lf ",i,temp);
  //printf("\n"); 
}

}








  // Loading the corresponding classifier
  printf("\nTotal interval entropy (H) maximizes at k=%lu\n",k);
  printf("Locading corresponding classifier at k=%lu....\n",k);
  char * inputfile, filenumber[10]; 
  inputfile = (char *)malloc(1000*sizeof(char));
  strcpy(inputfile,"Bagging_classifier_at_k=");
  sprintf(filenumber, "%lu", k);
  strcat(inputfile,filenumber);
  strcat(inputfile,".txt");
  fp = fopen(inputfile, "r");
  if(fp == NULL)  {
      perror ("Error opening file:\n");
      printf("No such trainning file: %s\n",inputfile);
      return(-1);
   }
   
  // Now the decision for TRIMEC6 or TRIMEC7 is available, and the corresponding training is loaded into the memory
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Reading the Mean TRIM vectors for all the taxa and for all instances
 //double myvariable=0.0;
 
 double **Mean_TRIM_Amphibian;
  Mean_TRIM_Amphibian = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Amphibian == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Amphibian[i] = (double *)malloc(pow(4,k)*sizeof(double));
     if (Mean_TRIM_Amphibian[i] == NULL) {
        printf("Error in allocating dist6[%lu] \n", i);
        //return 0;
      }
   }

  double **Mean_TRIM_Reptile;
  Mean_TRIM_Reptile = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Reptile == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Reptile[i] = (double *)malloc(pow(4,k)*sizeof(double));
     if (Mean_TRIM_Reptile[i] == NULL) {
        printf("Error in allocating dist6[%lu] \n", i);
        //return 0;
      }
   }
   
  
 double **Mean_TRIM_Bird;
  Mean_TRIM_Bird = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Bird == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Bird[i] = (double *)malloc(pow(4,k)*sizeof(double));
     if (Mean_TRIM_Bird[i] == NULL) {
        printf("Error in allocating dist6[%lu] \n", i);
        //return 0;
      }
   }


  double **Mean_TRIM_Insect;
  Mean_TRIM_Insect = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Insect == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Insect[i] = (double *)malloc(pow(4,k)*sizeof(double));
     if (Mean_TRIM_Insect[i] == NULL) {
        printf("Error in allocating dist6[%lu] \n", i);
        //return 0;
      }
   }


  double **Mean_TRIM_Fish;
  Mean_TRIM_Fish = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Fish == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Fish[i] = (double *)malloc(pow(4,k)*sizeof(double));
     if (Mean_TRIM_Fish[i] == NULL) {
        printf("Error in allocating dist6[%lu] \n", i);
        //return 0;
      }
   }
   
 double **Mean_TRIM_Mammal;
  Mean_TRIM_Mammal = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Mammal == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Mammal[i] = (double *)malloc(pow(4,k)*sizeof(double));
     if (Mean_TRIM_Mammal[i] == NULL) {
        printf("Error in allocating dist6[%lu] \n", i);
        //return 0;
      }
   }
  

 double **Mean_TRIM_Non_animal;
  Mean_TRIM_Non_animal = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Non_animal == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Non_animal[i] = (double *)malloc(pow(4,k)*sizeof(double));
     if (Mean_TRIM_Non_animal[i] == NULL) {
        printf("Error in allocating dist6[%lu] \n", i);
        //return 0;
      }
   }


 for (i =0;i<25;i++)
 {
  //printf("\nInstance %d ",i);
  for(j=0;j<pow(4,k);j++)
   {
    fscanf(fp,"%lf",&myvariable);
    Mean_TRIM_Amphibian[i][j] = myvariable;
   }
   //printf("  %d",j);
   for(j=0;j<pow(4,k);j++)
   {
    fscanf(fp,"%lf",&myvariable);
    Mean_TRIM_Reptile[i][j] = myvariable;
   }
   for(j=0;j<pow(4,k);j++)
   {
    fscanf(fp,"%lf",&myvariable);
    Mean_TRIM_Bird[i][j] = myvariable;
   }
   for(j=0;j<pow(4,k);j++)
   {
    fscanf(fp,"%lf",&myvariable);
    Mean_TRIM_Insect[i][j] = myvariable;
   }
   for(j=0;j<pow(4,k);j++)
   {
    fscanf(fp,"%lf",&myvariable);
    Mean_TRIM_Fish[i][j] = myvariable;
   }
   for(j=0;j<pow(4,k);j++)
   {
    fscanf(fp,"%lf",&myvariable);
    Mean_TRIM_Mammal[i][j] = myvariable;
   }
   for(j=0;j<pow(4,k);j++)
   {
    fscanf(fp,"%lf",&myvariable);
    Mean_TRIM_Non_animal[i][j] = myvariable;
   }
 
}
   
 fclose(fp);
 temp =0.0;

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Now identifyting the unknown genome

 // Voting system 
 int index =0;
 double min=0.0, part1 =0.0, part2=0.0;
 double JS_distance_A=0.0, JS_distance_R=0.0, JS_distance_B=0.0, JS_distance_I=0.0, JS_distance_F=0.0, JS_distance_M=0.0, JS_distance_NA=0.0;
 int A=0, R=0, B=0, I=0, F=0, M=0, NA=0; int taxa=-1, maxvote=0;
 double classifier[25][7];
 // if k=6, that is using TRIMEC-6
 if (k ==6)
 {
   double mid[4096];
   for (j=0;j<25;j++)
   {
      min =0; index =-1;
     printf("Output of classifier %lu\nTRIM distance of the given genme from Mean TRIM Vector of ",j);
    
     // initializing min JS distance and the index value
      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Amphibian
      /*for (i=0;i<4096;i++)
         mid[i] = (Mean_TRIM_Amphibian[j][i] + TRIM_Vector[0][i])/2;
      
      for (i=0;i<4096;i++)
         part1= part1 + TRIM_Vector[0][i] * log2(TRIM_Vector[0][i]/mid[i]); 
    
      for (i=0;i<4096;i++)
         part2= part2 + Mean_TRIM_Amphibian[j][i] * log2(Mean_TRIM_Amphibian[j][i]/mid[i]);
      
      JS_distance_A = (part1 + part2)/2;
      min = JS_distance_A; index =0; 
      
      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Reptile
      part1 =0.0; part2=0.0;
      for (i=0;i<4096;i++)
         mid[i] = (Mean_TRIM_Reptile[j][i] + TRIM_Vector[0][i])/2;
      
      for (i=0;i<4096;i++)
         part1= part1 + TRIM_Vector[0][i] * log2(TRIM_Vector[0][i]/mid[i]);
    
      for (i=0;i<4096;i++)
         part2= part2 + Mean_TRIM_Reptile[j][i] * log2(Mean_TRIM_Reptile[j][i]/mid[i]);
      
      JS_distance_R = (part1 + part2)/2;
      if ( JS_distance_R < min) {min = JS_distance_R; index =1;}
      
      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Bird
      part1 =0.0; part2=0.0;
      for (i=0;i<4096;i++)
         mid[i] = (Mean_TRIM_Bird[j][i] + TRIM_Vector[0][i])/2;
      
      for (i=0;i<4096;i++)
         part1= part1 + TRIM_Vector[0][i] * log2(TRIM_Vector[0][i]/mid[i]);
    
      for (i=0;i<4096;i++)
         part2= part2 + Mean_TRIM_Bird[j][i] * log2(Mean_TRIM_Bird[j][i]/mid[i]);
      
      JS_distance_B = (part1 + part2)/2;
      if ( JS_distance_B < min) {min = JS_distance_B; index =2;}
     */

     // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Insect
      part1 =0.0; part2=0.0;
      for (i=0;i<4096;i++)
         mid[i] = (Mean_TRIM_Insect[j][i] + TRIM_Vector[0][i])/2;
      
      for (i=0;i<4096;i++)
         if (TRIM_Vector[0][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[0][i] * log2(TRIM_Vector[0][i]/mid[i]);
         else part1 =part1;
    
      for (i=0;i<4096;i++)
         if (Mean_TRIM_Insect[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Insect[j][i] * log2(Mean_TRIM_Insect[j][i]/mid[i]);
         else part2=part2;
      
      JS_distance_I = (part1 + part2)/2; printf("\n\n Insect: %lf ",JS_distance_I);
      //if ( JS_distance_I < min) {min = JS_distance_I; index =3;}
      min = JS_distance_I; index =3;
      
     // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Fish
      part1 =0.0; part2=0.0;
      for (i=0;i<4096;i++)
         mid[i] = (Mean_TRIM_Fish[j][i] + TRIM_Vector[0][i])/2;
      
      for (i=0;i<4096;i++)
         if (TRIM_Vector[0][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[0][i] * log2(TRIM_Vector[0][i]/mid[i]);
         else part1=part1;
    
      for (i=0;i<4096;i++)
         if (Mean_TRIM_Fish[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Fish[j][i] * log2(Mean_TRIM_Fish[j][i]/mid[i]);
         else part2 =part2;
      
      JS_distance_F = (part1 + part2)/2; printf("\nFish: %lf ",JS_distance_F);
      if ( JS_distance_F < min) {min = JS_distance_F; index =4;}


      /*
      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Mammal
      part1 =0.0; part2=0.0;
      for (i=0;i<4096;i++)
         mid[i] = (Mean_TRIM_Mammal[j][i] + TRIM_Vector[0][i])/2;
      
      for (i=0;i<4096;i++)
         part1= part1 + TRIM_Vector[0][i] * log2(TRIM_Vector[0][i]/mid[i]);
    
      for (i=0;i<4096;i++)
         part2= part2 + Mean_TRIM_Mammal[j][i] * log2(Mean_TRIM_Mammal[j][i]/mid[i]);
      
      JS_distance_M = (part1 + part2)/2;
      if ( JS_distance_M < min) {min = JS_distance_M; index =5;}
     */

     // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Non_animal
      part1 =0.0; part2=0.0;
      for (i=0;i<4096;i++)
         mid[i] = (Mean_TRIM_Non_animal[j][i] + TRIM_Vector[0][i])/2;
      
      for (i=0;i<4096;i++)
         if (TRIM_Vector[0][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[0][i] * log2(TRIM_Vector[0][i]/mid[i]);
         else part1=part1;
    
      for (i=0;i<4096;i++)
         if (Mean_TRIM_Non_animal[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Non_animal[j][i] * log2(Mean_TRIM_Non_animal[j][i]/mid[i]);
         else part2 =part2;
      
      JS_distance_NA = (part1 + part2)/2; printf("\nNon-animal: %lf ",JS_distance_NA);
      if ( JS_distance_NA < min) {min = JS_distance_NA; index =6;}

       classifier[j][index] =1;
  
    }// end all instances of classifier at k=6  

    // collecting votes
  for (j=0;j<25;j++)
  {
    //A = A + classifier[j][0];
    //R = R + classifier[j][1];
    //B = B + classifier[j][2];
    I = I + classifier[j][3];
    F = F + classifier[j][4];
    //M = M + classifier[j][5];
    NA = NA + classifier[j][6];
  }

 }// end classifier at k=6
 else if (k==7) // that is using TRIMEC-7
  {

   double mid[16384];
   for (j=0;j<25;j++)
   {
     min =0; index =-1;// initializing min JS distance and the index value
      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Amphibian
     printf("Output of classifier %lu\nTRIM distance of the given genme from Mean TRIM Vector of ",j);
     part1 =0.0; part2=0.0; 
     for (i=0;i<16384;i++)
         mid[i] = (Mean_TRIM_Amphibian[j][i] + TRIM_Vector[1][i])/2;
      
     for (i=0;i<16384;i++)
         if (TRIM_Vector[1][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[1][i] * log2(TRIM_Vector[1][i]/mid[i]); 
         else part1 = part1;
          
     for (i=0;i<16384;i++)
         if (Mean_TRIM_Amphibian[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Amphibian[j][i] * log2(Mean_TRIM_Amphibian[j][i]/mid[i]); 
         else part2 =part2;
            
      JS_distance_A = (part1 + part2)/2; printf("\n\nAmphibian: %lf ",JS_distance_A);
      min = JS_distance_A; index =0; 
      
      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Reptile
      part1 =0.0; part2=0.0;
      for (i=0;i<16384;i++)
         mid[i] = (Mean_TRIM_Reptile[j][i] + TRIM_Vector[1][i])/2;
      
      for (i=0;i<16384;i++)
         if (TRIM_Vector[1][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[1][i] * log2(TRIM_Vector[1][i]/mid[i]);
         else part1 = part1;
    
      for (i=0;i<16384;i++)
         if (Mean_TRIM_Reptile[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Reptile[j][i] * log2(Mean_TRIM_Reptile[j][i]/mid[i]);
         else part2= part2;
      
      JS_distance_R = (part1 + part2)/2; printf("\nReptile: %lf ",JS_distance_R);
      if ( JS_distance_R < min) {min = JS_distance_R; index =1;}

      
      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Bird
      part1 =0.0; part2=0.0;
      for (i=0;i<16384;i++)
         mid[i] = (Mean_TRIM_Bird[j][i] + TRIM_Vector[1][i])/2;
      
      for (i=0;i<16384;i++)
         if (TRIM_Vector[1][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[1][i] * log2(TRIM_Vector[1][i]/mid[i]);
         else part1 = part1;
    
      for (i=0;i<16384;i++)
         if (Mean_TRIM_Bird[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Bird[j][i] * log2(Mean_TRIM_Bird[j][i]/mid[i]);
         else part2=part2;
      
      JS_distance_B = (part1 + part2)/2; printf("\nBird: %lf ",JS_distance_B);
      if ( JS_distance_B < min) {min = JS_distance_B; index =2;}



     // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Insect
      part1 =0.0; part2=0.0;
      for (i=0;i<16384;i++)
         mid[i] = (Mean_TRIM_Insect[j][i] + TRIM_Vector[1][i])/2;
      
      for (i=0;i<16384;i++)
        if (TRIM_Vector[1][i] !=0 && mid[i] !=0)       
         part1= part1 + TRIM_Vector[1][i] * log2(TRIM_Vector[1][i]/mid[i]);
        else part1 = part1;
    
      for (i=0;i<16384;i++)
         if (Mean_TRIM_Insect[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Insect[j][i] * log2(Mean_TRIM_Insect[j][i]/mid[i]);
         else part2 = part2;
      
      JS_distance_I = (part1 + part2)/2; printf("\nInsect: %lf ",JS_distance_I);
      if ( JS_distance_I < min) {min = JS_distance_I; index =3;}

      
     // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Fish
      part1 =0.0; part2=0.0;
      for (i=0;i<16384;i++)
         mid[i] = (Mean_TRIM_Fish[j][i] + TRIM_Vector[1][i])/2;
      
      for (i=0;i<16384;i++)
         if (TRIM_Vector[1][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[1][i] * log2(TRIM_Vector[1][i]/mid[i]);
         else part1= part1;
    
      for (i=0;i<16384;i++)
         if (Mean_TRIM_Fish[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Fish[j][i] * log2(Mean_TRIM_Fish[j][i]/mid[i]);
         else part2 =part2;
      
      JS_distance_F = (part1 + part2)/2; printf("\nFish: %lf ",JS_distance_F);
      if ( JS_distance_F < min) {min = JS_distance_F; index =4;}



      // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Mammal
      part1 =0.0; part2=0.0;
      for (i=0;i<16384;i++)
         mid[i] = (Mean_TRIM_Mammal[j][i] + TRIM_Vector[1][i])/2;
      
      for (i=0;i<16384;i++)
         if (TRIM_Vector[1][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[1][i] * log2(TRIM_Vector[1][i]/mid[i]);
         else part1 =part1;
    
      for (i=0;i<16384;i++)
         if (Mean_TRIM_Mammal[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Mammal[j][i] * log2(Mean_TRIM_Mammal[j][i]/mid[i]);
         else part2 = part2;
      
      JS_distance_M = (part1 + part2)/2; printf("\nMammal: %lf ",JS_distance_M);
      if ( JS_distance_M < min) {min = JS_distance_M; index =5;}



    // computing JS distance between the TRIM vector of unknown genome and Mean_TRIM_Non_animal
      part1 =0.0; part2=0.0;
      for (i=0;i<16384;i++)
         mid[i] = (Mean_TRIM_Non_animal[j][i] + TRIM_Vector[1][i])/2;
      
      for (i=0;i<16384;i++)
         if (TRIM_Vector[1][i] !=0 && mid[i] !=0)
         part1= part1 + TRIM_Vector[1][i] * log2(TRIM_Vector[1][i]/mid[i]);
         else part1 =part1;
    
      for (i=0;i<16384;i++)
         if (Mean_TRIM_Non_animal[j][i] !=0 && mid[i] !=0)
         part2= part2 + Mean_TRIM_Non_animal[j][i] * log2(Mean_TRIM_Non_animal[j][i]/mid[i]);
         else part2 = part2;
      
      JS_distance_NA = (part1 + part2)/2; printf("\nNon-animal: %lf \n\n",JS_distance_NA);
      if ( JS_distance_NA < min) {min = JS_distance_NA; index =6;}



     classifier[j][index] =1;
 
    } // end all instances of classifier at k=6
   
  // collecting votes
  for (j=0;j<25;j++)
  {
    A = A + classifier[j][0];
    R = R + classifier[j][1];
    B = B + classifier[j][2];
    I = I + classifier[j][3];
    F = F + classifier[j][4];
    M = M + classifier[j][5];
    NA = NA + classifier[j][6];
  }

  } // end classifier at k=7

  if (A > maxvote) {maxvote =A; taxa=0;}
  else if (R > maxvote) { maxvote =R; taxa =1;}
  else if (B > maxvote) { maxvote =B; taxa =2;}
  else if (I > maxvote) { maxvote =I; taxa =3;}
  else if (F > maxvote) { maxvote =F; taxa =4;}
  else if (M > maxvote) { maxvote =M; taxa =5;}
  else if (NA > maxvote) { maxvote =NA; taxa =6;}
  
  // Declaring the taxa 
  printf("\n\nBased on the maximum agreement among the 25 classifier under TRIMEC\n");
  if (taxa == 0)
   printf("This is a Amphibian\n");
  else if (taxa == 1)
   printf("This is a Reptile\n");
  else if (taxa == 2)
   printf("This is a Bird\n");
  else if (taxa == 3)
   printf("This is a Insect\n");
  else if (taxa == 4)
   printf("This is a Fish\n");
  else if (taxa == 5)
   printf("This is a Mammal\n");
  else if (taxa == 6)
   printf("This is a Non_animal\n");





 
  return(0);   
}



void TRIM(int k, char *fileName)
{
   
   // declaring local variables
   FILE * fp; 
   unsigned long long i=0, j=0; 
   unsigned long overall_STR_count=0, overall_Tandem_count=0, overall_Inter_STR_base =0, overall_STR_base=0; 
   int c;
   opterr=0; counter=0; atomic_counter=0;
   // reading the file containing genome sequence
   dna = (unsigned char *)malloc(len);
   fp = fopen(fileName,"rb");
   fread(dna, len , 1, fp);   
   dna_len = strlen(dna); 
   fclose(fp);
  // declaring variables for TRIM vector computation
   dist = (long int **)malloc(pow(4,k)*sizeof(long int *));
    if (dist == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<pow(4,k); i++) { 
   dist[i] = (long int *)malloc(9*sizeof(long int));
   if (dist[i] == NULL) {
      printf("Error in allocating dist6[%llu] \n", i);
      //return 0;
   }
 }
  
 for (i=0;i<pow(4,k);i++)
   {dist[i][0]=-1; dist[i][1]=-1; dist[i][2]=-1; dist[i][3]=0; dist[i][4]=-1; dist[i][5]=0; dist[i][6]=0; dist[i][7]=0; dist[i][8]=0;}
     
  

 interval = (double **)malloc(pow(4,k)*sizeof(double *));
    if (interval == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<pow(4,k); i++) { 
   interval[i] = (double *)malloc(5000*sizeof(double));
   if (interval[i] == NULL) {
      printf("Error in allocating dist[%llu] \n", i);
    //  return 0;
   }
 }



   repeats = (long int *)malloc(6*sizeof(long int));
 
   for (i=0;i<6;i++)
   repeats[i]=-1;
   

   /////////////////////////////////////////////////////////////////////////////////////////////

   strs(k); // mining TRs 
   length = (double)dna_len/(1024*1024);  
   ////////////////////////////////////////////////////////////////////////////////////////////
  
   // Knowledge extraction
   j=0; atomic_counter=0; 
   sum =0.0; total_entropy_density=0.0; total_entropy =0.0;

   for (i =0;i <pow(4,k); i++)
     {
        sum =0.0; entropy=0.0;//mean = 0.0; non_zero_interval_counter =0.0; std=0.0;
        for (j =0;j<5000;j++) // computing sum 
         {
            sum = sum + interval[i][j];
         }
        

        for (j =0;j<5000;j++) // computing pmf
         {
             if (interval[i][j] != 0 && sum !=0)  
               {
                  interval[i][j] = interval[i][j] / sum;
               }
             else 
               {
                  interval[i][j] =0.0;
               }
          }
       
     }

   // Computing entropy/seqence length
    for (i =0;i <pow(4,k); i++) 
    {    
      entropy=0.0;
      for (j =0;j<5000;j++) // computing entropy
          {
             if (interval[i][j] != 0) 
                {
                   entropy = entropy + interval[i][j]* log2(1/interval[i][j]);
                }
             else 
               {
                    entropy = entropy;
                }
          }

       if (entropy > 0.0)
        {
          //printf("%.14lf\n", entropy/length); 
          TRIM_Vector[k-6][i] = entropy;      
        }
      else 
        {
          //printf("0\n");
           TRIM_Vector[k-6][i] = 0; 
        }
                
      total_entropy_density += entropy;
     }
    printf(", H=%lf\n",total_entropy_density); 
      


  
}





unsigned int GetIdx(unsigned long i, int k)
{
  unsigned long j;
  static unsigned long idx = 0;
  if(i==0)
  {
    idx = 0;
    for(j=0; j<k; j++)
    {
      idx = (idx*4+lookup[dna[j]-'a']);
    }
    return idx;  
  }
   idx = (idx*4+lookup[dna[i+k-1]-'a']) % (int)pow(4,k);
   return idx;   
}


void strs(int k)
{
  int j, idx6, d=0; 
  //unsigned long length=0;
  long seconds, ns; 
  unsigned long i=0;
  struct timespec start, finish;
  counter=0;
  //printf("\nInside strs() : \n dna len = %d\n", dna_len);
  //printf("Location:   no. of repeats:    motif:\n" );
  //printf("%lu  %lf\n", dna_len,dna_len*0.4);
for(i=0; i< dna_len-k; i++) // scanning the sequence one character at a time
{
    // genertaing 4-base codes(index) for 6-mers starts at locaion i

    idx6 = GetIdx(i,k);
    
    if (dist[idx6][0] == -1) // 6-mer motif with index = idx6, appearing for first time 
      {
        dist[idx6][0] = i; // store current location the 6-mer motif with index = idx6
      } 
     else // not first occurrence, so check for repeats 
      { 
         d = i - dist[idx6][0];
         if (d < k) // overlapping occurrence of 6-mer
           {
             dist[idx6][0] = i; 
             dist[idx6][1] = -1; // advance the starting point of SSR will deal with atomicity or ; do nothing will report non-atomic SSRs 
           }
         else if (d > k) // not a tandem repeat
           {
              dist[idx6][0] = i; 
              dist[idx6][1] = -1; // terminating the current run of 6-mer SSR with index = idx6
           }  
          else // d == 6, tandem repeat
           {
               if (dist[idx6][0] < repeats[0] + repeats[5] - k && idx6 != repeats[2])// another 6-mer SSR is repeating, or its first time
                  { dist[idx6][0] = i; dist[idx6][1] =-1; } // so do not initiate new 6-mer SSR, just advanced its starting
               else
                  { // initiate new 6-mer SSR 
                    if ( dist[idx6][1]==-1) // initiation of 6-mer SSR with index = idx6
                    {
                        if (repeats[4] ==1 && repeats[5] >=cut_off) 
                        // as repeats6 are mutually exclusive, we can print the stored repeats6 
                        {
                           //if (repeats[2] == 32) {
                           //printf("%ld\t%ld\t%ld\t", repeats[0], repeats[1], repeats[3]);
                           //printf("%s\n",motif);}
                           counter++; // total number of STRs 
                           length = length + repeats[5]; // total STR length
                           dist[repeats[2]][2] +=1; // count the number of STRs of the motif with index = repeats6[2]
                           dist[repeats[2]][3] += repeats[1]; // count the number of tandem of the motif with index = repeats6[2]
                           
                           if (dist[repeats[2]][4] == -1) // STR of the motif with index = repeats[2] is appearing first time
                              { 
                                dist[repeats[2]][4] = repeats[0] + repeats[5] -1;  dist[repeats[2]][5] =0;
                                dist[repeats[2]][6] += repeats[5]; // count the number of bases of STR of motif with index = repeats6[2] 
                                dist[repeats[2]][7] = repeats[0]; //printf("\n Start location of motif(%ld)h id %ld\n",repeats[2],repeats[0]);
                              }
                           else // calculate inter-STR nucleotide bases
                              { 
                                dist[repeats[2]][5] += (repeats[0] - dist[repeats[2]][4]); 
                                if (dist[repeats[2]][2] < 5000)
                                interval[repeats[2]][dist[repeats[2]][2]] = abs(repeats[0] - dist[repeats[2]][4]);
                                dist[repeats[2]][4] = repeats[0] + repeats[5] -1; 
                                dist[repeats[2]][6] += repeats[5]; // count the number of bases of STR of motif with index = repeats6[2]
                                dist[repeats[2]][8] = repeats[0]; 
                              }
                           
                              
                         } 
                         repeats[0] = dist[idx6][0];// storing the start location of 6-mer repeats 
                         repeats[1] = 2; // number of repeats 
                         repeats[2] = idx6; // storing the 6-mer
                         repeats[3] = k; // storing the size of 6-mer, that is 6  
                         repeats[4] = 1; // setting valid bit
                         repeats[5] = repeats[1] * repeats[3];// computing cut_off 
                         motif[0] = dna[i]; motif[1] = dna[i+1];motif[2] = dna[i+2];motif[3] = dna[i+3];
                         motif[4] = dna[i+4]; motif[5] = dna[i+5]; motif[6]='\0'; 
                         dist[idx6][0] = i; 
                         dist[idx6][1]= 0; 
                     }
                    else 
                     {
                        repeats[1]++; repeats[5] = repeats[1] * repeats[3]; dist[idx6][0] = i;
                     } 
                  }// end of: initiate and increment 6-mer SSR     
         }// end of: d == 6, consecutive repeating occurrence

       }// end of: not first occurrence, so check for repeats 
   

}// end scanning 



if (repeats[4] ==1 && repeats[5] >=cut_off) // Printing the last stored repeats6 
  {
     //if (repeats[2] == 32) {
     //printf("%10ld%10ld%10ld", repeats[0], repeats[1], repeats[3]);
     //printf("%10s\n",motif);}
     counter++; // total number of STRs 
     length = length + repeats[5]; // total STR length
     dist[repeats[2]][2] +=1; // count the number of STRs of the motif with index = repeats6[2]
     dist[repeats[2]][3] += repeats[1]; // count the number of tandem of the motif with index = repeats6[2]
                          
     if (dist[repeats[2]][4] == -1) // STR of the motif with index = repeat6[2] is appearing first time
      { 
         dist[repeats[2]][4] = repeats[0] + repeats[5] -1;  dist[repeats[2]][5] =0;
         dist[repeats[2]][6] += repeats[5]; // count the number of bases of STR of motif with index = repeats6[2] 
         dist[repeats[2]][8] = repeats[0];//printf("\n End location of motif(%ld)h id %ld\n",repeats[2],repeats[0]);
      }
     else // calculate inter-STR nucleotide bases
      { 
         dist[repeats[2]][5] += (repeats[0] - dist[repeats[2]][4]); 
         if (dist[repeats[2]][2] < 5000)
         interval[repeats[2]][dist[repeats[2]][2]] = abs(repeats[0] - dist[repeats[2]][4]); 
         dist[repeats[2]][4] = repeats[0] + repeats[5] -1; 
         dist[repeats[2]][6] += repeats[5]; // count the number of bases of STR of motif with index = repeats6[2] 
         dist[repeats[2]][8] = repeats[0]+repeats[5];//printf("\n End location of motif(%ld) is %ld\n" ,repeats[2],repeats[0]);
     }
   } 

  //printf("Total STR length is %lf and Total number of STR is %lu \n", length,counter);

}
