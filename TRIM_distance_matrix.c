#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>
#include<unistd.h>
#include<getopt.h>


int main()
{
 
 int No_of_TRIM_Vectors; int motif_length;
 printf("Enter the number of TRIM vectors to load in-memory "); 
 scanf("%d", &No_of_TRIM_Vectors); 
 printf("\nEnter motif length used for constructing TRIM Vectors ");
 scanf("%d", &motif_length);
 int i=0,j=0,k=0; // local variables for indexing purpose
 double temp=0.0, sum=0.0;
 // Allocating memory to hold TRIM Vectors
 double **TRIM_Vectors; 
 TRIM_Vectors = (double **)malloc(No_of_TRIM_Vectors*sizeof(double *));
    if (TRIM_Vectors == NULL) {
     printf("Error in allocating TRIM_vectors \n");
    //return 0;
   }

   for(i=0; i < No_of_TRIM_Vectors; i++) 
   { 
     TRIM_Vectors[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
     if (TRIM_Vectors[i] == NULL) {
      printf("Error in allocating TRIM_vector[%d] \n", i);
      //return 0;
     }
   }
  
  // Allocating variables for file handing

 FILE *fp; double myvariable;   
 char * inputfile, filenumber[10]; 
 inputfile = (char *)malloc(1000*sizeof(char));
 //int fileNo=1; 
 
 
 for (i =0;i<No_of_TRIM_Vectors; i++)
 { 
   
   strcpy(inputfile,"T");
   strcat(inputfile,"RIM_vector_Genome");
   sprintf(filenumber, "%d", i+1);
   strcat(inputfile,filenumber);
   strcat(inputfile,".txt");      
   printf("Reading the TRIM vector from file: %s\n",inputfile);
    
   fp = fopen(inputfile, "r"); 
   if(fp == NULL)  {
      perror ("Error opening file");
      //printf("No such file: %s\n",inputfile);
      return(-1);
     }
   for (j=0;j<pow(4,motif_length);j++)
   { 
    fscanf(fp,"%lf",&myvariable);
    TRIM_Vectors[i][j] = myvariable;
    //printf("%.14lf \n", TRIM_Vectors[i][j]);
   } 

   fclose(fp);
 }
 
 

  
 // Normalizng the TRIM Vectors
 for (i=0;i<No_of_TRIM_Vectors;i++) 
 {
  temp =0.0; sum=0.0;
  for (j=0;j<pow(4,motif_length);j++) {
  //printf("%lf ", TRIM_Vectors[i][j]); 
  sum += TRIM_Vectors[i][j];
  }
  for (j=0;j<pow(4,motif_length);j++) {
  TRIM_Vectors[i][j] = TRIM_Vectors[i][j] /sum; 
  temp += TRIM_Vectors[i][j];
  }
  //printf("The sum of %dth TRIM vector is %lf ",i,temp);
  //printf("\n"); 
}
 // Now all TRIM Vectors are loaded into main memory.
 

 // Computing and writing the TRIM distance matrix into file
  printf("Computing TRIM distance matrix ...\n"); 
  double part1=0.0, part2=0.0;
  double *mid = (double *)malloc(pow(4,motif_length)*sizeof(double));  
  
  fp = fopen("TRIM_distance_matrix.txt", "wb");
  if(fp == NULL)  {
      perror ("Error opening file");
      //printf("No such file: %s\n",inputfile);
      return(-1);
     }
  
  
  for (i=0; i< No_of_TRIM_Vectors-1; i++)
  {
    for (j =i+1; j< No_of_TRIM_Vectors; j++)
     {
        part1 =0.0; part2=0.0; 
        for (k=0;k<pow(4,motif_length);k++)
         mid[k] = (TRIM_Vectors[i][k] + TRIM_Vectors[j][k])/2;
      
        for (k=0;k<pow(4,motif_length);k++)
         if (TRIM_Vectors[i][k] !=0 && mid[k] !=0)
          part1= part1 + TRIM_Vectors[i][k] * log2(TRIM_Vectors[i][k]/mid[k]); 
         else part1 = part1;
          
        for (k=0;k<pow(4,motif_length);k++)
         if (TRIM_Vectors[j][k] !=0 && mid[k] !=0)
          part2= part2 + TRIM_Vectors[j][k] * log2(TRIM_Vectors[j][k]/mid[k]); 
         else part2 =part2;
         fprintf(fp,"%.14f ",(part1 + part2)/2);    
              
     }
   }

   
 
 fclose(fp);  
 printf("TRIM distance matrix writen into file TRIM_distance_matrix.txt\n"); 
 
}
