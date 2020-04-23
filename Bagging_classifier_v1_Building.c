#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>
#include<unistd.h>
#include<getopt.h>

int printRandom(int lower, int upper)
{
  int num = (rand()%(upper -lower +1)) + lower;
  return num;
}

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
 
 for (i=0;i<No_of_TRIM_Vectors;i++) 
 {
  temp =0.0; sum=0.0;
  for (j=0;j<pow(4,motif_length);j++) {
  printf("%lf ", TRIM_Vectors[i][j]); 
  sum += TRIM_Vectors[i][j];
  }
  for (j=0;j<pow(4,motif_length);j++) {
  TRIM_Vectors[i][j] = TRIM_Vectors[i][j] /sum; 
  temp += TRIM_Vectors[i][j];
  }
  printf("The sum of %dth TRIM vector is %lf ",i,temp);
  printf("\n"); 
}
 // Now all TRIM Vectors are loaded into main memory.

  
 int lower = 0, upper =5, count =5;
 int * Amphibian_traning ;  int * Reptile_traning;  int * Bird_traning;
 int * Insect_traning;    int * Fish_traning;    int * Mammal_traning;
 int * Non_animal_traning;
   
  printf("Collecting prior knowledge about \"class-label\" of TRIM Vectors\n");
  
  ////////////////////////////////////////////// Amphibian //////////////////////////
  
  double **Mean_TRIM_Amphibian;
   Mean_TRIM_Amphibian = (double **)malloc(25*sizeof(double *));
     if (Mean_TRIM_Amphibian == NULL) {
      printf("Error in allocating dist \n");
     //return 0;
    }

    for(i=0; i<25; i++) { 
      Mean_TRIM_Amphibian[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
      if (Mean_TRIM_Amphibian[i] == NULL) {
         printf("Error in allocating dist6[%d] \n", i);
         //return 0;
       } 
    }

  printf("Provide the range of indices of TRIM Vectors for class Amphibian\n");
  printf("Enter lower and upper indices along with number of TRIM vectors for training: Put zeros to skip training"); 
  scanf("%d %d %d", &lower, &upper, &count);
  Amphibian_traning = (int * ) malloc(count*sizeof(int));
  if (lower >0 && upper >0 && count >0 && upper>lower)
  {
   
   // Randomizing for class Amphibian
   srand(time(0)); //fp = fopen("Mean_TRIM_Amphibian.txt", "wb");
   for (j =0;j< 25;j++) {
   printf("\nIndices of TRIM Vectors used in randomization instance %d for training class Amphibian\n",j);
   for (i =0; i< count; i++) {
   Amphibian_traning[i] = printRandom(lower-1, upper-1);
   printf("%d ", Amphibian_traning[i]);
   }
   // generating Mean TRIM vectors for class Amphibian
    printf("\n"); temp =0.0; 
    for (k =0;k < pow(4,motif_length); k++) {
      temp=0.0;
      for (i =0;i< count ; i++)
        temp += TRIM_Vectors[Amphibian_traning[i]][k];  
      Mean_TRIM_Amphibian[j][k] = temp/count;
    } 
    printf("Mean TRIM Vector of Amphibian for randomized instance %d is \n",j) ; temp =0.0;
   
    for (k =0;k < pow(4,motif_length); k++) {
     printf("%.14lf ", Mean_TRIM_Amphibian[j][k]); temp += Mean_TRIM_Amphibian[j][k];} 
     printf("temp %lf", temp);
 
   }// End generation of Mean TRIM Vectors for Amphibian
  }
  else { printf("Amphibian is class is not trained !\n"); }







  ////////////////////////////////////////////// Reptile //////////////////////////
  
  double **Mean_TRIM_Reptile;
  Mean_TRIM_Reptile = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Reptile == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Reptile[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
     if (Mean_TRIM_Reptile[i] == NULL) {
        printf("Error in allocating dist6[%d] \n", i);
        //return 0;
      }
   }
  printf("\n\nProvide the range of indices of TRIM Vectors for class Reptile\n");
  printf("Enter lower and upper indices along with number of TRIM vectors for training: Put zeros to skip training"); 
  scanf("%d %d %d", &lower, &upper, &count);
  Reptile_traning = (int * ) malloc(count*sizeof(int));
  if (lower >0 && upper >0 && count >0 && upper>lower)
  {
  

  // Randomizing for class Reptile
  srand(time(0));
  for (j =0;j< 25;j++) {
  printf("\nIndices of TRIM Vectors used in randomization instance %d for training class Reptile\n",j);
  for (i =0; i< count; i++) {
  Reptile_traning[i] = printRandom(lower-1, upper-1);
  printf("%d ", Reptile_traning[i]);
  }
  // generating Mean TRIM vectors for class Reptile
   printf("\n"); temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
     temp=0.0;
     for (i =0;i< count ; i++)
       temp += TRIM_Vectors[Reptile_traning[i]][k];  
     Mean_TRIM_Reptile[j][k] = temp/count;
   } 
   printf("Mean TRIM Vector of Reptile for randomized instance %d is \n",j) ; temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
    printf("%.14lf ", Mean_TRIM_Reptile[j][k]); temp += Mean_TRIM_Reptile[j][k]; } 
    printf("temp %lf", temp);
 
  }// End generation of Mean TRIM Vectors for Reptile
  }
  else { printf("Reptile class is not trained !");}





  ////////////////////////////////////////////// Bird //////////////////////////
  
  double **Mean_TRIM_Bird;
  Mean_TRIM_Bird = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Bird == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Bird[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
     if (Mean_TRIM_Bird[i] == NULL) {
        printf("Error in allocating dist6[%d] \n", i);
        //return 0;
      }
   }
  printf("\n\nProvide the range of indices of TRIM Vectors for class Bird\n");
  printf("Enter lower and upper indices along with number of TRIM vectors for training: Put zero to skip training"); 
  scanf("%d %d %d", &lower, &upper, &count);
  Bird_traning = (int * ) malloc(count*sizeof(int));
  if (lower >0 && upper >0 && count >0 && upper>lower)
  {
  

  // Randomizing for class Bird
  srand(time(0));
  for (j =0;j< 25;j++) {
  printf("\nIndices of TRIM Vectors used in randomization instance %d for training class Bird\n",j);
  for (i =0; i< count; i++) {
  Bird_traning[i] = printRandom(lower-1, upper-1);
  printf("%d ", Bird_traning[i]);
  }
  // generating Mean TRIM vectors for class Bird
   printf("\n"); temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
     temp=0.0;
     for (i =0;i< count ; i++)
       temp += TRIM_Vectors[Bird_traning[i]][k];  
     Mean_TRIM_Bird[j][k] = temp/count;
   } 
   printf("Mean TRIM Vector of Bird for randomized instance %d is \n",j) ; temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
    printf("%.14lf ", Mean_TRIM_Bird[j][k]); temp += Mean_TRIM_Bird[j][k]; } 
    printf("temp %lf", temp);
 
  }// End generation of Mean TRIM Vectors for Bird
  }
  else { printf("Bird class is not trained");}


  ////////////////////////////////////////////// Insect //////////////////////////

  double **Mean_TRIM_Insect;
  Mean_TRIM_Insect = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Insect == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Insect[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
     if (Mean_TRIM_Insect[i] == NULL) {
        printf("Error in allocating dist6[%d] \n", i);
        //return 0;
      }
   }

  printf("\n\nProvide the range of indices of TRIM Vectors for class Insect\n");
  printf("Enter lower and upper indices along with number of TRIM vectors for training: Put zeros to skip training"); 
  scanf("%d %d %d", &lower, &upper, &count);
  Insect_traning = (int * ) malloc(count*sizeof(int));
  
  if (lower >0 && upper >0 && count >0 && upper>lower)
  {
  
  // Randomizing for class Insect
  srand(time(0));
  for (j =0;j< 25;j++) {
  printf("\nIndices of TRIM Vectors used in randomization instance %d for training class Insect\n",j);
  for (i =0; i< count; i++) {
  Insect_traning[i] = printRandom(lower-1, upper-1);
  printf("%d ", Insect_traning[i]);
  }
  // generating Mean TRIM vectors for class Insect
   printf("\n"); temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
     temp=0.0;
     for (i =0;i< count ; i++)
       temp += TRIM_Vectors[Insect_traning[i]][k];  
     Mean_TRIM_Insect[j][k] = temp/count;
   } 
   printf("Mean TRIM Vector of Insect for randomized instance %d is \n",j) ; temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
    printf("%.14lf ", Mean_TRIM_Insect[j][k]); temp += Mean_TRIM_Insect[j][k]; } 
    printf("temp %lf", temp);
 
  }// End generation of Mean TRIM Vectors for Insect
  }
  else { printf("Insect class is not trained");}

 

  ////////////////////////////////////////////// Fish //////////////////////////

  
  double **Mean_TRIM_Fish;
  Mean_TRIM_Fish = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Fish == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Fish[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
     if (Mean_TRIM_Fish[i] == NULL) {
        printf("Error in allocating dist6[%d] \n", i);
        //return 0;
      }
   }

  printf("\n\nProvide the range of indices of TRIM Vectors for class Fish\n");
  printf("Enter lower and upper indices along with number of TRIM vectors for training: Put zeros to skip training"); 
  scanf("%d %d %d", &lower, &upper, &count);
  Fish_traning = (int * ) malloc(count*sizeof(int));
  if (lower >0 && upper >0 && count >0 && upper>lower)
  {
  

  // Randomizing for class Fish
  srand(time(0));
  for (j =0;j< 25;j++) {
  printf("\nIndices of TRIM Vectors used in randomization instance %d for training class Fish\n",j);
  for (i =0; i< count; i++) {
  Fish_traning[i] = printRandom(lower-1, upper-1);
  printf("%d ", Fish_traning[i]);
  }
  // generating Mean TRIM vectors for class Fish
   printf("\n"); temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
     temp=0.0;
     for (i =0;i< count ; i++)
       temp += TRIM_Vectors[Fish_traning[i]][k];  
     Mean_TRIM_Fish[j][k] = temp/count;
   } 
   printf("Mean TRIM Vector of Fish for randomized instance %d is \n",j) ; temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
    printf("%.14lf ", Mean_TRIM_Fish[j][k]); temp += Mean_TRIM_Fish[j][k]; } 
    printf("temp %lf", temp);
 
  }// End generation of Mean TRIM Vectors for Fish
  }
  else { printf("Fish class is not trained!\n");}



  ////////////////////////////////////////////// Mammal //////////////////////////

  double **Mean_TRIM_Mammal;
  Mean_TRIM_Mammal = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Mammal == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Mammal[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
     if (Mean_TRIM_Mammal[i] == NULL) {
        printf("Error in allocating dist6[%d] \n", i);
        //return 0;
      }
   }


  printf("\n\nProvide the range of indices of TRIM Vectors for class Mammal\n");
  printf("Enter lower and upper indices along with number of TRIM vectors for training: Put zeros to skip training"); 
  scanf("%d %d %d", &lower, &upper, &count);
  Mammal_traning = (int * ) malloc(count*sizeof(int));
  
  if (lower >0 && upper >0 && count >0 && upper>lower)
  {
    // Randomizing for class Mammal
  srand(time(0));
  for (j =0;j< 25;j++) {
  printf("\nIndices of TRIM Vectors used in randomization instance %d for training class Mammal\n",j);
  for (i =0; i< count; i++) {
  Mammal_traning[i] = printRandom(lower-1, upper-1);
  printf("%d ", Mammal_traning[i]);
  }
  // generating Mean TRIM vectors for class Mammal
   printf("\n"); temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
     temp=0.0;
     for (i =0;i< count ; i++)
       temp += TRIM_Vectors[Mammal_traning[i]][k];  
     Mean_TRIM_Mammal[j][k] = temp/count;
   } 
   printf("Mean TRIM Vector of Mammal for randomized instance %d is \n",j) ; temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) 
   {
    printf("%.14lf ", Mean_TRIM_Mammal[j][k]); 
    temp += Mean_TRIM_Mammal[j][k]; 
    } 
    printf("temp %lf", temp);
 
  }// End generation of Mean TRIM Vectors for Mammal
  }
  else { printf("Mammal class is not tranined !");}
  





////////////////////////////////////////////// Non_animal //////////////////////////

  double **Mean_TRIM_Non_animal;
  Mean_TRIM_Non_animal = (double **)malloc(25*sizeof(double *));
    if (Mean_TRIM_Non_animal == NULL) {
     printf("Error in allocating dist \n");
    //return 0;
   }

   for(i=0; i<25; i++) { 
     Mean_TRIM_Non_animal[i] = (double *)malloc(pow(4,motif_length)*sizeof(double));
     if (Mean_TRIM_Non_animal[i] == NULL) {
        printf("Error in allocating dist6[%d] \n", i);
        //return 0;
      }
   }


  printf("\n\nProvide the range of indices of TRIM Vectors for class Non_animal\n");
  printf("Enter lower and upper indices along with number of TRIM vectors for training: Put zeros to skip training"); 
  scanf("%d %d %d", &lower, &upper, &count);
  Non_animal_traning = (int * ) malloc(count*sizeof(int));
  
  if (lower >0 && upper >0 && count >0 && upper>lower)
  {
    // Randomizing for class Non_animal
  srand(time(0));
  for (j =0;j< 25;j++) {
  printf("\nIndices of TRIM Vectors used in randomization instance %d for training class Non_animal\n",j);
  for (i =0; i< count; i++) {
  Non_animal_traning[i] = printRandom(lower-1, upper-1);
  printf("%d ", Non_animal_traning[i]);
  }
  // generating Mean TRIM vectors for class Non_animal
   printf("\n"); temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) {
     temp=0.0;
     for (i =0;i< count ; i++)
       temp += TRIM_Vectors[Non_animal_traning[i]][k];  
     Mean_TRIM_Non_animal[j][k] = temp/count;
   } 
   printf("Mean TRIM Vector of Non_animal for randomized instance %d is \n",j) ; temp =0.0;
   for (k =0;k < pow(4,motif_length); k++) 
   {
    printf("%.14lf ", Mean_TRIM_Non_animal[j][k]); 
    temp += Mean_TRIM_Non_animal[j][k]; 
    } 
    printf("temp %lf", temp);
 
  }// End generation of Mean TRIM Vectors for Non_animal
  }
  else { printf("Non_animal class is not tranined !");}
  



  // writing the classifier in file 
  strcpy(inputfile,"Bagging_classifier_at_k=");
  sprintf(filenumber, "%d", motif_length);
  strcat(inputfile,filenumber);
  strcat(inputfile,".txt");
  fp = fopen(inputfile, "wb");
  for (j =0;j< 25 ; j++)
   {
    for (k=0; k < pow(4,motif_length); k++)
      fprintf(fp, "%.14lf ",Mean_TRIM_Amphibian[j][k]);
    fprintf(fp,"\n");
    
    for (k=0; k < pow(4,motif_length); k++)
      fprintf(fp, "%.14lf ",Mean_TRIM_Reptile[j][k]);
    fprintf(fp,"\n");  
    
    for (k=0; k < pow(4,motif_length); k++)
       fprintf(fp, "%.14lf ",Mean_TRIM_Bird[j][k]);
    fprintf(fp,"\n");
    
    for (k=0; k < pow(4,motif_length); k++)
      fprintf(fp, "%.14lf ",Mean_TRIM_Insect[j][k]);  
    fprintf(fp,"\n"); 
   
   for (k=0; k < pow(4,motif_length); k++)
      fprintf(fp, "%.14lf ",Mean_TRIM_Fish[j][k]);
   fprintf(fp,"\n");   
   
   for (k=0; k < pow(4,motif_length); k++)
      fprintf(fp, "%.14lf ",Mean_TRIM_Mammal[j][k]);  
   fprintf(fp,"\n");


   for (k=0; k < pow(4,motif_length); k++)
      fprintf(fp, "%.14lf ",Mean_TRIM_Non_animal[j][k]);  
   fprintf(fp,"\n");


   }

  fclose(fp); 
  printf("\nClassifier with motif length %d has been written to file %s\n",motif_length, inputfile);







}
