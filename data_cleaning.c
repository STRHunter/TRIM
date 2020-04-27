#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#define DMAX 

int main (int *argc, char *argv[]) 
{
   // declaring necessary variables to read data files
   FILE *fp1,  *fp2;
   unsigned long len, i=0, j=0; 
   int rand_arr[MAX];// MAX = number of places to translate
   unsigned char *buffer, *dna;
   time_t t;  
   float actual=0.0;
   // opening the data files

   fp1 = fopen(argv[1], "rb");// opening anonated file
   if( fp1 == NULL )  {
      perror ("Error opening file");
      return(-1);
   } 
   
   fp2 = fopen(argv[2], "wb");// opening file to write raw dna
   if( fp2 == NULL )  {
      perror ("Error opening file");
      return(-1);
   }
   
   // collecting length of the data in anonated file
   fseek(fp1, 0, SEEK_END);
   len = ftell(fp1);
   //printf("Length of the annotated file is %lu\n",len);    
   // allocating a memory space 'buffer' to hold anotated file content
   buffer = (unsigned char *) malloc(sizeof(unsigned char)*len +1);
   
   // reading the data file by single disk access
   fseek(fp1, 0, SEEK_SET);// going begin of the file
   fread(buffer, len , 1, fp1); // reading file content into buffer 
   fclose(fp1);// closing file
   //printf("Origical file contect is \n%s \n", buffer);
      
   // computing the raw dna size from annotated file 
   j=0;
   while (i <= len-2)
  {
     if (buffer[i]=='a' || buffer[i] =='c' || buffer[i]=='g' || buffer[i]=='t' || buffer[i]=='A' || buffer[i] =='C' || buffer[i]=='G' || buffer[i]=='T') 
        j++;
     i++;
  }// end while 

 printf("Raw DNA size from original file = %lu\n", j);
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
  
  // creating random locations
  /*
  srand((unsigned)time(&t));
  for (i =0; i < MAX;i++)
  { 
   rand_arr[i] = rand()%(int)strlen(dna);
  }
  // printing random locations
  //for (i =0;i < MAX;i++)
  //{ 
   //printf("%d\n",rand_arr[i]);
  //}
  // modifying the raw DNA
  for (i =0; i < MAX; i ++)
  {
   if (dna[rand_arr[i]] == 'a')
    dna[rand_arr[i]] = 'c';
   else if (dna[rand_arr[i]] == 'c')
    dna[rand_arr[i]] = 'g';
   else if (dna[rand_arr[i]] == 'g')
    dna[rand_arr[i]] = 't';
   else if (dna[rand_arr[i]] == 't')
    dna[rand_arr[i]] = 'a';
  }
  */  
  fwrite(dna,1,(unsigned long)strlen(dna),fp2);
  fclose(fp2);
  len = (unsigned long)strlen(dna); 
  //actual = ((double)(len - MAX)/ (double)len) * 100 ;  
  printf("Actual chromosome length %lu\n", len);
  //printf("%f\n",actual);
  free(buffer);
  free(dna);
  return(0);
}




