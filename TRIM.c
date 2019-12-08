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

void strs();

unsigned int lookup[20]={0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3};

void print_usage() {
    printf("Error in syntax ! Usage: STRS -l(cut_off)=num -k(motif length)=num\n");
}

int main (int argc, char *argv[]) 
{
   // declaring necessary variables to read raw dna file
   FILE *fp1;
   unsigned long long i=0, j=0; 
   unsigned long overall_STR_count=0, overall_Tandem_count=0, overall_Inter_STR_base =0, overall_STR_base=0; 
   int c;
   opterr=0; counter=0; atomic_counter=0;
   
   // opening the data files
   fp1 = fopen(argv[1], "rb");
   if( fp1 == NULL )  {
      perror ("Error opening file");
      return(-1);
   }
    
   //Specifying the expected options
   //The two options l, and k expect numbers as argument
    while ((c = getopt(argc, argv,"l:k:")) != -1) {
        switch (c) {
             case 'l' : cut_off= atoi(optarg);// minumum STR length(CUT_OFF)
                        break;
             case 'k' : k = atoi(optarg); // STR motif length
                        break;
             default: print_usage(); 
                 exit(EXIT_FAILURE);
        }
    }
    
   
   // checking argument validity

   if (cut_off < 0)
   {
    printf("CUT_OFF value must be positive: Using deafult CUT_OFF value 12\n");
    cut_off = 12;
   }
   if (k < 0)
   {
    printf("Motif length must be positive: Using deafult value 7\n");
    cut_off = 12;
   }
  
   
   dna = (unsigned char *)malloc(MAX);
   fread(dna, MAX , 1, fp1);   
   dna_len = strlen(dna); 
   //printf("%lu \n", dna_len);

   dist = (long int **)malloc(pow(4,k)*sizeof(long int *));
    if (dist == NULL) {
     printf("Error in allocating dist \n");
    return 0;
   }

   for(i=0; i<pow(4,k); i++) { 
   dist[i] = (long int *)malloc(9*sizeof(long int));
   if (dist[i] == NULL) {
      printf("Error in allocating dist6[%llu] \n", i);
      return 0;
   }
 }
  
 for (i=0;i<pow(4,k);i++)
   {dist[i][0]=-1; dist[i][1]=-1; dist[i][2]=-1; dist[i][3]=0; dist[i][4]=-1; dist[i][5]=0; dist[i][6]=0; dist[i][7]=0; dist[i][8]=0;}
     
  

 interval = (double **)malloc(pow(4,k)*sizeof(double *));
    if (interval == NULL) {
     printf("Error in allocating dist \n");
    return 0;
   }

   for(i=0; i<pow(4,k); i++) { 
   interval[i] = (double *)malloc(5000*sizeof(double));
   if (interval[i] == NULL) {
      printf("Error in allocating dist[%llu] \n", i);
      return 0;
   }
 }



   repeats = (long int *)malloc(6*sizeof(long int));
 
   for (i=0;i<6;i++)
   repeats[i]=-1;
   
   /////////////////////////////////////////////////////////////////////////////////////////////

   strs(); // calling short tandem repeat search(strs) 
   length = (double)dna_len/(1024*1024);  
   //printf("DNA length in MB %lf\n",length);
  
  ////////////////////////////////////////////////////////////////////////////////////////////

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
          printf("%.14lf\n", entropy/length);       
        }
      else printf("0\n");
                
      total_entropy_density += entropy/length;
     }
    //printf("%lf\n",total_entropy_density); 
     // Computing entropy without dividing by sequence length
   /*for (i =0;i <pow(4,k); i++)
   {  
     entropy=0.0;
     for (j=0;j<5000;j++) // computing entropy
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
      printf("%.14lf\n", entropy);       
      total_entropy += entropy;
      
   }
     printf("%.14lf\n",total_entropy);*/
     //printf("Total number of contributing motifs is %lu out of %lf\n", atomic_counter, pow(4,k));
     //printf("Percentage of contribution motifs is %lf ", (atomic_counter/pow(4,k))*100);
     //printf("%.14lf\n",total_entropy);
    return(0);   
}

unsigned int GetIdx(unsigned long i)
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


void strs()
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
for(i=dna_len*0.5; i< dna_len-k; i++) // scanning the sequence one character at a time
{
    // genertaing 4-base codes(index) for 6-mers starts at locaion i

    idx6 = GetIdx(i);
    
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


