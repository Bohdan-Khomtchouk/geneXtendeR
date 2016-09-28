/*
 EXPLANATION:
 This program takes in a variant list and goes through a single vcf file to find the matching variants within the vcf file.  It outputs the variants in the vcf file to a separate out vcf file that the user specifies.
 
 COMPILE CODE USING THIS COMMAND:
 gcc -o ex test
 
 USAGE assuming compiled with above command:
 ./ex <variant bed file> <gtf bed file>
 
 Example USAGE:
 ./ex test1.var test1.vcf
 
 ASSUMPTIONS:
 1.  Both the variant file and the gtf file are properly sorted in the same manner.
 2.  The line within the gtf file does not exceed 10000 characters
 3.  The user has write permissions to the executing directory
 4.  There is no header within either of the files.  The structure of the file is simply chr \t position \t position.
 
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int countlines(char ** fname){
    int ch =0;
    int lines = 0;
    FILE * fp =  fopen(*fname, "r");
    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            lines++;
        }
    }
    return lines;
}


void extractpeaks(char **f1, char **f2, char ** Rlist)
{
    char varstr[50000];
    char vcfstr[50000];
    char varstr1[50000];
    char vcfstr1[50000];
    char col1[100];
    char col2[100];
  
    
    char * pvarcol1;
    char * pvarcol2;
    char * pvarcol3;
   
    char * pvcfcol1;
    char * pvcfcol2;
    char * pvcfcol3;
    char * pvcfcol4;
    char * pvcfcol5;
    char * end_str;
    char * end_str2;
    FILE * vcffp=NULL;
    FILE * varfp=NULL;
    long i;
    long j;

    long k;
    long l;
    long l_1;
    long m;
    long n;
    
    long v1;
    long v2;
    long v3;
    long v4;
    int zero_count;
    int pflag;
    int wflag;
    
    wflag = 1;
    zero_count = 0;
    pflag = 0;
    /* open up all the files and make sure they are open, else return value 1 */
    varfp = fopen(*f1, "r");
    if (varfp == NULL)
    {
        Rprintf("Problem opening variant file\n");
    }
    vcffp = fopen(*f2, "r");
    if (vcffp == NULL)
    {
        Rprintf("Problem opening vcf file\n");
    }
    
    // Get the first lines of both files and see what's up.
    
    fgets(varstr, 50000, varfp);
    strcpy(varstr1,varstr);
    pvarcol1 = strtok_r(varstr,"\t\n",&end_str);
    k = strtol(pvarcol1,NULL,0);
    pvarcol2 = strtok_r(NULL,"\t\n",&end_str);
    i = strtol(pvarcol2,NULL,0);
    pvarcol3 = strtok_r(NULL,"\t\n",&end_str);
    m = strtol(pvarcol3,NULL,0);
    
    
    fgets(vcfstr, 10000, vcffp);
    strcpy(vcfstr1,vcfstr);
    pvcfcol1 = strtok_r(vcfstr,"\t",&end_str2);
    n = strtol(pvarcol1,NULL,0);
    pvcfcol2 = strtok_r(NULL,"\t",&end_str2);
    j = strtol(pvcfcol2,NULL,0);
    pvcfcol3 = strtok_r(NULL,"\t",&end_str2);
    l = strtol(pvcfcol3,NULL,0);
    pvcfcol4 = strtok_r(NULL,"\t",&end_str2);
    pvcfcol5 = strtok_r(NULL,"\t\n",&end_str2);
    

    strcpy(col1,pvcfcol4);
    strcpy(col2,pvcfcol5);
   

    l_1 = l;
    
    v1 = i - j;
    v2 = i - l;
    v3 = m - j;
    v4 = m - l;
    
    /*===============================================Mo'scode====================================================*/
    
        const int line_length = 100;
    
        const int arraySize = countlines(f2);
        char peaksArr[arraySize][line_length+1];
    
        int counter=0;
        char Buffer[100];
    /*============================================================================================================*/
    
    
    while(wflag > 0)
    {
        
        memset(Buffer,' ',100); //mos code
        
        // check chromosomes.  If gtf is on next chromosome then move var file appropriately.
        
        if (k < n) {
            
            

        //Rprintf("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%d\n", pvarcol1, pvarcol2, pvarcol3, n_1, j_1, l_1, col1, col2, abs((i - l_1))+1);
            pflag = 1;
        }
        else if (k > n) {
            pflag = 2;
            
        }
        else {
            if ((v1 > 0) && (v2 > 0) && (v3 > 0) && (v4 > 0)) {
                // if the peak in the gtf is all positive then we need to update the previous pointer to the current gtf peak
               
                
                strcpy(col1,pvcfcol4);
                strcpy(col2,pvcfcol5);
             
                l_1 = l;
                // set pflag to increase the gtf file
                pflag = 2;
            }
            else if ((v1 < 0) && (v2 < 0) && (v3 < 0) && (v4 < 0)) {
                
                // if all is negative then we need to which which is the closest peak and then move the peak file to the next flag.
                if ((j - m) < (i - l_1)) {
                    // the previous pointer is closre to the distance so we need to put the distance there
                    //Rprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n", pvarcol1, pvarcol2, pvarcol3, pvcfcol1, pvcfcol2, pvcfcol3, pvcfcol4, pvcfcol5, abs((j - m))+1);
                                   } else {
                    //Rprintf("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%d\n", pvarcol1, pvarcol2, pvarcol3, n_1, j_1, l_1, col1, col2, abs((i - l_1))+1);
                   
                    
                    
                }
                // set the pflag to represent a movement in the peak file
                pflag = 1;
            }
            else {
                /* add to the running total of the overlap and move ot the next query */
                zero_count = zero_count + 1;
               // Rprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\n", pvarcol1, pvarcol2, pvarcol3, pvcfcol1, pvcfcol2, pvcfcol3, pvcfcol4, pvcfcol5);
                
                
/*===============================================Mo's code====================================================*/
            sprintf(Buffer, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\n", pvarcol1, pvarcol2, pvarcol3,pvcfcol1,pvcfcol2,pvcfcol3,pvcfcol4,pvcfcol5 );
            strcpy(peaksArr[counter],Buffer);
            counter++;
  /*============================================================================================================*/
                
                
                pflag = 1;
            }
        }
        
        switch(pflag) {
            case 1 :
                // increase the gtf file
                if (fgets(varstr, 50000, varfp) != NULL) {
                    strcpy(varstr1,varstr);
                    pvarcol1 = strtok_r(varstr,"\t\n",&end_str);
                    k = strtol(pvarcol1,NULL,0);
                    pvarcol2 = strtok_r(NULL,"\t\n",&end_str);
                    i = strtol(pvarcol2,NULL,0);
                    pvarcol3 = strtok_r(NULL,"\t\n",&end_str);
                    m = strtol(pvarcol3,NULL,0);
                } else { wflag = 0; };
                break;
            case 2 :
                // increase the peak file
                if (fgets(vcfstr, 10000, vcffp) != NULL) {
                    strcpy(vcfstr1,vcfstr);
                    pvcfcol1 = strtok_r(vcfstr,"\t",&end_str2);
                    n = strtol(pvcfcol1,NULL,0);
                    pvcfcol2 = strtok_r(NULL,"\t",&end_str2);
                    j = strtol(pvcfcol2,NULL,0);
                    pvcfcol3 = strtok_r(NULL,"\t",&end_str2);
                    l = strtol(pvcfcol3,NULL,0);
                    pvcfcol4 = strtok_r(NULL,"\t",&end_str2);
                    pvcfcol5 = strtok_r(NULL,"\t\n",&end_str2);
                } else { wflag = 0; }
                break;
        }
        
        v1 = i - j;
        v2 = i - l;
        v3 = m - j;
        v4 = m - l;
    }
    fclose(varfp);
    fclose(vcffp);
        
        for(int x = 0; x<counter;x++){
            *(Rlist+x) = peaksArr[x];
        }
    
}

static const R_CMethodDef cMethods[] = {
    {"extractpeaks", (DL_FUNC) &extractpeaks, 3},
  	NULL
};

void R_init_extractpeaks(DllInfo *info)
{
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

