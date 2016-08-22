#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "auxiliar.h"


#define minValues 16 //mininum records to compute Higuchi Fractal Dimension Algorithm ! Impirical Value


int main(){

    FILE *fp; // general purpose file pointer

    int N = 0; // number of data records
    int i = 0; // general purpose
    int j = 0; // general purpose
    int k = 0; // general purpose
    int n = 0; // general purpose
    int limit = 0; // limit of the required calculations over data

    double *sequence; // input data
    double *Xk;
    double *Ik;
    double *P_origin;
    double *P;
    double *x;
    double *y;
    double *X;
    double *Y;
    double *Yfit;

    printf("Opening the raw data file...\n");
    fp = fopen("050.txt","r");

    printf("Counting number of records...\n");
    N = countRecordsInFile(fp);
    rewind(fp);

    printf( "\tNumber of records: %d\n",N);
    if(N <= minValues){
    	printf( "\t\tError: insufficient data\n");
    	exit(1);
    }

    limit = (int)N/7;
    printf("Limit of calculations: %d\n",limit);

    printf("Allocating memory...\n");
    sequence = (double *) calloc((N*2)+1, sizeof(double));
    Xk = (double *) calloc(N+1, sizeof(double));
    Ik = (double *) calloc(N+1, sizeof(double));
    P_origin = (double *) calloc(N+1, sizeof(double));
    P = (double *) calloc( (int)(N/2)+1, sizeof(double));
    x = (double *) calloc( (int)(N/2)+1, sizeof(double));
    y = (double *) calloc( (int)(N/2)+1, sizeof(double));
    X = (double *) calloc( (int)(N/2)+1, sizeof(double));
    Y = (double *) calloc( (int)(N/2)+1, sizeof(double));
    Yfit = (double *) calloc( (int)(N/2)+1, sizeof(double));

    if(sequence == NULL){
        printf( "\t\tError: allocating memory\n");
        exit(1);
    }

    printf("Reading raw data and closing the raw data file...\n");
    i = 0;
    while(i < N){
        fscanf(fp,"%lf\n",&sequence[i]);
        i++;
    }

    printf("Computing DFT...\n");
    dft(N,sequence,Xk,Ik,limit);

    printf("Computing P_origin...\n");

    for(i=0; i<limit+1; i++){
        Xk[i] = pow(Xk[i],2);
        Ik[i] = pow(Ik[i],2);

        P_origin[i] =  (Xk[i]+Ik[i]) / (PI2 * N); //2*PI*N
    }


    printf("Computing P...\n");
    for(i=0; i<limit; i++){
        P[i] = P_origin[i+1];  //first number is highly biased
    }

    printf("Computing x...\n");
    j=0;
    for(i=2;i<=limit;i++){
        x[j] = log((PI/N)*(i));
        j++;
    }

    printf("Computing y...\n");
    j=0;
    for(i=2;i<=limit;i++){
        if(P[i]!= 0.0) y[j] = log(P[i]);
        else y[j] = log(0.000001);
        j++;
    }

    printf("Linear Regression ...\n");
    reg_Linear r;
    r = reg_LeastSquareMeans(x,y,10,(limit));

    printf("\n_______________________________\n\n");
    printf( "\tHurst = %.4lf\n_______________________________\n\n",(1-r.b)/2);

    printf("Program ended...\n");

    fclose(fp);
    free(sequence);
    free(Xk);
    free(Ik);
    free(P_origin);
    free(P);
    free(x);
    free(y);
    free(X);
    free(Y);
    free(Yfit);

    return 0;
}
