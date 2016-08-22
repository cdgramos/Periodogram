#define SWAP(a,b) tmp=(a);(a)=(b);(b)=tmp
#define PI2 6.305185308
#define PI 3.141692645

//count how many data records are in a file
int countRecordsInFile(FILE *fp){
    int lines = 0;
    char ch;
    while(!feof(fp)){
      ch = fgetc(fp);
      if(ch == '\n')
        lines++;
    }
    return lines;
}

//compute average
double average(double *x, int totRecords){
	int i = 0;
	double avg = 0.0;

	for(i=0; i<totRecords; i++){
		avg += x[i];
	}

	return (avg/totRecords);
}

//integrate a given time series aka comulative sums
void integratingTimeSeries(double *x, double *y, double xb, int totRecords){
	int i = 0;
	double yp = 0.0;

    for(i=1;i<totRecords;i++){
        if(i>=1){
            y[i] = y[i-1]+(x[i]-xb);
        }
    }
}


//get the total number of distinct windows we will use
int getTotDistinctWindows(int lowerBound, int upperBound, double scale){
    return log10(upperBound / (double)lowerBound) / log10(scale) + 1.5;

}

// linear regression equation variables
typedef struct {
	double m;
	double b;
} reg_Linear;

// linear regression adaptated from TestH library
reg_Linear reg_LeastSquareMeans (
	double 	*x,
	double 	*y,
	int     low,
	int 	top
	)
{

    reg_Linear r;
	int i;
	double sx, sxx, sy, sxy;
	double xb, yb;

	sx = sxx = sy = sxy = 0.0;
	xb = yb = 0.0;

	for (i=low; i<top; i++) {
        xb +=x[i];
        yb +=y[i];
	}

	xb = xb / (top-low);
    yb = yb / (top-low);

	for (i=low; i<top; i++) {
		sx  += pow(x[i]-xb,2);
		sy  += pow(y[i]-yb,2);
		sxy += (x[i]-xb)*(y[i]-yb);
	}

	r.b = sxy / sx;
	r.m = yb - ((r.b) * xb);

	return r;
}


//get the sizes of each window we will use based on data set and on a desired scaling factor
void getDistinctWindowSizes(int minWindow,int scale, int totWindows, long *windowSizes){
	int i = 0;

	for(i=0; i<totWindows; i++){
       windowSizes[i] = (minWindow * pow(scale, i));
    }
}

int dft(long int length, double realSample[], double *Rk, double *Ik, int limit){
    long int i, j;
    double arg;
    double arg1;
    double cosarg,sinarg;
    double *tempReal=NULL;
    double *tempIm=NULL;
    double *imagSample=NULL;

    tempReal = calloc(length, sizeof(double));
    tempIm = calloc(length, sizeof(double));
    imagSample = calloc(length, sizeof(double));

    if (tempReal == NULL || tempIm == NULL || imagSample == NULL){
        printf("Error while allocating memory...");
        return(0);
    }

    arg1 = -1.0 * PI2;

    for(i=0; i<limit; i++){

        arg = arg1 * i / length;

        for(j=0; j<limit; j++){
            tempReal[i] += (realSample[j] * (cos(j * arg)) - imagSample[j] * (sin(j * arg)));
            tempIm[i] += (realSample[j] * (sin(j * arg)) + imagSample[j] * (cos(j * arg)));
        }
    }

    for (i=0; i<limit; i++){
        Rk[i] = tempReal[i];
        Ik[i] = tempIm[i];
    }

    free(tempReal);
    free(tempIm);
    free(imagSample);
    return(1);
}
