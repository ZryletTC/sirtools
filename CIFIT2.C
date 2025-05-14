#define VERSION "2.0"  /* "980209.1" */

/***********************************************************************/
/*  cifit.c - c version of ted mcclung's sifit program                 */
/*  see Muhandiram & McClung, J. Magn. Reson. 71, 187-192 (1987)       */
/*      but using Wilkinson's derivatives for eigenvectors             */
/*                                                                     */
/*  Written by Alex D. Bain, Dept of Chemistry, McMaster University    */
/*     (bain@mcmaster.ca)                                              */
/*    Copyright (c) Alex D. Bain, 1994                                 */
/*                                                                     */
/*  The program needs two files for input: a mechanism file and        */
/*  a data file.                                                       */
/*                                                                     */
/*  The data file just contains the observed magnetizations as a       */
/*  function of time and has the simpler format, as follows:           */
/*   1st line:  A title (<80 characters)                               */
/*   2nd line:  number of time values (number of spectra)              */
/*   3rd...     values of the time, with values of the magnetizations  */
/*       On the data lines, values should be separated by whitespace   */
/*       characters - tabs, spaces, carriage returns, otherwise        */
/*       format is free.                                               */
/*                                                                     */
/*  The mechanism file is as follows:                                  */
/*  The splitting into lines is not needed, except for the title       */
/*  Values only need be separated by whitespace, but pretty formatting */
/*        helps everyone.                                              */
/*   1st line:  A title (<80 characters)                               */
/*   2nd line:  Number of sites, number of rate processes              */
/*   3rd line:  Initial values of 1/T1 for each site                   */
/*   4th line:  Initial values of M(infinity) for each site            */
/*   5th line:  Initial values of M(zero) for each site                */
/*   And then, for each rate process                                   */
/*    Initial value for the rate                                       */
/*    Number of non-zero, off-diagonal elements of the exchange        */
/*        matrix for this process (an integer)                         */
/*     For each of these non-zero, off-diagonal elements,              */
/*      Row number and column number (integers), followed by the       */
/*          the value (a real number).                                 */
/*      Note that since this is C, row and column numbers start at 0!  */
/*                                                                     */
/*  For example, here is a mechanism file for a two-site, unequal      */
/*  populations case. Comments not allowed in the actual file (yet!)   */
/*     Mechanism for furfural at 170k       ;title                     */
/*     2      1          ;2 sites, 1 process                           */
/*    .99    .3          ;1/T1 for major site, minor site              */
/*     10    2.14        ;M(infinity) for major, minor                 */
/*    -9.64  2.23        ;M(0) for major, minor                        */
/*    .3                 ;Value for rate                               */
/*     2                 ;2 non-zero, off-diagonal elements of matrix  */
/*     0 1   1.0         ;element [0][1], rate of minor -> major       */
/*     1 0   .2          ;element [1][0], rate of major -> minor       */
/*                                                                     */
/*  Note that the actual rate major->minor is the element of the       */
/*    exchange matrix (0.2, in this case) times the rate (0.3)         */
/*                                                                     */
/***********************************************************************/

/***********************************************************************/
/*  Compiling the Program                                              */
/*   The program is written to accomodate ANSI and old-style (BOO!)    */
/*   function prototypes, based on predefined preprocessor directives  */
/*   If __STDC__ (note two underscores on each side) is true, we get   */
/*   ANSI style.  This is set in Turbo C by setting ANSI Keywords ON   */
/*   in the compiler option.  In Unix, this may or may not be pre-     */
/*   defined for any given ANSI compiler.  To make sure, compile with  */
/*     cc cifit.c -D__STDC__ -o cifit -lm                              */
/*   On the SunOS compiler, which does not like ANSI prototypes, use   */
/*     cc cifit.c -D__SUN__ -o cifit -lm                               */
/*   Note the two underscores on each side of __SUN__.                 */
/***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define MACH_EPS 2.220446049250313e-016 /* for eigenvalue routine */
#define BASIS   2
#define PI      3.141592653589793
#define TOLER   .0000001

#define TRUE      1
#define FALSE     0

#define sqr(X) ((X) * (X))

#define MAXIT 100
#define MAXEX 6
#define MAXWIDTH 3
#define OUTW 6

typedef struct
  {
    double **a;       /* exchange matrix */
    double **u;       /* matrix of eigenvectors */
    double **uinv;    /* inverse of eigenvectors */
    double **dmdp;    /* matrix of partial derivatives */
    double *elambd;   /* eigenvalues of exchange matrix */
    double *beta;     /* parameters: rates, T1's, initial & equil. magns. */
    double *magn;     /* calculated magnetizations  */
  }  ex_matrix;

typedef struct
 {              /* structure to define an exchange process */
   int nzel;    /* number of non-zero elements, since arrays are dynamic */
   int *row;           /* array with row numbers */
   int *col;           /* array with column numbers */
   double *element;    /* array with matrix elements */
 }  ex_process;

/*  function prototypes  */

#if __STDC__        /* ansi prototypes */
  double **dmatrix(int nrl,int nrh,int ncl,int nch);
  int *ivector(int nl,int nh);
  double *dvector(int nl,int nh);
  void free_ivector(int *v,int nl,int nh);
  void free_dvector(double *v,int nl,int nh);
  void free_dmatrix(double **m, int nrl,int nrh,int ncl,int nch);
  void gen_matrix(ex_matrix b,ex_process *primex);
  double dudk(int l,int i, ex_matrix g, double **pi);
  double dlamdk(int i,ex_matrix g,double **pi);
  double duindk(int i,int j,double **dudka, ex_matrix g);
  double dudt(int l,int i, int m, ex_matrix g);
  void msitef(ex_matrix,int,int,ex_process *,double);
  void get_eigenv(ex_matrix, int);
  void lubksbd(double **a,int n,int *indx,double *b);
  void ludcmpd(double **a,int n,int *indx,double *d);
  void printparms(FILE *fout,double chsq, ex_matrix g);
  void varyv(char *str, int *mf, int istart, int istop,
	     int *lista, ex_matrix g);
  int skipjunk(FILE *fp);
  int dofit(double *,int,int *,int,double **,double **,double *,
	   ex_matrix, ex_process *, FILE *);
  void endup(int *lista,int mfit,double **covar,double best_chisq,
	     ex_matrix g, ex_process *p,FILE *out);
  void mrqmin(double *sig,int ma,int *lista,int mfit,
	     double **covar,double **alpha,double *chisq,
	     double *alamda,ex_matrix g,ex_process *p, FILE *);
  void mrqcof(double *sig,int *lista,int mfit,
	     double **alpha,double *bb,double *chisq,int intl,
	     ex_matrix g, ex_process *p);
  void dtail_bal(int *ipiv, int *indxr, int *indxc,
		ex_matrix g, ex_process *primex, double *zmag);
  int eigen (int vec,int n,double **mat,double **eivec,double *valreal,
	   double *valim,int *cnt);
  int balance (int n,double **mat,double * skal,int *low,
	       int *high,int basis);
  int balback (int n,int low,int high,double *skal,double **eivec);
  int elmhes(int n,int low,int high,double **mat,int *perm);
  int elmtrans (int n,int low,int high,double **mat,int *perm,double **h);
  int hqr2 (int vec,int n,int low,int high,double **h,double *wr,
      double *wi,double **eivec,int *cnt);
  int hqrvec (int n,int low,int high,double **h,double *wr,double *wi,
	      double **eivec);
  int norm_1 (int n,double **v,double *wi);
  int get_file(int argc, char *rootname, char *prompt, char *ext,
	       char *opts, FILE **f);
  int comdiv (double ar,double ai,double br,double bi,double *cr,double *ci);
  double comabs (double ar ,double ai);
  void swap (double *x,double *y);

#elif __SUN__     /* old style prototypes */
  double **dmatrix();
  int *ivector();
  double *dvector();
  void free_ivector();
  void free_dvector();
  void free_dmatrix();
  void gen_matrix();
  double duindk();
  double dlamdk();
  double dudk();
  double dudt();
  void msitef();
  void get_eigenv();
  void lubksbd();
  void ludcmpd();
  void printparms();
  void varyv();
  int skipjunk();
  int dofit();
  void endup();
  void mrqmin();
  void mrqcof();
  void dtail_bal();
  int eigen ();
  int balance ();
  int balback ();
  int elmhes();
  int elmtrans ();
  int hqr2 ();
  int hqrvec ();
  int norm_1 ();
  int get_file();
  int comdiv ();
  double comabs ();
  void swap ();
#endif

/* global variables */
double *g_x;        /* times for observed points */
double *g_y;        /* observed magnetizations */
double *g_wksp;    /* work space for lu decomp */
int nsts;           /* number of sites */
int nrks;           /* number of rate processes */
int npoints;        /* number of time values */
int ndata;          /* total number of observed magns = npoints*nsites */

/*************************************************************************/
/*                  start of main program                                */
/*************************************************************************/

#if __STDC__
  main(int argc, char *argv[])

#elif __SUN__
  main(argc, argv)
    int argc;
    char *argv[];
#endif

{
	int i,j,k,l,m;
	int deriv,mfit,ma,intl,converged,*lista;
	int *ipiv, *indxr, *indxc;
	double best_chisq,chisq,sig2i,dy;
	double *sig,**covar,**alpha,*t1save,*equilz;

	ex_process primex[MAXEX];

	ex_matrix g;

	FILE *in,*data,*out;
	char ct[80];

     /* introduce ourselves and get file names */

	printf("cifit %s\n",VERSION);

	if (!get_file(argc, argv[1], "Mechanism file", ".mch", "r", &in))
	   exit(1);
	if (!get_file(argc, argv[1], "Data file", ".dat", "r", &data))
	   exit(1);
	if (!get_file(argc, argv[1], "Output file", ".out", "wt", &out))
	   exit(1);

	fprintf(out,"CIFIT %s\n",VERSION);

     /* get the labels from the files */
        skipjunk(in);
        fgets(ct,80,in);
	printf("Mechanism title: %s",ct);
	fprintf(out,"Mechanism title: %s",ct);
        skipjunk(data);
        fgets(ct,80,data);
	printf("Data title: %s",ct);
	fprintf(out,"Data title: %s",ct);

        skipjunk(in);
	fscanf(in,"%d %d",&nsts, &nrks); /* dim of exchange matrix, rates */
	if (nrks == 1)
	  fprintf(out,"Number of sites = %d, %d process\n",nsts,nrks);
	else
	  fprintf(out,"Number of sites = %d, %d processes\n",nsts,nrks);

	ma=3*nsts+nrks;                       /* size of parameter array */

     /* do memory allocation */
	g.a=dmatrix(0,nsts-1,0,nsts-1);       /* exchange matrix */
	g.u=dmatrix(0,nsts-1,0,nsts-1);       /* eigenvectors */
	g.uinv=dmatrix(0,nsts-1,0,nsts-1);    /* and inverse of u */
	g.elambd=dvector(0,nsts-1);           /* eigenvalues */
	g.beta= dvector(0, ma-1);             /* variable parameters */
	g.dmdp=dmatrix(0,nsts-1,0,ma-1);      /* partial derivs */
	t1save= dvector(0, nsts-1);           /* temporary storage for 1/T1 */
	lista=ivector(0,ma-1);                /* parameters to be varied */
	covar=dmatrix(0,ma-1,0,ma-1);
	alpha=dmatrix(0,ma-1,0,ma-1);
	if (nsts > ma) g_wksp=dvector(0,nsts-1);
	else g_wksp=dvector(0,ma-1);
	ipiv=(int *) malloc((unsigned) nsts*sizeof(int *));
	indxr=(int *) malloc((unsigned) nsts*sizeof(int *));
	indxc=(int *) malloc((unsigned) nsts*sizeof(int *));

     /* read in the guesses for parameters:
        all in one parameter array
	   0 - nsts-1:            1/T1 values
	   nsts - 2*nsts-1:       M(inf) values
	   2*nsts - 3*nsts-1:     M(0) values
     */

        skipjunk(in);
	for (i=0; i<nsts; i++) fscanf(in,"%lf ",&g.beta[i]);
        skipjunk(in);
	for (i=nsts; i<2*nsts; i++) fscanf(in,"%lf ",&g.beta[i]);
        skipjunk(in);
	for (i=2*nsts; i<3*nsts; i++) fscanf(in,"%lf ",&g.beta[i]);

	/* we actually use M(0)-M(inf) in the calculations,
	     so do subtraction */
	for (i=0;i<nsts;i++) g.beta[i+2*nsts] =
				g.beta[i+2*nsts]-g.beta[i+nsts];

     /* set up the exchange matrix */
	for (i=0; i<nsts; i++)  {    /* clear the exchange matrix */
	   for (j=0; j<nsts; j++)  g.a[i][j]=0.0;
	}

	for (i=0; i<nrks; i++)  {    /* allocate memory and set up arrays
					for the non-zero matrix elements */
	   skipjunk(in);				
	   fscanf(in,"%lf", &g.beta[3*nsts+i]);
	   skipjunk(in);
	   fscanf(in,"%d",  &primex[i].nzel);
	   primex[i].row = (int *) malloc(primex[i].nzel*(sizeof(int)));
	   primex[i].col = (int *) malloc(primex[i].nzel*(sizeof(int)));
	   primex[i].element =
		     (double *) malloc(primex[i].nzel*(sizeof(double)));
	   skipjunk(in);	     
	   for (j=0; j<primex[i].nzel; j++)   {
	      fscanf(in, "%d %d %lf", &primex[i].row[j], &primex[i].col[j],
				     &primex[i].element[j]);
	      /* change sign - off-diagonal elements are negative ! */
	      primex[i].element[j] = - primex[i].element[j];
	   }
	}

     /* read in observed values */
        skipjunk(data);
	fscanf(data,"%d",&npoints);         /* number of time points */
	fprintf(out,"%d Time points read in\n", npoints);

	ndata=npoints*nsts;                 /* total number of data */
	g_x=dvector(0,ndata-1);             /* space for input x */
	g_y=dvector(0,ndata-1);             /* and y  */
	sig=dvector(0,ndata-1);
	g.magn=dvector(0,nsts-1);           /* value of magnetizations */
	equilz=dvector(0,nsts-1);           /* z magnetizations at equil */

	/* read in the actual data */
	skipjunk(data);
	k=0;
	for (i=0; i<npoints; i++) {
	   fscanf(data, "%lf", &g_x[i]);
	   for (j=0;j<nsts;j++) fscanf(data,"%lf ", &g_y[k++]);
	   k -= nsts;
	   fprintf(out, "%6.3f:", g_x[i]);
	   for (j=0;j<nsts;j++) fprintf(out,"%8.4f ", g_y[k++]);
	   fprintf(out,"\n");
	}
	for (i=0; i<ndata; i++)  sig[i]=1.0;  /* all the same weights */

     /* set up initial values */

	gen_matrix(g, primex);
	deriv=0;  intl=1;  chisq=0.0;  m=0;
	get_eigenv(g, intl);
	for (i=0;i<npoints;i++) {
	    if (i>0) intl=0;
	    msitef(g,deriv,intl,primex,g_x[i]);
	    for (l=0;l<nsts;l++)  {
		sig2i=1.0/(sig[m]*sig[m]);
		dy=g_y[m]-g.magn[l];
		chisq += dy*dy*sig2i;
		m++;
	    }
	}

	fprintf(out,"\nInitial values of parameters:\n");
	printparms(out,chisq,g);
	fprintf(out,"\nIteration # 0\n");

     /* do detailed balance calculation */

	for (i=0; i<nsts; i++)
          {
             t1save[i] = g.beta[i];
             g.beta[i] = 0.0;  /* no relax'n */
          }
	dtail_bal(ipiv,indxr,indxc,g,primex,equilz);
	printf("\nCalculated values for z magns are:\n");
	for (i=0;i<nsts;i++) printf(" No.%2d=%8.4f",i,equilz[i]);
	printf("\n");
	printf("Normalized values for z magns are:\n");
	for (i=0;i<nsts;i++)
	    printf(" No.%2d=%8.4f",i,g.beta[i+nsts]/g.beta[2*nsts-1]);
	printf("\n\n");
	for (i=0; i<nsts; i++) g.beta[i] = t1save[i]; /* restore relax'n */

     /* decide which parameters we want to vary */

	mfit=0;
	varyv("1/T1", &mfit, 0, nsts-1, lista,g);
	varyv("M(inf)", &mfit, nsts, 2*nsts-1, lista,g);
	varyv("M(0)-M(inf)", &mfit, 2*nsts, 3*nsts-1, lista,g);
	varyv("Rates", &mfit, 3*nsts, 3*nsts+nrks-1, lista,g);
	printf("Varying %2d parameters, with numbers ",mfit);
	for (i=0;i<mfit;i++) printf("%2d ",lista[i]);

     /* get started in earnest */

	converged =
	   dofit(sig,ma,lista,mfit,covar,alpha,&best_chisq,g,primex,out);

	if (converged)
	   endup(lista, mfit, covar, best_chisq, g, primex, out);
	else printf("\nNot converged!\n");

	intl=-1;   /* release the memory msitef and get_eigenv grabbed */
	msitef(g,deriv,intl,primex,g_x[i]);
	get_eigenv(g,intl);

	fclose(in);
	fclose(out);
	fclose(data);

}   /* end of main program */

/********************************************************************/
#if __STDC__
  int dofit(double *sig,int ma,int *lista,int mfit,
	     double **covar,double **alpha,double *pchisq,
	     ex_matrix g, ex_process *p, FILE *out)

#elif __SUN__
  int dofit(sig, ma, lista, mfit, covar, alpha, pchisq, g, p, out)
    double *sig;
    int ma, *lista, mfit;
    double **covar, **alpha, *pchisq;
    ex_matrix g;
    ex_process *p;
    FILE *out;
#endif

/* driver for the marquardt fit routine */

{
	int i,j,k,itst;
	double alamda,ochisq;

	printf("\nStarting the calculation\n");

	/*  initialize the marquardt routine */
	alamda = -1;
	mrqmin(sig,ma,lista,mfit,covar,alpha,
	       pchisq,&alamda,g,p,out);

	/*  set up the iteration  */
	k=1;
	itst=0;
	while (itst < 2) {
		printf("\n%s %2d %17s %10.7f %10s %9.2e\n","Iteration #",k,
			"chi-squared:",*pchisq,"alamda:",alamda);
		fprintf(out,"\n%s %2d %10s %9.2e %10s %12.8f\n",
			 "Iteration #",k,"alamda:",alamda,"chi^2 =",*pchisq);
		k++;
		ochisq=*pchisq;
		mrqmin(sig,ma,lista,mfit,covar,alpha,pchisq,&alamda,g,p,out);
		if ((*pchisq - ochisq) > TOLER)
			itst=0;
		else if (fabs(ochisq-*pchisq) / *pchisq < 0.1)
			itst++;

	}

	/* check whether converged, finish up */
	if ((fabs(*pchisq-ochisq)/(*pchisq) < 0.1) && (alamda < 0.001))  {
	   alamda=0.0;
	   mrqmin(sig,ma,lista,mfit,covar,alpha,
	       pchisq,&alamda,g,p,out);
	   ochisq=*pchisq/(ndata-mfit);
	   printf("\nFinal Values and Uncertainties:\n");
	   for (j=0; j<mfit; j++)   {
	      printf("# %2d = %12.6f +/- %12.6f\n",lista[j],
			      g.beta[lista[j]], sqrt(ochisq*covar[j][j]));
	   }

	   printf("Normalized Covariance matrix\n");
	   for (j=0; j<mfit; j++)   {
	      for (i=0;i<mfit;i++)   {    /* print out normalized chi^2 */
		printf("%8.4f",(covar[i][j]/sqrt(covar[i][i]*covar[j][j])));
	      }
	   printf("\n");
	   }

	   for (j=0; j<mfit; j++)
	     {
	       for (i=0;i<mfit;i++)  covar[i][j] *=ochisq;
	     }  /* true covariance matrix to pass back to main routine  */

	   return(1);    /* converged */
	}
	else return(0);  /* not converged */

}

/********************************************************************/
#if __STDC__
  void endup(int *lista, int mfit, double **covar,double best_chisq,
	     ex_matrix g, ex_process *p, FILE *out)

#elif __SUN__
  void endup(lista, mfit, covar, best_chisq, g, p, out)
    int  *lista, mfit;
    double **covar, best_chisq;
    ex_matrix g;
    ex_process *p;
    FILE *out;
#endif

{
     int i,j,k,deriv,intl;
     double nchisq,chisq,denom,yc;
     double time,sttime,fintime,timincr;

     FILE *quattro;
     char outfile[80], answ[10];

     /* print out values */
	fprintf(out,"\nFinal Values of Fitted Parameters and Uncertainties:\n");
	for (j=0; j<mfit; j++)
	   fprintf(out,"# %2d = %12.6f +/- %12.6f\n",lista[j],g.beta[lista[j]],
		     sqrt(covar[j][j]));
	fprintf(out,"\nFinal values of all parameters:\n");
	printparms(out,best_chisq,g);
	fprintf(out,"\nCovariance matrix\n");
	for (i=0; i<mfit; i++)  {
	   for (j=0; j<mfit; j++)  fprintf(out,"%10.6f ",covar[i][j]);
	   fprintf(out,"\n");
	}

        /* generate the matrix with the final parameters */
	gen_matrix(g,p);

	intl=0;
	get_eigenv(g, intl);

	fprintf(out,"\nObserved and Calculated Values\n");
	k=0;   deriv=FALSE;
	chisq=0.0;    intl=0;
	for (i=0;i<npoints;i++)  {
	   fprintf(out,"%6.2f  Calcd:",g_x[i]);
	   msitef(g,deriv,intl,p,g_x[i]);
	   for (j=0;j<nsts;j++) fprintf(out," %8.4f",g.magn[j]);
	   fprintf(out,"\n         Obsd:");
	   for (j=0;j<nsts;j++) fprintf(out," %8.4f",g_y[k++]);
	   fprintf(out,"\n         Diff:");
	   k-=nsts;
	   for (j=0;j<nsts;j++)  {      /* while we are at it */
	      yc = g.magn[j]-g_y[k++];  /* calculate a normalized */
	      fprintf(out," %8.4f",yc); /* value for chi squared */
	      chisq += yc*yc;
	   }

	   fprintf(out,"\n\n");
	}  /* printing out values */

	denom=0;
	for (j=0; j<nsts; j++) denom += sqr(g.beta[j+nsts]);
	/* normalize to sum of squares of m(inf) */
	nchisq=chisq/denom;
	printf("Raw chi squared = %16.9f, Scaled by M(inf) = %16.9f\n",
		  chisq, nchisq);
	fprintf(out,"Raw chi squared = %16.9f, Scaled by M(inf) = %16.9f\n",
		  chisq, nchisq);
	printf("Percent Sqrt((scaled chisq)/(degr of freedom)) = %10.4f\n",
		100*(sqrt(nchisq/(ndata-mfit))));
	fprintf(out,"Percent Sqrt((scaled chisq)/(degr of freedom)) = %10.4f\n",
		100*(sqrt(nchisq/(ndata-mfit))));


	printf("Do you want to calculate a plot file? (y/n):");
	scanf("%s", answ);
	if ((answ[0] == 'y') || (answ[0] == 'Y'))
	     {
	       if (get_file(0, " ", "Plot file", " ", "wt", &quattro))
		 {
		   fprintf(quattro, "\"Observed and Calculated\"\n");
		   k=0;
		   for (i=0;i<npoints;i++)  {
		     fprintf(quattro,"%8.5f, ",g_x[i]);
		     msitef(g,deriv,intl,p,g_x[i]);
		     for (j=0;j<nsts;j++) fprintf(quattro," %8.4f, ",g.magn[j]);
		     for (j=0;j<nsts;j++) fprintf(quattro," %8.4f, ",g_y[k++]);
		     k-=nsts;
		     for (j=0;j<nsts;j++)  fprintf(quattro," %8.4f, ",g.magn[j]-g_y[k++]);
		     fprintf(quattro,"\n");
		   }  /* printing out values */

		   /* and a smooth curve for plotting */
		   fprintf(quattro, "\n\"Calculated Smooth Curve\"\n");
		   printf("Enter start time, final time, increment:");
		   scanf("%lf %lf %lf", &sttime, &fintime, &timincr);
		   time=sttime;
		   do
		     {
			msitef(g,deriv,intl,p,time);
			fprintf(quattro,"%8.5f,  " ,time);
			for (j=0;j<nsts;j++)
			    fprintf(quattro,"%8.4f, ",g.magn[j]);
			fprintf(quattro,"\n");
			time += timincr;
		     }  while (time < (fintime));

		   fclose(quattro);
		 }
	     }
}

/********************************************************************/
#if __STDC__
  void mrqmin(double *sig,int ma,int *lista,int mfit,
	     double **covar,double **alpha,double *chisq,
	     double *alamda,ex_matrix g,ex_process *p, FILE *out)

#elif __SUN__
  void mrqmin(sig, ma, lista, mfit, covar, alpha, chisq, alamda, g, p, out)
    double *sig;
    int ma, *lista, mfit;
    double **covar, **alpha, *chisq, *alamda;
    ex_matrix g;
    ex_process *p;
    FILE *out;
#endif

/*  set up and solve the equations for the increments in the
    parameters.  This is a general routine, and only uses
    the parameters that are varied, so matrix sizes are minimum */
{
	int i,k,kk,j,ihit,intl, done;
	static int *indx;
	static double *db,*bb,*col,**temp,ochisq,d;

	static ex_matrix tr;

	if (*alamda < 0.0) {    /* this is the initialization */
		intl=0;  /* flag for msitef, get_eigenv - they have
			    already been initalized */
		indx=ivector(0,mfit-1);
		col=dvector(0,mfit-1);
		db=dvector(0,mfit-1);
		bb=dvector(0,mfit-1);
		temp=dmatrix(0,mfit-1,0,mfit-1);
        	tr.a=dmatrix(0,nsts-1,0,nsts-1);
	        tr.u=dmatrix(0,nsts-1,0,nsts-1);
	        tr.uinv=dmatrix(0,nsts-1,0,nsts-1);
	        tr.elambd=dvector(0,nsts-1);
	        tr.beta= dvector(0, ma-1);
	        tr.dmdp=dmatrix(0,nsts-1,0,ma-1);
	        tr.magn=dvector(0,nsts-1);

		kk=mfit;
		for (j=0;j<ma;j++) {
			ihit=0;
			for (k=0;k<mfit;k++)
				if (lista[k] == j) ihit++;
			if (ihit == 0)
				lista[kk++]=j;
			else if (ihit > 1) printf("Bad LISTA permutation in MRQMIN-1");
		}
		if (kk != ma) printf("Bad LISTA permutation in MRQMIN-2");
		*alamda=0.001;
		mrqcof(sig,lista,mfit,alpha,bb,chisq,intl,g,p);
		ochisq=(*chisq);
	}

	/* normal case, do the real work */
	intl=0;
	for (j=0;j<mfit;j++) {
		for (k=0;k<mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		db[j]=bb[j];
	}

	/*  solve the equations for increments in parameters */
	ludcmpd(covar,mfit,indx,&d);
	lubksbd(covar,mfit,indx,db);

	if (*alamda == 0.0) {             /*  do we want to finish up?  */
		/* calculate inverse of alpha */
		for (j=0;j<mfit;j++) {
		   for (k=0;k<mfit;k++) temp[j][k]=alpha[j][k];
		}
		ludcmpd(temp,mfit,indx,&d);
		for (j=0; j<mfit; j++)  {
		   for (i=0; i<mfit; i++)  col[i]=0.0;
		   col[j]=1.0;
		   lubksbd(temp,mfit,indx,col);
		   for (i=0; i<mfit; i++) covar[i][j]=col[i];
		}

		free_ivector(indx,0,mfit-1);
		free_dvector(col,0,mfit-1);
		free_dvector(bb,0,mfit-1);
		free_dvector(db,0,mfit-1);
		free_dmatrix(temp,0,mfit-1,0,mfit-1);
        	free_dmatrix(tr.a,0,nsts-1,0,nsts-1);
	        free_dmatrix(tr.u,0,nsts-1,0,nsts-1);
	        free_dmatrix(tr.uinv,0,nsts-1,0,nsts-1);
	        free_dvector(tr.elambd,0,nsts-1);
	        free_dvector(tr.beta,0, ma-1);
	        free_dmatrix(tr.dmdp,0,nsts-1,0,ma-1);
	        free_dvector(tr.magn,0,nsts-1);

		return;
	}

	/* back to the normal work: add in increments and see if things
	   are any better  */
	for (j=0;j<ma;j++) tr.beta[j]=g.beta[j];
	for (j=0;j<mfit;j++)
		tr.beta[lista[j]] = g.beta[lista[j]]+db[j];

	if (mfit < OUTW)  {
	   for (j=0;j<mfit;j++)
	     {
	       printf("%9d", lista[j]);
	       fprintf(out,"%9d", lista[j]);
	     }
	   printf("\n");  fprintf(out,"\n");
	   for (j=0;j<mfit;j++)
	     {
	       printf(" %8.4f",g.beta[lista[j]]);
	       fprintf(out," %8.4f",g.beta[lista[j]]);
	     }
	   printf("\n");  fprintf(out,"\n");
	   for (j=0;j<mfit;j++)
	     {
	       printf(" %8.4f",db[j]);
	       fprintf(out," %8.4f",db[j]);
	     }
	   printf("\n");  fprintf(out,"\n");
	}
	else  {     /* too many parameters, split them up */
	   k=0;  done=0;
	   while (!done)  {
	      j=0;
	      printf("Number  ");   fprintf(out,"Number  ");
	      while ((k<mfit) && (j<OUTW))  {
		 printf("%9d", lista[k]);  fprintf(out,"%9d", lista[k++]);
		 j++;
	      }
	      printf("\n");  fprintf(out,"\n");
	      k=k-j;   j=0;
	      printf("Value     ");   fprintf(out,"Value     ");
	      while ((k<mfit) && (j<OUTW))  {
		 printf(" %8.4f",g.beta[lista[k]]);
		 fprintf(out, " %8.4f",g.beta[lista[k++]]);
		 j++;
	      }
	      printf("\n");  fprintf(out,"\n");
	      k=k-j;   j=0;
	      printf("Increment ");   fprintf(out,"Increment ");
	      while ((k<mfit) && (j<OUTW))  {
		 printf(" %8.4f",db[k]);  fprintf(out, " %8.4f",db[k++]);
		 j++;
	      }
	      printf("\n");  fprintf(out,"\n");

	      if (k==mfit) done=1;
	   }
	}

	/*  calculate fit with trial parameters  */
	mrqcof(sig,lista,mfit,covar,db,chisq,intl,tr,p);
	fprintf(out,"New chi^2 = %12.8f\n",*chisq);

	if (*chisq < ochisq)
          {    /* if better, decrease alamda  */
		*alamda *= 0.1;
		ochisq=(*chisq);  /* and remember new params  */
		for (j=0;j<mfit;j++) {
			for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
			bb[j]=db[j];
			g.beta[lista[j]]=tr.beta[lista[j]];
		}
	  }
	else
          {      /* otherwise, increase alamda and try again */
		*alamda *= 10.0;
		*chisq = ochisq;
	  }
	return;
}

/********************************************************************/
#if __STDC__
  void mrqcof(double *sig,int *lista,int mfit,
	     double **alpha,double *bb,double *chisq,int intl,
	     ex_matrix g, ex_process *p)

#elif __SUN__
  void mrqcof(sig, lista, mfit, alpha, bb, chisq, intl, g, p)
    double *sig;
    int *lista, mfit;
    double **alpha, *bb, *chisq;
    int intl;
    ex_matrix g;
    ex_process *p;
#endif

/*  calculate the coefficients for the non-linear least squares  */
{
	int k,j,i,l,m,deriv;
	double wt,sig2i,dy;

	for (j=0;j<mfit;j++) {
		for (k=0;k<=j;k++) alpha[j][k]=0.0;
		bb[j]=0.0;
	}
	*chisq=0.0;

	gen_matrix(g,p);

	get_eigenv(g, intl);

	m=0;
	deriv=1;
	for (i=0;i<npoints;i++) {
	    msitef(g,deriv,intl,p,g_x[i]);
	    for (l=0;l<nsts;l++)  {
		sig2i=1.0/(sig[m]*sig[m]);
		dy=g_y[m]-g.magn[l];
		for (j=0;j<mfit;j++) {
			wt=g.dmdp[l][lista[j]];
			for (k=0;k<=j;k++)
				alpha[j][k] += wt*g.dmdp[l][lista[k]];
			bb[j] += dy*wt;
		}
		(*chisq) += dy*dy*sig2i;
		m++;
	    }
	}

	for (j=1;j<mfit;j++)
		for (k=0;k<=j-1;k++) alpha[k][j]=alpha[j][k];
}

/********************************************************************/
#if __STDC__
  void dtail_bal(int *ipiv, int *indxr, int *indxc, ex_matrix g,
		ex_process *primex, double *zmag)

#elif __SUN__
  void dtail_bal(ipiv, indxr, indxc, g, primex, zmag)
    int *ipiv, *indxr, *indxc;
    ex_matrix g;
    ex_process *primex;
    double *zmag;
#endif

/*
   from the exchange matrix, calculate what the m(infinity) values should
   be from the principle of detailed balance.  This provides a cross check
   of the exchange mechanism for unequally populated cases, in case
   something is the wrong way around.  Not actually needed for the
   calculation, but this is a useful check.
*/

   {
	int i, *indx;
	double d;

        indx = ivector(0, nsts-2);

	/* make up matrix with only exchange terms */
	gen_matrix(g, primex);

	/* put last column in as constant terms */
	for (i=0; i<(nsts-1); i++) zmag[i]=-g.a[i][nsts-1];

	/*  solve the equations  */
	ludcmpd(g.a,nsts-1,indx,&d);
	lubksbd(g.a,nsts-1,indx,zmag);

        zmag[nsts-1] = 1.0;             /* set the last one to 1 */
        free_ivector(indx,0,nsts-2);    /* release memory */
   }                                    /* and return from dtail_bal */

/********************************************************************/
#if __STDC__
  void msitef(ex_matrix g,int deriv,int intl,ex_process *p,double ttime)

#elif __SUN__
  void msitef(g, deriv, intl, p, ttime)
    ex_matrix g;
    int deriv, intl;
    ex_process *p;
    double ttime;
#endif

/*
   this routine does the real work for the coupled exchange relaxation
   matrix.  It calculates the matrix itself, and the matrix of partial
   derivatives that is required for the non-linear least squares.  The
   memory for the matrices is all done dynamically, and uses a minimum
   amount of memory, since the order of the calculations allows us to
   use the same matrix over and over.
*/

   {
      int i,j,k,m;
      double sum;
      static double *tuin, *tuinv, *ex, *lamderiv;
      static double **rho, **uderiv, **uinvdv;

      if (intl == -1)       /* we are done, just release memory */
	{
	  free_dvector(tuin,0,nsts-1);
	  free_dvector(ex,0,nsts-1);
	  free_dmatrix(rho,0,nsts-1,0,nsts-1);
	  free_dmatrix(uderiv,0,nsts-1,0,nsts-1);
	  free_dmatrix(uinvdv,0,nsts-1,0,nsts-1);
	  free_dvector(lamderiv,0,nsts-1);
	  free_dvector(tuinv,0,nsts-1);
	  return;
	}

      /* otherwise, get started  */
      if (intl == 1) {    /* initialization */
	 tuin = dvector(0,nsts-1);
	 ex = dvector(0,nsts-1);
	 rho=dmatrix(0,nsts-1,0,nsts-1);
	 uderiv=dmatrix(0,nsts-1,0,nsts-1);
	 uinvdv=dmatrix(0,nsts-1,0,nsts-1);
	 lamderiv=dvector(0,nsts-1);
	 tuinv=dvector(0,nsts-1);
      }

      for (i=0; i<nsts; i++)  {
	 tuin[i] = 0.0;
	 ex[i] = exp(-g.elambd[i]*ttime);
	 for (j=0; j<nsts; j++)
	    tuin[i] += g.uinv[i][j]*g.beta[2*nsts+j];
      }

      for (i=0; i<nsts; i++)  {
	    g.magn[i] = 0.0;
	    for (j=0; j<nsts; j++)
	       g.magn[i] += g.u[i][j]*ex[j]*tuin[j];
	    g.magn[i] += g.beta[nsts+i];
      }

      if (deriv)  {       /* calculate partial derivatives */
      /* first do derivatives with respect to rates
	 - this is the worst  */

	 for (i=0; i<nsts; i++)  {
	    for (j=0; j<(3*nsts+nrks); j++)  g.dmdp[i][j]=0.0;
	 }

	 for (k=0; k<nrks; k++)  {
	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++)  rho[i][j]=0.0;
	    }
	    for (i=0; i<p[k].nzel; i++)
	       rho[p[k].row[i]][p[k].col[i]] = p[k].element[i];

	    /* fill in diagonal elements */
	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++)
		  if (j != i) rho[i][i] -= rho[j][i];
	    }

	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++)
		  uderiv[i][j] = dudk(i, j, g, rho);
	       lamderiv[i] = dlamdk(i, g, rho);
	    }

	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++)
		  uinvdv[i][j] = duindk(i, j, uderiv, g);
	    }

	    for (i=0; i<nsts; i++)  {
	       tuinv[i] = 0.0;
	       for (j=0; j<nsts; j++)
		  tuinv[i] += uinvdv[i][j]*g.beta[2*nsts+j];
	    }

	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++) g.dmdp[i][3*nsts+k] +=   /* first term */
		      ex[j]*(-ttime*lamderiv[j]*g.u[i][j]+uderiv[i][j])*tuin[j]
		     + ex[j]*g.u[i][j]*tuinv[j];
	    }

	 }    /* for k over all rates */

	 /* now do it for relaxation rates */
	 for (m=0;m<nsts;m++)  {
	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++)
		  uderiv[i][j] = dudt(i, j, m, g);
	       lamderiv[i] = g.uinv[i][m]*g.u[m][i];
	    }

	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++)
		  uinvdv[i][j] = duindk(i, j, uderiv, g);
	    }

	    for (i=0; i<nsts; i++)  {
	       tuinv[i] = 0.0;
	       for (j=0; j<nsts; j++)
		  tuinv[i] += uinvdv[i][j]*g.beta[2*nsts+j];
	    }

	    for (i=0; i<nsts; i++)  {
	       for (j=0; j<nsts; j++) g.dmdp[i][m] +=   /* first term */
		      ex[j]*(-ttime*lamderiv[j]*g.u[i][j]+uderiv[i][j])*tuin[j]
		     + ex[j]*g.u[i][j]*tuinv[j];
	    }

	 }  /* over all relaxation rates */

     /* now with respect m(infinity) */
	 for (i=0;i<nsts;i++)  g.dmdp[i][nsts+i] = 1.0;

     /* and w.r.t. m(0)-m(inf) */
	 for (k=0;k<nsts;k++)  {
	    for (j=0;j<nsts;j++)  {
	       sum=0.0;
	       for (i=0;i<nsts;i++)
		  sum += ex[i]*g.u[k][i]*g.uinv[i][j];
	       g.dmdp[k][2*nsts+j] = sum;
	    }
	 }


      }   /* if (deriv) */

      return;
   }

/********************************************************************/
#if __STDC__
  double dudk(int l,int i, ex_matrix g, double **pi)

#elif __SUN__
  double dudk(l, i, g, pi)
    int l, i;
    ex_matrix g;
    double **pi;
#endif


/*
    calculate partial of Uli w.r.t. K, as in McClung, eqn 9
    for each elementary process, put in exchange matrix pi
    Now replaced with Wilkinson's eigenvectors (jan 29, 98)
*/
   {
      int k,m,n;
      double du, duintem, dutem, dptem;

      du = 0.0;

      for (k=0; k<nsts; k++)   {
	 if ((k!=i)) {
	    dutem = (-g.u[l][k])/(g.elambd[k]-g.elambd[i]);
	    duintem = 0.0;
	    for (m=0; m<nsts; m++)  {
	       dptem = 0.0;
	       for (n=0; n<nsts; n++)  dptem += pi[m][n]*g.u[n][i];
	       duintem += g.uinv[k][m]*dptem;
	    }
	    du += dutem*duintem;
	 }
      }
      return(du);

   }

/********************************************************************/
#if __STDC__
  double dlamdk(int i,ex_matrix g,double **pi)

#elif __SUN__
  double dlamdk(i, g, pi)
    int i;
    ex_matrix g;
    double **pi;
#endif

/*
    calculate partial of lambda i w.r.t. K, as in McClung, eqn 11
    for each elementary process, put in exchange matrix pi
*/
   {
      int m,n;
      double dl, dltem;

      dl = 0.0;

      for (m=0; m<nsts; m++)   {
	 dltem = 0.0;
	 for (n=0; n<nsts; n++)  dltem += pi[m][n]*g.u[n][i];
	 dl += g.uinv[i][m]*dltem;
      }
      return(dl);

   }

/********************************************************************/
#if __STDC__
  double dudt(int l,int i, int m, ex_matrix g)

#elif __SUN__
  double dudt(l, i, m, g)
    int l, i, m;
    ex_matrix g;
#endif

/*
    calculate partial of Uli w.r.t. T1, as in McClung, eqn 6
    change for Wilkinson's eigenvectors
*/
   {
      int k;
      double du;

      du = 0.0;

      for (k=0; k<nsts; k++)   {
	 if ((k!=i))  {
	    du += (-g.u[l][k])*g.uinv[k][m]*g.u[m][i]/
			(g.elambd[k]-g.elambd[i]);
	 }
      }
      return(du);

   }

/********************************************************************/
#if __STDC__
  double duindk(int i,int j,double **dudka, ex_matrix g)

#elif __SUN__
  double duindk(i, j, dudka, g)
    int i, j;
    double **dudka;
    ex_matrix g;
#endif

/*
    calculate partial of U(inverse)ij w.r.t. K, as in McClung, eqn 10
    requires complete matrix of du/dk before we can start.
*/
   {
      int k,n;
      double du, dutem;

      du = 0.0;

      for (k=0; k<nsts; k++)   {
	 dutem = 0.0;
	    for (n=0; n<nsts; n++)  dutem += g.uinv[n][j]*dudka[k][n];
	 du = du - g.uinv[i][k]*dutem;  /* deriv is negative of sum */
      }
      return(du);

   }

/********************************************************************/
#if __STDC__
  void get_eigenv(ex_matrix g, int intl)

#elif __SUN__
  void get_eigenv(g, intl)
    ex_matrix g;
    int intl;
#endif

   {
      int i, j, vec, res;
      static int *indx, *cnt;
      static double d, *lami, *col, **yy;

       if (intl == -1)  /* we are done, simply release memory
			   and return */
	 {
	   free_dmatrix(yy,0,nsts-1,0,nsts-1);
	   free_dvector(lami,0,nsts-1);
	   free_dvector(col,0,nsts-1);
	   free_ivector(indx,0,nsts-1);
	   free_ivector(cnt,0,nsts-1);
	   return;
	 }

       if (intl == 1)  /* initial memory allocation */
	 {
	   yy=dmatrix(0,nsts-1,0,nsts-1);
	   lami= dvector(0,nsts-1);    /* imag part of eigenvectors
					  not actually needed */
	   col=dvector(0,nsts-1);
	   indx=ivector(0,nsts-1);
	   cnt=ivector(0,nsts-1);
	 }

       /* start the calculation */
       vec=1;

       for (i=0; i<nsts; i++)  {
	  for (j=0; j<nsts; j++)  yy[i][j]=g.a[i][j];
       }

       res=eigen(vec, nsts, yy, g.u, g.elambd, lami, cnt);
       if (res) printf("error = %d\n", res);

       /*  note that eigenvectors are not orthogonal, so we must
	   calculate inverse explicitly ! */

       for (i=0; i<nsts; i++)  {
	  for (j=0; j<nsts; j++)  {
	     yy[i][j]=g.u[i][j];
	  }
       }

       ludcmpd(yy,nsts,indx,&d);
       for (i=0; i<nsts; i++)  {
	  for (j=0; j<nsts; j++)  col[j]=0.0;
	  col[i]=1.0;
	  lubksbd(yy,nsts,indx,col);
	  for (j=0; j<nsts; j++)  g.uinv[j][i]=col[j];
       }

       return;
   }

/********************************************************************/
#if __STDC__
  void gen_matrix(ex_matrix g, ex_process *primex)

#elif __SUN__
  void gen_matrix(g, primex)
    ex_matrix g;
    ex_process *primex;
#endif
   {
       int i,j;

       for (i=0;i<nsts;i++)  {
	  for (j=0;j<nsts;j++) g.a[i][j]=0;
       }

       for (i=0; i<nrks; i++)  {
	  for (j=0; j<primex[i].nzel; j++)  {
	     if ( (primex[i].row[j] >= nsts) || (primex[i].col[j] >= nsts)
		  || (primex[i].row[j] < 0) || (primex[i].col[j] < 0) )
		fprintf(stderr, "Out of range %d %d\n", i, primex[i].row[j]);
	     else
		g.a[primex[i].row[j]][primex[i].col[j]] +=
		    primex[i].element[j]*g.beta[3*nsts+i];
	  }
       }

       /* fill in diagonal elements */
       for (i=0; i<nsts; i++)  {
	  for (j=0; j<nsts; j++)
	     if (j != i) g.a[i][i] -= g.a[j][i];
	  g.a[i][i] += g.beta[i];  /* add T1 contribution */
       }

   }

/********************************************************************/
#if __STDC__
  void printparms(FILE *fout,double chsq, ex_matrix g)

#elif __SUN__
  void printparms(fout, chsq, g)
    FILE *fout;
    double chsq;
    ex_matrix g;
#endif

{
    int i,j;

	fprintf(fout,"1/T1's\n");
	j=0;
	for (i=0;i<nsts;i++)  {
	   fprintf(fout," No.%2d=%8.4f",i,g.beta[i]);
	   j++;
	   if (j>MAXWIDTH) {
	      fprintf(fout,"\n");
	      j=0;
	   }
	}

	fprintf(fout,"\nM(inf)'s\n");
	j=0;
	for (i=0;i<nsts;i++)  {
	   fprintf(fout," No.%2d=%8.4f",i+nsts,g.beta[i+nsts]);
	   j++;
	   if (j>MAXWIDTH) {
	      fprintf(fout,"\n");
	      j=0;
	   }
	}

	fprintf(fout,"\nM(0)-M(inf)'s\n");
	j=0;
	for (i=0;i<nsts;i++)  {
	   fprintf(fout," No.%2d=%8.4f",i+2*nsts,g.beta[i+2*nsts]);
	   j++;
	   if (j>MAXWIDTH) {
	      fprintf(fout,"\n");
	      j=0;
	   }
	}

	fprintf(fout,"\nand M(0)'s for reference\n");
	j=0;
	for (i=0;i<nsts;i++)  {
	   fprintf(fout," No.%2d=%8.4f",i+2*nsts,
			 g.beta[i+2*nsts]+g.beta[i+nsts]);
	   j++;
	   if (j>MAXWIDTH) {
	      fprintf(fout,"\n");
	      j=0;
	   }
	}

	fprintf(fout,"\nRates\n");
	j=0;
	for (i=0;i<nrks;i++)  {
	   fprintf(fout," No.%2d=%8.4f",i+3*nsts,g.beta[i+3*nsts]);
	   j++;
	   if (j>MAXWIDTH) {
	      fprintf(fout,"\n");
	      j=0;
	   }
	}
	fprintf(fout,"\n");
	fprintf(fout,"Chi squared value = %12.8f\n",chsq);

}

/********************************************************************/
#if __STDC__
  void varyv(char *str, int *mf, int istart, int istop,
	     int *lista, ex_matrix g)

#elif __SUN__
  void varyv(str, mf, istart, istop, lista, g)
    char *str;
    int *mf, istart, istop, *lista;
    ex_matrix g;
#endif

{
   int i,itst,nvar,nmax,legal;

	nmax=istop-istart+1;

	printf("Values for %s are:\n",str);
	for (i=0;i<nmax;i++) printf(" No.%2d=%8.4f",i+istart,g.beta[i+istart]);
	do  {
	   legal=1;
	   printf("\nHow many do you want to vary?:");
	   scanf("%d",&nvar);
	   if ((nvar<0) || (nvar>nmax)) legal=0;
	}  while (!legal);

	if ((nvar>0) && (nvar != nmax)) {
	   do {
	      legal=1;
	      if (nvar == 1)
		 printf("Enter the parameter number\n");
	      else
		 printf("Enter %2d parameter numbers, separated by spaces\n",nvar);
	      for (i=0;i<nvar;i++)  {
		 scanf("%d",&itst);
		 if ((itst<istart) || (itst>istop)) legal=0;
		 lista[(*mf)++]=itst;
	      }
	      if (!legal)  {
		 *mf -= nvar;
		 printf("Numbers not in legal range!\n");
	      }
	   }   while (!legal);
	}

	if (nvar == nmax) {
	   printf("Varying all %2d values for %s\n",nvar,str);
	   for (i=0; i<nmax; i++) lista[(*mf)++]=i+istart;
	}

	if (nvar == 0) printf("All values of %s fixed\n",str);
}

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/* skipjunk -- skips white spaces and strings of the form #....\n
   Here .... is a comment string */
#if __STDC__
int     skipjunk(FILE *fp)

#elif __SUN__
int     skipjunk(fp)
FILE    *fp;
#endif

{
     int        c;
     
     for ( ; ; )        /* forever do... */
     {
	  /* skip blanks */
	  do
	       c = getc(fp);
	  while ( isspace(c) );
	  
	  /* skip comments (if any) */
	  if ( c == '#' )
	       /* yes it is a comment (line) */
	       while ( (c=getc(fp)) != '\n' )
		    ;
	  else
	  {
	       ungetc(c,fp);
	       break;
	  }
     }
     return 0;
}


/********************************************************************/
#if __STDC__
  int eigen (int vec,int n,double **mat,double **eivec,double *valreal,
	   double *valim,int *cnt)

#elif __SUN__
  int eigen (vec, n, mat, eivec, valreal, valim, cnt)
    int vec, n;
    double **mat, **eivec, *valreal, *valim;
    int *cnt;
#endif

   /*************************************************
   *                                                *
   *  eigen is a series of routines to generate the *
   *  eigenvalues and eigenvectors of a real,       *
   *  non-symmetric matrix.  In principle, both     *
   *  the eigenvalues and eigenvectors may be       *
   *  complex, but the eigenvectors will occur      *
   *  as complex conjugate pairs.  In this case,    *
   *  they are still packed into one matrix.        *
   *  This routine was written auf deutsch und auf  *
   *  C by Harald Schwischow, Dept. of Physics,     *
   *  McMaster University                           *
   *                                                *
   *************************************************/


{

 int    low, high, res;
 double *skal;


  if ( n < 2 ) return (1) ;
     skal = (double *) malloc(n* sizeof(double));
  if ( skal == NULL ) return(2);


res = balance (n, mat, skal, &low, &high, BASIS) ;

  if ( res != 0 )
    {
       if (skal) free (skal); return (100 + res);
    }

res = elmhes (n, low, high, mat, cnt);

  if ( res != 0 )
    {
       if (skal) free (skal) ; return (200 + res);
    }

  if ( vec )
    { res = elmtrans (n, low, high, mat, cnt, eivec);
      if ( res != 0 )
	{
	   if (skal) free (skal) ; return (300 + res);
	}
    }

res = hqr2 (vec, n, low, high, mat, valreal, valim, eivec, cnt);

  if ( res != 0 )
    {
      if ( skal ) free ( skal); return (400 + res);
    }

  if ( vec )
    { res = balback (n, low, high, skal, eivec );

      if ( res !=0 )
	{
	  if (skal) free (skal) ; return(500 + res);
	}
	  res = norm_1 (n, eivec, valim);
	  if ( res !=0 )
	    {
	      if (skal) free (skal) ; return ( 600 + res);
	    }
    }


	if (skal) free (skal) ;
	return (0);
}

/********************************************************************/
#if __STDC__
  int balance (int n,double **mat,double * skal,int *low,int *high,int basis)

#elif __SUN__
  int balance (n, mat, skal, low, high, basis)
    int n;
    double **mat, *skal;
    int *low, *high, basis;
#endif

{
 register  i,j,k,m;
 int       iter;
 double    b2,r,c,f,g,s;

 b2 = (double) ( basis * basis ); m = 0; k = n-1;

    do
      { iter = FALSE ;

	for (j=k; j >= 0; j--)
	  {  for (r = 0.0, i = 0; i <= k; i++)
	       if ( i != j ) r += fabs(mat[j][i]);
	     if ( r == 0.0 )
	       { skal[k] = (double) j;
		  if ( j != k )
		    { for (i = 0; i <= k; i++) swap(&mat[i][j],&mat[i][k]);
		      for (i = m; i < n;i++)   swap(&mat[j][i],&mat[k][i]);
		    }
		  k--; iter = TRUE;
	       }

	  }
      }
    while (iter);

    do
      { iter = FALSE;
	for (j = m; j <= k; j++)
	   { for ( c = 0.0, i = m; i <= k; i++)
		if ( i != j ) c += fabs(mat[i][j]);
		  if ( c == 0.0 )
		    { skal[m] = (double) j;
		       if ( j != m )
		    { for ( i = 0; i <= k; i++) swap(&mat[i][j],&mat[i][m]);
		      for ( i =m ; i < n; i++) swap(&mat[j][i],&mat[m][i]);
		    }
			 m++; iter = TRUE;

		    }
	   }
      }

while ( iter ) ;

 *low = m; *high = k;
  for( i = m; i <= k; i++) skal[i] = 1.0;

    do
      { iter = FALSE;
	for ( i = m ; i <= k; i++)
	   { for ( c = r = 0.0, j = m; j <= k; j++ )
	     if ( j != i )
	       { c += fabs(mat[j][i]);
		 r += fabs(mat[i][j]);
	       }
	g = r / basis; f = 1.0; s = c + r;
	while ( c < g )  { f *= basis; c *= b2; }
	g = r * basis;
	while ( c >= g)  { f /= basis; c /= b2; }
	     if ( (c+r) / f < 0.95 * s )
	       { g = 1.0 / f; skal[i] *= f; iter = TRUE;
		 for ( j = m; j < n; j++) mat[i][j] *= g;
		 for ( j = 0; j <= k; j++) mat[j][i] *= f;
	       }
	   }
      }

while (iter);
return (0);
}

/********************************************************************/
#if __STDC__
  int balback (int n,int low,int high,double *skal,double **eivec)

#elif __SUN__
  int balback (n, low, high, skal, eivec)
    int n, low, high;
    double *skal, **eivec;
#endif

 {
  register i, j, k;
  double s;

   for (i = low; i <= high; i++)
      {
       s = skal[i];
	for (j = 0; j < n; j++) eivec[i][j] *= s;
      }
   for (i = low-1; i >= 0; i--)
      {
       k = (int) skal[i];
	if ( k != i )
	  for (j = 0; j < n; j++) swap (&eivec[i][j], &eivec[k][j]);
      }
   for (i = high+1; i < n; i++)
      {
       k = (int) skal[i];
	if ( k != i )
	   for ( j = 0; j < n; j++) swap (&eivec[i][j], &eivec[k][j]);
      }
  return (0);
 }

/********************************************************************/
#if __STDC__
  int elmhes(int n,int low,int high,double **mat,int *perm)

#elif __SUN__
  int elmhes(n, low, high, mat, perm)
    int n, low, high;
    double **mat;
    int *perm;
#endif

 {
  register i, j, m;
  double   x, y;

  for (m = low+1; m < high; m++)
    {
      i = m;
      x = 0.0;
      for ( j = m; j <= high; j++)
	if ( fabs(mat[j][m-1]) > fabs(x) )
	  { x = mat[j][m-1]; i = j; }
      perm[m] = i;
      if ( i != m )
	{
	  for ( j = m-1; j < n; j++) swap (&mat[i][j], &mat[m][j]);
	  for ( j = 0; j <= high; j++) swap (&mat[j][i], &mat[j][m]);
	}
      if ( x != 0.0 )
	for ( i = m+1; i <= high; i++)
	  {
	     y = mat[i][m-1];
	     if ( y != 0.0 )
	       { y = mat[i][m-1] = y / x;
		 for ( j = m; j < n; j++) mat[i][j] -= y * mat[m][j];
		 for ( j = 0; j <= high; j++) mat[j][m] += y * mat[j][i];
	       }
	  }
    }


  return (0);
 }


/********************************************************************/
#if __STDC__
  int elmtrans (int n,int low,int high,double **mat,int *perm,double **h)

#elif __SUN__
  int elmtrans (n, low, high, mat, perm, h)
    int n, low, high;
    double **mat;
    int *perm;
    double **h;
#endif

{
  int k, i;
  int    j;

     for (i = 0; i < n; i++)
	{ for ( k = 0; k < n ; k++) h[i][k] = 0.0;
	  h[i][i] = 1.0;
	}
     for ( i = high-1; i > low; i--)
	{
	 j = perm[i];
	 for (k = i +1; k <= high; k++) h[k][i] = mat[k][i-1];
	 if ( i != j )
	   {
	    for ( k = i; k <= high; k++)
	       { h[i][k] = h[j][k] ; h[j][k] = 0.0;
	       }
	    h[j][i] = 1.0;
	   }
	}
      return(0);
}


/********************************************************************/
#if __STDC__
  int hqr2 (int vec,int n,int low,int high,double **h,double *wr,
	    double *wi,double **eivec,int *cnt)

#elif __SUN__
  int hqr2 (vec, n, low, high, h, wr, wi, eivec, cnt)
    int vec, n, low, high;
    double **h, *wr, *wi, **eivec;
    int *cnt;
#endif

{
 int      na, en, iter, i, j, k, l, m;
 double   p, q, r, s, t, w, x, y, z;

  for ( i = 0; i < n; i++)
    if ( i < low || i > high)
      { wr[i] = h[i][i]; wi[i] = 0.0; cnt[i] = 0;
      }
  en = high; t= 0.0;

  while ( en >= low )
       { iter = 0; na = en -1;

	     for (;;)
		{ for ( l = en; l > low; l--)
		    if ( fabs(h[l][l-1]) <=
			    MACH_EPS * (fabs(h[l-1][l-1]) + fabs(h[l][l])) )
		       break;

	     x = h[en][en];

	     if ( l == en )
	       { wr[en] = h[en][en] = x + t; wi[en] = 0.0;
		 cnt[en] = iter;
		 en--; break;
	       }

	     y = h[na][na]; w = h[en][na] * h[na][en];

	     if ( l == na)
	       { p = ( y - x ) * 0.5; q = p * p + w;
		 z = sqrt(fabs(q));
		 x = h[en][en] = x + t; h[na][na] = y + t;
		 cnt[en] = -iter; cnt[na] = iter;
		 if (q >= 0.0)
		   {
		     z = ( p < 0.0 ) ? ( p - z ) : ( p + z );
		     wr[na] = x + z; wr[en] = s = x - w / z;
		     wi[na] = wi[en] = 0.0;
		     x = h[en][na]; r = sqrt(x*x + z*z);
		     { if (vec)
			 { p = x / r; q = z / r;
			    for ( j = na; j < n; j++)
			       { z = h[na][j];
				 h[na][j] = q * z + p * h[en][j];
				 h[en][j] = q * h[en][j] - p * z;
			       }
			    for ( i = 0; i <= en; i++)
			       { z = h[i][na];
				 h[i][na] = q * z + p * h[i][en];
				 h[i][en] = q * h[i][en] - p * z;
			       }
			    for ( i = low; i <= high; i++)
			       { z = eivec[i][na];
				 eivec[i][na] = q * z + p * eivec[i][en];
				 eivec[i][en] = q * eivec[i][en] - p * z;
			       }
			 }
		     }
		   }
		 else
		   { wr[na] = wr[en] = x + p;
		     wi[na] = z; wi[en] = - z;
		   }


	       en -= 2;
	       break;
	       }

	     if ( iter >= MAXIT )
	       { cnt[en] = MAXIT + 1; return(en);
	       }

	     if ( ( iter != 0 ) && ( iter % 10 == 0 ) )
	       { t += x;
		 for ( i = low; i <= en; i++) h[i][i] -= x;
		    s = fabs(h[en][na]) + fabs(h[na][en-2]);
		    x = y = 0.75 * s; w = - 0.4375 * s * s;
	       }

	     iter ++;


	     for ( m = en-2; m >= l; m--)
		{ z = h[m][m]; r = x - z; s = y - z;
		  p = ( r * s - w ) / h[m+1][m] + h[m][m+1];
		  q = h[m+1][m+1] - z - r - s;
		  r = h[m+2][m+1];
		  s = fabs(p) + fabs(q) + fabs(r);
		  p /= s; q /= s; r /= s;
		  if ( m == l ) break;
		  if ( fabs(h[m][m-1]) * ( fabs(q) + fabs(r) ) <=
		       MACH_EPS * fabs(p)
		       * ( fabs(h[m-1][m-1]) + fabs(z) + fabs(h[m+1][m+1]) ) )
		     break;
		}

	     for ( i = m+2; i <= en; i++ ) h[i][i-2] = 0.0;
	     for ( i = m+3; i <= en; i++ ) h[i][i-3] = 0.0;

	     for ( k = m; k <= na; k++ )
		{
		 if ( k != m  )
		   {
		     p = h[k][k-1]; q = h[k+1][k-1];
		     r = ( k != na ) ? h[k+2][k-1] : 0.0;
		     x = fabs(p) + fabs(q) + fabs(r);
		     if ( x == 0.0 ) continue ;
		     p /= x; q /= x; r /= x;
		   }
		 s = sqrt ( p*p + q*q + r*r );
		 if ( p < 0.0 ) s = -s;
		 if ( k != m ) h[k][k-1] = -s * x;
		    else if ( l != m ) h[k][k-1] = -h[k][k-1];
		 p += s; x = p / s; y = q / s; z = r / s; q /= p;
		 r /= p;

	     for ( j = k; j < n; j++ )
	       { p = h[k][j] + q * h[k+1][j];
		 if ( k != na )
		   { p += r * h[k+2][j]; h[k+2][j] -= p * z;
		   }
		 h[k+1][j] -= p * y; h[k][j] -= p * x;
	       }
	     j = ( k+3 < en ) ?  ( k+3 ) : en ;
	     for ( i = 0; i <= j; i++ )
	       { p = x * h[i][k] + y * h[i][k+1];
		 if ( k != na )
		   { p += z * h[i][k+2]; h[i][k+2] -= p * r;
		   }
		 h[i][k+1] -= p * q; h[i][k] -= p;
	       }
	     if ( vec )
		  for ( i = low; i <= high; i++ )
		     { p = x * eivec[i][k] + y * eivec[i][k+1];
		       if ( k != na )
			 { p += z * eivec[i][k+2]; eivec[i][k+2] -= p * r;
			 }
		       eivec[i][k+1] -= p * q;
		       eivec[i][k] -= p;
		     }
		}
		}
       }

    if ( vec )
      if ( hqrvec ( n, low, high,h, wr, wi, eivec) ) return(99);
    return(0);
}

/********************************************************************/
#if __STDC__
  int hqrvec (int n,int low,int high,double **h,double *wr,double *wi,double **eivec)

#elif __SUN__
  int hqrvec (n, low, high, h, wr, wi, eivec)
    int n, low, high;
    double **h, *wr, *wi, **eivec;
#endif

 {

  register    i, k, j, l, m, en, na;
  double      p, q , r, s, t, w, x, y, z, ra, sa, vr, vi, norm;

    for ( norm = 0.0, k = 0, i = 0; i < n; i++ )
      {
       for ( j = k; j < n; j++ ) norm += fabs(h[i][j]);
       k = i;
      }
    if ( norm == 0.0 ) return(1);

    for ( en = n-1; en >= 0; en--)
      {
	p = wr[en]; q = wi[en]; na = en - 1;
	if ( q == 0.0 )
	  {
	    m = en; h[en][en] = 1.0;
	    for  ( i = na; i >= 0; i--)
	      {
	       w = h[i][i] - p; r = h[i][en];
	       for ( j = m; j <= na; j++ ) r += h[i][j] * h[j][en];
	       if ( wi[i] < 0.0 )
		 { z = w; s = r; }
	       else
		 {
		  m = i;
		  if ( wi[i] == 0.0 )
		     h[i][en] = -r / ( ( w != 0.0 ) ? (w) : ( MACH_EPS * norm));
	       else
		 {
		   x = h[i][i+1]; y = h[i+1][i];
		   q = sqr(wr[i] - p) + sqr(wi[i]);
		   h[i][en] = t = (x * s - z * r) / q;
		   h[i+1][en] = ( ( fabs(x) > fabs(z) ) ?
				  (-r -w * t) / x : (-s -y * t) / z);
		 }
		 }
	      }
	  }
     else
	if ( q < 0.0 )
	  {
	    m = na;
	    if ( fabs(h[en][na]) > fabs(h[na][en]) )
	      {
		h[na][na] = - ( h[en][en] - p ) / h[en][na];
		h[na][en] = - q / h[en][na];
	      }
	    else
	      comdiv(-h[na][en], 0.0, h[na][na] - p, q, &h[na][na], &h[na][en]);

	    h[en][na] = 1.0; h[en][en] = 0.0;
	    for ( i = na - 1; i >= 0; i--)
	      {
		w = h[i][i] - p; ra = h[i][en]; sa = 0.0;
		for ( j = m; j <= na; j++ )
		   {
		     ra += h[i][j] * h[j][na];
		     sa += h[i][j] * h[j][en];
		   }
		if ( wi[i] < 0.0 )
		   {
		    z = w; r = ra; s = sa;
		   }
		else
		   {
		     m = i;
		     if ( wi[i] == 0.0 )
			comdiv(-ra, -sa, w, q, &h[i][na], &h[i][en]);
		     else
		       {
			 x = h[i][i+1]; y = h[i+1][i];
			 vr = sqr(wr[i] - p) + sqr(wi[i]) - sqr(q);
			 vi = 2.0 * q * (wr[i] - p);
			 if ( vr == 0.0 && vi == 0.0 )
			    vr = MACH_EPS * norm *
				 ( fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(z) );
			 comdiv (x*r-z*ra+q*sa, x*s-z*sa-q*ra, vr, vi, &h[i][na], &h[i][en]);
			 if ( fabs(x) > fabs(z) + fabs(q) )
			   {
			     h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
			     h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
			   }
			 else
			   comdiv(-r -y * h[i][na], -s - y * h[i][en], z, q,
						   &h[i+1][na], &h[i+1][en]);
		       }
		   }
	      }
	  }
      }


  for ( i = 0; i < n; i++)
    if ( i < low || i > high )
      for ( k = i+1; k < n; k++ ) eivec[i][k] = h[i][k];

  for ( j = n-1; j >= low; j-- )
    {
     m = ( j <= high ) ? j : high;
     if ( wi[j] < 0.0 )
       for ( l = j-1, i = low; i <= high; i++ )
	 {
	   for ( y = z = 0.0, k = low; k <= m; k++ )
	      { y += eivec[i][k] * h [k][l];
		z += eivec[i][k] * h [k][j];
	      }
	   eivec[i][l] = y; eivec[i][j] = z;
	 }
       else
	 if ( wi[j] == 0.0 )
	   for ( i = low; i <= high; i++ )
	      {
		for ( z = 0.0, k = low; k <= m; k++ )
		    z += eivec[i][k] * h [k][j];
		eivec[i][j] = z;
	      }
    }
   return(0);
 }

/********************************************************************/

#if __STDC__
  int norm_1 (int n,double **v,double *wi)

#elif __SUN__
  int norm_1 (n, v, wi)
    int n;
    double **v, *wi;
#endif

 {
  register i, j;
  double   maxi, tr, ti;

  for ( j = 0; j < n; j++ )
    { if ( wi[j] == 0.0 )
	{ maxi = v[0][j];
	  for ( i = 1; i < n; i++)
	     if ( fabs(v[i][j]) > fabs(maxi) ) maxi = v[i][j];
	  if ( maxi != 0.0 )
	    { maxi = 1.0 / maxi;
	      for ( i = 0; i < n; i++ ) v[i][j] *= maxi;
	    }
	}
      else
	{ tr = v[0][j]; ti = v[0][j+1];
	  for ( i = 1; i < n; i++ )
	    if ( comabs (v[i][j], v[i][j+1]) > comabs (tr,ti) )
	      { tr = v[i][j]; ti = v[i][j+1];
	      }
	  if ( tr != 0.0 || ti != 0.0 )
	    for ( i = 0; i < n; i++ )
	      comdiv (v[i][j], v[i][j+1], tr, ti, &v[i][j], &v[i][j+1]);
	  j++;
	}
    }
 return(0);
 }

/********************************************************************/
#if __STDC__
  int get_file(int argc, char *rootname, char *prompt, char *ext,
	       char *opts, FILE **f)

#elif __SUN__
  get_file(argc, rootname, prompt, ext, opts, f)
    int argc;
    char *rootname, *prompt, *ext, *opts;
    FILE **f;
#endif

    {
	char infile[80];

	/* get filenames and open files */
	if (argc == 2)   /* exactly one command line argument */
	  {
	    strcpy(infile, rootname);
	    strcat(infile, ext);
	  }
	else
	  {
	    printf("%s: ", prompt);
	    scanf("%s",infile);
	  }
	if ((*f=fopen(infile,opts)) == NULL) {
	     fprintf(stderr, "cannot open file %s .\n",prompt);
	     return(0);
	}
	else return(1);
   }

/********************************************************************/
#if __STDC__
  int comdiv (double ar,double ai,double br,double bi,double *cr,double *ci)

#elif __SUN__
  int comdiv (ar, ai, br, bi, cr, ci)
    double ar,ai,br,bi,*cr,*ci;
#endif

  {
   double temp;

   if ( br == 0.0 && bi == 0.0 ) return (1);
   if ( fabs(br) > fabs(bi) )
     { temp = bi / br; br = temp * bi + br;
       *cr = (ar + temp * ai) / br;
       *ci = (ai - temp * ar) / br;
     }
   else
     { temp = br / bi; bi = temp * br + bi;
       *cr = ( temp * ar + ai ) / bi;
       *ci = ( temp * ai - ar ) / bi;
     }
   return (0);
  }

/********************************************************************/
#if __STDC__
  double comabs (double ar ,double ai)

#elif __SUN__
  double comabs(ar, ai)
    double ar, ai;
#endif

  {
   double temp;

   if ( ar == 0.0 && ai == 0.0 ) return (0.0);
   ar = fabs(ar); ai = fabs(ai);

   if ( ai > ar )
     { temp = ai; ai = ar; ar = temp; }

   return ( (ai == 0.0) ? (ar) : (ar * sqrt(1.0 + ai/ar * ai/ar)) );
  }

/********************************************************************/
#if __STDC__
  void swap (double *x,double *y)

#elif __SUN__
  void swap(x,y)
    double *x, *y;
#endif

  {
    double temp;
    temp = *x; *x = *y; *y = temp;

  }

/********************************************************************/

   /*************************************************
   *                                                *
   *     ludcmpd and lubksb are modified versions   *
   *     of the numerical recipes ludcmp and lubksb *
   *     for LU decomposition and backsubstitution  *
   *     for solution of linear equations.          *
   *     modified for double precision and zero-    *
   *     based arrays.                              *
   *                                                *
   *************************************************/

#define TINY 1.0e-20;
#if __STDC__
  void ludcmpd(double **a,int n,int *indx,double *d)

#elif __SUN__
  void ludcmpd(a, n, indx, d)
    double **a;
    int n, *indx;
    double *d;
#endif

{
	int i,imax,j,k;
	double big,dum,sum,temp;

	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) printf("Singular matrix in routine LUDCMP\n");
		g_wksp[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=g_wksp[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			g_wksp[imax]=g_wksp[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}

#undef TINY

/********************************************************************/
#if __STDC__
  void lubksbd(double **a,int n,int *indx,double *b)

#elif __SUN__
  void lubksbd(a, n, indx, b)
    double **a;
    int n, *indx;
    double *b;
#endif

{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii>-1)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

/********************************************************************/
#if __STDC__
  double **dmatrix(int nrl,int nrh,int ncl,int nch)

#elif __SUN__
  double **dmatrix(nrl, nrh, ncl, nch)
    int nrl, nrh, ncl, nch;
#endif

   /*************************************************
   *                                                *
   *  Allocate memory for a double precision matrix *
   *                                                *
   *************************************************/

   {
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) fprintf(stderr,"allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) fprintf(stderr,"allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
   }

/********************************************************************/
#if __STDC__
  int *ivector(int nl,int nh)

#elif __SUN__
  int *ivector(nl, nh)
    int nl, nh;
#endif

{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) fprintf(stderr,"allocation failure in ivector()");
	return v-nl;
}

/********************************************************************/
#if __STDC__
  double *dvector(int nl,int nh)

#elif __SUN__
  double *dvector(nl, nh)
    int nl, nh;
#endif

{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) fprintf(stderr,"allocation failure in dvector()");
	return v-nl;
}

/********************************************************************/
#if __STDC__
  void free_dvector(double *v,int nl,int nh)

#elif __SUN__
  void free_dvector(v, nl, nh)
    double *v;
    int nl, nh;
#endif

{
	free((char*) (v+nl));
}

/********************************************************************/

#if __STDC__
  void free_ivector(int *v,int nl,int nh)

#elif __SUN__
  void free_ivector(v, nl, nh)
    int *v;
    int nl, nh;
#endif

{
	free((char*) (v+nl));
}

/********************************************************************/
#if __STDC__
  void free_dmatrix(double **m, int nrl,int nrh,int ncl,int nch)

#elif __SUN__
  void free_dmatrix(m, nrl, nrh, ncl, nch)
    double **m;
    int nrl, nrh, ncl, nch;
#endif

{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

/* end of program cifit */
