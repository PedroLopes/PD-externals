
/*
DTW prototype version 0.1 by PedroLopes (plopes)
-------------------------
This external object was written to puredata. It computes the similarity value between to series of numbers, received on the inlets as a pd-list, using the Dynamic Time Warping algorithm (DTW).
-------------------------
The DTW implementation is based on Andrew Slater and John Coleman's DTW code written in C, and available here: 
- http://www.phon.ox.ac.uk/files/slp/Extras/dtw.html
-------------------------
The array size inlets are based on the arraysize code available on puredata/SVN which in its turn is based on arraysize from pixlib, this here:
- http://pix.test.at/pd/pixlib
*/

/*
	arraysize -- report the size of an array
	
	usage: |arraysize <array name>|

	methods: bang, set <array name>
*/

#include <m_pd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define VERY_BIG  (1e30)


char* filestr1 = "myarray0.data";
char* filestr2 = "myarray2.data";

clock_t c1;
clock_t c2;

time_t start,end;


//---------------------------------------

static t_class *arraysize_class;

typedef struct _arraysize {
  t_object  x_obj;
  t_symbol *array_name;
  int size;
} t_arraysize;

void arraysize_set(t_arraysize *x, t_symbol *s)
{
  x->array_name = s;
}

void arraysize_size(t_arraysize *x, t_float *s)
{
  float var = atof("0.1");;
  //not working
  //x->size = atoi(*s);
  post("Set size not working in this version:%f",var);
}


void arraysize_bang(t_arraysize *pdx)
{

double dif;

c1=clock();
time(&start);

  //t_garray *garray;
double **globdist;
double **Dist;

double top, mid, bot, cheapest, total;
unsigned short int **move;
unsigned short int **warp;
unsigned short int **temp;

unsigned int I, X, Y, n, i, j, k;
//array size is not working in this version (should be fixed soon)
unsigned int xsize = 10;
unsigned int ysize = 10;
//params is bound to 1 for now
unsigned int params = 1;

unsigned int debug; /* debug flag */

float **x, **y; /*now 2 dimensional*/

FILE *file1, *file2, *glob, *debug_file, *output_file;

 /* open x-parameter file */

if ((file1=fopen(filestr1,"rb"))==NULL)
{fprintf(stderr,"File %s cannot be opened\n",filestr1);
exit(1);
}

/* open y-parameter file */

if ((file2=fopen(filestr2,"rb"))==NULL)
{fprintf(stderr,"File %s cannot be opened\n",filestr2);
exit(1);
}


printf("xsize:%d,ysize:%d,params:%d\n",xsize,ysize,params);

     if ((debug_file = fopen("debugDTW.data","wb")) == NULL)
       {fprintf(stderr,"Cannot open debug file\n");
       exit(1);
       }

     debug = 1;


if (debug==1) fprintf(debug_file,"xsize %d ysize %d params %d\n",xsize,ysize,params);

/* allocate memory for x and y matrices */

if ((x = malloc(xsize * sizeof(float *))) == NULL)
     fprintf(stderr,"Memory allocation error (x)\n");

for (i=0; i < xsize; i++)
     if ((x[i] = malloc(params * sizeof(float))) == NULL)
     fprintf(stderr,"Memory allocation error (x)\n");

if ((y = malloc(ysize * sizeof(float *))) == NULL)
     fprintf(stderr,"Memory allocation error (y)\n");

for (i=0; i < ysize; i++)
     if ((y[i] = malloc(params * sizeof(float))) == NULL)
     fprintf(stderr,"Memory allocation error (y)\n");

     /* allocate memory for Dist */

if ((Dist = malloc(xsize * sizeof(double *))) == NULL)
     fprintf(stderr,"Memory allocation error (Dist)\n");

for (i=0; i < xsize; i++)
if ((Dist[i] = malloc(ysize * sizeof(double))) == NULL)
     fprintf(stderr,"Memory allocation error (Dist)\n");

     /* allocate memory for globdist */

if ((globdist = malloc(xsize * sizeof(double *))) == NULL)
     fprintf(stderr,"Memory allocation error (globdist)\n");

for (i=0; i < xsize; i++)
if ((globdist[i] = malloc(ysize * sizeof(double))) == NULL)
     fprintf(stderr,"Memory allocation error (globdist)\n");

     /* allocate memory for move */

if ((move = malloc(xsize * sizeof(short *))) == NULL)
     fprintf(stderr,"Memory allocation error (move)\n");

for (i=0; i < xsize; i++)
if ((move[i] = malloc(ysize * sizeof(short))) == NULL)
     fprintf(stderr,"Memory allocation error (move)\n");

     /* allocate memory for temp */

if ((temp = malloc(xsize * 2 * sizeof(short *))) == NULL)
     fprintf(stderr,"Memory allocation error (temp)\n");

for (i=0; i < xsize*2; i++)
if ((temp[i] = malloc(2 * sizeof(short))) == NULL)
     fprintf(stderr,"Memory allocation error (temp)\n");

     /* allocate memory for warp */

if ((warp = malloc(xsize * 2 * sizeof(short *))) == NULL)
     fprintf(stderr,"Memory allocation error (warp)\n");

for (i=0; i < xsize*2; i++)
if ((warp[i] = malloc(2 * sizeof(short))) == NULL)
     fprintf(stderr,"Memory allocation error (warp)\n");

for (i=0; i < xsize; i++)
{
  for (k=0; k < params; k++)
    {
if (feof(file1))
  {fprintf(stderr,"Premature EOF in %s\n",filestr1);
  exit(1);
  }

  int retf = fscanf(file1,"%f ",&x[i][k]);
  post("value1:%f",x[i][k]);
  if (retf == -1)
  	post("scan error");
  

if (debug == 1)
  fprintf(debug_file,"float_x[%d %d] = %f\n",i,k,x[i][k]);
    }
}

/*read y parameter in y[]*/

for (i=0; i < ysize; i++)
{
  for (k=0; k < params; k++)
    {

if (feof(file2))
  {fprintf(stderr,"Premature EOF in %s\n",filestr2);
  exit(1);
  }

int retf= fscanf(file2,"%f ",&y[i][k]);
if (retf == -1) post("scan error");
post("value2:%f",x[i][k]);

 if (debug == 1)
   fprintf(debug_file,"float_y[%d %d] = %f\n",i,k,y[i][k]);
    }
}

post("Computing distance matrix ...\n");

//current

/*Compute distance matrix*/

for(i=0;i<xsize;i++) {
  for(j=0;j<ysize;j++) {
    total = 0;
    for (k=0;k<params;k++) {
      total = total + ((x[i][k] - y[j][k]) * (x[i][k] - y[j][k]));
    }

    Dist[i][j] = total;

    if (debug == 1)
      fprintf(debug_file,"Dist: %d %d %.0f %.0f\n",i,j,total,Dist[i][j]);
  }
}

post("Warping in progress ...\n");

/*% for first frame, only possible match is at (0,0)*/

globdist[0][0] = Dist[0][0];
for (j=1; j<xsize; j++)
	globdist[j][0] = VERY_BIG;

globdist[0][1] = VERY_BIG;
globdist[1][1] = globdist[0][0] + Dist[1][1];
move[1][1] = 2;

for(j=2;j<xsize;j++)
	globdist[j][1] = VERY_BIG;

for(i=2;i<ysize;i++) {
	globdist[0][i] = VERY_BIG;
	globdist[1][i] = globdist[0][i-1] + Dist[1][i];
	if (debug == 1)
	  fprintf(debug_file,"globdist[2][%d] = %.2e\n",i,globdist[2][i]);

	for(j=2;j<xsize;j++) {
		top = globdist[j-1][i-2] + Dist[j][i-1] + Dist[j][i];
		mid = globdist[j-1][i-1] + Dist[j][i];
		bot = globdist[j-2][i-1] + Dist[j-1][i] + Dist[j][i];
		if( (top < mid) && (top < bot))
		{
		cheapest = top;
		I = 1;
		}
	else if (mid < bot)
		{
		cheapest = mid;
		I = 2;
		}
	else {cheapest = bot;
		I = 3;
		}

/*if all costs are equal, pick middle path*/
       if( ( top == mid) && (mid == bot))
	 I = 2;

	globdist[j][i] = cheapest;
	move[j][i] = I;
	if (debug == 1) {
	  fprintf(debug_file,"globdist[%d][%d] = %.2e\n",j,i,globdist[j][i]);
	  fprintf(debug_file,"move j:%d:i:%d=%d\n",j,i,move[j][i]);
	      }
      }
}

if (debug == 1) {
  for (j=0; j<xsize; j++) {
    for (i=0; i<ysize; i++) {
      fprintf(debug_file,"[%d %d] globdist: %.2e    move: %d    \n",j,i,globdist[j][i],move[j][i]);
    }
  }
}



X = ysize-1; Y = xsize-1; n=0;
warp[n][0] = X; warp[n][1] = Y;


while (X > 0 && Y > 0) {
n=n+1;


if (n>ysize *2) {fprintf (stderr,"Warning: warp matrix too large!");
exit(1);
} 

 if (debug == 1)
   fprintf(debug_file,"Move %d %d %d\n", Y, X, move[Y][X]);

if (move[Y] [X] == 1 )
	{
	warp[n][0] = X-1; warp[n][1] = Y;
	n=n+1;
	X=X-2; Y = Y-1;
	}
else if (move[Y] [X] == 2)
	{
	X=X-1; Y = Y-1;
	}
else if (move[Y] [X] == 3 )
	{
	warp[n] [0] = X;
	warp[n] [1] = Y-1; 
	n=n+1;
	X=X-1; Y = Y-2;
      }
else {fprintf(stderr,"Error: move not defined for X = %d Y = %d\n",X,Y); 
}
warp[n] [0] =X;
warp[n] [1] =Y;

}


/*flip warp*/
for (i=0;i<=n;i++) {
  temp[i][0] = warp[n-i][0];
  temp[i][1] = warp[n-i][1];

}

for (i=0;i<=n;i++) {
  warp[i][0] = temp[i][0];
  warp[i][1] = temp[i][1];

}

post("Warping complete. Writing output file.\n");

//current2

post("%f\n",globdist[xsize-1][ysize-1]);
fprintf(stdout,"%s     %s     %.3f\n",filestr1,filestr2,globdist[xsize-1][ysize-1]);

 outlet_float(pdx->x_obj.ob_outlet, globdist[xsize-1][ysize-1]);

c2 = clock();
time(&end);


post("timecpu:%d", (c2-c1)/CLOCKS_PER_SEC);

dif = difftime (end,start);
post("Took %.2lf seconds in diff mode", dif );


//end of func
} 

/*void arraysize_bang(t_arraysize *x)
{
  t_garray *garray;

  if(!(garray = (t_garray *)pd_findbyclass(x->array_name,garray_class))) {
    pd_error(x, "%s: no such table", x->array_name->s_name);
  } else {
    //outlet_float(x->x_obj.ob_outlet, garray_npoints(garray));
    post("SIZE: %d", garray_npoints(garray));
  }
}*/



void my_list_method(t_arraysize *x, t_symbol *s, int argc, t_atom *argv)
{
  int i = 0;
  t_float f1=0, f2=0;

  while(i<argc)
  {	
	f1=atom_getfloat(argv+i);
	post("atom:%f",atom_getfloat(argv+i));
	i++;
  }

} 

void *arraysize_new(t_symbol *s)
{
  t_arraysize *x = (t_arraysize *)pd_new(arraysize_class);

  symbolinlet_new(&x->x_obj, &x->array_name);
  outlet_new(&x->x_obj, gensym("float"));

  x->array_name = s;

  return (void *)x;
}

void arraysize_setup(void)
{
  arraysize_class = class_new(gensym("arraysize"), (t_newmethod)arraysize_new, 0, sizeof(t_arraysize), CLASS_DEFAULT, A_DEFSYMBOL, 0);
  
  class_addmethod(arraysize_class,(t_method)arraysize_set,gensym("set"), A_DEFSYMBOL, 0);
  class_addbang(arraysize_class,arraysize_bang);
  class_addlist(arraysize_class, my_list_method);
  class_addmethod(arraysize_class,(t_method)arraysize_size,gensym("size"), A_FLOAT, 0);

}
