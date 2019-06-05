/* analyze_dx
 * -----------------------
 * Program that reads a potential and solvent accessibility, etc. files 
 * in dx format (such as output by APBS) 
 * Various options are implemented:
 *
 * WORK IN PROGRESS
 *
 * Author: Jan Lipfert
 * Date: September 2005
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Some constants */
#define TINY 10e-20
#define ANG3_TO_LITER 1e-27
#define N_AVOGADRO 6.02214199e23

/* Function declarations */
void exit_error();
void read_dx_file(char *fname, int *nx, int *ny, int *nz, int *npoints, double *dx, double *dy, double *dz, double **arr);
void count_charge(int nspecies, double *val, double *cons, double *pot, double *sol, int npoints, double dx, double dy, double dz);

/* Global variables*/
int verbose;  /* 0 = off, 1 = on */
int excel;    /* 0 = off, 1 = on */ 


int main(int argc, char *argv[]) {

  int run_type; /* 1 = count charges */
  int first_param_arg;
  
  int i, j, k;
  int nspecies, nx, ny, nz, npoints, nxs, nys, nzs;  
  
  FILE *pot_fp, *sol_fp;
  double *val, *cons, dx, dy, dz, *pot, dxs, dys, dzs, *sol;
    
  /* Set defaults */
  verbose = 0;
  excel = 0;
  first_param_arg = 0;

  /* Read in command line arguments */
  if(argc < 2) exit_error();
  else if(argv[1][0]!='-') exit_error();

  /* Read in all option arguments */
  for(i=1; i<argc; ++i) {
    if(argv[i][0]=='-') 
      {
	first_param_arg = i+1; 
	if(strcmp(&argv[i][1],"c")==0) 
	  {
	    run_type = 1;
	    fprintf(stderr, "Run mode: Counting charges \n");
	  } 
	else if(strcmp(&argv[i][1],"v")==0) 
	  {
	    verbose = 1;
	    fprintf(stderr, "Verbose ON \n");
	  } 
	else if(strcmp(&argv[i][1],"N")==0) 
	  {	    
	    excel = 1;
	    fprintf(stderr, "NUCLEAR (excel output) ON \n");
	  } 
	else if(strcmp(&argv[i][1],"help")==0) 
	  {
	    exit_error();
	  }
      }
    else
      {
	break;
      }
  }
  

  /* CHARGE COUNTING */
  if(run_type == 1)
    {
      if(argc-first_param_arg < 4) exit_error();
	if((pot_fp=fopen(argv[first_param_arg],"r"))==NULL) 
	  {
	    fprintf(stderr,"Error reading potential file: %s\n",argv[first_param_arg]);
	    exit(-1);
	}
	if((sol_fp=fopen(argv[first_param_arg+1],"r"))==NULL) 
	  {
	    fprintf(stderr,"Error reading solvent accessibility file: %s\n",argv[first_param_arg+1]);
	    exit(-1);
	}
	fclose(pot_fp);
	fclose(sol_fp);

	nspecies = (argc - first_param_arg - 2)/2;
	val  = (double *) malloc(nspecies * sizeof(double));
	cons = (double *) malloc(nspecies * sizeof(double));
	j=0;
	for(i=first_param_arg+2; i<argc; i += 2)
	  {	    
	    sscanf(argv[i],"%lf",&val[j]);
	    sscanf(argv[i+1],"%lf",&cons[j]);	    
	    j++;
	  }	
	fprintf(stderr, "Found %d ion specie(s) \n", nspecies);
	for(i=0; i<nspecies; i++)
	fprintf(stderr, "Ion %d has valence %f (e), present at concentration %f (M) \n", i+1, val[i], cons[i]);

	fprintf(stderr, "Reading potential file: %s\n", argv[first_param_arg]);
	read_dx_file(argv[first_param_arg], &nx, &ny, &nz, &npoints, &dx, &dy, &dz, &pot);
	fprintf(stderr, "Reading solvent accessibility file: %s\n", argv[first_param_arg+1]);
	read_dx_file(argv[first_param_arg+1], &nxs, &nys, &nzs, &npoints, &dxs, &dys, &dzs, &sol);
	
	if(nx != nxs || ny != nys || nz != nzs || fabs(dx - dxs) > TINY || fabs(dy - dys) > TINY || fabs(dz - dzs) > TINY)
	  {
	    fprintf(stderr, "Solvent accessibility and potential have differing ni or dr! \n");
	    return -1;
	  }
	
	//	fprintf(stderr, "pot[0] = %f  pot[1] = %f  pot[2] = %f\n", pot[0], pot[1], pot[2]);
	//      fprintf(stderr, "sol[0] = %f  sol[1] = %f  sol[2] = %f\n", sol[0], sol[1], sol[2]);

	count_charge(nspecies, val, cons, pot, sol, npoints, dx, dy, dz);
	
    }
  
  
  return 0;

}

/* Function: count_charge
 * ----------------------
 * Given solvent accessibility and potential arrays,
 * this function computes the 
 * - total number of ions
 * - bulk number of ions
 * - excess number of ions (= total - bulk)
 * in the box 
 * for all nspecies ions with valence = val[i] and concentrations cons[i]
 */
void count_charge(int nspecies, double *val, double *cons, double *pot, double *sol, int npoints, double dx, double dy, double dz)
{
  int i, sp;
  double ave_pot;
  double total[nspecies];
  double bulk[nspecies];
  double excess[nspecies];
  
  double dV, V;

  dV = dx*dy*dz;
  /* Compute average potential */ 
  ave_pot = 0.0;
  V = 0.0;
  for(i=0; i<npoints; i++)
    {  
      if(fabs(sol[i]) > TINY)
	{
	  ave_pot += pot[i];
	  V += dV;
	}
    }
  ave_pot /= (double) npoints;
  fprintf(stderr, "Average potential = %f\n", ave_pot);

  /* Take out the average from pot */
  for(i=0; i<npoints; i++)
    {  
      pot[i] = pot[i] - ave_pot;
    }

  if(excel)
    {
      for(sp=0; sp<nspecies; sp++)
	{
	  printf("Val %d   Total %d    Bulk %d     Excess %d    ", sp+1, sp+1, sp+1, sp+1);
	}
      printf("\n");
    }

  for(sp=0; sp<nspecies; sp++)
    {
      total[sp] = 0.0;
      for(i=0; i<npoints; i++)
	{  
	  if(fabs(sol[i]) > TINY)
	    {
	      total[sp] += exp(-val[sp]*pot[i]); 
	    }
	}
      
      total[sp] = total[sp] * exp(-val[sp]*ave_pot) * cons[sp] * N_AVOGADRO * ANG3_TO_LITER * dV; 
      bulk[sp]  = cons[sp] * N_AVOGADRO * ANG3_TO_LITER * V;
      excess[sp]= total[sp] - bulk[sp];
      
      printf("%4.2f  %9.4f  %9.4f  %9.4f    ", val[sp], total[sp], bulk[sp], excess[sp]);
    }
  printf("\n");

  return;
}


/* Function: read_dx_file
 * ---------------------------
 * This function reads a dx file, like the ones outputted from APBS
 * First, the header is read, and the parameters dx, dy, dz and nx, ny, nz are extracted,
 * making fairly specific assumptions on the syntax
 * 
 * Memory for the the grid point array arr is allocated and the array is read in.
 * Note the somewhat crazy pointer arithmetic to allow this function to allocate memory for
 * the array
 */
void read_dx_file(char *fname, int *nx, int *ny, int *nz, int *npoints, double *dx, double *dy, double *dz, double **arr)
{
  char buf[99999];
  FILE *fp;
  double temp;
  int counter, tempint;

  if((fp=fopen(fname,"r"))==NULL) 
    {
      fprintf(stderr,"Error reading file: %s\n",fname);
      exit(-1);
    }

  *nx = *ny = *nz = *npoints = 0;
  *dx = *dy = *dz = 0.0;

  /* Read the dx file header info */
  while(fgets(buf,99999,fp)!=NULL) 
    {
      if(buf[0]!='#') 
	{
	  if(strncmp(buf,"object 1 class gridpositions counts", 34)==0)
	    {
	      sscanf(buf, "object 1 class gridpositions counts %d %d %d", nx, ny, nz);
	    }
	  else if(strncmp(buf,"delta", 5)==0)
	    {
	      sscanf(buf, "delta %lf %lf %lf", dx, &temp, &temp);
	      fgets(buf,99999,fp);
	      sscanf(buf, "delta %lf %lf %lf", &temp, dy, &temp);
	      fgets(buf,99999,fp);
	      sscanf(buf, "delta %lf %lf %lf", &temp, &temp, dz);
	    }
	  else if(strncmp(buf,"object 3 class array type double rank 0 items", 40)==0)
	    {
	      sscanf(buf, "object 3 class array type double rank 0 items %d", npoints);
	      break;
	    }
	}
    }
  
  /* Check the header info */
  if(!*nx || !*ny || !*nz || !*npoints)
    {
      fprintf(stderr, "(nx, ny, nz) values not found in header of file %s\n", fname);
      exit(-1);
    }
  if(*dx < TINY || *dy < TINY || *dz < TINY )
    {
      fprintf(stderr, "(dx, dy, dz) values not found in header of file %s\n", fname);
      exit(-1);
    }
  if(*npoints != *nx * *ny * *nz)
    {
      fprintf(stderr, "npoints (%d) not equal product of ni (%d, %d, %d) \n", *npoints, *nx, *ny, *nz);
    }

  // fprintf(stderr, "nx = %d  ny = %d  nz = %d \n", *nx, *ny, *nz);
  // fprintf(stderr, "dx = %f  dy = %f  dz = %f \n", *dx, *dy, *dz);
  
  /* Allocate memory for the big array */
  *arr = (double *) calloc(*npoints , sizeof(double));
  if(*arr == NULL) 
    {
      fprintf(stderr, "Memory allocation failed! \n");
      exit(-1);
    }
  
  /* Read the grid array */
  counter = 0;
  while(fgets(buf,99999,fp)!=NULL) 
    {      
      tempint = sscanf(buf, "%lf %lf %lf", &(*arr)[counter], &(*arr)[counter+1], &(*arr)[counter+2]);
      counter += tempint;
    }
  if(counter != *npoints)
    {
      fprintf(stderr, "Read wrong number of points: npoints (%d), read (%d) \n", *npoints, counter);
      exit(-1);
    }
  //  fprintf(stderr, "Read points: npoints (%d), read (%d) \n", *npoints, counter);


  fclose(fp);
  return;
}


void exit_error() {
	fprintf(stderr,"Usage: analyze_dx <option>  [input parameters]  \n");
	fprintf(stderr,"       -c                   -> Count up charges\n");
	fprintf(stderr,"          [input parameters] = <pot file> <sol file> <valence 1> <concentration 1> [more ions] \n");	
	fprintf(stderr,"       -v                   -> verbose ON  \n");	
	fprintf(stderr,"       -help                -> print this message\n");
	fprintf(stderr,"       -N                   -> The NUCLEAR option (excel output) \n");
	fprintf(stderr,"\n");
        fprintf(stderr," Copyright (C) 2005 Jan Lipfert\n");
        fprintf(stderr," analyze_dx comes with ABSOLUTELY NO WARRANTY; for details see the file LICENSE\n");
        fprintf(stderr," This is free software, and you are welcome to redistribute it\n");
        fprintf(stderr," under certain conditions; see LICENSE for details\n");
	exit(-1); 
}
