#include <stdio.h>
#include <vector>
#include "read_athena_header.hpp"
#include "load_tracers.hpp"

using namespace std;

long load_tracers(char fdir[], char filebase[], char fsuffix[], char fdir_out[], int isnap, int isub, vector<tracer> *t, float dthresh)
{
  FILE *fp;
  char filename[200];   /*name of file containing the tracers*/
  long n_tracers;       /*number of tracers in the file*/
  float *d;             /*tracer densities*/
  float *x;             /*tracer x positions*/
  float *y;             /*tracer y positions*/
  float *z;             /*tracer z positions*/
  float *vx;            /*tracer x velocities*/
  float *vy;            /*tracer y velocities*/
  float *vz;            /*tracer z velocities*/
  long  *l;             /*tracer ids*/
  long ntd = 0;         /*number of tracers above the density threshold*/
  AthenaHeader *h;

  //corners of the bounding box
  float t_min[3] = {1e9,1e9,1e9};
  float t_max[3] = {-1e9,-1e9,-1e9};


  //buffer for storing tracers into the tracer vector *t
  tracer tin;
  
  if(isub==0)
  {
    /*create a new filename*/
    sprintf(filename,"%s/%s.%04d.%s",fdir,filebase,isnap,fsuffix);
  }else{
    /*create a new filename*/
    sprintf(filename,"%s/%s-id%d.%04d.%s",fdir,filebase,isub,isnap,fsuffix);
  }

  /*open tracer file*/
  if(!(fp = fopen(filename,"r")))
  {
    printf("Error opening %s in load tracers (fdir=%s, filebase=%s, fdout=%s.\n",filename,fdir,filebase,fdir_out);
    fflush(stdout); 
  }else{
    printf("Opening %s.\n",filename);
    fflush(stdout);
  }

  /* Read header */
  h = ReadAthenaHeader(fp);

  ShowAthenaHeader(h);


  /* read the number of tracer in this file */ 
  fread(&n_tracers,1,sizeof(long),fp);

  printf("n_tracers = %ld\n",n_tracers);
  fflush(stdout);

  /* Allocate buffer */
  if(!(d = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property d (n_tracers = %ld).\n",n_tracers);
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(x = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property buf.\n");
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(y = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property y.\n");
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(z = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property z.\n");
    fflush(stdout);
    exit(-1);
  }

  /* Allocate buffer */
  if(!(vx = (float *)malloc(n_tracers*sizeof(float))))
  {
    printf("Error allocating tracer property z.\n");
    fflush(stdout);
    exit(-1);
  }

    /* Allocate buffer */
    if(!(vy = (float *)malloc(n_tracers*sizeof(float))))
    {
      printf("Error allocating tracer property z.\n");
      fflush(stdout);
      exit(-1);
    }


    /* Allocate buffer */
    if(!(vz = (float *)malloc(n_tracers*sizeof(float))))
    {
      printf("Error allocating tracer property z.\n");
      fflush(stdout);
      exit(-1);
    }


    /* read density */
    fread(d,n_tracers,sizeof(float),fp);
    //printf("d[0] %e\n",d[0]);
    //printf("d[%ld] %e\n",n_tracers-1,d[n_tracers-1]);

    /* read M1 */
    fread(x,n_tracers,sizeof(float),fp);
    for(long tt=0;tt<n_tracers;tt++)
      vx[tt] = x[tt]/d[tt];
    //printf("vx[0] %e\n",vx[0]);
    //printf("vx[%ld] %e\n",n_tracers-1,vx[n_tracers-1]);

    /* read M2 */
    fread(x,n_tracers,sizeof(float),fp);
    for(long tt=0;tt<n_tracers;tt++)
      vy[tt] = x[tt]/d[tt];
    //printf("vy[0] %e\n",vy[0]);
    //printf("vy[%ld] %e\n",n_tracers-1,vy[n_tracers-1]);
    
    /* read M3 */
    fread(x,n_tracers,sizeof(float),fp);
    for(long tt=0;tt<n_tracers;tt++)
      vz[tt] = x[tt]/d[tt];
    //printf("vz[0] %e\n",vz[0]);
    //printf("vz[%ld] %e\n",n_tracers-1,vz[n_tracers-1]);

    /* read x1 */
    fread(x,n_tracers,sizeof(float),fp);
    //printf("x[0] %e\n",x[0]);
    //printf("x[%ld] %e\n",n_tracers-1,x[n_tracers-1]);
    
    /* read x2 */
    fread(y,n_tracers,sizeof(float),fp);
    //printf("y[0] %e\n",y[0]);
    //printf("y[%ld] %e\n",n_tracers-1,y[n_tracers-1]);
    
    /* read x3 */
    fread(z,n_tracers,sizeof(float),fp);  
    //printf("z[0] %e\n",z[0]);
    //printf("z[%ld] %e\n",n_tracers-1,z[n_tracers-1]);
    
    //* Allocate buffer */
    if(!(l = (long *)malloc(n_tracers*sizeof(long))))
    {
      printf("Error allocating tracer property buf.\n");
      fflush(stdout);
    }

    /* read id */
    fread(l,n_tracers,sizeof(long),fp);

    /*close tracer file*/
    fclose(fp);


    /*keep only particles with densities above or = threshold*/
    ntd = 0;
    long id_min = 100000000000;
    long id_max = -1;
    for(long tt=0;tt<n_tracers;tt++)
    {
      if(l[tt]<id_min)
        id_min = l[tt];
      if(l[tt]>id_max)
        id_max = l[tt];


      if( x[tt] < t_min[0])
        t_min[0] = x[tt];
      if( y[tt] < t_min[1])
        t_min[1] = y[tt];
      if( z[tt] < t_min[2])
        t_min[2] = z[tt];

      if( x[tt] > t_max[0])
        t_max[0] = x[tt];
      if( y[tt] > t_max[1])
        t_max[1] = y[tt];
      if( z[tt] > t_max[2])
        t_max[2] = z[tt];

      //if particle is above the threshold, keep it
      if(d[tt]>=dthresh)
      {  
        tin.id = l[tt];
        tin.d = d[tt];
        tin.x[0] = x[tt];
        tin.x[1] = y[tt];
        tin.x[2] = z[tt];
        tin.v[0] = vx[tt];
        tin.v[1] = vy[tt];
        tin.v[2] = vz[tt];

        //add to tracer list
        (*t).push_back(tin);

        //remember that we've kept a particle
        ntd++;
      }
    }
    //printf("t_min %e %e %e t_max %e %e %e\n",t_min[0],t_min[1],t_min[2],t_max[0],t_max[1],t_max[2]);
    //printf("id_min %ld id_max %ld\n",id_min,id_max);

    /*free buffer memory*/
    free(d);
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(l);


  /*free header*/
  free(h);


  //return number of tracers > density
  return ntd;

}
