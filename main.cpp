#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "load_tracers.hpp" //correct catalogues

//#define SKIP_SEARCH
//#define SKIP_MERGE


//global variable definition

std::vector<tracer> tv;	/* tracers read in from file */
std::vector<tracer> td;	/* tracers sorted by density */
tracer tin;		          /* buffer for adding tracers to tracer vectors */

//main program

int main(int argc, char **argv)
{


  float dthresh = 0.0; //density threshold

  int n_peak_all;
  long ntd;

  char	filebase[200];	//base filename (e.g., "turbulence")
  char	fdirbase[200];	//base directory name containing snaps (e.g., "mcd.M5.512")
  char	fdir_out[200];  //output directory
  char  fdir[200];		//id* directory
  char  fsuffix[200];   //suffix on the tracer file name
  int	isnap;          //snapshot number
  int	nfiles;         //number of subvolumes in snapshot



  //MPI wallclock timers
  double t_start;
  double t_end;


  //MPI rank, size, comm
  int rank;
  int np;

  //mpi world communicator
  MPI_Comm world = MPI_COMM_WORLD;


  int nx = 100;
  int ny = 100;
  int ix;
  int iy;
  double x;
  double y;
  double *image;
  double *image_sum;

  //x-y projection
  int a0 = 0;
  int a1 = 1;

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////	Begin program
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////



  //initialize MPI
  MPI_Init(&argc,&argv);  
  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&np);

  // check number of arguments
  if(!(argc>=4))
  {
    if(rank==0)
    {
      //printf("./athena4_tracers_regions filebase snapshot nfiles [fdirbase] [dthresh] [N_grid]\n");
      printf("./find-dr filebase snapshot nfiles [fdirbase] [fdir_out] [dthresh] [b]\n");
      fflush(stdout);
    }
    MPI_Finalize();
    return 0;
  }

  //suffix on the filename
  if(argc>=9)
  {
    sprintf(fsuffix,"%s.tra",argv[8]);
  }else{
    sprintf(fsuffix,"tra");
  }

  if(argc>=10)
    nx = atoi(argv[9]);

  if(argc>=11)
    ny = atoi(argv[10]);


  //start timing
  t_start = MPI_Wtime();



  if(rank==0)
  {
    printf("Reading files...\n");
    fflush(stdout);
  }

  //get filename base
  sprintf(filebase,"%s",argv[1]);

  //get snapshot number
  isnap = atoi(argv[2]);

  /*get number of files in snapshot*/
  nfiles = atoi(argv[3]);

  //define the base working directory
  if(argc>=5)
  {
    sprintf(fdirbase,"%s",argv[4]);
  }else{
    sprintf(fdirbase,"./");
  }

  //print to screen info about the run
  if(rank==0)
  {
    printf("filename base     = %s\n",filebase);
    printf("fdir     base     = %s\n",fdirbase);
    printf("number of files   = %d\n",nfiles);
    printf("density threshold = %e\n",dthresh);
    //printf("N_grid            = %d\n",N_grid);
    fflush(stdout);
  }

  //allocate the image
  image = (double *) calloc(nx*ny,sizeof(double));
  image_sum = (double *) calloc(nx*ny,sizeof(double));

  //loop over mpi ranks
  for(int isub=rank;isub<nfiles;isub+=np)
  {

  	//create a new directory name
    sprintf(fdir,"%s/id%d",fdirbase,isub);

  	//print progress
  	printf("rank %d isub %d nfiles %d fdir %s\n",rank,isub,nfiles,fdir);



  	//load the tracers
  	ntd = load_tracers(fdir, filebase, fsuffix, fdir_out, isnap, isub, &tv, dthresh);

    printf("rank %d tv.size %ld\n",rank,tv.size());

    //add tracer to the image
    for(long tt=0;tt<tv.size();tt++)
    {
      //spatial location of tracer
      x = tv[tt].x[a0];

      if(x>=1.0)
        x-=1.0;
      if(x<0.0)
        x+=1.0;

      y = tv[tt].x[a1];
      if(y>=1.0)
        y-=1.0;
      if(y<0.0)
        y+=1.0;

      //2-d integer locations
      ix = (int) (x*((double) nx));
      iy = (int) (y*((double) ny));

      //increase the local count
      image[nx*iy + ix] += 1.0;

    }//end loop over tracer

  	//destroy the tracers
  	vector<tracer>().swap(tv);
  }

  //perform a reduction
  MPI_Reduce(image,image_sum,nx*ny, MPI_DOUBLE, MPI_SUM, 0, world);

  //make a barrier
  MPI_Barrier(world);

  //output the data to file
  if(rank==0)
  {
    //for(int i=0;i<nx;i++)
      //for(int j=0;j<ny;j++)
        //printf("i %d j %d image %e\n",i,j,image_sum[nx*j+i]);



    FILE *fp;
    char fname[200];
    sprintf(fname,"images/image.%04d.dat",isnap);
    fp = fopen(fname,"w");
    fwrite(&nx,sizeof(int),1,fp);
    fwrite(&ny,sizeof(int),1,fp);
    fwrite(image_sum,sizeof(double),nx*ny,fp);
    fclose(fp);

  }


  //free image
  free(image);


  //free image
  free(image_sum);

  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  //		    DONE!
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  if(rank==0)
  {
    t_end = MPI_Wtime();
    printf("Total time = %fs.\n",t_end-t_start);
    printf("Done!\n");
    fflush(stdout);
  }

  MPI_Finalize();
  return 0;
 
}

