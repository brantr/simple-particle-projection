#include "search_for_groups.hpp"
#include "shock_data_types.hpp"
#include "assign_peak.hpp"	          /*routines for assigning tracers to a density peak*/
#include "comparison_functions.hpp"   /*comparison functions for sorting*/
#include "write_shock_catalogues.hpp" /*read and write catalogues*/
#include "construct_fof_catalogue.hpp"
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int search_for_groups(float dthresh, float bsq, char filebase[], char fdirbase[], char fsuffix[], char fdir_out[], int isnap, int nfiles, int rank, int np, MPI_Comm world)
{

  int dim = 3;		        //number of physical dimensions for tree searching
  int n_peak_total = 0;   //total number of peaks found by a process
  int n_peak_all = 0;     //total peaks found by all processes across all subvolumes

  int isub;               //snapshot subregion

  int flag_files = 1;     //Are there still subvolume files left to process?

  double t_start = 0;	    //starting time of group search
  double t_end   = 0;	    //ending time of group search
  double t_assign  = 0;	  //time for assigning groups
  double t_io      = 0;	  //time for reading and writing group files
  double t_io_start;      //start time for i/o
  double t_io_end;        //end time for i/o
  double t_io_max;        //max time for i/o
  double t_assign_start;  //start time for assigning tracers to peaks
  double t_assign_end;    //end time for assigning tracers to peaks
  double t_assign_max;    //max time for assigning tracers to peaks


  char fdir[200];	        //directory for writing group information

  long ntd       = 0;     //number of tracers in subvolume on this process
  long ntd_total = 0;	    //total number of tracers handled on this process
  long ntd_all   = 0;     //total number of tracers handled on all processes

  MPI_Status status;      //MPI status for communication

  int sb[2];              /*snapshot buffer*/


  //start MPI timer
  t_start = MPI_Wtime();

  //everyone starts on their own snapshot
  isub = rank-1;
  if(rank==0)
    isub = np-1;

  //loop until isnap >= nfiles
  while(flag_files)
  {
    if(rank!=0)
    {
      //create a new directory name
      sprintf(fdir,"%s/id%d",fdirbase,isub);

      printf("rank %3d: %s\n",rank,fdir);
      fflush(stdout);

      //time io
      t_io_start = MPI_Wtime(); 

      //keep only particles with densities above or = threshold
      //checked load_tracers() Aug 25, 2014
      ntd = load_tracers(fdir, filebase, fsuffix, fdir_out, isnap, isub, &tv, dthresh);

#if DEBUG >= 3
      printf("tracers loaded %d %d ntd = %ld\n",isnap,isub,ntd);
      fflush(stdout);
#endif //DEBUG

      //add to total number of tracers
      ntd_total += ntd;

      //time io
      t_io_end = MPI_Wtime(); 
      t_io += (t_io_end-t_io_start);

      //only proceed if there are particles
      //above the density threshold
      if(ntd>0)
      {

        /*sort tracer particles by density*/ 
        std::sort( tv.begin(), tv.end(), tracer_density_comparison);

#if DEBUG >= 3
        printf("sorting done %d %d ntd = %ld\n",isnap,isub,ntd);
        fflush(stdout);
#endif //DEBUG

        // store the tree data 
        data.resize(extents[ntd][dim]);

        //loop over tracers
        for(long tt=ntd-1;tt>=0;tt--)
        {
          //store tracer positions in data for tree building
          for(int k=0;k<3;k++)
            data[tt][k] = tv[tt].x[k];

          //initialize peak index
          //that will be set in refine_fof_groups()
          tv[tt].peak_index = -1;
        }


        //build the kdtree
        tree = new kdtree2(data,true);

#if DEBUG >= 3
        printf("tree made %d %d ntd = %ld\n",isnap,isub,ntd);
        fflush(stdout);
#endif //DEBUG

        //time assign 
        t_assign_start = MPI_Wtime(); 


#if DEBUG >= 3
        printf("FOF search for sub = %d, n_tracers = %ld (peak index [0] %ld)\n",isub,tv.size(),tv[0].peak_index);
        fflush(stdout);
        if(ntd<4)
          for(long tt=0;tt<ntd;tt++)
            printf("Entering find_fof_mst tt %ld peak_index %ld\n",tt,tv[tt].peak_index);
#endif //DEBUG

        //perform FOF search, must assign correct peak_index
        find_fof_mst(data, bsq, dthresh, fdir_out, isnap, isub);  //FOF SEARCH


#if DEBUG >= 3
        printf("fof mst done %d %d ntd = %ld (peak index [0] %ld)\n",isnap,isub,ntd,tv[0].peak_index);
        fflush(stdout);
#endif

        //time assign 
        t_assign_end = MPI_Wtime(); 
        t_assign += (t_assign_end-t_assign_start);


        //Build the shocks and write peak data to file
        //The correct peak index must be set before
        //entering this subroutine.

        n_peak_total += construct_fof_catalogue(ntd, fdir_out, isnap, isub);

#if DEBUG >= 3
        printf("construct fof done %d %d ntd = %ld\n",isnap,isub,ntd);
        fflush(stdout);
#endif //DEBUG

        //time io
        t_io_end = MPI_Wtime();
        t_io += (t_io_end-t_io_start);

        /*delete the tree*/
        delete tree;

        //free tracer particle vector
        data.resize(extents[1][1]);

      }else{//ntd>0

        //this subvolume has no tracers above the density
        //threshold, so we write out null shock catalogues
        write_null_shock_list_isub(fdir_out, isnap, isub);
        write_null_shock_data_isub(fdir_out, isnap, isub);
      }
    }

    //Here we communicate with rank 0 to get assigned the
    //next subvolume to operate on
    if(rank!=0)
    {
      sb[0] = isub;
      sb[1] = rank;
      MPI_Send(&sb[0], 2, MPI_INT, 0, rank, world);
      MPI_Recv(&sb[0], 2, MPI_INT, 0, rank, world, &status);
      isub = sb[0];
      if(isub>=nfiles)
      {
        flag_files = 0;
        printf("Processor %d leaving...\n",rank);
        fflush(stdout);
      }
    }else{
      MPI_Recv(&sb[0], 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, world, &status);
      sb[0] = isub;
      MPI_Send(&sb[0], 2, MPI_INT, sb[1], sb[1], world);

      //increment snapshot subregion
      isub++;

      if(isub>=nfiles+np-1)
      {
        flag_files = 0;
        printf("Processor %d leaving...\n",rank);
        fflush(stdout);
      }
    } // rank!=0
  } //while_flags

  //sum number of peaks across all processes
  MPI_Allreduce(&n_peak_total,&n_peak_all,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&ntd_total,&ntd_all,1,MPI_LONG,MPI_SUM,world);


  if(rank==0)
  {
    printf("Total number of peaks     = %d.\n",n_peak_all);
    printf("Total number of particles = %ld.\n",ntd_all);
    fflush(stdout);
  }

  MPI_Barrier(world);
  t_end = MPI_Wtime();

  MPI_Reduce(&t_assign,&t_assign_max,1,MPI_DOUBLE,MPI_MAX,0,world);
  MPI_Reduce(&t_io,&t_io_max,1,MPI_DOUBLE,MPI_MAX,0,world);

  if(rank==0)
  {
    printf("Total shock identification time = %fs.\n",t_end-t_start);
    printf("Max peak assignment time        = %fs.\n",t_assign_max);
    printf("Max file i/o time               = %fs.\n",t_io_max);

    fflush(stdout);
  }

  //return total number of peaks found by all processes
  //across all subvolumes
  return n_peak_all;
}
