#ifndef  LOAD_TRACERS_H
#define  LOAD_TRACERS_H

using namespace std;

struct tracer
{
  long id;
  float d;
  float x[3];
  float v[3];
};

long load_tracers(char fdir[], char filebase[], char fsuffix[], char fdir_out[], int isnap, int isub, vector<tracer> *t, float dthresh);

#endif /*LOAD_TRACERS_H*/
