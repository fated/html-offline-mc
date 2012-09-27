
/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mira.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.11 $ */
/* $State: Exp $ */

#include <time.h>
#include "mira.h"
#include "utilities.h"
#include "mucutils.h"
#include "redopt.h"

#define ROUND_PER 10

static long class_errors=0;
static double *kernel_values;

static MCSolution mc_sol;
static RedOptDef redopt_def; /*initizlize*/
static KF kernel_fun;
static REDOPT_FUN redopt_fun;
static MCStateSparse mcsp;
static MCDataDef mc_datadef;
static MIRAParamDef mira_pd;

long mira_example(long i);
void mira_update_state(long i, long mistake_k);
void mira_construct(long n_iteration);
void mira_destruct();


MCSolution mira(MCDataDef _mc_datadef, MIRAParamDef _mira_pd) {
  long mistake_k;
  long iteration=0;
  long n_iteration=0;
  long i; /* current example */
  long round;
  
  double mira_time0;/* , mira_time1; */

  mira_time0 = (double)(clock())/(double)CLOCKS_PER_SEC;

  printf("\nLearner (MIRA)  version 1.0\n");
  printf("Initializing ... start\n"); fflush(stdout);

  mc_datadef = _mc_datadef;
  mira_pd = _mira_pd;
  n_iteration = (long) (mira_pd.epochs * (double)mc_datadef.m);
  round = n_iteration / ROUND_PER;
  printf("Requested margin (beta) %f\n", mira_pd.beta);
  printf("Total iterations %ld ( %f epochs x sample size %ld )\n", n_iteration, mira_pd.epochs, mc_datadef.m);
  printf("Initializing ... done\n"); 
  mira_construct(n_iteration);
  printf("\n");
  printf("Iteration    Example    Class errors    No. SPS\n");
  printf("---------    -------    ------------    -------\n");

  fflush(stdout);
  
  for (iteration=0, i=0; iteration<n_iteration; iteration++, i=(++i)%mc_datadef.m) {
    mistake_k = mira_example(i);
    mira_update_state(i, mistake_k);
    
/*      mira_pd.beta *= .95; */
    
    if ((iteration % round) == (round - 1)) {
      printf("%8ld%10ld%14ld%13ld\n", iteration+1, i, class_errors, mcsp.n_supp_pattern);
      fflush(stdout);
    }
  }

  /*    printf("%8ld%10ld%14ld%13ld\n", iteration, i, class_errors, mcsp.n_supp_pattern); */
  /*   mira_time1 = (double)(clock())/(double)CLOCKS_PER_SEC; */
  
  /*   printf("\nMIRA Run time: %20.6f(sec)\n", mira_time1-mira_time0); */
  
  /* update matrix_f */
  
  mc_sol.size              = mc_datadef.m;
  mc_sol.k                 = mcsp.k;
  mc_sol.n_supp_pattern    = mcsp.n_supp_pattern;
  mc_sol.is_voted          = mcsp.is_voted;
  mc_sol.supp_pattern_list = mcsp.supp_pattern_list;
  mc_sol.votes_weight      = mcsp.votes_weight;
  mc_sol.tau               = mcsp.tau;

  mira_destruct();
  return (mc_sol);  
}


void mira_construct(long n_iteration) {
  kernel_construct(mc_datadef.l);
  kernel_fun = kernel_get_function(mira_pd.kernel_def);
  kernel_values = ut_calloc(n_iteration, sizeof(double));

  redopt_construct(mc_datadef.k);
  redopt_fun = redopt_get_function(mira_pd.redopt_type);
  redopt_def.delta = mira_pd.delta;
  redopt_def.b = (double*) ut_calloc(mc_datadef.k, sizeof(double));
  
  mc_statesparse_construct(&mcsp, n_iteration, mc_datadef.k, mira_pd.is_voted);
}


void mira_destruct() {
  redopt_destruct();
  kernel_destruct();
  free(kernel_values);
  free(redopt_def.b);
  free(*mcsp.matrix_f);
  free(mcsp.matrix_f);
  mc_statesparse_clear(&mcsp);
}


long mira_example(long i) {
  long r;
  long sps_i;
  long mistake_k;
  
  for (sps_i=0; sps_i<mcsp.n_supp_pattern; sps_i++)
    kernel_values[sps_i] = kernel_fun(mc_datadef.x[i], mc_datadef.x[mcsp.supp_pattern_list[sps_i]]);
  
  if (mcsp.is_voted == 1) {
    for (r=0; r<mc_datadef.k; r++) {
      redopt_def.b[r]=0;
      for (sps_i=0; sps_i<mcsp.n_supp_pattern; sps_i++)
	redopt_def.b[r] += mcsp.votes_weight[sps_i] * mcsp.tau[sps_i][r] * kernel_values[sps_i] ;
    }
  }
  else {
    for (r=0; r<mc_datadef.k; r++) {
      redopt_def.b[r]=0;
      for (sps_i=0; sps_i<mcsp.n_supp_pattern; sps_i++)
	redopt_def.b[r] += mcsp.tau[sps_i][r] * kernel_values[sps_i] ;
    }
  }
  
/*    if (iteration > 0)  assumtion : if iteration == 0 redopt_def.b[r]=0 */
/*      for (r=0; r<mc_datadef.k; r++) */
/*        redopt_def.b[r] /= iteration; */
  

  redopt_def.a = (kernel_fun)(mc_datadef.x[i], mc_datadef.x[i]);
  redopt_def.y = mc_datadef.y[i];
  redopt_def.alpha = mcsp.tau[mcsp.n_supp_pattern]; /* candidate for next pattern */
  if (redopt_margin_error(&redopt_def, 0) > 0)
    class_errors ++;
  redopt_def.b[mc_datadef.y[i]] -= mira_pd.beta;
  mistake_k = redopt_fun(&redopt_def); 

  return (mistake_k);
}


void mira_update_state(long i, long mistake_k) {

  if (mistake_k > 0) { 
    kernel_values[mcsp.n_supp_pattern] = redopt_def.a;
    mcsp.supp_pattern_list[mcsp.n_supp_pattern] = i;
    if (mcsp.is_voted == 1)
      mcsp.votes_weight[mcsp.n_supp_pattern]=1;
    mcsp.n_supp_pattern++;
  }
  else {
    if (mcsp.is_voted == 1)
      mcsp.votes_weight[mcsp.n_supp_pattern-1]++;
  } 

} 
