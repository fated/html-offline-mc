
/* $Author: kobics $ */
/* $Date: 2003/06/26 07:21:18 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mconline.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.6 $ */
/* $State: Exp $ */

#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#include "mconline.h"
#include "utilities.h"
#include "mucutils.h"
#include "redopt.h"

#define ROUND_PER 20

static long n_inserts=0;
static long class_errors=0;
static double *kernel_values;
static double *kernel_values_i;

static MCSolution mc_sol;
static RedOptDef redopt_def; /*initizlize*/
static KF kernel_fun;
static REDOPT_FUN redopt_fun;
static MCStateFull mcfl;
static MCDataDef mc_datadef;
static MCOnlineParamDef mconline_pd;

static double *delta_tau           = NULL;
static double *delta_tau1           = NULL;


static long     spp_pattern_size =0;
static long*    spp_pattern_indicators =NULL;
static double** spp_pattern_scores =NULL;
static double*  spp_pattern_weights =NULL;
static double*  spp_pattern_norm =NULL;
static double   norm2 = 0;
static double   sum_tau = 0;
static double*  example_score=NULL;
static double*  example_score1=NULL;
static double*  example_norms=NULL;

static long current_example; 

typedef struct {
  long index_abs;
  long index_spp;
  double val;
  char sgn;
} ExampleFind;



int   longcmp(const void *n1, const void *n2);
long  mconline_example();
void  mconline_update_stage(long mistake_k);
void  mconline_construct(long n_iteration, MCClassifier* mc_cls_init);
void  mconline_destruct();
char* mconline_update_type_text(const enum MCOnlineUpdateType ut);
char* mconline_stage_type_text(const enum MCOnlineStageType st);
char* mconline_find_type_text(const enum MCOnlineFindType ft);
double compute_margin(double* scores, long y);

ExampleFind find_example();
void  example_update(const long i, const double* diff_tau);


MCSolution mconline(MCDataDef _mc_datadef, MCOnlineParamDef _mconline_pd, MCClassifier* mc_cls_init) {
  long mistake_k;
  long iteration=0;
  long n_iteration=0;


  long round;
  
  double mconline_time0;/* , mconline_time1; */

  mconline_time0 = (double)(clock())/(double)CLOCKS_PER_SEC;

  printf("\nLearner (MCONLINE)  version 1.1\n");
  printf("Initializing ... start\n"); fflush(stdout);
  
  mc_datadef = _mc_datadef;
  mconline_pd = _mconline_pd;
  n_iteration = (long) (mconline_pd.epochs * (double)mc_datadef.m);
  round = n_iteration / ROUND_PER;

  mconline_construct(n_iteration, mc_cls_init);

  printf("Requested margin (beta) %f\n", mconline_pd.beta);
  printf("Total iterations %ld ( %f epochs x sample size %ld )\n", n_iteration, mconline_pd.epochs, mc_datadef.m);
  printf("Example set size %ld\n", mconline_pd.spp_pattern_size);
  printf("Stage  type %s\n", mconline_stage_type_text(mconline_pd.stage_type));
  printf("Find   type %s\n", mconline_find_type_text(mconline_pd.find_type));
  printf("Update type %s\n", mconline_update_type_text(mconline_pd.update_type));
  printf("Through example parameter (gamma) %f\n", mconline_pd.gamma);
  printf("Initializing ... done\n"); 

  printf("\n");
  printf("Iteration    Example    Class errors    No. SPS    No. Inserts\n");
  printf("---------    -------    ------------    -------    -----------\n");

  fflush(stdout);
  
  for (iteration=0, current_example=0; 
       iteration<n_iteration; 
       iteration++, current_example=(++current_example)%mc_datadef.m) {
    
    /* check if error +  get new weights */    
    mistake_k = mconline_example();     
    
    
    /* if need to replace - replace */
    mconline_update_stage(mistake_k);
    
    if ((iteration % round) == (round - 1)) {
      printf("%8ld%10ld%14ld%13ld%13ld\n", iteration+1, current_example, class_errors, mcfl.n_supp_pattern,n_inserts);
      fflush(stdout);
    }
  }
  
  /*   printf("Cumulative training errors : %ld\n", class_errors); */
  /*   fflush(stdout); */
  
  
  {
    long i,r;
    qsort(mcfl.supp_pattern_list, mcfl.n_supp_pattern, sizeof(long), &longcmp);
    for (i=0; i<mcfl.n_supp_pattern; i++)
      for (r=0; r<mcfl.k; r++)
	mcfl.tau[i][r] = mcfl.tau[mcfl.supp_pattern_list[i]][r];
    for (i=mcfl.n_supp_pattern; i<mc_datadef.m; i++)
      for (r=0; r<mcfl.k; r++)
	mcfl.tau[i][r]=0;
  }
  
  
  mc_sol.size              = mc_datadef.m;
  mc_sol.k                 = mcfl.k;
  mc_sol.n_supp_pattern    = mcfl.n_supp_pattern;
  mc_sol.is_voted          = 0;
  mc_sol.supp_pattern_list = mcfl.supp_pattern_list;
  mc_sol.votes_weight      = NULL;
  mc_sol.tau               = mcfl.tau;
  

  mconline_destruct();
  
  return (mc_sol);  
}


void mconline_construct(long n_iteration, MCClassifier* mc_cls_init) {
  long i, r;
  
  if (mc_cls_init == NULL) {
    fprintf(stdout, "Zero initialization\n");fflush(stdout);
  }
  else {
    fprintf(stdout, "Using init classifier\n");fflush(stdout);
  }
  kernel_construct(mc_datadef.l);
  kernel_fun = kernel_get_function(mconline_pd.kernel_def);
  kernel_values = ut_calloc(mc_datadef.m, sizeof(double));
  kernel_values_i = ut_calloc(mc_datadef.m, sizeof(double));

  redopt_construct(mc_datadef.k);
  redopt_fun = redopt_get_function(mconline_pd.redopt_type);
  redopt_def.delta = mconline_pd.delta;
  redopt_def.b = (double*) ut_calloc(mc_datadef.k, sizeof(double));
  
  delta_tau = (double *) ut_calloc(mc_datadef.k, sizeof(double));
  delta_tau1 = (double *) ut_calloc(mc_datadef.k, sizeof(double));

  spp_pattern_size = mconline_pd.spp_pattern_size;

  spp_pattern_indicators = (long *) ut_calloc(mc_datadef.m, sizeof(long));
  
  spp_pattern_scores = (double **) ut_calloc(mc_datadef.m, sizeof(double*));
  *spp_pattern_scores = (double *) ut_calloc(mc_datadef.m*mc_datadef.k, sizeof(double));
  for (i=1; i<mc_datadef.m; i++) 
    spp_pattern_scores[i] = spp_pattern_scores[i-1] + mc_datadef.k;
  
  spp_pattern_weights = (double *) ut_calloc(mc_datadef.m, sizeof(double));

  spp_pattern_norm = (double *) ut_calloc(mc_datadef.m, sizeof(double));

  example_score = (double*) ut_calloc(mc_datadef.k, sizeof(double));
  example_score1 = (double*) ut_calloc(mc_datadef.k, sizeof(double));
  example_norms = (double*) ut_calloc(mc_datadef.m, sizeof(double));



  mc_statefull_construct(&mcfl, n_iteration, mc_datadef.k);

  for (i=0; i<mc_datadef.m; i++)
    spp_pattern_indicators[i]=-1;
  for (i=0; i<mc_datadef.m; i++)
    for (r=0; r<mc_datadef.k; r++)
      spp_pattern_scores[i][r] =0;
  for (i=0; i<mc_datadef.m; i++)
    spp_pattern_weights[i] =0;
  for (i=0; i<mc_datadef.m; i++)
    spp_pattern_norm[i] =0;
  for (i=0; i<mc_datadef.m; i++)
    example_norms[i] =0;

  
  if (mc_cls_init != NULL) {
    fprintf(stderr, "Can't init a classifer now\n");
    exit(0);


    mcfl.n_supp_pattern = mc_cls_init->solution.n_supp_pattern;

    for (i=0; i<mcfl.n_supp_pattern; i++) {
      mcfl.supp_pattern_list[i] = mc_cls_init->solution.supp_pattern_list[i]; 
    
      for (r=0; r<mcfl.k; r++)
	mcfl.tau[mcfl.supp_pattern_list[i]][r] = mc_cls_init->solution.tau[i][r];
      
      spp_pattern_indicators[mcfl.supp_pattern_list[i]] = i;
      spp_pattern_weights[mcfl.supp_pattern_list[i]] = 
	2 * mcfl.tau[mcfl.supp_pattern_list[i]][mc_datadef.y[mcfl.supp_pattern_list[i]]];
    }

  }
}



void mconline_destruct() {


  redopt_destruct();
  kernel_destruct();
  free(kernel_values);
  free(kernel_values_i);
  free(redopt_def.b);
  free(*mcfl.matrix_f);
  free(mcfl.matrix_f);
  mc_statefull_clear(&mcfl);

  if (delta_tau != NULL)
    free(delta_tau);

  if (delta_tau1 != NULL)
    free(delta_tau1);

  if (spp_pattern_indicators != NULL)
    free(spp_pattern_indicators);

  if (*spp_pattern_scores != NULL)
    free(*spp_pattern_scores);
  if (spp_pattern_scores != NULL)
    free(spp_pattern_scores);

  if (spp_pattern_weights != NULL)
    free(spp_pattern_weights);
  if (spp_pattern_norm != NULL)
    free(spp_pattern_norm);
  if (example_score1 != NULL)
    free(example_score1);

  if (example_norms != NULL)
    free(example_norms);
}


long mconline_example() {
  long r;
  long sps_i;
  long mistake_k;
  
  for (sps_i=0; sps_i<mcfl.n_supp_pattern; sps_i++)
    kernel_values[mcfl.supp_pattern_list[sps_i]] = 
      kernel_fun(mc_datadef.x[current_example], mc_datadef.x[mcfl.supp_pattern_list[sps_i]]);
  kernel_values[current_example] = kernel_fun(mc_datadef.x[current_example],mc_datadef.x[current_example]);
  example_norms[current_example] = kernel_values[current_example];
  
  for (r=0; r<mc_datadef.k; r++) {
    redopt_def.b[r]=0;
    for (sps_i=0; sps_i<mcfl.n_supp_pattern; sps_i++)
      redopt_def.b[r] += 
	mcfl.tau[mcfl.supp_pattern_list[sps_i]][r] * kernel_values[mcfl.supp_pattern_list[sps_i]] ;
  }
  
  redopt_def.a = kernel_values[current_example];
  redopt_def.y = mc_datadef.y[current_example];
  redopt_def.alpha = delta_tau; 
  if (redopt_margin_error(&redopt_def, 0) > 0)
    class_errors ++;
  redopt_def.b[mc_datadef.y[current_example]] -= mconline_pd.beta;
  mistake_k = redopt_fun(&redopt_def); 
  redopt_def.b[mc_datadef.y[current_example]] += mconline_pd.beta;
  return (mistake_k);
}


char* mconline_update_type_text(const enum MCOnlineUpdateType ut) {
  static char *mconline_update_type_names[] = {
    "ILLEGAL",
    "update_uniform (0)",
    "update_max (1)",
    "update_prop (2)",
    "update_mira (3)",
    "update_perceptron (4)",
    "update_ovr_mira (5)",
    "update_mira (6)",
    "update_rand (7)"
  };

  return ((ut>7) ? mconline_update_type_names[0] : mconline_update_type_names[ut+1]); 
}


char* mconline_stage_type_text(const enum MCOnlineStageType ut) {
  static char *mconline_stage_type_names[] = {
    "ILLEGAL",
    "bound-find-update (0)",
    "update-bound-find (1)",
    "update-find-all   (2)",
    "update-find-one   (3)"
  };

  return ((ut>3) ? mconline_stage_type_names[0] : mconline_stage_type_names[ut+1]); 
}



char* mconline_find_type_text(const enum MCOnlineFindType ut) {
  static char *mconline_find_type_names[] = {
    "ILLEGAL",
    "maximal_margin (0)",
    "minimal_weight (1)",
    "maximal_margin_wo_example (CC-I) (2)",
    "maximal_diff_norm (CC-II) (3)",
    "maximal_diff_norm_normalized (CC-II) (4)"
  };

  return ((ut>4) ? mconline_find_type_names[0] : mconline_find_type_names[ut+1]); 
}


int   longcmp(const void *n1, const void *n2) {
  if (*(int*)(n1)>*(int*)(n2))
    return 1;
  else if (*(int*)(n1)<*(int*)(n2))
    return -1;
  else return 0;
}



ExampleFind find_example() {

  long j, r, abs_i;
  double temp;
 
  ExampleFind example_find;
  example_find.index_abs = -1;
  example_find.index_spp = -1;



  switch (mconline_pd.find_type) {
  case FIND_MAXIMAL_MARGIN:
    example_find.val = -DBL_MAX;
    for (j=0; j<mcfl.n_supp_pattern; j++) {
      temp = compute_margin(spp_pattern_scores[mcfl.supp_pattern_list[j]], 
			    mc_datadef.y[mcfl.supp_pattern_list[j]]);
      
      if (temp > example_find.val) {
	example_find.index_spp = j;
	example_find.index_abs = mcfl.supp_pattern_list[j];
	example_find.val = temp;
      }
    }
    example_find.sgn = +1;
    break;
    
  case FIND_MIMIMAL_WEIGHT:
    example_find.val = DBL_MAX;
    for (j=0; j<mcfl.n_supp_pattern; j++) {
      if (spp_pattern_weights[mcfl.supp_pattern_list[j]]<example_find.val) {
	example_find.index_spp = j;
	example_find.index_abs = mcfl.supp_pattern_list[j];
	example_find.val = spp_pattern_weights[example_find.index_abs];
      }
    }
    example_find.sgn = -1;
    break;
    
  case FIND_MAXIMAL_MARGIN_WO_EXAMPLE: /* CleanCache 1 */
    example_find.val = -DBL_MAX;
    for (j=0; j<mcfl.n_supp_pattern; j++) {
      for (r=0; r<mc_datadef.k; r++)
	example_score1[r] = spp_pattern_scores[mcfl.supp_pattern_list[j]][r] 
	  - mcfl.tau[mcfl.supp_pattern_list[j]][r] * example_norms[mcfl.supp_pattern_list[j]];
      temp = compute_margin(example_score1,mc_datadef.y[mcfl.supp_pattern_list[j]]);
      
      if (temp > example_find.val) {
	example_find.index_spp = j;
	example_find.index_abs = mcfl.supp_pattern_list[j];
	example_find.val = temp;
      }
    }
    example_find.sgn = +1;
    break;
	
	
  case FIND_MAXIMAL_DIFF_NORM: /* CleanCache 2 */
  case FIND_MAXIMAL_DIFF_NORM_NORMALIZED: /* CleanCache 2 */

    example_find.val = DBL_MAX;
    if (mcfl.n_supp_pattern > 1) {
      for (j=0; j<mcfl.n_supp_pattern; j++) {
	abs_i = mcfl.supp_pattern_list[j];
	temp = spp_pattern_norm[abs_i]/ 
	  (sum_tau - mcfl.tau[abs_i][mc_datadef.y[abs_i]]);
	
	if (temp < example_find.val) {
	  example_find.index_spp = j;
	  example_find.index_abs = mcfl.supp_pattern_list[j];
	  example_find.val = temp;
	}
      }
    }
    example_find.sgn = +1;
    example_find.val = norm2/sum_tau - example_find.val;
    if (mconline_pd.find_type == FIND_MAXIMAL_DIFF_NORM)
      example_find.val *= (sum_tau*sum_tau);
    fprintf(stderr, "%20.10f\n", example_find.val); 
    break;
      
  }

  return (example_find);
}






void example_update(const long update_example, const double* diff_tau) {
  /* Assume kernel_val is not need to be used in this round */
  double* kv;
  long sps_i, abs_i, r;
  
  if (update_example != current_example) {
    kv = kernel_values_i;
    for (sps_i=0; sps_i<mcfl.n_supp_pattern; sps_i++)
      kv[mcfl.supp_pattern_list[sps_i]] = 
	kernel_fun(mc_datadef.x[update_example], mc_datadef.x[mcfl.supp_pattern_list[sps_i]]);
    kernel_values[update_example] = kernel_fun(mc_datadef.x[update_example],mc_datadef.x[update_example]);
  }
  else {
    kv = kernel_values;
  }
  
  if (spp_pattern_indicators[update_example] == -1) {
    /* a new spp, always there is a place for it !! */
    mcfl.supp_pattern_list[mcfl.n_supp_pattern] = update_example;
    spp_pattern_indicators[update_example] = mcfl.n_supp_pattern;
    spp_pattern_norm[update_example] = norm2;
    for (r=0; r<mc_datadef.k; r++)
      spp_pattern_scores[update_example][r] = redopt_def.b[r];
    mcfl.n_supp_pattern++;
    n_inserts++;

  }

  /* update weights */
  spp_pattern_weights[update_example] += diff_tau[mc_datadef.y[update_example]];

  /* compute scores for example i */
  for (r=0; r<mc_datadef.k; r++) {
    example_score[r]=0;
    for (sps_i=0; sps_i<mcfl.n_supp_pattern; sps_i++)
      example_score[r] += 
	mcfl.tau[mcfl.supp_pattern_list[sps_i]][r] * kv[mcfl.supp_pattern_list[sps_i]];
  }
  
  /* update norms of all examples */
  for (sps_i=0; sps_i<mcfl.n_supp_pattern; sps_i++) {
    abs_i = mcfl.supp_pattern_list[sps_i];
    for (r=0; r<mc_datadef.k; r++)
      example_score1[r] = example_score[r] ;
    if (update_example != abs_i) {
      for (r=0; r<mc_datadef.k; r++)
	example_score1[r] -= mcfl.tau[abs_i][r] * kv[abs_i];
      
      for (r=0; r<mc_datadef.k; r++)
	spp_pattern_norm[abs_i] +=  
	  2*diff_tau[r]*example_score1[r] + 
	  diff_tau[r]*diff_tau[r]*kv[update_example];
    }
  }
  
  
  /* update scores */
  for (sps_i=0;  sps_i<mcfl.n_supp_pattern; sps_i++) {
    abs_i = mcfl.supp_pattern_list[sps_i];
    for (r=0; r<mc_datadef.k; r++) {
      spp_pattern_scores[abs_i][r] += diff_tau[r] * kv[abs_i];
    }
  }
  
  
  /* update norm2 */
  for (r=0; r<mc_datadef.k; r++)
    norm2 += 2 * diff_tau[r] * example_score[r] + diff_tau[r]*diff_tau[r]*kv[update_example];


  /* update tau */
  for (r=0; r<mc_datadef.k; r++)
    mcfl.tau[update_example][r] += diff_tau[r];
  sum_tau += diff_tau[mc_datadef.y[update_example]];



  /* remove example if needed */
  if (mcfl.tau[update_example][mc_datadef.y[update_example]]<ACTUAL_ZERO) {
    sps_i = spp_pattern_indicators[update_example];
    
    for (r=0; r<mcfl.k; r++) 
      mcfl.tau[update_example][r] =0;

    spp_pattern_weights[update_example] = 0;
    for (r=0; r<mcfl.k; r++) 
      spp_pattern_scores[update_example][r] = 0;
    spp_pattern_norm[update_example] = 0;

    mcfl.n_supp_pattern--;
    mcfl.supp_pattern_list[sps_i] = mcfl.supp_pattern_list[mcfl.n_supp_pattern];
    spp_pattern_indicators[mcfl.supp_pattern_list[mcfl.n_supp_pattern]] = sps_i;
    spp_pattern_indicators[update_example] = -1;
  }

}


void mconline_update_stage(long mistake_k) {
  
  long r;
  ExampleFind example_find;
  
  switch (mconline_pd.stage_type) {
  case STAGE_BOUND_FIND_UPDATE:
    if (mistake_k > 0) {   
      if (spp_pattern_indicators[current_example]==-1 
	  && mcfl.n_supp_pattern  >= spp_pattern_size) {
	/* need to throw */
	example_find = find_example();
	if (example_find.index_abs>=0) {
	  for (r=0; r<mc_datadef.k; r++)
	    delta_tau1[r] = -mcfl.tau[example_find.index_abs][r];
	  example_update(example_find.index_abs, delta_tau1);
	}
      }
      example_update(current_example, delta_tau);
    }
    break;
    
  case STAGE_UPDATE_BOUND_FIND:
    if (mistake_k > 0) {   
      example_update(current_example, delta_tau);
      if (mcfl.n_supp_pattern > spp_pattern_size) {
	/* need to throw */
	example_find = find_example();
	if (example_find.index_abs>=0) {
	  for (r=0; r<mc_datadef.k; r++)
	    delta_tau1[r] = -mcfl.tau[example_find.index_abs][r];
	  example_update(example_find.index_abs, delta_tau1);
	}
      }
    }
    break;
    
  case STAGE_UPDATE_FIND_ALL:
    if (mistake_k > 0) {   
      example_update(current_example, delta_tau);
      while (1) {
	example_find = find_example();
	if (example_find.index_abs<0) 
	  break;
	if (example_find.val * example_find.sgn < 
	    mconline_pd.gamma * example_find.sgn) 
	  break;
	
	for (r=0; r<mc_datadef.k; r++)
	  delta_tau1[r] = -mcfl.tau[example_find.index_abs][r];
	example_update(example_find.index_abs, delta_tau1);
      }
    }
    
    break;

  case STAGE_UPDATE_FIND_ONE:
    if (mistake_k > 0) {   
      example_update(current_example, delta_tau);
      example_find = find_example();
      if (example_find.index_abs<0) 
	break;
      if (example_find.val * example_find.sgn < 
	  mconline_pd.gamma * example_find.sgn) 
	break;
      
      for (r=0; r<mc_datadef.k; r++)
	delta_tau1[r] = -mcfl.tau[example_find.index_abs][r];
      example_update(example_find.index_abs, delta_tau1);
    }
    break;
    
  }


}


double compute_margin(double* scores, long y) {
  
  long r;
  double margin=-DBL_MAX;
    
  for (r=0; r<y; r++) 
    if (scores[r] > margin)
      margin = scores[r];
   
  for (r++; r<mc_datadef.k; r++) 
    if (scores[r] > margin)
      margin = scores[r];
  
  margin = scores[y] - margin;

  return (margin);
}

