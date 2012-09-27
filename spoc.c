/* $Author: kobics $ */
/* $Date: 2003/08/06 11:19:40 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/spoc.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.6 $ */
/* $State: Exp $ */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>

#include "spoc.h"
#include "cachelru.h"
#include "redopt.h"
#include "utilities.h"

//#define DEBUG_SPOC 1
//#define DEBUG_SPOC1 1

/* pointers */
static MCDataUnit **x   = NULL;
static long *y       = NULL;
static long m        = 0;
static long k        = 0;
static long l        = 0;
static double **tau = NULL;
static double beta  = 0;

static SPOCParamDef spoc_pd;
static MCDataDef mc_datadef;
static MCSolution mc_sol;
static RedOptDef redopt_def;
static REDOPT_FUN redopt_fun;
static KF kernel_function;

static double *vector_a            = NULL;
static double **matrix_f           = NULL;
static double *row_matrix_f_next_p = NULL;
static long   **matrix_eye         = NULL;
static double *delta_tau           = NULL;
static double *old_tau             = NULL;
static double *vector_b            = NULL; 

static double  max_psi             =0;
static long    next_p;
static long    next_p_list;

static long    *supp_pattern_list = NULL;
static long    *zero_pattern_list = NULL;
static long    n_supp_pattern     = 0;
static long    n_zero_pattern     = 0;


int    longcmp(const void *n1, const void *n2);
void   spoc_construct();
long   spoc_epsilon(double epsilon);
void   choose_next_pattern(long *pattern_list, long n_pattern);
void   update_matrix_f(double *kernel_pattern_p);
long   allocate_memory();
void   free_auxilary_memory();
void   free_external_memory();
void   spoc_initialize();
void   update_matrix_f(double *kernel_next_p);
double next_epsilon1(double epsilon_cur, double epsilon);
double next_epsilon2(double epsilon_cur, double epsilon);
long   get_no_supp1();
double get_train_error(double b);
void   dump_matrix_f();

int   intcmp(const void *n1, const void *n2) {
  return (*(int*)(n1) - *(int*)(n2));
}
  


/* -------------------------------------------------------- */


MCSolution spoc(MCDataDef _mc_datadef, SPOCParamDef _spoc_pd) {

  double epsilon_current =1;
  clock_t spoc_time0, spoc_time1;

  mc_datadef = _mc_datadef;
  spoc_pd = _spoc_pd;
  spoc_construct();
  
  spoc_time0 = clock();
  
#ifdef MONITOR
  fprintf(spoc_pd.monitor_file, "run def\n");
  fprintf(spoc_pd.monitor_file, "\tepsilon :%10.5e\n",spoc_pd.epsilon);
  fprintf(spoc_pd.monitor_file, "\tepsilon0 :%10.5e\n",spoc_pd.epsilon0);
  fprintf(spoc_pd.monitor_file, "\tdelta :%10.5e\n",spoc_pd.delta);
  
  fprintf(spoc_pd.monitor_file, "\n*\n");

  fprintf(spoc_pd.monitor_file, "epsilon\n");
  fprintf(spoc_pd.monitor_file, "iterations\n");
  fprintf(spoc_pd.monitor_file, "iterations_opt\n");
  fprintf(spoc_pd.monitor_file, "cum_supp_only\n");
  fprintf(spoc_pd.monitor_file, "cum_mistake_k\n");
  fprintf(spoc_pd.monitor_file, "cum_new_supp\n");
  fprintf(spoc_pd.monitor_file, "cum_cache_miss\n");
  fprintf(spoc_pd.monitor_file, "cum_delta_q\n");
  fprintf(spoc_pd.monitor_file, "cum_max_psi\n");
  fprintf(spoc_pd.monitor_file, "\n*\n");

#endif

  //printf("Epsilon decreasing from %e to %e\n", spoc_pd.epsilon0, spoc_pd.epsilon);
  redopt_construct(k);
  redopt_fun = redopt_get_function(spoc_pd.redopt_type);
  //if (spoc_pd.redopt_type == REDOPT_TYPE_APPROX)
    //printf("Delta %e\n", spoc_pd.delta);

  redopt_def.delta = spoc_pd.delta;
  redopt_def.b = (double*) ut_calloc(k, sizeof(double));

  epsilon_current = spoc_pd.epsilon0;
  
  //printf("\n"); 
  //printf("New Epsilon   No. SPS      Max Psi   Train Error   Margin Error\n"); 
  //printf("-----------   -------      -------   -----------   ------------\n"); 
  //fflush(stdout);
  
  while (max_psi > spoc_pd.epsilon * beta) {
    
    //printf("%11.5e   %7ld   %10.3e   %7.2f%%      %7.2f%%\n",
	//   epsilon_current, n_supp_pattern, max_psi/beta, get_train_error(beta), get_train_error(0));fflush(stdout);
    spoc_epsilon(epsilon_current);
    epsilon_current  = next_epsilon2(epsilon_current , spoc_pd.epsilon);
  }   
  //printf("%11.5e   %7ld   %10.3e   %7.2f%%      %7.2f%%\n",
  //	 spoc_pd.epsilon, n_supp_pattern, max_psi/beta, get_train_error(beta), get_train_error(0));fflush(stdout);


  free (redopt_def.b); 

  

  if ( 0 ) {
    long i;
    int* a = (int*) ut_calloc((size_t)n_supp_pattern, sizeof(int));
    for (i=0; i<n_supp_pattern; i++)
      a[i]=(int)supp_pattern_list[i];
    
    
    qsort(a, n_supp_pattern, sizeof(int), &intcmp); 
    
    
    fprintf(stderr, "%ld\n", n_supp_pattern);
    for (i=0; i<n_supp_pattern; i++)
      fprintf(stderr, "%ld\t", supp_pattern_list[i]);
    fprintf(stderr, "done");
  }

  {
    long i,r;

    qsort(supp_pattern_list, n_supp_pattern, sizeof(long), &longcmp); 

    for (i=0; i<n_supp_pattern; i++)
      for (r=0; r<k; r++)
	tau[i][r] = tau[supp_pattern_list[i]][r];
    for (i=n_supp_pattern; i<m; i++)
      for (r=0; r<k; r++)
	tau[i][r]=0;
    
  }


  mc_sol.size              = m;
  mc_sol.k                 = k;
  mc_sol.l                 = l;
  mc_sol.n_supp_pattern    = n_supp_pattern;
  mc_sol.is_voted          = 0;
  mc_sol.supp_pattern_list = supp_pattern_list;
  mc_sol.votes_weight      = NULL;
  mc_sol.tau               = tau;
  


  spoc_time1 = clock();
  
  //printf("\nNo. support pattern %ld ( %ld at bound )\n", n_supp_pattern, get_no_supp1());

  /*     dump_matrix_f(); */

  return (mc_sol);
}



long spoc_epsilon(double epsilon) {

  long supp_only =1;
  long cont = 1;
  long mistake_k;
  double *kernel_next_p;
  long r;
  long i;

#ifdef MONITOR
  long   iterations     =0;
  long   iterations_opt =0;
  long   cum_supp_only  =0;
  long   cum_mistake_k  =0;
  long   cum_new_supp   =0;
  long   cum_cache_miss =0;
  double cum_delta_q    =0;
  double cum_max_psi    =0;
#endif

  
  while (cont) {
#ifdef DEBUG_SPOC1
    fprintf(stderr,"spoc : %10.4f %10ld %10.4f %10ld\n",epsilon, n_supp_pattern, max_psi / beta, supp_only);
#endif
    max_psi = 0;
    if (supp_only) 
      choose_next_pattern(supp_pattern_list, n_supp_pattern);
    else
      choose_next_pattern(zero_pattern_list, n_zero_pattern);
    
#ifdef MONITOR1
    fprintf(spoc_pd.file_next_p, "%d\t", next_p);
#endif
    
#ifdef MONITOR
    iterations ++;
    cum_supp_only += supp_only;
    cum_max_psi   += max_psi;
#endif
    
    if (max_psi > epsilon * beta) {
#ifdef DEBUG_SPOC1
    fprintf(stderr,"max_psi > epsilon * beta\n");
#endif	
      redopt_def.a = vector_a[next_p];
      for (r=0; r<k; r++)
	redopt_def.b[r] = matrix_f[next_p][r] - redopt_def.a * tau[next_p][r];
      redopt_def.y = y[next_p];

      for (r=0; r<k; r++)
	old_tau[r] = tau[next_p][r];
      redopt_def.alpha = tau[next_p];

      mistake_k = (*redopt_fun)(&redopt_def);

      for (r=0; r<k; r++)
	delta_tau[r] = tau[next_p][r]- old_tau[r];
      
      if (!cachelru_retrive(next_p, (void*)(&kernel_next_p))) {
	for (i=0; i<m; i++)
	  kernel_next_p[i] = kernel_function(x[i], x[next_p]);

#ifdef MONITOR
	cum_cache_miss++;
#endif
      }
#ifdef MONITOR
      iterations_opt ++;
      cum_mistake_k += mistake_k;
      {
	double delta_q = 0;
	for (r=0; r<k; r++)
	  delta_q += delta_tau[r] * delta_tau[r];
	delta_q *= (0.5 * vector_a[next_p]);

	for (r=0; r<k; r++)
	  delta_q += delta_tau[r] * matrix_f[next_p][r];

	cum_delta_q += delta_q;

	fprintf(spoc_pd.q_monitor_file, "\t%20.10e",  delta_q);

      }
#endif
            
      update_matrix_f(kernel_next_p);
      
      if (supp_only) {
	for (r=0; r<k; r++)
	  if (tau[next_p][r] != 0) break;
	if (r == k) {
	  zero_pattern_list[n_zero_pattern++] = next_p;
	  supp_pattern_list[next_p_list] = supp_pattern_list[--n_supp_pattern];
	}
      }
      else {
	supp_pattern_list[n_supp_pattern++] = next_p;
	zero_pattern_list[next_p_list] = zero_pattern_list[--n_zero_pattern]; 
#ifdef MONITOR
	cum_new_supp++;
#endif
	supp_only =1;
      }

    }
    else {
      if (supp_only)
	supp_only =0;
      else
	cont =0;
    }
  }

#ifdef MONITOR
  fprintf(spoc_pd.monitor_file, "\t%20.10e",  epsilon);
  fprintf(spoc_pd.monitor_file, "\t%ld", iterations);
  fprintf(spoc_pd.monitor_file, "\t%ld", iterations_opt);
  fprintf(spoc_pd.monitor_file, "\t%ld", cum_supp_only);
  fprintf(spoc_pd.monitor_file, "\t%ld", cum_mistake_k);
  fprintf(spoc_pd.monitor_file, "\t%ld", cum_new_supp);
  fprintf(spoc_pd.monitor_file, "\t%ld", cum_cache_miss);
  fprintf(spoc_pd.monitor_file, "\t%20.10e",  cum_delta_q);
  fprintf(spoc_pd.monitor_file, "\t%20.10e",  cum_max_psi);
  fprintf(spoc_pd.monitor_file, "\n");
#endif

  return (1);
}


void spoc_construct() {
  
  long req_n_blocks_lrucache;
  long n_blocks_lrucache, tmp;
  
  //printf("\nOptimizer (SPOC)  (version 1.0)\n");
  //printf("Initializing ... start \n"); fflush(stdout);

  x = mc_datadef.x;
  y = mc_datadef.y;
  m = mc_datadef.m;
  k = mc_datadef.k;
  l = mc_datadef.l;

  /* scale beta with m */
  //fflush(stdout);
  beta = spoc_pd.beta;

  //printf("Requested margin (beta) %e\n", spoc_pd.beta);
  kernel_function = kernel_get_function(spoc_pd.kernel_def);

  allocate_memory();
  redopt_construct(k);
  kernel_construct(l);
  
  n_supp_pattern = 0;
  n_zero_pattern = 0;
  
  req_n_blocks_lrucache = MIN(m, ((long)(((double)spoc_pd.cache_size) / ((((double)sizeof(double)) * ((double)m)) / ((double)MB))))); 
  //printf("Requesting %ld blocks of cache (External bound %ldMb)\n",
//	 req_n_blocks_lrucache, spoc_pd.cache_size);
  n_blocks_lrucache = cachelru_construct(m, req_n_blocks_lrucache, sizeof(double)*m);
  
#ifdef DEBUG_SPOC
  fprintf(stderr,"requested cache size = \t%ld\n",req_n_blocks_lrucache);
  fprintf(stderr,"final cache size = \t%ld\n",n_blocks_lrucache);
#endif

  spoc_initialize();

#ifdef MONITOR

  {
    fprintf(spoc_pd.monitor_file, "problem def\n");
    fprintf(spoc_pd.monitor_file, "\tm :%10d\n",m);
    fprintf(spoc_pd.monitor_file, "\tk :%10d\n",k);
    fprintf(spoc_pd.monitor_file, "\tl :%10d\n",l);
    
    
    fprintf(spoc_pd.monitor_file, "spoc_pd\n");
    fprintf(spoc_pd.monitor_file, "\tbeta :%10.5e\n",spoc_pd.beta);
    fprintf(spoc_pd.monitor_file, "\tcache_size :%ld\n",spoc_pd.cache_size);

    fprintf(spoc_pd.monitor_file, "kernel\n");
    fprintf(spoc_pd.monitor_file, "\tkernel_type :\t%s\n",kernel_get_type_name(spoc_pd.kernel_def.kernel_type));
    fprintf(spoc_pd.monitor_file, "\tpolynom_degree :\t%d\n",spoc_pd.kernel_def.polynom_degree);
    fprintf(spoc_pd.monitor_file, "\tpolynom_a0 :\t%f\n",spoc_pd.kernel_def.polynom_a0);
    fprintf(spoc_pd.monitor_file, "\texponent_sigma :\t%f\n",spoc_pd.kernel_def.exponent_sigma);

    fprintf(spoc_pd.monitor_file, "cachelru\n");
    fprintf(spoc_pd.monitor_file, "\treq_n_blocks_lrucache :\t%ld\n",req_n_blocks_lrucache);
    fprintf(spoc_pd.monitor_file, "\tn_blocks_lrucache :\t%ld\n",n_blocks_lrucache);
    fprintf(spoc_pd.monitor_file, "\tblocks_size :\t%dB\n",sizeof(double)*m);

  }
#endif
  
  //printf("Initializing ... done\n"); fflush(stdout);
}


void spoc_initialize() {
  
  long i;
  long s,r;

#ifdef DEBUG_SPOC
  fprintf(stderr,"spoc_initialize ....\t");
#endif

  /* vector_a */
  for (i=0; i<m; i++) {
    vector_a[i] = kernel_function(x[i], x[i]);
#ifdef DEBUG_SPOC1
    fprintf(stderr,"vector_a[i] = %10.4f\n", vector_a[i]);
#endif
    
  }

  /* matrix_eye */
  for (r=0; r<k; r++)
    for (s=0; s<k; s++)
      if (r != s)
	matrix_eye[r][s] = 0;
      else
	matrix_eye[r][s] = 1;
  
  /* matrix_f */
  for (i=0; i<m; i++) {
    for (r=0; r<k; r++) {
      if (y[i] != r)
	matrix_f[i][r] = 0;
      else
	matrix_f[i][r] = -beta;
    }
  }

  /* tau */
  for (i=0; i<m; i++)
    for (r=0 ;r<k; r++)
      tau[i][r] = 0;

/*  assume tau=0 */
/*  ============ */  
/*    for (j=0; j<m; j++) { */
/*      for (r=0 ;r<k; r++) { */
/*        if (tau[j][r] != 0) { */
/*  	for (i=0; i<m; i++) */
/*  	  matrix_f[i][r] += kernel_function(x[i], x[j]) * tau[j][r]; */
/*        } */
/*      } */
/*    } */
  
/*    for (i=0; i<m; i++) { */
/*      is_tau_r_zero = 1; */
    
/*      for (r=0; r<k; r++) { */
/*        if (tau[i][r] != 0) { */
/*  	is_tau_r_zero = 0; */
/*  	break; */
/*        } */
/*      } */
/*      if (is_tau_r_zero) */
/*        zero_pattern_list[n_zero_pattern++] =i; */
/*      else */
/*        supp_pattern_list[n_supp_pattern++] =i; */
/*    } */

  supp_pattern_list[0] =0;
  n_supp_pattern =1;

  for (i=1; i<m; i++)
    zero_pattern_list[i-1] =i;
  n_zero_pattern = m-1;
  choose_next_pattern(supp_pattern_list, n_supp_pattern);

#ifdef DEBUG_SPOC
  fprintf(stderr,"done\n");
#endif
 
}




void spoc_destruct() {

  free_auxilary_memory();
  //free_external_memory();
  redopt_destruct();
  kernel_destruct();
  cachelru_destruct(); 

}


/* -----------------------------------------------------------  */

long allocate_memory() {
  long i;
  
#ifdef DEBUG_SPOC
  fprintf(stderr,"allocate_memory ....\t");
#endif

  /* tau */
  tau = (double **) ut_calloc(m, sizeof(double*));
  *tau = (double *) ut_calloc(m*k, sizeof(double));
  for (i=1; i<m; i++) 
    tau[i] = tau[i-1] + k;
  
  /* vector_a */
  vector_a = (double *) ut_calloc(m, sizeof(double));

  /* matrix_f */
  matrix_f = (double **) ut_calloc(m, sizeof(double*));
  *matrix_f = (double *) ut_calloc(m*k, sizeof(double));
  for (i=1; i<m; i++) 
    matrix_f[i] = matrix_f[i-1] + k;

  /* matrix_eye */
  matrix_eye = (long **) ut_calloc(k, sizeof(long*));
  *matrix_eye = (long *) ut_calloc(k*k, sizeof(long));
  for (i=1; i<k; i++) 
    matrix_eye[i] = matrix_eye[i-1] + k;
  
  /* delta_tau */
  delta_tau = (double *) ut_calloc(k, sizeof(double));

  /* old_tau */
  old_tau = (double *) ut_calloc(k, sizeof(double));

  /* vector_b */
  vector_b = (double *) ut_calloc(k, sizeof(double));
 
  /* supp_pattern_list */
  supp_pattern_list = (long *) ut_calloc(m, sizeof(long));

  /* zero_pattern_list */
  zero_pattern_list = (long *) ut_calloc(m, sizeof(long));
  
#ifdef DEBUG_SPOC
 fprintf(stderr,"done\n");
#endif
 
 return (1);
 
}



void free_auxilary_memory() {
  
#ifdef DEBUG_SPOC
  fprintf(stderr,"free_auxilary_memory ....\t");
#endif

  /* vector_a */
  if (vector_a != NULL)
    free(vector_a);
      
  /* matrix_f */
  if (matrix_f != NULL) {
    if (*matrix_f != NULL)
      free(*matrix_f);
    free(matrix_f);
  }

  /* matrix_eye */
  if (matrix_eye != NULL) {
    if (matrix_eye != NULL)
      free(*matrix_eye);
    free(matrix_eye);
  }
  
  /* delta_tau */
  if (delta_tau != NULL)
    free(delta_tau);

  /* old_tau */
  if (old_tau != NULL)
    free(old_tau);

  /* vector_b */
  if (vector_b != NULL)
    free(vector_b);

  /* zero_pattern_list */
  if (zero_pattern_list != NULL)
    free(zero_pattern_list);
  
#ifdef DEBUG_SPOC
 fprintf(stderr,"done\n");
#endif

}



void free_external_memory() {
  
#ifdef DEBUG_SPOC
  fprintf(stderr,"free_external_memory ....\t");
#endif
  
  
  /* supp_pattern_list */
  if (supp_pattern_list != NULL)
    free(supp_pattern_list);

  /* tau */
  if (*tau != NULL)
    free(*tau);
  
  if (tau != NULL)
    free(tau);

#ifdef DEBUG_SPOC
  fprintf(stderr,"done\n");
#endif
 
}



void  choose_next_pattern(long *pattern_list, long n_pattern) {

  double psi;    /* KKT value of example */
  double psi1;   /* max_r matrix_f[i][r] */
  double psi0;   /* min_{r, tau[i][r]<delta[yi][r]}  matrix_f[i][r] */
  
  long i;
  long r;
  long p=0;

  double *matrix_f_ptr;

#ifdef DEBUG_SPOC1
  fprintf(stderr,"in choose_pattern_next\n");
#endif
  
  for (i=0; i<n_pattern; i++) {
    psi1 = -DBL_MAX;
    psi0 =  DBL_MAX;
    
    p=pattern_list[i];
    matrix_f_ptr = matrix_f[p];

    for (r=0; r<k; r++, matrix_f_ptr++) {
      if (*matrix_f_ptr > psi1)
	psi1 = *matrix_f_ptr;
      
      if (*matrix_f_ptr < psi0)
	if (tau[p][r] < matrix_eye[y[p]][r])
	  psi0 = *matrix_f_ptr;
    }
    
    psi = psi1 - psi0;
   
    if (psi > max_psi) {
      next_p_list = i;
      max_psi = psi;
    }
  }
  next_p = pattern_list[next_p_list];
  row_matrix_f_next_p = matrix_f[p];

#ifdef DEBUG_SPOC1
  fprintf(stderr,"next_p = %ld %f. done choose_pattern_next\n", next_p, max_psi);
#endif

}


void update_matrix_f(double *kernel_next_p) {

  long i;
  long r;
  
  double *delta_tau_ptr = delta_tau;
  double *kernel_next_p_ptr;
   
  for (r=0; r<k; r++, delta_tau_ptr++)
    if (*delta_tau_ptr != 0)
      for (i=0, kernel_next_p_ptr = kernel_next_p ; i<m; i++, kernel_next_p_ptr++)
	matrix_f[i][r] += (*delta_tau_ptr) * (*kernel_next_p_ptr);
}



double next_epsilon1(double epsilon_cur, double epsilon) {

  double e;
  e= epsilon_cur; /* just for the compiler */
  e=  ((max_psi / beta) * .95);
  
#ifdef DEBUG_SPOC1
  fprintf(stderr,"next_epsilon1 : max(%f , %f)\n", e,  epsilon);
#endif  
  return (MAX(e, epsilon));
}


double next_epsilon2(double epsilon_cur, double epsilon) {
  static double iteration =12;  
  double e = epsilon_cur / log10(iteration);
  
  iteration+=2;

#ifdef DEBUG_SPOC1
  fprintf(stderr,"next_epsilon2 : iteration = %f\n", iteration);
  fprintf(stderr,"next_epsilon2 : max(%f , %f)\n", e ,  epsilon);
#endif  
  return (MAX( e , epsilon));
}





void dump_matrix_f() {
  long i, r;
  printf("\n");
  for (i=0; i<MIN(10,m); i++) {
    for (r=0; r<k; r++) {
      //fprintf(stderr, "%15.9f", matrix_f[i][r]);
    }
    //fprintf(stderr, "\n");
  }
  printf("\n");
}


double get_train_error(double b) {
  long i, r;
  double max;
  double errors = 0;
  for (i=0; i<m; i++) {
    max = -DBL_MAX;
    for (r=0; r<y[i]; r++) {
      if (matrix_f[i][r]>max) {
	max = matrix_f[i][r];
      }
    }
    for (r++; r<k; r++) {
      if (matrix_f[i][r]>max) {
	max = matrix_f[i][r];
      }
    }
    if ((max-b) >= matrix_f[i][y[i]])
      errors += 1;
  }
  return (100*errors/((double)m));
/*   printf("Training Error = %7.2f\n", 100*errors/((double)m)); */
/*    fflush(stdout); */
}

long get_no_supp1() {
  long n =0;
  long i;

  for (i=0; i<m; i++)
    if (tau[i][y[i]] == 1)
      n++;
  return (n);
}


int   longcmp(const void *n1, const void *n2) {
  return (*(long*)(n1) - *(long*)(n2));
  /*  
      if (*(long*)(n1)>*(long*)(n2))
      return 1;
      else if (*(long*)(n1)<*(long*)(n2))
      return -1;
      else return 0;
  */
}

