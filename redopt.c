/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/redopt.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 5.3 $ */
/* $State: Exp $ */



#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "redopt.h"
#include "utilities.h"

/* solve reduced problem with k=2 */
#define REDOPT_TWO(x0,x1,i0,i1) { \
double temp = 0.5*(x0-x1); \
temp = (temp < 1) ? temp : 1; \
(rod->alpha)[i0] = -temp; \
(rod->alpha)[i1] =  temp; }


struct IndexedDouble {
  double value;
  long    index;
};

int    indexed_double_cmp(const void *id1, const void *id2);
long   redopt_exact(const RedOptDef *rod);
long   redopt_approx(const RedOptDef *rod);
long   redopt_exact_long(const RedOptDef *rod);
long   redopt_analytic_binary(const RedOptDef *rod);

void  dump_redopt_def(const RedOptDef *rod);

/* --------------------------------------------------------------------------- */

/* no. of classes */
static long k = 0;
static struct IndexedDouble *vector_d = NULL;

/* --------------------------------------------------------------------------- */

long redopt_construct(long _k) {
  k=_k;
  vector_d = (struct IndexedDouble *) ut_calloc(k, sizeof(struct IndexedDouble));
  return (1);
}


REDOPT_FUN  redopt_get_function(enum RedOptType red_opt_type) {
  //printf("Reduced optimization algorithm : ");
  switch (red_opt_type) {
  case  REDOPT_TYPE_EXACT:
    //printf("exact\n"); fflush(stdout);
    return (redopt_exact);
    break;
  case REDOPT_TYPE_APPROX:
    //printf("approximate\n"); fflush(stdout);
    return (redopt_approx);
    break;
  case REDOPT_TYPE_ANALYTIC_BINARY:
    //printf("binary\n"); fflush(stdout);
    return (redopt_analytic_binary);
    break;
  default:
    break;
  }
  return (NULL);
}


/* solve reduced exactly, use sort */
long  redopt_exact(const RedOptDef *rod) {
  
  double phi0=0;       /* potenial functions phi(t)   */
  double phi1;       /* potenial functions phi(t+1) */
  double sum_d =0;
  double theta;      /* threshold */
  long mistake_k =0;  /* no. of labels with score greater than the correct label */ 
  long r;
  long r1;

#ifdef DEBUG_REDOPT1
  dump_redopt_def(rod);
#endif
  
  /* pick only problematic labels */
  for (r=0; r<k; r++) {
    if ((rod->b)[r] > (rod->b)[rod->y]) {
      vector_d[mistake_k].index = r;
      vector_d[mistake_k].value = ((rod->b[r]) / rod->a);
      sum_d += vector_d[mistake_k].value; 
      mistake_k++;
    }
    /* for other labels, alpha=0 */
    else {
      (rod->alpha)[r]=0;
    }
  }

  /* if no mistake labels return */
  if (mistake_k == 0) {
    return (0);
  }
  /* add correct label to list (up to constant one) */
  vector_d[mistake_k].index = rod->y;
  vector_d[mistake_k].value = ((rod->b[rod->y]) / rod->a);
  
  /* if there are only two bad labels, solve for it */
  if (mistake_k == 1) {
    REDOPT_TWO(vector_d[0].value, vector_d[1].value, vector_d[0].index, vector_d[1].index);
    return (2);
  }

  /* finish calculations */
  sum_d += vector_d[mistake_k].value;
  vector_d[mistake_k].value++;
  mistake_k++;
  
  /* sort vector_d by values */
  qsort(vector_d, mistake_k, sizeof(struct IndexedDouble), indexed_double_cmp);
  
  /* go down the potential until sign reversal */
  for (r=1, phi1=1; phi1>0 && r<mistake_k; r++) {
    phi0 = phi1;
    phi1 = phi0 - r * (vector_d[r-1].value - vector_d[r].value);
  }

  /* theta < min vector_d.value */
  /* nu[r] = theta */
  if (phi1 > 0) {
    sum_d /= mistake_k;
    for (r=0; r<mistake_k; r++) {
      (rod->alpha)[vector_d[r].index] = sum_d - vector_d[r].value;
    }
    (rod->alpha)[rod->y]++;
  
  }
  /* theta > min vector_d.value */
  else {
    theta = - phi0 / (--r);
    theta += vector_d[--r].value;
    /* update tau[r] with nu[r]=theta */
    for (r1=0; r1 <= r; r1++) {
      (rod->alpha)[vector_d[r1].index] = theta - vector_d[r1].value;
    }
    /* update tau[r]=0, nu[r]=vector[d].r */
    for (; r1 < mistake_k; r1++) {
      (rod->alpha)[vector_d[r1].index] = 0;
    }
    (rod->alpha)[rod->y]++;
  }

  
#ifdef DEBUG_REDOPT1
  dump_redopt_def(rod);
  fprintf(stderr,"resopt_exact : mistake_k = %d\n\n",mistake_k);
#endif
    
  return (mistake_k);
}

long  redopt_approx(const RedOptDef *rod) {

  double old_theta = DBL_MAX;  /* threshold */
  double theta = DBL_MAX;      /* threshold */
  double temp;
  long mistake_k =0; /* no. of labels with score greater than the correct label */ 
  long r;

  /* pick only problematic labels */
  for (r=0; r<k; r++) {
    if ((rod->b)[r] > (rod->b)[rod->y]) {
      vector_d[mistake_k].index = r;
      vector_d[mistake_k].value = ((rod->b[r]) / rod->a);
      mistake_k++;
    }
    /* for other labels, alpha=0 */
    else {
      (rod->alpha)[r]=0;
    }
  }

  /* if no mistake labels return */
  if (mistake_k == 0) {
    return (0);
  }
  
  /* add correct label to list (up to constant one) */
  vector_d[mistake_k].index = rod->y;
  vector_d[mistake_k].value = ((rod->b[rod->y]) / rod->a);

  /* if there are only two bad labels, solve for it */
  if (mistake_k == 1) {
    REDOPT_TWO(vector_d[0].value, vector_d[1].value, vector_d[0].index, vector_d[1].index);
    return (2);
	       
  }

  /* finish calculations */
  vector_d[mistake_k].value++;
  mistake_k++;
  
  /* initialize theta to be min D_r */
  for (r=0; r<mistake_k; r++) {
    if (vector_d[r].value < theta)
      theta = vector_d[r].value;
  }

  /* loop until convergence of theta */
  while (1) {
    old_theta = theta;
    
    /* calculate new value of theta */
    theta = -1;
    for (r=0; r<mistake_k; r++) {
      if (old_theta > vector_d[r].value)
	theta += old_theta;
      else
	theta += vector_d[r].value;
    }
    theta /=mistake_k;
    
    if (FABS((old_theta - theta)/theta) < rod->delta)
      break;
  }
  
  /* update alpha using threshold */
  for (r=0; r<mistake_k; r++) {
    temp = theta - vector_d[r].value;
    if (temp < 0)
      (rod->alpha)[vector_d[r].index]=temp;
    else
      (rod->alpha)[vector_d[r].index]=0;
  }
  (rod->alpha)[rod->y]++;
  
  return(mistake_k);
}

/* same as redopt_exact except : */
/* no seeking for bad labels, use all labels (miatake_k = k) */
/* no speial case of mistake_k = 2 */
/* here for validation and time comparsion */

long  redopt_exact_long(const RedOptDef *rod) {

  double phi0=0, phi1;
  double theta;
  long r, r1;
  double sum_d;

  sum_d = 0;
  for (r=0; r<k; r++) { 
    vector_d[r].value = ((rod->b[r]) / rod->a);
    vector_d[r].index = r;
    sum_d += vector_d[r].value;
  }
  vector_d[rod->y].value += 1;

  qsort(vector_d, k, sizeof(struct IndexedDouble), indexed_double_cmp);
  
  for (r=1, phi1=1; phi1>0 && r<k; r++) {
    phi0 = phi1;
    phi1 = phi0 - r * (vector_d[r-1].value - vector_d[r].value);
  }

  if (phi1 > 0) {
    sum_d /= k;
    for (r=0; r<k; r++) {
      (rod->alpha)[vector_d[r].index] = sum_d - vector_d[r].value;
    }
    (rod->alpha)[rod->y]++;
  }
  else {
    theta = - phi0 / --r;
    theta += vector_d[--r].value;
    for (r1=0; r1<=r; r1++) {
      (rod->alpha)[vector_d[r1].index] = theta - vector_d[r1].value;
    }
    for (; r1<k; r1++) {
      (rod->alpha)[vector_d[r1].index] =  0;
    }
    (rod->alpha)[rod->y]++;
  }

  return(0);
}

/* solve for k=2 */
long redopt_analytic_binary(const RedOptDef *rod) {

  long y0 = 1- rod->y; /* other label */
  long y1 = rod->y;    /* currect label */
  
  if ( (rod->b)[y0] > (rod->b)[y1] ) {
    vector_d[y0].value = ((rod->b[y0]) / rod->a);
    vector_d[y1].value = ((rod->b[y1]) / rod->a);

    REDOPT_TWO(vector_d[y0].value, vector_d[y1].value, y0, y1);
    return (2);
  }
  else {
    (rod->alpha)[0] = (rod->alpha[1]) = 0;
    return (0);
  }
}


void redopt_destruct() {
  free(vector_d);
}


int indexed_double_cmp(const void *id1, const void *id2) {
  if (((struct IndexedDouble *)id1)->value > ((struct IndexedDouble *)id2)->value)
    return (-1);
  else if (((struct IndexedDouble *)id1)->value < ((struct IndexedDouble *)id2)->value)
    return (1);
  else
    return (0);
}




void  dump_redopt_def(const RedOptDef *rod) {
  long r;

  fprintf(stderr, "dump_redopt_def\n");
  fprintf(stderr, "a= %10.4f\n",rod->a);
  fprintf(stderr, "b=\n");
  for (r=0; r<k; r++)
    fprintf(stderr, "%10.4f",rod->b[r]);
  fprintf(stderr, "\n");
  fprintf(stderr, "y= %ld\n",rod->y);
  fprintf(stderr, "alpha=\n");
  for (r=0; r<k; r++)
    fprintf(stderr, "%10.4f",rod->alpha[r]);
  fprintf(stderr, "\n");
  fprintf(stderr, "delta= %10.4f\n",rod->delta);
  
}


char*       redopt_get_type_name(enum RedOptType rot) {
  static char *redopt_type_names[] = {
    "ILLEGAL",
    "exact",
    "approc",
    "binary"
  };

  return ( (rot<1 || rot>3) ? redopt_type_names[0] : redopt_type_names[rot]); 
}



long redopt_margin_error(const RedOptDef *rod, const double beta) {
  long errors =0;
  long r;
    
  for (r=0; r<rod->y; r++) {
    if (rod->b[r] >= rod->b[rod->y]-beta) {
      errors++;
    }
  }
   
  for (r++; r<k; r++) {
    if (rod->b[r] >= rod->b[rod->y]-beta) {
      errors++;
    }
  }
  return (errors);
}
