/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/kernel.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.6 $ */
/* $State: Exp $ */



#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kernel.h"
#include "utilities.h"


#define N_BUILT_IN_POLY_KERNELS 9

/* innder product */
#define KERNEL_INNER(x,y) { \
register long   tx =0; \
register long   ty =0; \
 b=0; \
 tx = x->ind; \
 ty = y->ind; \
 while (tx>=0 && ty>=0) { \
  if (tx == ty) { \
   b += x->val * y->val; \
   x++; tx = x->ind; \
   y++; ty = y->ind; \
  } \
  else if (tx > ty) {\
   y++; ty = y->ind; \
  } \
  else { \
   x++; tx = x->ind; \
  } \
 } \
}

#define KERNEL_INNER_ID(x) { \
register long   tx =0; \
 b=0; \
 tx = x->ind; \
 while (tx>=0) { \
  b += x->val * x->val; \
  x++; tx = x->ind; \
 } \
}

static long l;
static KernelDef kernel_def;

/* -------------------------------------------------- */

double kernel_exponent(MCDataUnit *x, MCDataUnit *y);
double kernel_exponent_np(MCDataUnit *x, MCDataUnit *y);
double kernel_polynom_homo(MCDataUnit *x, MCDataUnit *y);
double kernel_polynom_non_homo(MCDataUnit *x, MCDataUnit *y);
double kernel_polynom_non_homo_np(MCDataUnit *x, MCDataUnit *y);

double kernel_inner(MCDataUnit *x, MCDataUnit *y); 

/* ========================================================== */


void kernel_construct(long _l) {
#ifdef DEBUG_KERNEL
  fprintf(stdout,"kernel_construct\n");
#endif
  l=_l;
}

KF kernel_get_function(KernelDef _kernel_def) {

#ifdef DEBUG_KERNEL
  fprintf(stdout,"kernel_get_function\n");
#endif
  kernel_def = _kernel_def;
  
  switch (kernel_def.kernel_type) {
  case KERNEL_EXPONENT:
    return (kernel_exponent);
    break;
  case KERNEL_EXPONENT_NP:
    return (kernel_exponent_np);
    break;
  case KERNEL_POLYNOMIAL_HOMO:
    if (kernel_def.polynom_degree < 1) {
      fprintf(stdout,"kernel_get_function :  Error - polynom degree must be positive\n");
      return (NULL);
    }
    return (kernel_polynom_homo);
    break;

  case KERNEL_POLYNOMIAL_NON_HOMO_NP: 
    if (kernel_def.polynom_degree < 1) {
      fprintf(stdout,"kernel_get_function :  Error - polynom degree must be positive\n");
      return (NULL);
    }
    return (kernel_polynom_non_homo_np);
    break;
  
  case KERNEL_POLYNOMIAL_NON_HOMO:
    if (kernel_def.polynom_degree < 1) {
      fprintf(stdout,"kernel_get_function :  Error - polynom degree must be positive\n");
      return (NULL);
    }
        
    return (kernel_polynom_non_homo);
    break;
  default:
    fprintf(stdout,"kernel_get_function :  Error - no such kernel type\n");
      
  }
  return (NULL);
}

long kernel_destruct() {
#ifdef DEBUG_KERNEL
  fprintf(stdout,"kernel_destruct\n");
#endif
  return (0);
}



/* ================================================================================== */
/* kernels */
/* --------*/


/* KERNEL_EXPONENT */
double kernel_exponent(register MCDataUnit *x, register MCDataUnit *y) {

  register double a, c=0;
  register double b=0;
  MCDataUnit *temp_x = x; 
  MCDataUnit *temp_y = y;

#ifdef DEBUG_KERNEL1
  fprintf(stdout,"kernel_exponent : sigma=%f\n",kernel_def.exponent_sigma);
#endif

  KERNEL_INNER_ID(x);
  c += b;
  KERNEL_INNER_ID(y); 
  c += b;
  
  x = temp_x;
  y = temp_y;
  KERNEL_INNER(x,y);
  c -= 2*b;
  a = 0.5*(c/(kernel_def.exponent_sigma * kernel_def.exponent_sigma));

  return exp(-a);
}


/* KERNEL_EXPONENT_NP */
double kernel_exponent_np(register MCDataUnit *x, register MCDataUnit *y) {
  register double c=0;
  register double b=0;
  double d=0;
  MCDataUnit *temp_x = x; 
  MCDataUnit *temp_y = y;

  KERNEL_INNER_ID(x);
  c += b;

  KERNEL_INNER_ID(y); 
  c += b;

  x = temp_x;
  y = temp_y;
  KERNEL_INNER(x,y);
  c -= 2*b;
  d = exp(-c);
/*    fprintf(stderr, "%lf %lf\n", c, d); */

  return (d);
}


/* KERNEL_POLYNOMIAL_HOMO */
double kernel_polynom_homo(register MCDataUnit *x, register MCDataUnit *y) {
  register double a=1;
  register double b=0;
  long p = kernel_def.polynom_degree;

#ifdef DEBUG_KERNEL1
  fprintf(stdout,"kernel_polynom_homo : degree=%ld\n", kernel_def.polynom_degree);
#endif

  KERNEL_INNER(x,y);

  while (p) {
    while (!(p&1)) {
      p>>=1;
      b = b*b;
    }
    p>>=1;
    a=a*b;
    b=b*b;
  }
  return (a);
}


/* KERNEL_POLYNOMIAL_HOMO_NON */
double kernel_polynom_non_homo(register MCDataUnit *x, register MCDataUnit *y) {
  register double b = 0, a = 1;
  long p = kernel_def.polynom_degree;

#ifdef DEBUG_KERNEL1
  fprintf(stdout,"kernel_polynom_homo_non : degree=%ld\ta0=%f\n", kernel_def.polynom_degree, kernel_def.polynom_a0);
#endif

  KERNEL_INNER(x,y);

  b += kernel_def.polynom_a0;
  while (p) {
    while (!(p&1)) {
      p>>=1;
      b = b*b;
    }
    p>>=1;
    a=a*b;
    b=b*b;
  }

  return (a);
}
  
  
/* KERNEL_POLYNOMIAL_HOMO_NON_NP */
double kernel_polynom_non_homo_np(register MCDataUnit *x, register MCDataUnit *y) {
  register double b = 0, a = 1;
  long p = kernel_def.polynom_degree;

#ifdef DEBUG_KERNEL1
  fprintf(stdout,"kernel_polynom_homo_non_np : degree=%ld\n", kernel_def.polynom_degree);
#endif
  KERNEL_INNER(x,y);
  b++;

  while (p) {
    while (!(p&1)) {
      p>>=1;
      b = b*b;
    }
    p>>=1;
    a=a*b;
    b=b*b;
  }

  return (a);
}


double kernel_inner(MCDataUnit *x, MCDataUnit *y) {

  register double b  =0;
  register long   tx =0;
  register long   ty =0;

  tx = x->ind;
  ty = y->ind;
  while (tx>=0 && ty>=0) {
   
    if (tx == ty) {
      b += x->val * y->val;
      x++; tx = x->ind;
      y++; ty = y->ind;
    }
    else if (tx > ty) {
      y++; ty = y->ind;
    }
    else {
      x++; tx = x->ind;
    }
  }
  return (b);
}


/* ------------------------------------------------------------- */


char*  kernel_get_type_name(enum KernelType kt) {

  static char *kernel_type_names[] = {
    "ILLEGAL",
    "kernel_exponent (0)",
    "kernel_exponent_np (1)",
    "kernel_polynomial_homo (2)",
    "kernel_polynomial_non_homo (3)",
    "kernel_polynomial_non_homo_np (4)"
  };

  return ((kt>4) ? kernel_type_names[0] : kernel_type_names[kt+1]); 
}


void   kernel_read(KernelDef *kernel_def, FILE *file){
  ut_fread(&kernel_def->kernel_type,              sizeof(int), 1, file);
  ut_fread(&kernel_def->polynom_degree,           sizeof(int), 1, file);
  ut_fread(&kernel_def->polynom_a0,               sizeof(double), 1, file);
  ut_fread(&kernel_def->exponent_sigma,           sizeof(double), 1, file);
}

void   kernel_write(const KernelDef *kernel_def, FILE *file) { 
  ut_fwrite(&kernel_def->kernel_type,    sizeof(int), 1, file);
  ut_fwrite(&kernel_def->polynom_degree, sizeof(int), 1, file);
  ut_fwrite(&kernel_def->polynom_a0,     sizeof(double), 1, file);
  ut_fwrite(&kernel_def->exponent_sigma, sizeof(double), 1, file);
}

void   kernel_text(const KernelDef *kernel_def, FILE *file) {
  fprintf(file, "Kernel : %s \n( polynom_degree= %ld, polynom_a0= %.10f, exponent_sigma= %.10f)\n" ,
	 kernel_get_type_name(kernel_def->kernel_type),
	 kernel_def->polynom_degree,
	 kernel_def->polynom_a0,
	 kernel_def->exponent_sigma);  fflush(stdout);

}
