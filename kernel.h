/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/kernel.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.1 $ */
/* $State: Exp $ */




#ifndef __KERNEL_H
#define __KERNEL_H

#include <stdio.h>
#include "mcdata.h"

enum KernelType 
{ KERNEL_EXPONENT=0, 
  KERNEL_EXPONENT_NP, 
  KERNEL_POLYNOMIAL_HOMO,
  KERNEL_POLYNOMIAL_NON_HOMO, 
  KERNEL_POLYNOMIAL_NON_HOMO_NP };

typedef double (*KF)(MCDataUnit*, MCDataUnit*);

typedef struct {
  enum KernelType kernel_type;
  long polynom_degree;
  double polynom_a0;
  double exponent_sigma;
} KernelDef;

void   kernel_construct(long _l);
KF     kernel_get_function(KernelDef _kernel_def);
long   kernel_destruct();
char*  kernel_get_type_name(enum KernelType kt);
void   kernel_read(KernelDef *kernel_def, FILE *file);
void   kernel_write(const KernelDef *kernel_def, FILE *file); 
void   kernel_text(const KernelDef *kernel_def, FILE *file); 

#endif
