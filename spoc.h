/* $Author: kobics $ */
/* $Date: 2003/08/06 11:19:40 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/spoc.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.2 $ */
/* $State: Exp $ */



#ifndef __SPOC_H
#define __SPOC_H

#include "kernel.h"
#include "redopt.h"
#include "mucutils.h"


typedef struct {
  double beta;
  long   cache_size; /* Mb */
  double epsilon;
  double epsilon0;
  double delta;
  KernelDef kernel_def;
  enum RedOptType redopt_type;

} SPOCParamDef;


MCSolution spoc(MCDataDef _mc_datadef, SPOCParamDef _spoc_pd); 
void spoc_destruct();

#endif
