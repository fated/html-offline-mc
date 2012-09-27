/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mira.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 5.1 $ */
/* $State: Exp $ */



#ifndef __MIRA_H
#define __MIRA_H

#include "kernel.h"
#include "redopt.h"
#include "mucutils.h"

typedef struct {
  double beta;
  double epochs;
  long   is_voted;
  KernelDef kernel_def;
  enum RedOptType redopt_type;
  double delta;
} MIRAParamDef;

MCSolution mira(MCDataDef mc_datadef, MIRAParamDef mira_pd);

#endif
