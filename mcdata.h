/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mcdata.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.3 $ */
/* $State: Exp $ */



#ifndef __MCDATA_H
#define __MCDATA_H

#include <stdio.h>

typedef struct {
  double val;
  long   ind;
} MCDataUnit;


typedef struct {
  MCDataUnit **x;
  long *y;
  long m;
  long k;
  long l;
  long max_l;
} MCDataDef;

void mc_datadef_read(MCDataDef *mc_datadef, FILE *file);
void mc_datadef_destruct(MCDataDef *mc_datadef);
void mc_datadef_initialize(MCDataDef *mc_datadef);
void mc_datadef_text(const MCDataDef *mc_datadef, FILE *file);


#endif
