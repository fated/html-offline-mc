/* $Author: kobics $ */
/* $Date: 2003/06/02 08:35:48 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/redopt.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 5.2 $ */
/* $State: Exp $ */



#ifndef __REDOPT__H
#define __REDOPT__H

typedef struct {
  double a;
  double *b;
  long    y;
  double *alpha;
  double delta;
  /*  
  To move later beta decrease into redopt
  double beta; 
  */
} RedOptDef;


enum RedOptType 
{ REDOPT_TYPE_EXACT=0,
  REDOPT_TYPE_APPROX,
  REDOPT_TYPE_ANALYTIC_BINARY};


typedef long (*REDOPT_FUN)(const RedOptDef*);


long        redopt_construct(long _k);
REDOPT_FUN  redopt_get_function(enum RedOptType red_opt_type);
void        redopt_destruct();
char*       redopt_get_type_name(enum RedOptType rot);

long        redopt_margin_error(const RedOptDef *rod, const double beta);
#endif
