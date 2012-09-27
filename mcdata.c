/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mcdata.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.5 $ */
/* $State: Exp $ */



#include <stdio.h>
#include "mcdata.h"
#include "utilities.h"

long scan_for_max_l(FILE *file, long m, long l);

void mc_datadef_read(MCDataDef *mc_datadef, FILE *file) {
  long i, t;
  double temp;
  long cur_l;
  
  mc_datadef->max_l = scan_for_max_l(file, mc_datadef->m, mc_datadef->l)+1;
  rewind(file);
  
  mc_datadef->x = (MCDataUnit **) ut_calloc(mc_datadef->m, sizeof(MCDataUnit*));
  *(mc_datadef->x) = (MCDataUnit *) ut_calloc(mc_datadef->m * mc_datadef->max_l, sizeof(MCDataUnit));
  
  for (i=1; i<mc_datadef->m; i++) {
    mc_datadef->x[i] = mc_datadef->x[i-1] + mc_datadef->max_l;
  }

  mc_datadef->y = (long *) ut_calloc(mc_datadef->m, sizeof(long));
     
  for (i=0; i<mc_datadef->m; i++) {
    fscanf(file, "%ld", mc_datadef->y+i);
    cur_l=0;
    for (t=0; t<mc_datadef->l; t++) {
      fscanf(file, "%lf", &temp);
      if (temp != 0) {
	mc_datadef->x[i][cur_l].ind = t;
	mc_datadef->x[i][cur_l].val = temp;
	cur_l++;
      }
    }
    mc_datadef->x[i][cur_l].ind = -1;
    mc_datadef->x[i][cur_l].val = 0;
  }
}

void mc_datadef_destruct(MCDataDef *mc_datadef) {
  if (*mc_datadef->x != NULL)
    free(*mc_datadef->x);
  if (mc_datadef->x != NULL)
    free(mc_datadef->x);
  if (mc_datadef->y != NULL)
    free(mc_datadef->y);
  mc_datadef->x = NULL;
  mc_datadef->y = NULL;
}

void mc_datadef_initialize(MCDataDef *mc_datadef) { 
  mc_datadef->m = 0;
  mc_datadef->k = 0;
  mc_datadef->l = 0;
  mc_datadef->max_l = 0;
  mc_datadef->x = NULL;
  mc_datadef->y = NULL;
}


void mc_datadef_text(const MCDataDef *mc_datadef, FILE *file) {
  
  long i, t;
  
  for (i=0; i<mc_datadef->m; i++) {
    fprintf(file, "%ld ", mc_datadef->y[i]);
    for (t=0; mc_datadef->x[i][t].ind>=0; t++) {
      fprintf(file, "%ld:%f ", mc_datadef->x[i][t].ind, mc_datadef->x[i][t].val);
    } 
    fprintf(file, "%ld:%f ", mc_datadef->x[i][t].ind, mc_datadef->x[i][t].val);
    fprintf(file, "\n");
  }
}

long scan_for_max_l(FILE *file, long m, long l) {
  long i,t;
  long max_l =0;
  long cur_l;
  double temp;
  
  for (i=0; i<m; i++) {
    cur_l=0;
    fscanf(file, "%lf", &temp);
    for (t=0; t<l; t++) {
      fscanf(file, "%lf", &temp);
      if (temp !=0)
	cur_l++;
    }
    if (cur_l > max_l)
      max_l = cur_l;
  }
  return (max_l);

}


/* --------------------------------------------------------------------------------------- */
