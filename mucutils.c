/* $Author: kobics $ */
/* $Date: 2003/06/02 08:35:48 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mucutils.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.3 $ */
/* $State: Exp $ */



#include <stdio.h>
#include "mucutils.h"
#include "utilities.h"


void mc_scaledef_read(MCScaleDef *scale_def, FILE *file) {
  ut_fread(scale_def->comment,              sizeof(char),   COMMENT_LEN, file); 
  ut_fread(&scale_def->to_scale,             sizeof(long),    1, file);
  ut_fread(&scale_def->l,                    sizeof(long),    1, file);
  ut_fread(&scale_def->scale_factor,         sizeof(double),  1, file);
  ut_fread(&scale_def->to_zero_data_mean,    sizeof(long),    1, file);
  if (scale_def->to_zero_data_mean == 1) {
    scale_def->data_mean = (double *) ut_calloc(scale_def->l, sizeof(double));
    ut_fread(scale_def->data_mean,           sizeof(double), scale_def->l, file);
  }
}

void mc_scaledef_write(const MCScaleDef *scale_def, FILE *file) {
  ut_fwrite(scale_def->comment,              sizeof(char),   COMMENT_LEN, file); 
  ut_fwrite(&scale_def->to_scale,             sizeof(long),    1, file);
  ut_fwrite(&scale_def->l,                    sizeof(long),    1, file);
  ut_fwrite(&scale_def->scale_factor,         sizeof(double),  1, file);
  ut_fwrite(&scale_def->to_zero_data_mean,    sizeof(long),    1, file);
  if (scale_def->to_zero_data_mean == 1) {
    ut_fwrite(scale_def->data_mean,           sizeof(double), scale_def->l, file);
  }
}


void mc_scaledef_text(const MCScaleDef *scale_def, FILE *file) {
  fprintf(file, "%s\n",  scale_def->comment);
  fprintf(file, "%ld\n", scale_def->to_scale); 
  fprintf(file, "%ld\n", scale_def->l);
  fprintf(file, "%20.10f\n", scale_def->scale_factor);
  fprintf(file, "%ld\n", scale_def->to_zero_data_mean);
  if (scale_def->to_zero_data_mean == 1) {
    long t;
    for (t=0; t<scale_def->l; t++) {
      if (t%5==4) fprintf(file, "\n");
      fprintf(file, "%20.10f ", scale_def->data_mean[t]);
    }
  }
  fprintf(file, "\n");
  fflush(stdout);
}


void mc_scaledef_scale(const MCScaleDef *scale_def, MCDataUnit **data, long m) {

  if (scale_def->to_scale ==1) {

    long i, t;
    
        
    if (scale_def->to_zero_data_mean) {
    printf("Forcing zero mean ... PROBLEM ..."); fflush(stdout);
    exit(1);
      for (i=0; i<m; i++)
	for (t=0; t<scale_def->l; t++) 
	  data[i][t].val -= scale_def->data_mean[t];
      printf("done\n"); fflush(stdout);
    }
  
    printf("Scaling by factor %f ... ", scale_def->scale_factor); fflush(stdout);
    for (i=0; i<m; i++)
      for (t=0; data[i][t].ind>=0; t++) 
	data[i][t].val *= scale_def->scale_factor;
    printf("done\n"); fflush(stdout);
  }
}


void mc_scaledef_initialize(MCScaleDef *scale_def) {
  *scale_def->comment='\0';
  scale_def->to_scale = 0;
  scale_def->l = 0; 
  scale_def->scale_factor = 0;
  scale_def->to_zero_data_mean = 0;
  scale_def->data_mean = NULL;
  
}

void mc_scaledef_destruct(MCScaleDef *scale_def) {
  if (scale_def->data_mean != NULL)
    free(scale_def->data_mean);
}

/* --------------------------------------------------------------------------------------- */

void mc_solution_construct(MCSolution *sol, long size, long k, long l, long n_supp_pattern, long is_voted) {
  long i, r;
  sol->size = size;
  sol->k = k;
  sol->l = l;
  sol->n_supp_pattern = n_supp_pattern;
  sol->is_voted = is_voted;
  
  sol->supp_pattern_list = (long *) ut_calloc(n_supp_pattern, sizeof(long));
  for (i=0; i<n_supp_pattern; i++)
    sol->supp_pattern_list[i] =0;
  
  if (is_voted==1) {
    sol->votes_weight = (long *) ut_calloc(n_supp_pattern, sizeof(long)); 
    for (i=0; i<n_supp_pattern; i++)
      sol->votes_weight[i] =0;
  }
  else
    sol->votes_weight = NULL;

  sol->tau = (double **) ut_calloc(n_supp_pattern, sizeof(double*));
  *(sol->tau) = (double *) ut_calloc(n_supp_pattern*k, sizeof(double));
  for (i=1; i<n_supp_pattern; i++) 
    sol->tau[i] = sol->tau[i-1] + k;

  for (i=0; i<n_supp_pattern; i++)
    for (r=0; r<k; r++)
      sol->tau[i][r]=0;
  
}

void mc_solution_destruct(MCSolution *sol) {
  free(sol->supp_pattern_list);
  sol->supp_pattern_list = NULL;
  if (sol->votes_weight != NULL)
    free(sol->votes_weight);
  sol->votes_weight = NULL;
  free(*sol->tau);
  free(sol->tau);
  sol->tau=NULL;
}


void mc_solution_clear(MCSolution *sol) {
  sol->supp_pattern_list = NULL;
  sol->votes_weight = NULL;
  sol->tau = NULL;
}

void mc_solution_read(MCSolution *sol, FILE *file) {
  long i;
  ut_fread(&sol->size     , sizeof(long), 1, file);
  ut_fread(&sol->k     , sizeof(long), 1, file);
  ut_fread(&sol->l     , sizeof(long), 1, file);
  ut_fread(&sol->n_supp_pattern, sizeof(long), 1, file);
  ut_fread(&sol->is_voted, sizeof(long), 1, file);
  mc_solution_construct(sol, sol->size, sol->k, sol->l, sol->n_supp_pattern, sol->is_voted);
  ut_fread(sol->supp_pattern_list, sizeof(long), sol->n_supp_pattern, file);
  if (sol->is_voted == 1)
    ut_fread(sol->votes_weight, sizeof(long), sol->n_supp_pattern, file);

  for (i=0; i<sol->n_supp_pattern; i++)
    ut_fread(sol->tau[i]   , sizeof(double), sol->k, file);
  ut_fread(sol->comment,              sizeof(char),   COMMENT_LEN, file); 
}


void mc_solution_write(const MCSolution *sol, FILE *file) {
  long i;
  ut_fwrite(&sol->size     , sizeof(long), 1, file);
  ut_fwrite(&sol->k     , sizeof(long), 1, file);
  ut_fwrite(&sol->l     , sizeof(long), 1, file);
  ut_fwrite(&sol->n_supp_pattern, sizeof(long), 1, file);
  ut_fwrite(&sol->is_voted, sizeof(long), 1, file);
  ut_fwrite(sol->supp_pattern_list, sizeof(long), sol->n_supp_pattern, file);
  if (sol->is_voted == 1)
     ut_fwrite(sol->votes_weight, sizeof(long), sol->n_supp_pattern, file);
  for (i=0; i<sol->n_supp_pattern; i++)
    ut_fwrite(sol->tau[i]   , sizeof(double), sol->k, file);
  ut_fwrite(sol->comment,              sizeof(char),   COMMENT_LEN, file);
}

  
void mc_solution_text(const MCSolution *sol, FILE *file) {
  long i, r;
  
  fprintf(file, "%ld\n", sol->size);
  fprintf(file, "%ld\n", sol->k);
  fprintf(file, "%ld\n", sol->l);
  fprintf(file, "%ld\n", sol->n_supp_pattern);
  fprintf(file, "%ld\n", sol->is_voted);
  
  for (i=0; i<sol->n_supp_pattern; i++)
    fprintf(file, "%ld ", sol->supp_pattern_list[i]);
  fprintf(file, "\n");
  
  if (sol->is_voted == 1) {
    for (i=0; i<sol->n_supp_pattern; i++)
      fprintf(file, "%ld ", sol->votes_weight[i]);
    fprintf(file, "\n");
  }
  for (i=0; i<sol->n_supp_pattern; i++) {
    for (r=0; r<sol->k; r++) {
      fprintf(file, "%15.10f", sol->tau[i][r]);
    }
    fprintf(file, "\n");
  }
  
  fprintf(file, "\n%s\n", sol->comment);
  fflush(file);
}


/* --------------------------------------------------------------------------------------- */

void mc_statesparse_construct(MCStateSparse *stsp, long size, long k, long is_voted) {

  long i, r;
  
  stsp->size = size;
  stsp->k = k;
  stsp->n_supp_pattern = 0;
  stsp->is_voted = is_voted;

  stsp->supp_pattern_list = (long *) ut_calloc(size, sizeof(long));
  for (i=0; i<size; i++)
    stsp->supp_pattern_list[i]=0;

  if (is_voted==1) {
    stsp->votes_weight = (long *) ut_calloc(size, sizeof(long)); 
    for (i=0; i<size; i++)
      stsp->votes_weight[i] =0;
  }
  else
    stsp->votes_weight = NULL;
  
  stsp->tau = (double **) ut_calloc(size, sizeof(double*));
  *(stsp->tau) = (double *) ut_calloc(size*k, sizeof(double));
  for (i=1; i<size; i++) 
    stsp->tau[i] = stsp->tau[i-1] + k;
  for (i=0; i<size; i++)
    for (r=0; r<k; r++) 
      stsp->tau[i][r] = 0;

  stsp->matrix_f = (double **) ut_calloc(size, sizeof(double*));
  *(stsp->matrix_f) = (double *) ut_calloc(size*k, sizeof(double));
  for (i=1; i<size; i++) 
    stsp->matrix_f[i] = stsp->matrix_f[i-1] + k;
  for (i=0; i<size; i++)
    for (r=0; r<k; r++) 
      stsp->matrix_f[i][r] = 0;
}


void mc_statesprase_destruct(MCStateSparse *stsp) {

  free(stsp->supp_pattern_list);
  stsp->supp_pattern_list = NULL;
  if (stsp->is_voted)
    free(stsp->votes_weight);
  stsp->votes_weight = NULL;
  free(*stsp->tau);
  free(stsp->tau);
  stsp->tau=NULL;
  free(*stsp->matrix_f);
  free(stsp->matrix_f);
  stsp->matrix_f=NULL;
}




void mc_statesparse_clear(MCStateSparse *stsp) {

  stsp->supp_pattern_list = NULL;
  stsp->votes_weight = NULL;
  stsp->tau=NULL;
  stsp->matrix_f=NULL;
}

void mc_statesparse_read(MCStateSparse *stsp, FILE *file) {
  long i;
  ut_fread(&stsp->size, sizeof(long), 1, file);
  ut_fread(&stsp->k, sizeof(long), 1, file);
  ut_fread(&stsp->n_supp_pattern, sizeof(long), 1, file);
  ut_fread(&stsp->is_voted, sizeof(long), 1, file);
  mc_statesparse_construct(stsp, stsp->size, stsp->k, stsp->is_voted);
  ut_fread(stsp->supp_pattern_list, sizeof(long), stsp->size, file);
  if (stsp->is_voted)
    ut_fread(stsp->votes_weight, sizeof(long), stsp->size, file);
  for (i=0; i<stsp->size; i++)
    ut_fread(stsp->tau[i]   , sizeof(double), stsp->k, file);
  for (i=0; i<stsp->size; i++)
    ut_fread(stsp->matrix_f[i]   , sizeof(double), stsp->k, file);
}


void mc_statesparse_write(const MCStateSparse *stsp, FILE *file) {
 long i;
  ut_fwrite(&stsp->size, sizeof(long), 1, file);
  ut_fwrite(&stsp->k, sizeof(long), 1, file);
  ut_fwrite(&stsp->n_supp_pattern, sizeof(long), 1, file);
  ut_fwrite(&stsp->is_voted, sizeof(long), 1, file);
  ut_fwrite(stsp->supp_pattern_list, sizeof(long), stsp->size, file);
  if (stsp->is_voted)
    ut_fwrite(stsp->votes_weight, sizeof(long), stsp->size, file);
  for (i=0; i<stsp->size; i++)
    ut_fwrite(stsp->tau[i]   , sizeof(double), stsp->k, file);
  for (i=0; i<stsp->size; i++)
    ut_fwrite(stsp->matrix_f[i]   , sizeof(double), stsp->k, file);
}

void mc_statesparse_text(const MCStateSparse *stsp, FILE *file) {

  long i, r;
  
  fprintf(file, "%ld\n", stsp->size);
  fprintf(file, "%ld\n", stsp->k);
  fprintf(file, "%ld\n", stsp->n_supp_pattern);
  fprintf(file, "%ld\n", stsp->is_voted);

  for (i=0; i<stsp->n_supp_pattern; i++)
    fprintf(file, "%ld ", stsp->supp_pattern_list[i]);
  fprintf(file, "\n");
  
  if (stsp->is_voted == 1) {
    for (i=0; i<stsp->size; i++)
      fprintf(file, "%ld ", stsp->votes_weight[i]);
    fprintf(file, "\n");
  }

  for (i=0; i<stsp->size; i++) {
    for (r=0; r<stsp->k; r++) {
      fprintf(file, "%15.10f", stsp->tau[i][r]);
    }
    fprintf(file, "\n");
  }
  
  fprintf(file, "\n");
  for (i=0; i<stsp->size; i++) {
    for (r=0; r<stsp->k; r++) {
      fprintf(file, "%15.10f", stsp->matrix_f[i][r]);
    }
    fprintf(file, "\n");
  }
  
  fprintf(file, "\n%s\n", stsp->comment);
  fflush(stdout);

}

/* --------------------------------------------------------------------------------------- */


void mc_statefull_construct(MCStateFull *stfl, long size, long k){

  long i, r;
  stfl->size = size;
  stfl->k = k;
  stfl->n_supp_pattern =0;
  stfl->n_zero_pattern =0;

  stfl->supp_pattern_list = (long *) ut_calloc(size, sizeof(long));
  for (i=0; i<size; i++)
    stfl->supp_pattern_list[i]=0;
  stfl->zero_pattern_list = (long *) ut_calloc(size, sizeof(long));
  for (i=0; i<size; i++)
    stfl->zero_pattern_list[i]=0;
  
  stfl->tau = (double **) ut_calloc(size, sizeof(double*));
  *(stfl->tau) = (double *) ut_calloc(size*k, sizeof(double));
  for (i=1; i<size; i++) 
    stfl->tau[i] = stfl->tau[i-1] + k;
  for (i=0; i<size; i++)
    for (r=0; r<k; r++) 
      stfl->tau[i][r] = 0;

  stfl->matrix_f = (double **) ut_calloc(size, sizeof(double*));
  *(stfl->matrix_f) = (double *) ut_calloc(size*k, sizeof(double));
  for (i=1; i<size; i++) 
    stfl->matrix_f[i] = stfl->matrix_f[i-1] + k;
  for (i=0; i<size; i++)
    for (r=0; r<k; r++) 
      stfl->matrix_f[i][r] = 0;
  
}




void mc_statefull_destruct(MCStateFull *stfl) {

  free(stfl->supp_pattern_list);
  stfl->supp_pattern_list = NULL;
  free(stfl->zero_pattern_list);
  stfl->zero_pattern_list = NULL;
  free(*stfl->tau);
  free(stfl->tau);
  stfl->tau=NULL;
  free(*stfl->matrix_f);
  free(stfl->matrix_f);
  stfl->matrix_f=NULL;


}

void mc_statefull_clear(MCStateFull *stfl) {

  stfl->supp_pattern_list = NULL;
  stfl->zero_pattern_list = NULL;
  stfl->tau=NULL;
  stfl->matrix_f=NULL;
}


void mc_statefull_read(MCStateFull *stfl, FILE *file) {

  long i;
  ut_fread(&stfl->size, sizeof(long), 1, file);
  ut_fread(&stfl->k, sizeof(long), 1, file);
  ut_fread(&stfl->n_supp_pattern, sizeof(long), 1, file);
  mc_statefull_construct(stfl, stfl->size, stfl->k);
  ut_fread(stfl->supp_pattern_list, sizeof(long), stfl->size, file);
  ut_fread(stfl->zero_pattern_list, sizeof(long), stfl->size, file);
  for (i=0; i<stfl->size; i++)
    ut_fread(stfl->tau[i]   , sizeof(double), stfl->k, file);
  for (i=0; i<stfl->size; i++)
    ut_fread(stfl->matrix_f[i]   , sizeof(double), stfl->k, file);
}

void mc_statefull_write(const MCStateFull *stfl, FILE *file) {

  long i;
  ut_fwrite(&stfl->size, sizeof(long), 1, file);
  ut_fwrite(&stfl->k, sizeof(long), 1, file);
  ut_fwrite(&stfl->n_supp_pattern, sizeof(long), 1, file);
  ut_fwrite(stfl->supp_pattern_list, sizeof(long), stfl->size, file);
  ut_fwrite(stfl->zero_pattern_list, sizeof(long), stfl->size, file);
  for (i=0; i<stfl->size; i++)
    ut_fwrite(stfl->tau[i]   , sizeof(double), stfl->k, file);
  for (i=0; i<stfl->size; i++)
    ut_fwrite(stfl->matrix_f[i]   , sizeof(double), stfl->k, file);

}


void mc_statefull_text(const MCStateFull *stfl, FILE *file) {
  long i, r;
  
  fprintf(file, "%ld\n", stfl->size);
  fprintf(file, "%ld\n", stfl->k);
  fprintf(file, "%ld\n", stfl->n_supp_pattern);
  fprintf(file, "%ld\n", stfl->n_zero_pattern);

  for (i=0; i<stfl->n_supp_pattern; i++)
    fprintf(file, "%ld ", stfl->supp_pattern_list[i]);
  fprintf(file, "\n");
  
  for (i=0; i<stfl->n_zero_pattern; i++)
    fprintf(file, "%ld ", stfl->zero_pattern_list[i]);
  fprintf(file, "\n");
  
  for (i=0; i<stfl->size; i++) {
    for (r=0; r<stfl->k; r++) {
      fprintf(file, "%15.10f", stfl->tau[i][r]);
    }
    fprintf(file, "\n");
  }
  
  fprintf(file, "\n");
  for (i=0; i<stfl->size; i++) {
    for (r=0; r<stfl->k; r++) {
      fprintf(file, "%15.10f", stfl->matrix_f[i][r]);
    }
    fprintf(file, "\n");
  }
  
  fprintf(file, "\n%s\n", stfl->comment);
  fflush(file);


}


void mc_classifier_read(MCClassifier *cls, FILE *file) {
  ut_fread(cls->comment, sizeof(char), COMMENT_LEN, file); 
  mc_scaledef_read(&cls->scale_def, file);
  kernel_read(&cls->kernel_def, file);
  mc_solution_read(&cls->solution, file);
}

void mc_classifier_write(const MCClassifier *cls, FILE *file) {

  ut_fwrite(cls->comment, sizeof(char), COMMENT_LEN, file); 
  mc_scaledef_write(&cls->scale_def, file);
  kernel_write(&cls->kernel_def, file);
  mc_solution_write(&cls->solution, file);

}

void mc_classifier_text(const MCClassifier *cls, FILE *file) {
  fprintf(file, "\n%s\n", cls->comment);
  mc_scaledef_text(&cls->scale_def, file);
  kernel_text(&cls->kernel_def, file);
  mc_solution_text(&cls->solution, file);

}

void mc_classifier_destruct(MCClassifier *cls) {
  mc_scaledef_destruct(&cls->scale_def);
  mc_solution_destruct(&cls->solution);

}
