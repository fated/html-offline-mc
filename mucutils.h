/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:35 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mucutils.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 6.2 $ */
/* $State: Exp $ */



#ifndef __MUCUTILS_H
#define __MUCUTILS_H

#include "mcdata.h"
#include "kernel.h"


#define COMMENT_LEN 512

typedef struct {
  char   comment[COMMENT_LEN];
  long   to_scale;
  long   l;
  double scale_factor;
  long   to_zero_data_mean;
  double *data_mean;
} MCScaleDef;


typedef struct {
  long size, k, l;
  long n_supp_pattern;
  long is_voted;
  long *supp_pattern_list;
  long *votes_weight;
  double **tau;
  char comment[COMMENT_LEN];
} MCSolution;


typedef struct {
  long size, k;
  long n_supp_pattern;
  long is_voted;
  long *supp_pattern_list;
  long *votes_weight;
  double **tau;
  double **matrix_f;
  char comment[COMMENT_LEN];
} MCStateSparse;

typedef struct {
  long size, k;
  long n_supp_pattern;
  long n_zero_pattern;
  long *supp_pattern_list;
  long *zero_pattern_list;
  double **tau;
  double **matrix_f;
  char comment[COMMENT_LEN];
} MCStateFull;

typedef struct {
  char       comment[COMMENT_LEN];
  MCScaleDef scale_def;
  MCSolution solution;
  KernelDef  kernel_def;
} MCClassifier;

void mc_scaledef_read(MCScaleDef *scale_def, FILE *file);
void mc_scaledef_write(const MCScaleDef *scale_def, FILE *file);
void mc_scaledef_text(const MCScaleDef *scale_def, FILE *file);
void mc_scaledef_scale(const MCScaleDef *scale_def, MCDataUnit **data, long m);
void mc_scaledef_initialize(MCScaleDef *scale_def);
void mc_scaledef_destruct(MCScaleDef *scale_def);

void mc_solution_construct(MCSolution *solution, long size, long k, long l, long n_supp_pattern, long is_voted);
void mc_solution_destruct(MCSolution *solution);
void mc_solution_clear(MCSolution *solution);
void mc_solution_read(MCSolution *solution, FILE *file);
void mc_solution_write(const MCSolution *solution, FILE *file);
void mc_solution_text(const MCSolution *solution, FILE *file);

void mc_statesparse_construct(MCStateSparse *stsp, long size, long k, long is_voted);
void mc_statesparse_destruct(MCStateSparse *stsp);
void mc_statesparse_clear(MCStateSparse *stsp);
void mc_statesparse_read(MCStateSparse *stsp, FILE *file);
void mc_statesparse_write(const MCStateSparse *stsp, FILE *file);
void mc_statesparse_text(const MCStateSparse *stsp, FILE *file);

void mc_statefull_construct(MCStateFull *stfl, long size, long k);
void mc_statefull_destruct(MCStateFull *stfl);
void mc_statefull_clear(MCStateFull *stfl);
void mc_statefull_read(MCStateFull *stfl, FILE *file);
void mc_statefull_write(const MCStateFull *stfl, FILE *file);
void mc_statefull_text(const MCStateFull *stfl, FILE *file);

void mc_classifier_read(MCClassifier *cls, FILE *file);
void mc_classifier_write(const MCClassifier *cls, FILE *file);
void mc_classifier_text(const MCClassifier *cls, FILE *file);
void mc_classifier_destruct(MCClassifier *cls);


#endif
