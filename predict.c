#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <errno.h>
#include <math.h>
#include "cgic.h"
#include "spoc.h"
#include "utilities.h"
#include "mucutils.h"

#define BufferLen 1024

////////////////

typedef struct {
  long n_errors;
  long **error_statistics;
} ErrStatistics;

static MCScaleDef  mc_scaledef;
static MCDataDef   datadef_train;
static MCDataDef   datadef_test;
static SPOCParamDef spoc_pd; 
static MCClassifier mc_cls;

static char train_name[BufferLen], test_name[BufferLen]; 

static long **supp_tau_lists;
static long *n_supp_tau;

static long *predict_label;
static double **similarity_score;

static ErrStatistics error_test;
static KF kernel_fun;

void read_data(char *file_name, MCDataDef *pd);
void initialize_train();
void initialize_test();
void data_errors(ErrStatistics *es);
long mc_score();
void venn_predict(long pl);

void free_memory();

////////////////

void HandleSubmit();
void File();
void ShowForm();

int cgiMain() {
  if ((cgiFormSubmitClicked("predict") == cgiFormSuccess)) {
    HandleSubmit();
  } else {
    ShowForm();
  }
  return 0;
}

void ShowForm() {
  cgiHeaderContentType("text/html");
  fprintf(cgiOut, "<html><head>\n");
  fprintf(cgiOut, "<title>Upload Data Set</title></head>\n");
  fprintf(cgiOut, "<body><h1>Upload Data Set</h1>\n");
  fprintf(cgiOut, "<p>Data set should be prepared in a single file (training set and a test set should be saved in different files).</p>\n");
  fprintf(cgiOut, "<p>The first three lines should be :<br>\n no._of_instances<br />\n no._of_classes<br />\n no._of_attributes<br /></p>\n");
  fprintf(cgiOut, "<p>Each instance is placed in a line in the following format :<br>\n label feat_1 ... feat_l<br /></p>\n");
  fprintf(cgiOut, "<p>A label is an integer between 0 and k-1.<br />\n");
  fprintf(cgiOut, "A feat is any number in a floating point format.</p>\n");

  fprintf(cgiOut, "<!-- 2.0: multipart/form-data is required for file uploads. -->");
  fprintf(cgiOut, "<form method=\"POST\" enctype=\"multipart/form-data\" ");
  fprintf(cgiOut, "	action=\"");
  cgiValueEscape(cgiScriptName);
  fprintf(cgiOut, "\">\n");
  fprintf(cgiOut, "<p>Training Data Upload:\n");
  fprintf(cgiOut, "<input type=\"file\" name=\"train\" value=\"\"> (Select A Local File)<br />\n");
  fprintf(cgiOut, "Testing Data Upload:\n");
  fprintf(cgiOut, "<input type=\"file\" name=\"test\" value=\"\"> (Select A Local File)<br />\n");
  fprintf(cgiOut, "<input type=\"submit\" name=\"predict\" value=\"Submit\">\n");
  fprintf(cgiOut, "<input type=\"reset\" value=\"Reset\"></p>\n");
  fprintf(cgiOut, "</form></body></html>\n");
}

void HandleSubmit() {
  int i;
  double sum, min;

  cgiHeaderContentType("text/html");
  fprintf(cgiOut, "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\">\n\n");
  fprintf(cgiOut, "<title>Offline Venn Machine Prediction</title></head><body background=\"../gifs/rhultile.gif\">\n");
  fprintf(cgiOut, "<center><h2>Offline Venn Machine Prediction <font size=\"-1\">[go to <a href=\"../index.html\">Main Menu</a>]</font></h2></center>\n");

  File();
  
  mc_score();

  /*
  min = similarity_score[0][0];
  for (i=1; i<10; i++)
    if (similarity_score[0][i]<min)
      min = similarity_score[0][i];
  sum = similarity_score[0][predict_label]-min;
  fprintf(cgiOut, "<center><table border=\"1\">");
  for (i=0; i<10; i++) {
    fprintf(cgiOut, "<tr><td>");
    if (i==predict_label)
      fprintf(cgiOut, "<font color=\"red\">");
    fprintf(cgiOut, "%d", i);
    if (i==predict_label)
      fprintf(cgiOut, "</font>");
    fprintf(cgiOut, "</td><td>");
    if (i==predict_label)
      fprintf(cgiOut, "<font color=\"red\">");
    fprintf(cgiOut, "%.2f%%", ((similarity_score[0][i]-min) * 100 / sum));
    if (i==predict_label)
      fprintf(cgiOut, "</font>");
    fprintf(cgiOut, "</td></tr>");

  }

  
  //venn_predict(predict_label);
  */
  for (i=0; i<datadef_test.m; i++)
    fprintf(cgiOut, "%ld ", predict_label[i]);

  fprintf(cgiOut, "<br />Test Error : %4.2f%% (%ld / %ld) <br />\n", 100*((double)error_test.n_errors / (double)datadef_test.m), error_test.n_errors, datadef_test.m);
  fprintf(cgiOut, "<center><a href=\"");
  cgiValueEscape(cgiScriptName);
  fprintf(cgiOut, "\">Back</a></center></body></html>\n");

  free_memory();

}

void venn_predict(long pl) {
  FILE *file;
  long int i, j, k;
  int st[5][10], idx, pre, label, tlabel;
  float score, t[10][10], q[10], sum, min, max, lb, ub;

  for (i=0; i<5; i++)
    for (j=0; j<10; j++)
      st[i][j] = 0;
  file = fopen("usps1.rpt", "r");
  for (i=0; i<7291; i++) {
    fscanf(file, "%d", &label);
    fscanf(file, "%f", &score);
    st[(int)(score/0.1)][label]+=1;
  }
  fclose(file);
  score = similarity_score[0][pl];
  for (j=0; j<10; j++) {
    if (score<0.1) 
      idx=0;
    else if (score<0.2) 
      idx=1;
    else if (score<0.3) 
      idx=2;
    else if (score<0.4) 
      idx=3;
    else idx=4;
    st[idx][j] = st[idx][j]+1;
    sum = 0;
    for (k=0; k<10; k++) 
      sum = sum + st[idx][k];
    min = 1;
    for (k=0; k<10; k++) {
      t[j][k]=((float)st[idx][k])/sum;
      if (t[j][k]<min) {
	min = t[j][k];
      }
    }
    q[j] = min;
    st[idx][j] = st[idx][j]-1;
  }
  // for(i=0; i<10; i++) {
    max = 0;
    for (k=0; k<10; k++)
      if (q[k]>max) {
	max = q[k];
	pre = k;
      }
    // if (pre != -1) {
    min = q[k]; max = 0;
    for (k=0; k<10; k++)
      if (t[pre][k] > max)
	max = t[pre][k];
      // }
      //      q[pre] = -1;
      //pre = -1;
      //}
    
fprintf(cgiOut, "Venn Machine predict %d with confidence [%.2f%%, %.2f%%]", pre, (1-max)*100, (1-min)*100);
}

long mc_score() {

  read_data(train_name, &datadef_train);
  read_data(test_name, &datadef_test);

  initialize_train();
  mc_cls.solution= spoc(datadef_train, spoc_pd);
  mc_cls.scale_def = mc_scaledef;
  mc_cls.kernel_def = spoc_pd.kernel_def;

  initialize_test();
  kernel_construct(mc_cls.solution.l);
  kernel_fun = kernel_get_function(mc_cls.kernel_def);

  data_errors(&error_test);

  return (1);
}

void read_data(char *file_name, MCDataDef *pd) {
  FILE *file;

  file = ut_fopen(file_name, "r");
  mc_datadef_read(pd, file);
  fclose(file);
  mc_scaledef_scale(&mc_cls.scale_def, pd->x, pd->m);
  
}

void initialize_train() {
  spoc_pd.beta = 1e-4;
  spoc_pd.cache_size = 4096;
  
  spoc_pd.kernel_def.kernel_type = KERNEL_EXPONENT_NP;
  spoc_pd.kernel_def.polynom_degree = 1;
  spoc_pd.kernel_def.polynom_a0 = 1;
  spoc_pd.kernel_def.exponent_sigma = .5;
  
  spoc_pd.epsilon =  1e-3;
  spoc_pd.epsilon0 = (1-1e-6);
  spoc_pd.delta =    1e-4;
  spoc_pd.redopt_type = REDOPT_TYPE_EXACT;
  mc_scaledef_initialize(&mc_scaledef);
}


void initialize_test() {
  long r;
  long si;
   
  n_supp_tau = (long *) ut_calloc(mc_cls.solution.k, sizeof(long));
  predict_label = (long *) ut_calloc(datadef_test.m, sizeof(long));
  similarity_score = (double **) ut_calloc(datadef_test.m, sizeof(double*));
  *similarity_score = (double *) ut_calloc(datadef_test.k * datadef_test.m, sizeof(double));
  
  for (r=0; r<mc_cls.solution.k; r++)
    n_supp_tau[r] = 0;

  for (r=0; r<datadef_test.m; r++)
    predict_label[r] = 0;

  for (r=1; r<datadef_test.m; r++) 
    similarity_score[r] = similarity_score[r-1] + datadef_test.k;

  for (r=0; r<datadef_test.m; r++)
    for (si=0; si<datadef_test.k; si++)
      similarity_score[r][si] =0;

  supp_tau_lists = (long **) ut_calloc(mc_cls.solution.k, sizeof(long*));
  *supp_tau_lists = (long *) ut_calloc(mc_cls.solution.n_supp_pattern * mc_cls.solution.k, sizeof(long));
  for (r=1; r<mc_cls.solution.k; r++) 
    supp_tau_lists[r] = supp_tau_lists[r-1] + mc_cls.solution.n_supp_pattern;
  
  for (r=0; r<mc_cls.solution.k; r++) {
    for (si=0; si<mc_cls.solution.n_supp_pattern; si++) {
      if (mc_cls.solution.tau[si][r] != 0) {
	supp_tau_lists[r][n_supp_tau[r]++] = si;
      }
    }
  }

}

void data_errors(ErrStatistics *es) {

  long i, si;
  long r, s;
  double *kernel_values = NULL;
  double sim_score;
  double max_sim_score;
  long n_max_sim_score;
  long best_y;
  long supp_pattern_index;

  es->n_errors = 0;
  kernel_values = (double *) ut_calloc(mc_cls.solution.n_supp_pattern, sizeof(double));
  es->error_statistics = (long **) ut_calloc(mc_cls.solution.k, sizeof(long*));
  *es->error_statistics = (long *) ut_calloc(mc_cls.solution.k * mc_cls.solution.k, sizeof(long));

  for (r=1; r<mc_cls.solution.k; r++) 
    es->error_statistics[r] = es->error_statistics[r-1] + datadef_test.k;
  
  for (r=0; r<mc_cls.solution.k; r++)
    for (s=0; s<mc_cls.solution.k; s++)
      es->error_statistics[r][s] =0;

  for (i=0; i<datadef_test.m; i++) {
    for (si=0; si<mc_cls.solution.n_supp_pattern; si++)
      kernel_values[si] = kernel_fun(datadef_test.x[i], datadef_train.x[mc_cls.solution.supp_pattern_list[si]]);
    
    n_max_sim_score =0;
    max_sim_score = -DBL_MAX;
    best_y = -1;
    
    for (r=0; r<datadef_test.k; r++) {
      sim_score =0;
      for (si=0; si<n_supp_tau[r]; si++) {
	supp_pattern_index = supp_tau_lists[r][si];
	sim_score += mc_cls.solution.tau[supp_pattern_index][r] * kernel_values[supp_pattern_index];     
      }
      similarity_score[i][r] = sim_score;
      
      if (sim_score > max_sim_score) {
	max_sim_score = sim_score;
	n_max_sim_score =1;
	best_y = r;
	predict_label[i] = r;
      }
      else if (sim_score == max_sim_score)
	n_max_sim_score++;
    }
       
    es->error_statistics[datadef_test.y[i]][best_y]++;
    if ((n_max_sim_score>1) || (best_y != datadef_test.y[i])) {
      es->n_errors++;
    } 
  }
  free(kernel_values);
  
}


void free_memory() {
  mc_datadef_destruct(&datadef_train);
  mc_datadef_destruct(&datadef_test);
  mc_classifier_destruct(&mc_cls);
  spoc_destruct();

  free(n_supp_tau);
  free(*supp_tau_lists);
  free(supp_tau_lists);

  free(predict_label);
  free(*similarity_score);
  free(similarity_score);

  free(*error_test.error_statistics);
  free(error_test.error_statistics);
}

void File() {   
  cgiFilePtr file;
  FILE *targetFile;   
  char buffer[BufferLen]; 
  int size, got, i, t;

  if ((cgiFormFileName("train", train_name, sizeof(train_name)) != cgiFormSuccess) || (cgiFormFileName("test", test_name, sizeof(test_name)) != cgiFormSuccess)) {
    fprintf(cgiOut, "<p>Please select both files.</p>\n");
    return;
  }
  if (cgiFormFileOpen("train", &file) != cgiFormSuccess) {
    fprintf(cgiOut, "Could not open the file.<p>\n");
    return;
  }
  if (access(train_name, 0) == 0) {
    fprintf(cgiOut, "File already exists on server!<p>\n");
    return;
  }
  targetFile=fopen(train_name, "w");
  if (cgiFormFileRead(file, buffer, sizeof(buffer), &got) == cgiFormSuccess) {
    sscanf(buffer, "%ld %ld %ld", &datadef_train.m, &datadef_train.k, &datadef_train.l);
    t = 0;
    for (i=0; i<got; i++)
      if (buffer[i] == '\n') {
	t++;
	if (t == 3) break;
      }
    got = got - i - 1;
    strcpy(buffer, buffer+i+1);  
    fwrite(buffer, 1, got, targetFile);
  }
  while (cgiFormFileRead(file, buffer, sizeof(buffer), &got) ==
	 cgiFormSuccess)
    {
      fwrite(buffer, 1, got, targetFile);
    }
  cgiFormFileClose(file);
  fclose(targetFile);
  if (cgiFormFileOpen("test", &file) != cgiFormSuccess) {
    fprintf(cgiOut, "Could not open the file.<p>\n");
    return;
  }
  if (access(test_name, 0) == 0) {
    fprintf(cgiOut, "File already exists on server!<p>\n");
    return;
  }
  targetFile=fopen(test_name, "w");
  if (cgiFormFileRead(file, buffer, sizeof(buffer), &got) == cgiFormSuccess) {
    sscanf(buffer, "%ld %ld %ld", &datadef_test.m, &datadef_test.k, &datadef_test.l);
    t = 0;
    for (i=0; i<got; i++)
      if (buffer[i] == '\n') {
	t++;
	if (t == 3) break;
      }
    got = got - i - 1;
    strcpy(buffer, buffer+i+1);  
    fwrite(buffer, 1, got, targetFile);
  }
  while (cgiFormFileRead(file, buffer, sizeof(buffer), &got) ==
	 cgiFormSuccess)
    {
      fwrite(buffer, 1, got, targetFile);
    }
  cgiFormFileClose(file);
  fclose(targetFile);
}

