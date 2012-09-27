/* $Author: kobics $ */
/* $Date: 2003/06/02 08:35:48 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/utilities.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 5.2 $ */
/* $State: Exp $ */



#ifndef __UTILITIES
#define __UTILITIES

#include <stdlib.h>
#include <stdio.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define FABS(x) ( (x > 0) ? (x) : (-x)) 

#define MB (1048576)
#define KB (1024)

#define ACTUAL_ZERO (1e-15)

void *ut_calloc(size_t nobj, size_t size);

FILE  *ut_fopen(const char *filename, const char *mode);
long   ut_fclose(FILE *stream);
size_t ut_fread(void *ptr, size_t size, size_t nobj, FILE *stream);
size_t ut_fwrite(const void *ptr, size_t size, size_t nobj, FILE *stream);

#endif
