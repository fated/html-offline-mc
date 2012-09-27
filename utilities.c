/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:36 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/utilities.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 5.1 $ */
/* $State: Exp $ */



#include "utilities.h"
#include <errno.h>

void *ut_calloc(size_t nobj, size_t size) {
  void *p = NULL;
  p = calloc(nobj, size);
  if (p==NULL) {
    perror("can't allocate memory");
    exit(1);
  }
  return (p);
}


FILE *ut_fopen(const char *filename, const char *mode) {
  FILE *f = NULL;
  f = fopen(filename, mode);
  if (f==NULL) {
    perror("can't open file");
    exit(1);
  }
  return (f);
}



long ut_fclose(FILE *stream) {
  if (fclose(stream)==EOF) {
    perror("can't close file");
    exit(1);
  }
  return (0);
}

size_t ut_fread(void *ptr, size_t size, size_t nobj, FILE *stream) {
  if (fread(ptr, size, nobj, stream) < nobj) {
    perror("can't read file");
    exit(1);
  }
  return (nobj);
}

size_t ut_fwrite(const void *ptr, size_t size, size_t nobj, FILE *stream) {
   if (fwrite(ptr, size, nobj, stream) < nobj) {
    perror("can't write file");
    exit(1);
  }
  return (nobj);
}

