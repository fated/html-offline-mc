/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:34 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/cachelru.h,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 5.1 $ */
/* $State: Exp $ */



#ifndef __CACHELRU_H
#define __CACHELRU_H

/* MEMORY SIZE = 2GB = 2048 MB */
#define MEMORY_SIZE (1024*2)

long cachelru_construct(long _n_data, long _n_blocks, size_t block_size);
long cachelru_retrive(long i, void** pdata);
void cachelru_destruct();

#endif
