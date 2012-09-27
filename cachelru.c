/* $Author: kobics $ */
/* $Date: 2003/01/22 15:17:34 $ */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/cachelru.c,v $ */
/* $Name:  $ */
/* $Locker:  $ */
/* $Revision: 5.1 $ */
/* $State: Exp $ */



#include <stdlib.h>
#include <stdio.h>
#include "cachelru.h"
#include "utilities.h"


/* data */
/* ---- */

/* double linked list */
typedef struct node *Node_p;

typedef struct node{
  Node_p prev, next;
  void *data;
  long index;
} Node;

/* variables */
static long n_data;            /* size of whole memory */
static long n_blocks;          /* size of cahce */

static Node** data_index;      /* data_index of blocks */
static Node* tail;             /* pointer to circle double linked list */


static void *memory_block =NULL;
static void *start_memory_block =NULL;

Node* allocate_node(size_t block_size);

void dump_cachelru();             /* for debugging */

/* functions */
/* --------- */

/* return the number of blocks or zero if not enough memory */
long cachelru_construct(long _n_data, long _n_blocks, size_t block_size) {
  
  long i;
  Node *node;
  Node *current;
  long memory_block_size =0;

  n_data = _n_data;
  n_blocks = _n_blocks;

#ifdef DEBUG_CACHELRU
  fprintf(stderr, "cachelru_construct : n_data\t = %ld\n", n_data);
  fprintf(stderr, "cachelru_construct : n_blocks\t = %ld\n", n_blocks);
  fprintf(stderr, "cachelru_construct : block_size\t = %d\n", block_size);
#endif 
  
  /* allocate and initialize data_index */
  data_index = (Node **) ut_calloc (n_data, sizeof(Node *));
  for (i=0; i<n_data; i++)
    data_index[i]=NULL;
  
#ifdef DEBUG_CACHELRU1
  fprintf(stderr, "cachelru_construct : data_index initilized %ld units of size %d\n", n_data,sizeof(Node*));
#endif 
  

#ifdef DEBUG_CACHELRU
  fprintf(stderr, "cachelru_construct : n_blocks = min(%ld,%ld)\n",n_blocks, 
(long)(((((double)MB) * ((double)MEMORY_SIZE)) / (((double)block_size) + ((double)(sizeof(Node)))))));
#endif 
  

  memory_block_size = block_size + sizeof(Node);
  n_blocks = MIN(n_blocks, ((long) ((((double)MB) * ((double)MEMORY_SIZE)) / ((double) memory_block_size))));

  memory_block = calloc(n_blocks, memory_block_size);

#ifdef DEBUG_CACHELRU
  fprintf(stderr, "cachelru_construct : n_blocks\t = %ld\n", n_blocks);
  fprintf(stderr, "cachelru_construct : memory_block_size\t= %ld B\n", memory_block_size);
#endif 
  
  while (memory_block == NULL) {
    n_blocks--;
    
#ifdef DEBUG_CACHELRU1
  fprintf(stderr, "cachelru_construct : n_blocks\t = %ld\n", n_blocks);
#endif 
  
  memory_block = calloc(n_blocks, memory_block_size);
  }
  
  start_memory_block = memory_block;

  //printf("CacheLRU allocated %ld blocks of size %10.3fKb\n",n_blocks,(double)memory_block_size/(double)KB); fflush(stdout);
  //printf("Allocated total %7.3fMb (Internal bound %dMb)\n"
	 //,((double)n_blocks * (double)memory_block_size)/(double)MB, 
	 //MEMORY_SIZE); fflush(stdout);
  

#ifdef DEBUG_CACHELRU
  fprintf(stderr, "cachelru_construct : allocated %ld blocks of size %ldB at address %p\n", 
	  n_blocks,memory_block_size, memory_block);
#endif


#ifdef DEBUG_CACHELRU1
  {
    long i,j;
    fprintf(stderr, "cachelru_construct : scanning memory\n");
      
      for (i=0; i<n_blocks; i++)
	for (j=0; j<memory_block_size; j++)
	  *((char*)(memory_block+j+i*memory_block_size)) = 0;
  }
#endif 
  

  /* allocate and initialize tail */
  tail = allocate_node(block_size);
  if (tail == NULL) {
    free(data_index);
    return (0);
  }
  current = tail;
    
#ifdef DEBUG_CACHELRU1
  fprintf(stderr, "cachelru_construct : initilized tail in pointer %p\n", tail);
#endif 
  
  /* allocate and initialize rest of list */
  for (i=1; i<n_blocks; i++) {
    node = allocate_node(block_size);
  
    if (node == NULL) {          /* if not enough memory */
      
#ifdef DEBUG_CACHELRU1
      fprintf(stderr, "cachelru_construct : can't allocate more blocks\n");
#endif 
      n_blocks = i;            /* and auto size        */
      break;                   /* resize cache         */
    }
    else {                       /* node allocated add it to end of list */
            
      current->next = node;
      node->prev = current;
      current = node;
    }
  }

  /* circle the list */
  current->next = tail;
  tail->prev = current;
  
#ifdef DEBUG_CACHELRU1
  fprintf(stderr, "cachelru_construct : n_blocks = %ld\n",n_blocks);
  fprintf(stderr, "cachelru_construct : total allocated memory > %ld Kb\n",((n_blocks*(block_size+sizeof(Node)))+n_data*sizeof(Node*))/1024);
#endif 
  return (n_blocks);
}


long cachelru_retrive(long i, void** pdata) {
  Node *node;

  if (data_index[i] != NULL) {         /* data in cache */

#ifdef DEBUG_CACHELRU1
    fprintf(stderr,"cachelru_retrive : in cache for %ld\n",i);
#endif

    node = data_index[i];
    *pdata = node->data;
    
    /* move node to end of list */
    if (node == tail)
      tail = tail->prev;
    else {
      (node->next)->prev = node->prev;
      (node->prev)->next = node->next;
      
      node->next = tail->next;
      node->prev = tail;
      (tail->next)->prev = node;
      tail->next = node;

    }
    return (1);
  }
  else {                               /* data not in cahce */

#ifdef DEBUG_CACHELRU1
    fprintf(stderr,"cachelru_retrive : out of cache for %ld throwing %ld\n",i,tail->index);
#endif
    if (tail->index != -1)
      data_index[tail->index] = NULL;
    data_index[i] = tail;
    tail->index=i;
    *pdata = tail->data;
    tail=tail->prev;
    return (0);
  }
}

/* destroy all alocated memory */
void cachelru_destruct() {

/*    long i; */
/*    Node *node; */
  

#ifdef DEBUG_CACHELRU
  fprintf(stderr, "cachelru_destruct : freeing memory_block\n");
#endif 

  free(start_memory_block);
  
/*  #ifdef DEBUG_CACHELRU */
/*    fprintf(stderr, "cachelru_destruct : freeing linked list\n"); */
/*  #endif  */

/*    for (i=0; i<n_blocks; i++) { */
/*      node = tail->prev; */
/*      free(tail); */
/*      tail=node; */
/*    } */

#ifdef DEBUG_CACHELRU1
  fprintf(stderr, "cachelru_destruct : freeing data_index\n");
#endif 
  free(data_index);
}


/* local function */
/* -------------- */


/* return node and its data memory */
/* return NULL otherwise           */
Node* allocate_node(size_t block_size) {
    Node* node;

    node = (Node *) memory_block;
    memory_block += sizeof(Node);
  
    node->data = memory_block;
    memory_block += block_size;
    
    node->index = -1;
    return (node);
}


void dump_cache() {
  long i;
  Node* current;

  current=tail;
  for (i=0; i<n_blocks; i++) {
    fprintf(stderr,"%5ld",current->index);
    current=current->prev;
  }
  fprintf(stderr,"\n");

  current=tail;
  for (i=0; i<n_blocks; i++) {
    fprintf(stderr,"%10ld",*((long*)current->data));
    current=current->prev;
  }
  fprintf(stderr,"\n");
  
  
  for (i=0; i<n_data; i++) {
    fprintf(stderr,"%p ",data_index[i]);
  }
  fprintf(stderr,"\n");
}






