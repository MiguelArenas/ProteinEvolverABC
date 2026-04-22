/*  
   Copyright (c)2005, Trevor Bruen
   All rights reserved.                          

   Any feedback is very welcome.
   email: trevor@mcb.mcgill.ca
*/

#include <stdlib.h>
#include <stdio.h>
#include "queue.h"
#include "mem.h"

void init_queue(queue *q)
{
  (*q).elements=(queue_elem *)mcalloc( (q->max_size),sizeof(queue_elem));
  q->front=0;
  q->back=0;
  q->cur_size=0;
}

void destroy(queue *q)
{
  free((*q).elements);
}

queue *enqueue(queue *q, queue_elem elem)
{
  if(q->cur_size == q->max_size)
    {
      return NULL;
    }
  else
    {
      (*q).elements[q->back]=elem;
      (q->cur_size)=(q->cur_size)+1;
      /* Loop around storage space */
      q->back=(q->back+1) % (q->max_size);
    }
  return q;
}

queue *element_at(queue *q, queue_elem *elem, int index)
{
  int mapped_index=0;
  if(index > q->cur_size)
    return NULL;
  else
    {
      mapped_index=((q->front)+(index)) % (q->max_size);
      *elem=(*q).elements[mapped_index];
      
    }
  return q;
}

void clear_queue(queue *q)
{
  q->front=0;
  q->back=0;
  q->cur_size=0;
}
void print_queue(queue *q)
{
  int cur=q->cur_size,i;
  queue_elem elem=0;
  fprintf(stdout,"Queue: ");
  for(i=0;i<cur;i++)
    {
      element_at(q,&elem,i);
      fprintf(stdout,"%d ",elem);

    }
  fprintf(stdout,"\n");
}
queue *dequeue_front(queue *q, queue_elem *elem)
{
  if(q->cur_size == 0)
    {
      //fprintf(stderr,"Cannot dequeue\n");
      return NULL;
    }
  else
    {
      *elem=(*q).elements[q->front];
      (q->cur_size)=(q->cur_size)-1;
      /* Loop around storage space */
      q->front=(q->front+1) % (q->max_size);
    
    }
  return q;
}
