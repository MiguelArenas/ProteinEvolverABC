#ifndef QUEUE
#define QUEUE 1

typedef int queue_elem;

typedef struct {
  int max_size;
  int front;
  int back;
  int cur_size;
  queue_elem *elements;
} queue;



void init_queue(queue *q);

void destroy(queue *q);

queue *enqueue(queue *q, queue_elem elem);

queue *element_at(queue *q, queue_elem *elem, int index);

void clear_queue(queue *q);

void print_queue(queue *q);

queue *dequeue_front(queue *q, queue_elem *elem);

#endif 
