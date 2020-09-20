
#ifndef MEM
#define MEM


size_t alloc_size();

void *mcalloc(size_t num, size_t size);

void *mmalloc(size_t size);


void *mrealloc(void *ptr, size_t old_size,size_t new_size);

void mfree(void *ptr, size_t size);

void error(const char * message); 


char ffclose(FILE *handle);


FILE *ffopen(char *name,char *mod);

#endif
