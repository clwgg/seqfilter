#include <zlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Heng Li klib
#include "klib/kseq.h"
#include "klib/kbtree.h"

// kbtree struct
typedef struct {
    char *key;
    int count;
} elem_t;

// initialize kbtree
#define elem_cmp(a, b) (strcmp((a).key, (b).key))
KBTREE_INIT(str, elem_t, elem_cmp)

// initialize kseq
KSEQ_INIT(gzFile, gzread)

// prototype output function
void make_output(kseq_t *seq, int *l, int *n, long *nseq, int *minseq, int *maxseq, FILE *pass);

int main(int argc, char *argv[])
{

  int c; 
  int max = 0;
  int min = 0;
  char *in = 0;
  char *out = 0;
  char *id = 0;
  bool flag = false;
  while ((c = getopt(argc, argv, "i:o:l:x:m:n")) >= 0) {
    switch (c) {
            case 'i': in = optarg; break;
            case 'o': out = optarg; break;
            case 'l': id = optarg; break;
            case 'x': max = atoi(optarg); break;
            case 'm': min = atoi(optarg); break;
            case 'n': flag = true; break;
    }
  }

  if (!in || !out) {
    fprintf(stderr, "Usage: %s -i <in.fa/q> -l <ID.list> -o <out.fa/q>\n\n", argv[0]);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "        -i   Input file (fasta/fastq)\n");
    fprintf(stderr, "        -o   Output file\n");
    fprintf(stderr, "        -l   ID list file (ID\\nID)\n");
    fprintf(stderr, "        -x   Maximum sequence length to put out\n");
    fprintf(stderr, "        -m   Minimum sequence length to put out\n");
    fprintf(stderr, "        -n   Set for negative filtering instead of positive\n");
    return 1;
  }

  // Initialize kbtree
  kbtree_t(str) *b;
  elem_t *p, t;
  b = kb_init(str, KB_DEFAULT_SIZE);

  // Initialize string stream
  kstream_t *ks;
  kstring_t s = {0,0,0};

  // Fill kbtree if id-file is defined
  if (id) {

    // open ID file
    gzFile f;
    f = gzopen(id, "r");
    ks = ks_init(f);

    while (ks_getuntil(ks, '\n', &s, 0) >= 0) {

      char *tmp = malloc( s.m );
      strcpy(tmp, s.s);

      t.key = tmp, t.count = 1;
      p = kb_getp(str, b, &t);
      if (!p) kb_putp(str, b, &t);
      else ++p->count;

    }

    ks_destroy(ks);
    gzclose(f);
    free(s.s);
  }

  // open output file
  FILE *pass;
  pass = fopen(out,"w+");
  
  // open sequence file
  gzFile fp;
  fp = gzopen(in, "r");

  // initialize sequence variables and btree elements for checking existence of IDs
  elem_t check;

  kseq_t *seq;
  seq = kseq_init(fp);
  int l;

  // initialize variables for filter statistics
  int n = 0;
  long nseq = 0;
  int minseq = 0;
  int maxseq = 0;

  // loop through input sequence file, apply filters and produce output
  while ((l = kseq_read(seq)) >= 0) {
    check.key = seq->name.s, check.count = 1;
   
    p = kb_getp(str, b, &check);

    if (min && l < min) {
      continue;
    }
    if (max && l > max) {
      continue;
    }

    if(p && !flag) {
      make_output(seq,&l,&n,&nseq,&minseq,&maxseq,pass);
    }
    else if(!p && flag) {
      make_output(seq,&l,&n,&nseq,&minseq,&maxseq,pass);
    }
  }

  fprintf(stderr, "After Filtering:\n");
  fprintf(stderr, "Total n:\t%i\n", n);
  fprintf(stderr, "Total seq:\t%ld bp\n", nseq);
  fprintf(stderr, "Min seq:\t%i bp\n", minseq);
  fprintf(stderr, "Max seq:\t%i bp\n", maxseq);

  // close and free everything thats left
  fclose(pass);
  kseq_destroy(seq);
  gzclose(fp);

  // iterate through kbtree and free all keys
  kbitr_t itr;
  kb_itr_first(str, b, &itr);
  for (; kb_itr_valid(&itr); kb_itr_next(str, b, &itr)) { 
    p = &kb_itr_key(elem_t, &itr);
    free(p->key);
  }
  kb_destroy(str, b);
  return 0;

}

// function to produce output
// WARNING: works on pointers to variables in main()
void make_output(kseq_t *seq, int *l, int *n, long *nseq, int *minseq, int *maxseq, FILE *pass)
{ 
  if(seq->qual.s) {
    fprintf(pass, "@%s\n%s\n+\n%s\n", seq->name.s, seq->seq.s, seq->qual.s);
  }
  else {
    fprintf(pass, ">%s\n%s\n", seq->name.s, seq->seq.s);
  }
  *n = *n + 1;
  *nseq = *nseq + *l;
  if (*minseq) {
    if (*minseq > *l) {
      *minseq = *l;
    }
  }
  else {
    *minseq = *l;
  }
  if (*maxseq) {
    if (*maxseq < *l) {
      *maxseq = *l;
    }
  }
  else {
    *maxseq = *l;
  }
}

