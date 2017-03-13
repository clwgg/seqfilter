seqfilter: seqfilter.c klib/kseq.h klib/kbtree.h
	gcc $< -O3 -lz -o $@
