#!/bin/bash
esearch -db sra -query PRJNA1071720 \
  | efetch -format runinfo \
  | cut -d',' -f1 \
  | grep SRR \
  | shuf \
  | head -n  10 >  sra_ids_10.txt

