#!/bin/bash
while read SRR; do
    fasterq-dump --split-files --threads 8 $SRR -O data/
done < sra_ids_10.txt
