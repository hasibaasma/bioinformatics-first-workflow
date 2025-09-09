#!/bin/bash

for f in data/*.fastq; do
    base=$(basename $f .fastq)
    seqtk sample -s100 $f 0.1 > data/${base}_sub.fastq
done
