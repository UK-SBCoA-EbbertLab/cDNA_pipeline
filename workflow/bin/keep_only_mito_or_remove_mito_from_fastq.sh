#!/usr/bin/env bash

original_fastq=$1



## Suppose you keep sequence names that you want to exclude in ids.txt and sequences in seq.fa:
awk 'BEGIN{while((getline<"ids.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' seq.fa

