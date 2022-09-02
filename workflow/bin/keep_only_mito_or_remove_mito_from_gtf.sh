#!/usr/bin/env bash

original_gtf=$1


awk -F '\t' '$1=="chrM"' $original_gtf >> only_mito.gtf
awk -F '\t' '$1=="chrM"' $original_gtf >> only_nuclear.gtf
