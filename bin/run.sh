#!/bin/sh

SUFFIX=cmo
LIB_SUFFIX=cma

BASE_OBJECTS="str.$LIB_SUFFIX seq.$SUFFIX lazy_io.$SUFFIX misc_lib.$SUFFIX string_ext.$SUFFIX fasta.$SUFFIX genbank.$SUFFIX hmm.$SUFFIX orf_bounding.$SUFFIX gene_predictor.$SUFFIX gene_prediction_2.$SUFFIX gene_prediction_4.$SUFFIX gene_prediction_7.$SUFFIX gene_prediction_10.$SUFFIX gene_prediction_3_4.$SUFFIX"

ocaml $BASE_OBJECTS
