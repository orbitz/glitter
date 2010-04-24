#!/bin/sh

OCAMLC="ocamlc"

if [ $OCAMLC = "ocamlc" ]
then
    SUFFIX=cmo
    LIB_SUFFIX=cma
else
    SUFFIX=cmx
    LIB_SUFFIX=cmxa
fi

BASE_OBJECTS="str.$LIB_SUFFIX seq.$SUFFIX lazy_io.$SUFFIX misc_lib.$SUFFIX string_ext.$SUFFIX fasta.$SUFFIX genbank.$SUFFIX hmm.$SUFFIX gene_predictor.$SUFFIX gene_prediction_2.$SUFFIX gene_prediction_4.$SUFFIX gene_prediction_7.$SUFFIX gene_prediction_10.$SUFFIX gene_prediction_3_4.$SUFFIX"

mkdir -p ../output

for i in *.ml; do $OCAMLC -pp "camlp4o pa_extend.cmo" -I +camlp4 -g -c $i; done

# $OCAMLC -g -o ../output/gb_to_genelist $BASE_OBJECTS gb_to_genelist.$SUFFIX

# $OCAMLC -g -o ../output/glitter_2 $BASE_OBJECTS glitter_2.$SUFFIX
# $OCAMLC -g -o ../output/glitter_4 $BASE_OBJECTS glitter_4.$SUFFIX
# $OCAMLC -g -o ../output/glitter_7 $BASE_OBJECTS glitter_7.$SUFFIX
# $OCAMLC -g -o ../output/glitter_10 $BASE_OBJECTS glitter_10.$SUFFIX
# $OCAMLC -g -o ../output/glitter_3_4 $BASE_OBJECTS glitter_3_4.$SUFFIX

