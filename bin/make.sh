#!/bin/sh

OCAMLC=ocamlc

for i in *.ml; do $OCAMLC -pp "camlp4o pa_extend.cmo" -I +camlp4 -g -c $i; done
