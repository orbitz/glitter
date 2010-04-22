#!/bin/sh
for i in *.ml; do ocamlc -pp "camlp4o pa_extend.cmo" -I +camlp4 -g -c $i; done
