(* 4 state hmm *)

type coding_state = Q0 | NotGene | C1 | C2 | C3


(*
 * This creates our training.  The algorithm here uses
 * the trainer for gene_prediction_2 and then takes the
 * results and cuts them into 3 state pieces representing
 * the codons
 *)
let create_training_data gene_boundaries fin =
  let split_gene gd =
    let rec sg acc gd = 
      let l = String_ext.length gd in
      if l < 3 then
	(Seq.of_list (List.rev acc), gd)
      else
	sg ((C3, String_ext.make 1 gd.[2])::(C2, String_ext.make 1 gd.[1])::(C1, String_ext.make 1 gd.[0])::acc) (String_ext.sub gd 3 (l - 3))
    in
    sg [] gd
  in
  let rec ctd d sin =
    if d = "" then
      begin
	match Seq.next sin with
	    Some (Gene_prediction_2.NotGene, d) ->
	      [< '(NotGene, d); ctd "" sin >]
	  | Some (Gene_prediction_2.Gene, d) ->
	      let (c, d') = split_gene d in
	      [< c; ctd d' sin >]
	  | Some (Gene_prediction_2.Q0, _) ->
	      raise (Failure "Not supposed to happen")
	  | None ->
	      [< >]
      end
    else
      begin
	match Seq.next sin with
	    Some (Gene_prediction_2.NotGene, _) ->
	      raise (Failure ("Expecting a gene: d: " ^ d))
	  | Some (Gene_prediction_2.Gene, v) ->
	      let (c, d') = split_gene (d ^ v) in
	      [< c; ctd d' sin >]
	  | Some (Gene_prediction_2.Q0, _) ->
	      raise (Failure "Not supposed to happen")
	  | None ->
	      raise (Failure "Expecting more gene data but got nothing")
      end
  in
  let sin = Gene_prediction_2.create_training_data gene_boundaries fin in
  ctd "" sin
