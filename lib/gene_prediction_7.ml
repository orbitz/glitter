(* 4 state hmm *)

type coding_state = Q0 | NotGene | A | T | G | C1 | C2 | C3

module H = Hmm.Make(struct type s = coding_state type a = char let compare = compare end)

(*
 * This creates our training.  The algorithm here uses
 * the trainer for gene_prediction_2 and then takes the
 * results and cuts them into 3 state pieces representing
 * the codons
 *)
let create_training_data gene_boundaries fin =
  let rec gene_start d sin =
    if String_ext.length d < 3 then
      match Seq.next sin with
	  Some (Gene_prediction_2.Gene, d') ->
	    gene_start (d ^ d') sin
	| Some (Gene_prediction_2.NotGene, _) ->
	    raise (Failure "At start of gene and gene end, wuuut?")
	| _ ->
	    raise (Failure "At start of gene and sequence end, wuuuut?")
    else
    let (l, d) = Hmm.labeled_of_string_all [A; T; G;] d in
    (Seq.of_list (List.map (fun (label, c) -> (label, String_ext.make 1 c)) l), d)
  in
  let split_gene gd =
    let (l, d) = Hmm.labeled_of_string_all [C1; C2; C3] gd in
    (Seq.of_list (List.map (fun (label, c) -> (label, String_ext.make 1 c)) l), d)
  in
  let rec ctd d sin =
    if d = "" then
      begin
	match Seq.next sin with
	    Some (Gene_prediction_2.NotGene, d) ->
	      [< '(NotGene, d); ctd "" sin >]
	  | Some (Gene_prediction_2.Gene, d) ->
	      let (c, d') = gene_start d sin in
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

let predict training_fname fasta_fname =
  Gene_predictor.predict training_fname fasta_fname Q0 H.train H.forward_viterbi (Gene_prediction_2.map_td_to_list create_training_data) Misc_lib.identity NotGene A
