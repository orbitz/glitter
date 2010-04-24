(* 4 state hmm *)

type coding_state = 
    Q0 | NotGene 
  | Start1 | Start2 | Start3 | C1 | C2 | C3 | Stop1 | Stop2 | Stop3

module H = Hmm.Make(struct type s = coding_state type a = char let compare = compare end)

let to_stream l = Seq.of_list (List.map (fun (label, c) -> (label, String_ext.make 1 c)) l)

let rec label_start d fin =
  if String_ext.length d >= 3 then
    match Hmm.labeled_of_string [Start1; Start2; Start3] d with
	Some (l, d') ->
	  (to_stream l, d')
      | _ -> raise (Failure "Labeling start condons failed")
  else
    match Seq.next fin with
	Some (Gene_prediction_2.Gene, d') ->
	  label_start (d ^ d') fin
      |	Some (Gene_prediction_2.NotGene, _) ->
	  raise (Failure "Looking for start codon but got NotGene")
      | Some (Gene_prediction_2.Q0, _) | None ->
	  raise (Failure "Looking for start codon but got EOF")

let rec consume_notgenes f fin =
  match Seq.next fin with
      Some (Gene_prediction_2.NotGene, d) ->
	[< '(NotGene, d); consume_notgenes f fin>]
    | Some (Gene_prediction_2.Gene, d) ->
	let (l, d') = label_start d fin in
	[< l; f d' fin >]
    | Some (Gene_prediction_2.Q0, _) | None ->
	[< >]

let label_codons d fin =
  let (l, d') = Hmm.labeled_of_string_all [C1; C2; C3] d in
  (to_stream l, d')
  
let rec consume_gene d fin =
  let dl = String_ext.length d in
  if dl > 6 then
    let start = String_ext.sub d 0 (dl - 3) in
    let last_3 = String_ext.sub d (dl - 3) 3 in
    let (l, d') = label_codons start fin in
    [< l; consume_gene (d' ^ last_3) fin >]
  else
    match Seq.next fin with
	Some (Gene_prediction_2.NotGene, d') when dl = 3 ->
	  let (l, _) = Hmm.labeled_of_string_all [Stop1; Stop2; Stop3] d in
	  [< to_stream l; '(NotGene, d'); consume_notgenes consume_gene fin >]
      | Some (Gene_prediction_2.NotGene, _) ->
	  raise (Failure "This should never happen")
      | Some (Gene_prediction_2.Q0, _) | None ->
	  let (l, _) = Hmm.labeled_of_string_all [Stop1; Stop2; Stop3] d in
	  [< to_stream l >]
      | Some (Gene_prediction_2.Gene, d') ->
	  consume_gene (d ^ d') fin

let create_training_data gene_boundaries fin =
  let sin = Gene_prediction_2.create_training_data gene_boundaries fin in
  consume_notgenes consume_gene sin


let predict training_fname fasta_fname =
  Gene_predictor.predict training_fname fasta_fname Q0 H.train H.forward_viterbi create_training_data NotGene Start1

