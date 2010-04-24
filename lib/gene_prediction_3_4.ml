(* A 3rd order HMM with 4 states *)

type coding_state = Q0 | NotGene | Start | Codon | Stop

module H = Hmm.Make(struct type s = coding_state type a = string let compare = compare end)

let join (c1, c2, c3) = String_ext.concat "" (List.map snd [c1; c2; c3])

let consume_next f (c1, c2, c3) sin = 
  match Seq.next sin with
      Some c4 ->
	f (c2, c3, c4) sin
    | None ->
	[< >]


let rec consume_intragene cont ((c1, _, _) as cs) sin =
  let joined = join cs in
  match fst c1 with
      Gene_prediction_10.Start2 
    | Gene_prediction_10.Start3 
    | Gene_prediction_10.C1
    | Gene_prediction_10.C2
    | Gene_prediction_10.C3 ->
	[< '(Codon, joined); consume_next (consume_intragene cont) cs sin >]
    | Gene_prediction_10.Stop1 ->
	[< '(Stop, joined); consume_next cont cs sin >]
    | _ ->
	raise (Failure "Not sure what I got in intragene...")

let rec consume_notgene ((c1, _, _) as cs) sin =
  match fst c1 with
      Gene_prediction_10.NotGene 
    | Gene_prediction_10.Stop2
    | Gene_prediction_10.Stop3 ->
	let joined = join cs in
	[< '(NotGene, joined); consume_next consume_notgene cs sin >]
    | Gene_prediction_10.Start1 ->
	let joined = join cs in
	[< '(Start, joined); consume_next (consume_intragene consume_notgene) cs sin >]
    | _ ->
	raise (Failure "Waiting for NotGene | Start1 | Start2, dunno what I got...")


let create_training_data gene_boundaries fin =
  let sin = Gene_prediction_10.create_training_data gene_boundaries fin in
  match Seq.take 3 sin with
      [c1; c2; c3] ->
	consume_notgene (c1, c2, c3) sin
    | _ ->
	raise (Failure "Expecting at least 3 values")


let map_td_to_list f gb fasta = Seq.map (fun (s, v) -> (s, [v])) (f gb fasta)

let convert fin =
  let rec f =
    function
	[] ->
	  raise (Failure "An empty entry in convert should never happen")
      | x::xs as xss -> begin
	  let str = String_ext.concat "" (List.map (String_ext.make 1) xss) in
	  match Seq.next fin with
	      Some v ->
		[< 'str; f (xs @ [v]) >]
	    | None ->
		[< 'str >]
	end
  in
  f (Seq.take 3 fin)

let predict training_fname fasta_fname =
  Gene_predictor.predict training_fname fasta_fname Q0 H.train H.forward_viterbi (map_td_to_list create_training_data) convert NotGene Start



