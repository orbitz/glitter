(* 2 state HMM *)

type coding_state = Q0 | NotGene | Gene

(*
 * If the next to read in the fasta file is a header it returns
 * Some ((0, 0), header value)
 * If a sequence
 * Some ((s, e), sequence data)
 * If the sequence is done
 * None
 *)
let read_sequence_with_boundaries oe fin =
  match Seq.next fin with
      Some (Fasta.Sequence d) -> Some ((oe + 1, oe + String_ext.length d), d)
    | Some (Fasta.Header h) -> Some ((0, 0), h)
    | None -> None

let remove_overlap l =
  let rec ro' acc =
    function
	[] ->
	  acc
      |	(_gs1, ge1)::(gs2, _ge2)::gbs when ge1 >= gs2 ->
	  ro' acc gbs
      | xs::gbs ->
	  ro' (xs::acc) gbs
  in
  List.rev (ro' [] l)


(*
 * Returns a stream of training data
 *)
let create_training_data gene_boundaries fin =
  let rec read_next f oe =
    match read_sequence_with_boundaries oe fin with
	Some ((0, 0), _h) -> read_next f 0
      | Some ((s, e), d) -> f (s, e) d
      | None -> [< >]
  in
  let rec ctd gb (s, e) d =
    try
      match gb with
	  [] ->
	    [< '(NotGene, d); read_next (ctd []) e >]
	| (gs, _ge)::_gbs when s < gs && gs <= e ->
	    let subd = String_ext.sub d 0 (gs - s) in
	    let restd = String_ext.sub d (gs - s + 1) (String_ext.length d - (gs - s + 1)) in
	    [< '(NotGene, subd); ctd gb (gs + 1, e) restd >]
	| (gs, ge)::gbs when gs <= s && ge < e ->
	    let subd = String_ext.sub d 0 (ge - s) in
	    let restd = String_ext.sub d (ge - s + 1) (String_ext.length d - (ge - s + 1)) in
	    [< '(Gene, subd); ctd gbs (ge + 1, e) restd >]
	| (gs, ge)::gbs when s = gs && ge = e ->
	    [< '(Gene, d); read_next (ctd gbs) e >]
	| (gs, ge)::_gbs when s = gs && e < ge ->
	    [< '(Gene, d); read_next (ctd gb) e >]
	| (gs, _ge)::_gbs when e < gs ->
	    [< '(NotGene, d); read_next (ctd gb) e >]
	| (gs, ge)::_gbs when gs < s && e < ge ->
	    [< '(Gene, d); read_next (ctd gb) e >]
	| (gs, ge)::_ -> raise (Failure (Printf.sprintf "Should not be here AT ALL: s: %d e: %d gs: %d ge: %d" s e gs ge))
    with
	Invalid_argument _ ->
	  raise (Failure (Printf.sprintf "foo - s: %d e: %d gs: %d ge: %d d.length: %d d: %s" s e (fst (List.hd gb)) (snd (List.hd gb)) (String_ext.length d) d))
  in
  let gene_boundaries = remove_overlap (List.sort compare gene_boundaries) in
  read_next (ctd gene_boundaries) 0

  
