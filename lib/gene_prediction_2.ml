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
      Some (Fasta.Sequence d) -> 
	Some ((oe, oe + String_ext.length d), d)
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
	    (* If there are no gene boundaries left, all of the remaining values aren't genes *)
	    [< '(NotGene, d); read_next (ctd []) e >]
	| (gs, ge)::_gbs when gs <= s && ge >= e ->
	    (* If the gene spands all of the range we have, return the gene *)
	    [< '(Gene, d); read_next (ctd gb) e >]
	| (gs, ge)::gbs when gs <= s && ge < e ->
	    (* If the gene starts at s or before s it means we are currently inside
	       the gene.  If the gene ends before the ending here we pull out the 
	       beginning and return that as a gene and recurse weith the rest *)
	    let subd = String_ext.sub d 0 (ge - s) in
	    let restd = String_ext.sub d (ge - s) (String_ext.length d - (ge - s)) in
	    [< '(Gene, subd); ctd gbs (ge, e) restd >]
	| (gs, _ge)::_gbs when e < gs ->
	    (* If teh start of the next gene is after the end of this, return as not a gene *)
	    [< '(NotGene, d); read_next (ctd gb) e >]
	| (gs, _ge)::_gbs when s <= gs && gs <= e ->
	    (* If the gene boundaries is in side the boundary that we currently have, 
	       return the first as NotGene and recursively call the next piece *)
	    let subd = String_ext.sub d 0 (gs - s) in
	    let restd = String_ext.sub d (gs - s) (String_ext.length d - (gs - s)) in
	    [< '(NotGene, subd); ctd gb (gs, e) restd >]
	| (gs, ge)::_ -> raise (Failure (Printf.sprintf "Should not be here AT ALL: s: %d e: %d gs: %d ge: %d" s e gs ge))
    with
	Invalid_argument _ ->
	  raise (Failure (Printf.sprintf "foo - s: %d e: %d gs: %d ge: %d d.length: %d d: %s" s e (fst (List.hd gb)) (snd (List.hd gb)) (String_ext.length d) d))
  in
  let gene_boundaries = remove_overlap (List.sort compare gene_boundaries) in
  read_next (ctd gene_boundaries) 0

  

(*
 * We verify that the training data looks like we expect it.  In this case,
 * the length of genes are a multiple of 3 and start with a start codon and
 * end with stop codon
 *)
let verify_training_data td =
  let match_start_codon s = List.mem s ["ATG"; "GTG"; "TTG"]  in
  let match_stop_codon s = List.mem s ["TAA"; "TAG"; "TGA"] in
  let verify_gene start_pos gene =
    let start_codon = String_ext.sub gene 0 3 in
    let stop_codon = String_ext.sub gene (String_ext.length gene - 3) 3 in
    if match_start_codon start_codon && match_stop_codon stop_codon && String_ext.length gene mod 3 = 0 then
      ()
    else
      Printf.printf "Unknown codons, start_pos: %d start: %s stop: %s length: %d sequence \"%s\"" start_pos start_codon stop_codon (String_ext.length gene) gene
  in
  let rec vtd count acc =
    let sl = String_ext.length in
    match Seq.next td with
	Some (NotGene, d) ->
	  if acc <> "" then begin
	    verify_gene (1 + count - sl acc) acc;
	    vtd (count + sl d) ""
	  end
	  else
	    vtd (count + sl d) ""
      | Some (Gene, d) ->
	  vtd (count + sl d) (acc ^ d)
      | Some (Q0, d) ->
	  vtd (count + sl d) ""
      | None ->
	  ()
  in
  vtd 0 ""
