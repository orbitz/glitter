(*
 * This implements ORF-bounding, a term I made up.  It means taking the boundaries of a
 * predicted gene and looking for ORFs in an n-nucleotide surrounding area of it for exactly
 * matching start/stop codons.
 *)


let common_start_codons = ["ATG"; "GTG"; "TTG"]
let common_stop_codons = ["TAA"; "TAG"; "TGA"]

let read_by_codons fin =
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


let orf_bounds_from_stream starts stops sin =
  let f (start, stop) (idx, codon) =
    if List.mem codon starts then
      (idx::start, stop)
    else if List.mem codon stops then
      (start, idx::stop)
    else
      (start, stop)
  in
  let (start, stop) = Seq.fold_left f ([], []) (Seq.enumerate sin) in
  (* This is not necessary but just doing it for visual debugging purposes *)
  (List.rev start, List.rev stop)


let create_orf_bounds starts stops fasta_fname =
  let sin = read_by_codons (Fasta.to_seq (Fasta.read_file ~chunk_size:10000 fasta_fname)) in
  orf_bounds_from_stream starts stops sin


(*
 * This takes a tuple containg the start/stop of a predicted gene
 * and a tuple from the output of create_orf_bounds and a cut off for
 * how far in both directions of the gene to look for bounds.
 * If no acceptable bound can be found in the vicinity, the gene
 * is returned as is
 * 
 * An acceptable bounds is defined as being close to the gene start
 * and stop given (on the gene start and stop is preferred) and
 * the length must be a multiple of 3
 *)
let find_nearest_bound cutoff (starts, stops) (gs, ge) =
  let pstarts = List.filter (fun v -> v > gs - cutoff && v < gs + cutoff) starts in
  let pstops = List.filter (fun v -> v > ge - cutoff && v < ge + cutoff) stops in
  let local_orfs = 
    List.fold_left 
      (fun acc s -> 
	 List.fold_left (fun acc e -> 
			   if e - s mod 3 == 0 then
			     (s, e)::acc
			   else
			     acc) acc pstops)
      []
      pstarts
  in
  let compare (s1, e1) (s2, e2) = 
    let ldiff = abs (s1 - gs) + abs (e1 - ge) in
    let rdiff =  abs (s2 - gs) + abs (e2 - ge) in
    if ldiff < rdiff then
      -1
    else if ldiff = rdiff then
      0
    else
      1
  in
  if local_orfs = [] then
    (gs, ge)
  else
    List.hd (List.sort compare local_orfs) 
    

let test fname = 
  let _ = create_orf_bounds common_start_codons common_stop_codons fname in
  if find_nearest_bound 20 ([1; 100], [10; 200]) (0, 9) <> (0, 9) then
    raise (Failure "Found wrong nearest neighbor...")
  else
    ()
