(*
 * Functions for prediction genes
 *)

let min_gene_length = 500

let genes ng sg l = 
  Seq.to_list (Seq.map 
		 (fun s -> 
		    if s <> sg && s <> ng then
		      sg
		    else
		      s)
		 (Seq.of_list l))

let predict training_fname fasta_fname start_end train viterbi create_training_data fasta_converter ng sg =
  let fasta = Fasta.read_file ~chunk_size:10000 fasta_fname in
  (* Limit out training data to those ORFs that are min_gene_length+ nt in length *)
  let gene_boundaries = List.filter (fun (s, e) -> e - s > min_gene_length) (Genbank.read training_fname) in
  let training_data = create_training_data gene_boundaries fasta in
  let hmm = train start_end training_data in
  let (total, path, prob) = viterbi (fasta_converter (Fasta.to_seq (Fasta.read_file ~chunk_size:10000 fasta_fname))) hmm in
  let gene_list = List.map (fun (_, b) -> b) (List.filter (fun (s, _) -> s <> ng) (Misc_lib.boundaries (genes ng sg path))) in
  let count = List.length gene_list in
  let training_count = List.length gene_boundaries in
  ((total, path, prob), hmm, gene_list, count, training_count)


let write_gene_list gene_list fout =
  List.iter (fun (s, e) -> Printf.fprintf fout "%d\t%d\n" s e) gene_list

