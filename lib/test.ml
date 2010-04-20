let list_of_string s = Seq.to_list (Seq.of_string s)

(* type states = Q0 | Q1 | Q2 *)

(* let transitions = [ (Q0, [ (Q1, 1.0) ]) *)
(* 		  ; (Q1, [ (Q1, 0.8); (Q2, 0.15); (Q0, 0.05) ]) *)
(* 		  ; (Q2, [ (Q2, 0.7); (Q1, 0.3) ]) *)
(* 		  ] *)


(* let emissions = [ (Q1, [ ('Y', 1.0) ]) *)
(* 		; (Q2, [ ('R', 1.0) ]) *)
(* 		] *)

(* type states = Q0 | Rainy | Sunny *)

(* type alphabet = Walk | Shop | Clean *)

(* let transitions = [ (Q0, [(Rainy, 0.6); (Sunny, 0.4)]) *)
(* 		  ; (Rainy, [(Rainy, 0.7); (Sunny, 0.3)]) *)
(* 		  ; (Sunny, [(Rainy, 0.4); (Sunny, 0.6)]) *)
(* 		  ] *)


(* let emissions = [ (Rainy, [(Walk, 0.1); (Shop, 0.4); (Clean, 0.5)]) *)
(* 		; (Sunny, [(Walk, 0.6); (Shop, 0.3); (Clean, 0.1)]) *)
(* 		] *)


(* let observations = [ Walk; Shop; Clean ] *)


(* module H = Hmm.Make(struct type s = states type a = alphabet let compare = compare end) *)

(* let hmm = H.make transitions emissions Q0 *)

(* let fv = H.forward_viterbi (Seq.of_list observations) hmm *)

(* let () = let _ = fv in () *)


(* type states = Q0 | One | Two *)

(* module H = Hmm.Make(struct type s = states type a = char let compare = compare end) *)

(* let training_data = [ (One, list_of_string "CGATATT") *)
(* 		    ; (Two, list_of_string "CGATTCT") *)
(* 		    ; (One, list_of_string "ACGCGC") *)
(* 		    ; (Two, list_of_string "GTAT") *)
(* 		    ; (One, list_of_string "ACTAGCT") *)
(* 		    ; (Two, list_of_string "TATC") *)
(* 		    ; (One, list_of_string "TGATC") *)
(* 		    ] *)


(* let training_data = [ (One, list_of_string "CGAT") *)
(* 		    ; (One, list_of_string "ATT") *)
(* 		    ; (Two, list_of_string "CGATTCT") *)
(* 		    ; (One, list_of_string "ACGCGC") *)
(* 		    ; (Two, list_of_string "GTAT") *)
(* 		    ; (One, list_of_string "ACTAGCT") *)
(* 		    ; (Two, list_of_string "TATC") *)
(* 		    ; (One, list_of_string "TGATC") *)
(* 		    ] *)


module H = Hmm.Make(struct type s = Gene_prediction_4.coding_state type a = char let compare = compare end)

let fasta = Fasta.read_file ~chunk_size:10000 "../datasets/E.coli.O103.H2_str.12009.fasta"
let gene_boundaries = Genbank.read "../datasets/E.coli.O103.H2_str.12009.gb"
let td = Seq.to_list (Gene_prediction_4.create_training_data gene_boundaries fasta)
let training_data = Seq.map (fun (s, v) -> (s, list_of_string v)) (Seq.of_list td)
let hmm = H.train Gene_prediction_4.Q0 training_data
let viterbi_results = H.forward_viterbi (Fasta.to_seq (Fasta.read_file ~chunk_size:10000 "../datasets/E.coli.O103.H2_str.12009.fasta")) hmm


