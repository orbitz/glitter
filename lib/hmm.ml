(*
 * Tools to construct and HMM and perform calculations on it
 *)


module type Ordered_type =
  sig
    type s
    type a
    val compare: s -> s -> int
  end

module Make =
  functor (Elt : Ordered_type) ->
    struct
      type state = Elt.s
      type alphabet = Elt.a
      type probability = float
      type transition = state * (state * probability) list
      type emission = state * (alphabet * probability) list

      module StateMap = Map.Make(struct type t = state let compare = Elt.compare end)
	  
      type hmm_state = { transitions : transition list
		       ; emissions : emission list
		       ; start_end : state
		       }

      (* misc function for calculating a random weight *)
      let random_by_weight l =
	let rec get_random a v =
	  function
	      [] -> raise (Failure "You didn't represent all your base....")
	    | (s, x)::xs -> 
		let a = a + int_of_float (x *. 1000.0) in
		if v <= a then
		  s
		else
		  get_random a v xs
	in
	get_random 0 (Random.int 1001) l


      (* a misc funciton to work on lists *)
      let rec drop_while f =
	function
	    [] -> []
	  | x::xs when f x -> x::xs
	  | _::xs -> drop_while f xs
	      
      let get_transitions hmm s =
	List.filter (fun (s, p) -> p > 0.0) (List.assoc s hmm.transitions)

      let get_transition hmm s1 s2 =
	List.assoc s2 (get_transitions hmm s1)

      let get_transition_def hmm s1 s2 def =
	try
	  get_transition hmm s1 s2
	with
	    Not_found -> def

      let get_states_by_emission hmm e =
	List.map fst (List.filter (fun (s, es) -> List.mem_assoc e es && List.assoc e es > 0.0) hmm.emissions)

      let get_emissions hmm s =
	List.filter (fun (s, p) -> p > 0.0) (List.assoc s hmm.emissions)

      let get_emission hmm s1 s2 =
	List.assoc s2 (get_emissions hmm s1)

      let get_emission_def hmm s1 s2 def =
	List.assoc s2 (get_emissions hmm s1)
	  
      let get_random_transition hmm s =
	random_by_weight (get_transitions hmm s)
	  
      let get_random_emission hmm s =
	random_by_weight (get_emissions hmm s)

      let get_states hmm =
	List.filter ((<>) hmm.start_end) (List.map fst hmm.transitions)
	  
      let make t e se =
	{ transitions = t
	; emissions = e
	; start_end = se
	}
	  
	  
      let sample_emission hmm =
	let rec samp acc s =
	  let ns = get_random_transition hmm s in
	  if ns = hmm.start_end then
	    acc
	  else
	    samp ((get_random_emission hmm ns)::acc) ns
	in
	samp [] hmm.start_end
	  
      let sample_emission_seq hmm =
	let rec samp s =
	  let ns = get_random_transition hmm s in
	  if ns = hmm.start_end then
	    [< >]
	  else
	    let o = get_random_emission hmm ns in
	    [< 'o; samp ns >]
	in
	samp hmm.start_end
	  
      let prob_of_sequence hmm seq =
	let rec p_of_s s =
	  function
	      [] ->
		snd (List.hd (List.filter (fun (s, _) -> s = hmm.start_end) (get_transitions hmm s)))
	    | x::xs ->
		let trans = get_transitions hmm s in
		let (ts, tp) =
		  List.hd
		    (drop_while
		       (fun s -> List.mem_assoc x
			  (get_emissions hmm (fst s)))
		       trans)
		in
		let ep = List.assoc x (get_emissions hmm ts) in
		tp *. ep *. p_of_s ts xs
	in
	p_of_s hmm.start_end seq

      (* 
       * Calculates viterbi path + probabilities
       * returns (total probability, path, probability of path)
       *)
      let forward_viterbi obs hmm =
	let logsum x y = x +. log (1. +. exp (y -. x))
	in
	let calc_max t ob ns (total, argmax, valmax) s =
	  let (prob, v_path, v_prob) = StateMap.find s t in
	  let p = log (get_emission hmm s ob) +. log (get_transition hmm s ns) in
	  let prob = prob +. p in
	  let v_prob = v_prob +. p in
	  let total = if total = neg_infinity then prob else logsum total prob in
	  if v_prob > valmax then
	    (total, ns :: v_path, v_prob)
	  else
	    (total, argmax, valmax)
	in
	let accum_stats t ob u ns =
	  let states = get_states hmm in
	  StateMap.add ns (List.fold_left (calc_max t ob ns) (neg_infinity, [], neg_infinity) states) u
	in
	let walk_states t ob =
	  List.fold_left (accum_stats t ob) StateMap.empty (get_states hmm)
	in
	let walk_obs t =
	  Seq.fold_left walk_states t obs
	in
	let make_t =
	  List.fold_left
	    (fun t s ->
	       let start_p = log (get_transition_def hmm hmm.start_end s 0.0) in
	       StateMap.add s (start_p, [s], start_p) t)
	    StateMap.empty
	    (get_states hmm)
	in
	let find_max t (total, argmax, valmax) s =
	  let (prob, v_path, v_prob) = StateMap.find s t in
	  let total = logsum total prob in
	  if v_prob > valmax then
	    (total, v_path, v_prob)
	  else
	    (total, argmax, valmax)
	in
	let t = walk_obs make_t in
	let states = get_states hmm in
	let (total, v_path, v_prob) = 
	  List.fold_left 
	    (find_max t) 
	    (StateMap.find (List.hd states) t)
	    (List.tl states) 
	in
	(exp total, List.rev v_path, exp v_prob)

      (*
       * Train an HMM given an input sequence.  Sequence is lazy
       * this also needs to know the start_end state for training
       *)
      type hmm_builder = { trans : (state * int32) list StateMap.t
			 ; emiss : (alphabet * int32) list StateMap.t
			 ; se : state
			 }

      let update_m m ss vk =
	try
	  let v = StateMap.find ss m in
	  try
	    let c = List.assoc vk v in
	    StateMap.add ss ((vk, Int32.succ c)::(List.remove_assoc vk v)) m
	  with
	      Not_found ->
		StateMap.add ss ((vk, Int32.of_int 1)::v) m
	with
	    Not_found ->
	      StateMap.add ss [(vk, Int32.of_int 1)] m


      let train start_end seq =
	let rec train' ps cs hmm_b =
	  function
	      [] ->
		begin
		  match Seq.next seq with
		      Some (s, d) when s = cs ->
			train' ps cs hmm_b d
		    | Some (s, d) ->
			train' cs s hmm_b d
		    | None ->
			{hmm_b with trans = update_m hmm_b.trans cs start_end}
		end
	    | x::xs ->
		train'
		  cs
		  cs 
		  {hmm_b with
		     trans = update_m hmm_b.trans ps cs
		     ; emiss = update_m hmm_b.emiss cs x
		  }
		  xs
	in
	let make_probabilities trans =
	  StateMap.fold
	    (fun k v a ->
	       let sum = List.fold_left (Int32.add) (Int32.of_int 0) (List.map snd v) in
	       (k, (List.map (fun (s, v) -> (s, Int32.to_float v /. Int32.to_float sum)) v))::a)
	    trans
	    []
	in
	let hmm_b = train' start_end start_end {trans = StateMap.empty; emiss = StateMap.empty; se = start_end} [] in
	make (make_probabilities hmm_b.trans) (make_probabilities hmm_b.emiss) start_end
	  

    end
	    