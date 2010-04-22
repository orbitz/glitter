
(* Very expensive function for what it does, but easy implementation *)
let list_of_string s = Seq.to_list (Seq.of_string s)

let uniq l =
  List.rev (List.fold_left 
	      (fun acc e -> 
		 match acc with
		     [] -> [e]
		   | x::_xs when e = x -> acc
		   | x::_xs -> e::acc)
	      []
	      l)
					  

let boundaries l =
  List.rev (Seq.fold_left
	      (fun acc (idx, e) ->
		 match acc with
		     [] ->
		       [(e, (idx, idx))]
		   | (s, (x, y))::accs when s = e ->
		       (s, (x, idx))::accs
		   | (s, (x, y))::accs ->
		       (e, (idx, idx))::(s, (x, idx))::accs)
	      []
	      (Seq.enumerate (Seq.of_list l)))


let rec zip l1 l2 =
  match (l1, l2) with
      (x1::x1s, x2::x2s) ->
	(x1, x2)::zip x1s x2s
    | _ ->
	[]

