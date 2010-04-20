
let uniq l =
  List.rev (List.fold_left 
	      (fun acc e -> 
		 match acc with
		     [] -> [e]
		   | x::_xs when e = x -> acc
		   | x::_xs -> e::acc)
	      []
	      l)
					  

