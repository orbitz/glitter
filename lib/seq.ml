
let next s = 
  try
    Some (Stream.next s)
  with 
      Stream.Failure -> None


let rec of_list =
  function
      [] -> 
	[< >]
    | x::xs -> 
	[< 'x; of_list xs >]
  
let of_string s =
  let rec of_string' idx =
    if idx < String.length s then
      let c = s.[idx] in
      [< 'c; of_string' (idx + 1) >]
    else
      [< >]
  in
  of_string' 0


let rec fold_left f a s =
  match next s with
      Some v ->
	fold_left f (f a v) s
    | None -> a


let to_list s = List.rev (fold_left (fun a e -> e::a) [] s)

(*
 * Takes a function and an initial value and returns a stream incrementally
 * calling the function
 *)
let rec iterate f v = 
  let fv = f v in
  [< 'fv; iterate f fv >]

let take n s =
  let rec take' acc = function
      0 -> List.rev acc
    | n -> 
	match next s with
	    Some e -> take' (e::acc) (n - 1)
	  | None -> List.rev acc
  in
  take' [] n


(*
 * Given a sequene will iterate that sequence and for each
 * element return a tuple of (intvalue, seq)
 *)
let rec enumerate ?(start = 0) ?(step = 1) seq =
  match next seq with
      Some v ->
	let d = (start, v) in
	[< 'd; enumerate ~start:(start + step) ~step:step seq >]
    | None ->
	[< >]

let rec map f s =
  match next s with
      Some v ->
	[< 'f v; map f s >]
    | None ->
	[< >]
