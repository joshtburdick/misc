(* Iterates over all terms of a certain size. *)

open Term
open Reduce

(* Computes a list of all lambda terms up to a given size.
  XXX uses much memory; using term_iter is preferred. *)
(*
let rec all_terms vars depth =
if depth <= 0
then []
else
	let var = map (fun x -> Var x) vars in
	let abs =
		let x = ref ("a" ^ (string_of_int depth)) in
		map (fun a -> Abstr (x, a)) (all_terms (x::vars) (depth-1)) in
	let app =
		let t = all_terms vars (depth-1) in
		concat (map (fun a -> (map (fun b -> App (a, b)) t)) t) in
	concat [var; abs; app]
*)

(* Iterates over all terms up to a given size. *)
let term_iter f depth =
let rec term_iter_1 vars depth f =
	if depth <= 0
	then ()
	else
		((List.iter (fun x -> f (Var x)) vars);
		(let x = ref "a" in
			(term_iter_1 (x::vars) (depth - 1) (fun a -> f (Abstr (x, a)))));
		(term_iter_1 vars (depth - 1)
			(fun a ->
				term_iter_1 vars (depth - 1)
					(fun b -> f (App (a, b)))))) in
term_iter_1 [] depth f

(** (approximate) normal form *)
let nf x = reduce 50 x

let is_nf x = alpha_equiv x (nf x)

(** "approximately normalizable" *)
let is_normy x = is_nf (reduce 50 x)

(** Tests whether x is in the "retract set" r *)
let is_in_retract r z =
let x1 = nf z in alpha_equiv x1 (nf (App (r, z)))
(*
alpha_equiv x1 (nf (App r x1))
*)

(** supply of fresh variables *)
let x1 = ref "x1"
let x2 = ref "x2"
let x3 = ref "x3"
let x4 = ref "x4"
let x5 = ref "x5"
let x6 = ref "x6"
let x7 = ref "x7"
let x8 = ref "x8"

(** useful function *)
let i = Abstr (x1, Var x1)

(** this seems like it would only "contain" i *)
let a = ref "a"
let si = Abstr(a, App (Var (a), i))

(** checks whether x = x;x *)
let b = ref "b"
let x1 = ref "x1"
let r x =
  let f = Abstr(b, App(Abstr(x1, App(x, App(x, Var(b)))), x)) in
	  ((is_nf x) && (alpha_equiv (nf x) (nf f)))

let search_for_inhabitants _ =
	let process_term x =
		(if (is_in_retract si x) then
		(print_term x; print_string "\n\n")) in
	term_iter process_term 4

(*
let search_for_inhabitants _ =
	let terms = all_terms [] 4 in
	let in_r = filter (is_in_retract si) terms in
  let f x = print_term x;
		print_string " --> ";
		print_term (nf x);
		print_string "\n\n" in *)
(*    var_count := 0;
    (print_string ((term_to_string x) ^ " --> " ^
		  (term_to_string (nf x)) ^ "\n\n")) in *)
    (*           List.iter f in_r *)
(*	List.iter (fun x -> print_term x; print_string "\n") terms *)

(*
let _ = search_for_inhabitants si
*)

(*
let _ = reduce_and_print 20 (App (si,si))
*)

(* test by printing all terms up to some size *)
(*
let _ =
	term_iter (fun x -> print_term x; print_string "\n") 3
*)
(*
let _ = (print_term si); print_string "\n"
*)
(* FIXME alpha-equiv is broken, apparently when comparing terms which
  contain the same actual reference. This may be fixed... *)
(*
let _ =
  let x1 = ref "x1" in
  let y1 = ref "y1" in
  let x2 = ref "x2" in
  let y2 = ref "y2" in 
  let a = Abstr (x1, Abstr (y1, Var y1)) in
  let b = Abstr (x2, Abstr (y2, Var y2)) in
  let r = alpha_equiv a b in
  print_string (if r then "true\n" else "false\n")
*)

(* main: various tests *)
(* prints terms *)
(* let _ = term_iter (fun x -> print_term x; print_string "\n") 5 *)

(* print terms, and reduce them *)
let _ = term_iter (fun x -> reduce_and_print 10 x; print_string "\n" ) 6

