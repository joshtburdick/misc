(* Reducing terms and/or normalizing them. *)

open Term

(* Reduce a term, possibly more than once. 
   Only does beta-conversion, not eta. *)
let rec reduce_step a =
  match a with
    (* XXX hopefully this is leftmost-innermost? *)
    | App (Abstr (x, body), rand) ->
      let b = reduce_step body in
      if (body = b)
      then subst x (copy_term rand) body
      else App (Abstr (x, b), rand)
    | App (rator, rand) ->
      App (reduce_step rator, reduce_step rand)
    | Abstr (x, body) -> Abstr (x, reduce_step body)
    | _ -> a

let is_nf x = alpha_equiv x (reduce_step x)

(* reduces some number of times, up to a limit.
   (reduce() can reduce more than one thing, so this is
   not an exact reduction count.  But increasing this
   increases the number of reductions which happen.) *)
let rec reduce n x = 
	if n <= 0
	then x
	else reduce (n-1) (reduce_step x)

(* same as above, but prints reducts as they happen. *)
let rec reduce_and_print n x =
	if (n <= 0 || is_nf x)
	then (print_term x; print_string "\n")
	else
		(print_term x;
		print_string "\n";
		(reduce_and_print (n-1) (reduce_step x)))

(* Check if something is in normal form. *)
let is_nf x = alpha_equiv x (reduce_step x)

(** (approximate) normal form *)
let nf x = reduce 50 x

(** "approximately normalizable" *)
let is_normalizable x = is_nf (nf x)

(** Tests whether x = rx (that is, whether
  "x is in the retract r" *)
let is_in_retract r x =
let x1 = nf x in alpha_equiv x1 (nf (App (r, x)))

(*
alpha_equiv x1 (nf (App r x1))
*)

