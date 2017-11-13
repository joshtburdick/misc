(* Lambda calculus, minimalistically.
20080308, jtb: started
20090111, jtb: am having doubts about this idea.
20090614, jtb: still have doubts about it.  However, iterating over
	terms seems to work.  Possibly just reduction is broken?
20090627, jtb: well, this code may or may not be broken.  However,
  my conjecture about \x.xI only being a retract containing I is counterexampled,
	since (\xyz.z)I = \yz.I = \yzw.w.   ??? wtf you sure?
(\xyz.z) I = \yz.z
	Um, d'oh...
20110103, jtb: alpha_equiv is definitely broken.
20110628, jtb: alpha_equiv bug appears fixed.
20161110, jtb: splitting this into *)

open List

type varname = string ref

(* for numbering variables *)
let var_count = ref 0

type term =
| Var of varname
| App of term * term
| Abstr of varname * term

(* Substitution. *)
let rec subst x s t =
  match t with
    | Var (v) -> if (x == v) then s else t
    | App (rator, rand) -> App (subst x s rator, subst x s rand)
(*    | Abstr (v, body) -> Abstr (v, subst x s body) *)
    | Abstr (v, body) ->
     let fresh = ref "" in
      Abstr (fresh, subst x s body)  *)

(* alpha equivalence *)
(* XXX now using a (hopefully) totally "fresh" variable *)
and alpha_equiv a b =
  match (a, b) with
    | (Var (x1), Var (x2)) -> x1 == x2
    | (App (a1,b1), App (a2,b2)) ->
      alpha_equiv a1 a2 && alpha_equiv b1 b2
    | (Abstr (x1,b1), Abstr (x2,b2)) ->
      let fresh = Var (ref "") in
			alpha_equiv (subst x1 fresh b1) (subst x2 fresh b2)
    | _ -> false

(* converting terms to strings, and printing them *)
let rec term_to_string x _ =
match x with
| Var v -> !v
| Abstr (x, a) ->
  var_count := !var_count + 1;
  let x1 = ref ("x" ^ (string_of_int (!var_count))) in
	let a1 = subst x (Var x1) a in
  "(\\" ^ !x1 ^ "." ^ (term_to_string a1 ()) ^ ")"
| App (a, b) -> "(" ^ (term_to_string a ()) ^ " " ^ (term_to_string b ()) ^ ")"

let print_term x =
	var_count := 0;
	print_string (term_to_string x ())


(* Reduce a term, generally more than once. 
   Only does beta-conversion, not eta. *)
let rec reduce_step a =
  match a with
    | App (Abstr (x, body), rand) -> subst x rand body
(*    | App (Abstr (x, body), rand) -> *)
    | App (rator, rand) -> App (reduce_step rator, reduce_step rand)
    | Abstr (x, body) -> Abstr (x, reduce_step body)
    | _ -> a

let is_nf x = alpha_equiv x (reduce_step x)

(* reduces some number of times, up to a limit.
   (reduce() can reduce more than one thing, so this is not
   an exact reduction count.  But increasing this strictly
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


(* Computes a list of all lambda terms up to a given size.
  XXX uses much memory; using term_iter is preferred. *)
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
	term_iter process_term 5

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
let _ =
	term_iter (fun x -> print_term x; print_string "\n") 1

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

