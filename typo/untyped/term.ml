(* Lambda calculus terms. *)

open List

type varname = string ref

(* for numbering variables. FIXME make non-global? *)
(* let var_count = ref 0 *)

(* The type of terms. *)
type term =
| Var of varname
| App of term * term
| Abstr of varname * term
(* | Eq of term * term   omitting this for now *)

(* Substitution. *)
let rec subst x s t =
  match t with
    | Var (v) -> if (x == v) then s else t
    | App (rator, rand) ->
      App (subst x s rator, subst x s rand)
    | Abstr (v, body) -> Abstr (v, subst x s body)

(* Create numbered variables, starting with "a","b","c"..."z",
  then switching to "x0", "x1", "x2", "x3"... *)
let number_to_var i =
  (* ??? why does numbering start at 1? *)
  if i <= 26
  then String.make 1 (char_of_int (i + 96))
  else "x" ^ (string_of_int (i-26))

(* Utility to rename all variables in a term to "fresh
  variable names". Akin to Prolog's copy_term/2.
  Kind of inefficient. *)
let copy_term a = 
let var_count = ref 0 in
let rec copy_term_1 a =
  match a with
  | Var (v) -> a
  | App (r, s) -> App (copy_term_1 r, copy_term_1 s)
  | Abstr (x, s) ->
    var_count := !var_count + 1;
    let x1 = ref (number_to_var !var_count) in
    Abstr (x1, subst x (Var x1) (copy_term_1 s))
in copy_term_1 a

(* Check for alpha-equivalence (slow but conceptually
  simple implementation.) *)
let alpha_equiv a b = ((copy_term a) = (copy_term b))

(* Alpha equivalence (deprecated).
let rec alpha_equiv_old a b =
  match (a, b) with
    | (Var (x1), Var (x2)) -> x1 == x2
    | (App (a1,b1), App (a2,b2)) ->
      alpha_equiv a1 a2 && alpha_equiv b1 b2
    | (Abstr (x1,b1), Abstr (x2,b2)) ->
      let fresh = Var (ref "") in
			alpha_equiv (subst x1 fresh b1) (subst x2 fresh b2)
    | _ -> false
*)

(* converting terms to strings, and printing them
let term_to_string x _ =
let var_count = ref 0 in
let rec term_to_string_orig x _ =
match x with
| Var v -> !v
| Abstr (x, a) ->
  var_count := !var_count + 1;
  let x1 = ref ("x" ^ (string_of_int (!var_count))) in
	let a1 = subst x (Var x1) a in
  "(\\" ^ !x1 ^ "." ^ (term_to_string_1 a1 ()) ^ ")"
| App (a, b) -> "(" ^ (term_to_string_1 a ()) ^
  " " ^ (term_to_string_1 b ()) ^ ")"
in term_to_string_1 x () *)

(* Converts a term to a string, without renaming variables
  (probably term_to_string() is more useful.)
  This doesn't stint on parentheses, in the name of
  clarity. *)
let rec term_to_string_1 x =
match x with
| Var v -> !v
| Abstr (x, a) ->
  "(\\" ^ !x ^ "." ^ (term_to_string_1 a) ^ ")"
| App (a, b) ->
  "(" ^ (term_to_string_1 a) ^ " " ^ (term_to_string_1 b) ^ ")"

let print_term x =
	print_string (term_to_string_1 (copy_term x))

