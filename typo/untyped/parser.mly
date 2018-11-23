/* Parser for untyped terms. */

%{
  open Hashtbl

  open Term

  let symbol = create 1000	(* hash table for variable names *)
%}

%token <int> INT
%token <string> ID
%token NUMBERSIGN BACKSLASH COLON DOT EXCLAMATION SEMICOLON
%token ARROW LPAREN RPAREN
%token LET
%token EQUAL
%token EOF
%token EQUALEQUAL
%right ARROW
%start main
%type <Term.term> main
/* %type <string * Term.term> let_binding */
%%
main:
    expr EOF			{ $1 }

expr:
    app				{ $1 }
  | BACKSLASH bind DOT expr	{ let x = find symbol $2
			 	  in remove symbol $2; Abstr (x, $4) }
/*
  | let_binding; expr           { remove symbol $1; $2 }

let_binding:
  LET bind EQUAL expr SEMICOLON  { ($2, $4) }
*/

bind: ID 	{ let x = ref $1
            in add symbol $1 x; $1 }

app:
    expr1			{ $1 }
  | app expr1			{ App ($1, $2) }

expr1:
  | ID				{ let x = find symbol $1 in Var (x) }
  | LPAREN expr RPAREN		{ $2 }

