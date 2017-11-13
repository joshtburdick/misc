
(* Lexer. *)

{
open Parser
exception Eof
}

rule token = parse
  [' ''\t''\n']	{ token lexbuf }
| "\\"		{ BACKSLASH }
| ";"           { SEMICOLON }
| "."		{ DOT }
| "=="		{ EQUALEQUAL }
| "("		{ LPAREN }
| ")"		{ RPAREN }
| "="           { EQUAL }
| "let"         { LET }
| ['0'-'9']+	{ INT (int_of_string (Lexing.lexeme lexbuf)) }
| ['A'-'Z''a'-'z']['A'-'Z''a'-'z''0'-'9''_']*	{ ID (Lexing.lexeme lexbuf) }
| eof		{ EOF }

