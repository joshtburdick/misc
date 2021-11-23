type token =
  | INT of (int)
  | ID of (string)
  | NUMBERSIGN
  | BACKSLASH
  | COLON
  | DOT
  | EXCLAMATION
  | SEMICOLON
  | ARROW
  | LPAREN
  | RPAREN
  | LET
  | EQUAL
  | EOF
  | EQUALEQUAL

val main :
  (Lexing.lexbuf  -> token) -> Lexing.lexbuf -> Term.term
