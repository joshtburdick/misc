
SRC = lexer.ml parser.ml term.ml reduce.ml term_iter.ml

# OBJS=term.cmx lexer.cmx parser.cmx reduce.cmx \
#       type.cmx pprint.cmx search.cmx prove.cmx main.cmx

lc : $(SRC)
	ocamlopt -o lc $(SRC)

term.cmx:
	ocamlopt -c term.ml
lexer.ml: lexer.mll parser.ml term.cmx
	ocamllex lexer.mll
parser.ml: parser.mly term.cmx
	ocamlyacc parser.mly;
	ocamlopt -c parser.mli

clean :
	rm lc lexer.ml parser.ml *.o *.cmi *.cmx

