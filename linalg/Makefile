LIBNAME = 	librrnum.a
INC =		-I../include
LIB_DIR =	../lib
RRLINALG =	$(LIB_DIR)/$(LIBNAME)
AR = 		ar
#CFLAGS =	-Wall -pg -g
#CC =		cc $(INC)
include ../compiler_rules.in
OBJ = 		tdiag.o sparse.o bicg.o

$(RRLINALG)($(OBJ)): $(OBJ)
	$(AR) rcs $(RRLINALG) $(OBJ)
	@echo "*** Completed compiling linear algebra routines."
	sh makeheader.sh

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.a $(LIBNAME:lib%.a=%.h.auto) ../lib/$(LIBNAME) ../include/$(LIBNAME:lib%.a=%.h)
