LIBNAME = 	librrnum.a
LIB =		-L../lib
INC =		-I../include
LIB_DIR = 	../lib
RRCALC =	$(LIB_DIR)/$(LIBNAME)
AR = 		ar
include ../include_dirs.in
include ../compiler_rules.in
OBJ = 		drawpscontours.o

$(RRCALC)($(OBJ)): $(OBJ)
	$(AR) -rcs $(RRCALC) $(OBJ)
	@echo "*** Completed compiling ps routines."
	sh makeheader.sh

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.a ../include/$(LIBNAME:lib%.a=%.h) $(LIBNAME:lib%.a=%.h.auto) ../lib/$(LIBNAME)
