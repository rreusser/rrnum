LIBNAME = 	librrnum.a
LIB =		-L../lib
INC =		-I../include
LIB_DIR =	../lib
RRCALC =	$(LIB_DIR)/$(LIBNAME)
AR = 		ar
include ../compiler_rules.in
OBJ = 		func.o \
		color.o \
		arrayops.o

$(RRCALC)($(OBJ)): $(OBJ)
	$(AR) -rcs $(RRCALC) $(OBJ)
	@echo "*** Completed compiling data structure routines."
	sh makeheader.sh

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.a ../include/$(LIBNAME:lib%.a=%.h) $(LIBNAME:lib%.a=%.h.auto) ../lib/$(LIBNAME)
