LIBNAME = 	librrnum.a
LIB =		-L../lib
INC =		-I../include
LIB_DIR = 	../lib
RRCALC =	$(LIB_DIR)/$(LIBNAME)
AR = 		ar
#CC =		cc
#CFLAGS =	-Wall -g $(INC) -pg
include ../compiler_rules.in
OBJ = 		rk2.o \
		rk4.o \
		ab2.o \
		euler.o \
		integrate.o \
		nr.o \
		nder.o \
		stencil.o

$(RRCALC)($(OBJ)): $(OBJ)
	$(AR) -rcs $(RRCALC) $(OBJ)
	@echo "*** Completed compiling calculus routines."
	sh makeheader.sh

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.a ../include/$(LIBNAME:lib%.a=%.h) $(LIBNAME:lib%.a=%.h.auto) ../lib/$(LIBNAME)
