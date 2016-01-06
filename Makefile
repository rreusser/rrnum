OBJ_DIR =	lib
CALC_DIR =	calc
ALG_DIR =	linalg
INC_DIR = 	include
DAT_DIR =	data
FIO_DIR =	fileio
GEN_DIR =	general
OGL_DIR =	opengl
GRID_DIR =	grid
PS_DIR =	postscript
LIB_INST_DIR =	/usr/local/lib
INC_INST_DIR =  /usr/local/include
SRCDIRS =	$(GEN_DIR) $(DAT_DIR) $(CALC_DIR) $(GRID_DIR) \
		$(FIO_DIR) $(ALG_DIR) $(OGL_DIR) $(PS_DIR)
OUTDIRS =	$(OBJ_DIR) $(INC_DIR)

# Default Target:

all: 
	@echo
	@echo "##### BUILDING RRNUM LIBRARY #####"
	@echo "##### Ricky Reusser, 2006"
	@echo "##### rreusser@iastate.edu"
	@echo
	BUILD=0 ; \
	for i in $(SRCDIRS); do \
	    $(MAKE) -w -C $$i ;  \
	    SUCCESS=$$? ; \
	    BUILD="$$BUILD && $$SUCCESS" ; \
	    if [ $$SUCCESS -eq 0 ]; then \
		echo "Successful build in $$i/" ; echo ; \
	    else \
		echo "Error compiling in $$i/" ; echo ; exit ; \
	    fi \
	done
	@echo
	@echo "##### DONE COMPILING..."
	@echo "##### type 'sudo make install' to install in /usr/local"
	@echo

# Rules to build individual targets

grid:
	@cd $(GRID_DIR) ; make

calc:
	@cd $(CALC_DIR) ; make

linalg:
	@cd $(ALG_DIR) ; make

header:
	@cd $(INC_DIR) ; make

fileio:
	@cd $(FIO_DIR) ; make

data:
	@cd $(DAT_DIR) ; make

general:
	@cd $(GEN_DIR) ; make

opengl:
	@cd $(OGL_DIR) ; make

postscript:
	@cd $(PS_DIR) ; make


clean:
	for i in $(SRCDIRS) ; do \
	( cd $$i ; make clean ) ; \
	done
	for i in $(OUTDIRS) ; do \
	( cd $$i ; make clean ) ; \
	done

install: 
	for i in $(OUTDIRS) ; do \
	( cd $$i ; make install ) ; \
	done

