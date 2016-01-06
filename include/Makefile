INST_DIR =	/usr/local/include
HEADERS =	rrcalc.h \
		rrfileio.h \
		rrgeneral.h \
		rrdata.h \
		rrgrid.h \
		rropengl.h \
		rrps.h \
		rrlinalg.h

all: $(HEADERS)
	@echo "Nothing to be done for include/"

install: $(LIBS)
	@echo "*** Installing headers in /usr/local/include"
	for i in $(HEADERS); do \
	( cp $$i $(INST_DIR)/$$i ) ; \
	done
	@echo "*** Successfully installed libraries"

clean:
	rm -f $(HEADERS) *.h
