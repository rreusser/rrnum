INST_DIR =	/usr/local/lib
LIBS = 		librrnum.a
all: $(LIBS)
	@echo "*** Libraries compiled.  Type 'sudo make install' to install in /usr/lib."

install: $(LIBS)
	@echo "*** Installing libraries in /usr/local/lib"
	for i in $(LIBS); do \
	( ranlib $$i; cp $$i $(INST_DIR)/$$i ) ; \
	done
	@echo "*** Successfully installed libraries"
clean:
	rm -f $(LIBS) *.a
