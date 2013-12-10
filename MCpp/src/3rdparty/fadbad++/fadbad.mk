# This makefile creates symbolic links to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

incobjs = fadbad.h fadiff.h badiff.h tadiff.h

#####

install:
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(fbpath)/$$INC $$INC; \
		fi; \
	done
	@echo

#####

cleaninstall:
	cd $(incpath) ; rm $(incobjs)
