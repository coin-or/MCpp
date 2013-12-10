# This makefile creates symbolic links to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

incobjs = mcfunc.hpp mcop.hpp interval.hpp mccormick.hpp tmodel.hpp \
          specbnd.hpp mcprofil.hpp mcfilib.hpp mcfadbad.hpp mclapack.hpp

#####

install: dispInstall
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(mcpath)/$$INC $$INC; \
		fi; \
	done
	@echo

dispInstall:
	@echo
	@(echo '***Installing MC++ library (ver.' $(version)')***')
	@echo

#####

cleaninstall: dispCleanInstall
	cd $(incpath) ; rm $(incobjs)

dispCleanInstall:
	@echo
	@(echo '***Uninstalling MC++ library (ver.' $(version)')***')
	@echo
