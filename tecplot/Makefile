# set paths and compiler flags
TECLIB_VERSION = 12.0

all: teclib

teclib: 
	@# Tecplot Binary library  
	@ if [ -n "$(TECLIB_VERSION)" ] && [ ! -f tecio-$(TECLIB_VERSION)/tecio.a ] ; then \
		echo "===============================================================";\
		echo "= M A K E   T E C P L O T                                     =";\
		echo "===============================================================";\
		if [ ! -d tecio-$(TECLIB_VERSION) ] ; then \
			tar -xvzf tecio-$(TECLIB_VERSION).tar.gz 1>/dev/null && cd .. ; \
		fi ; \
		cd tecio-$(TECLIB_VERSION) && ./Runmake linuxg27x64.24 -tecio && rm pltview pltview.o && ln -s tecio.a libtecio.a; \
		echo "===============================================================";\
		echo "= T E C P L O T   D O N E                                     =";\
		echo "===============================================================";\
	fi


# utility targets
.PHONY: clean
clean: cleanteclib 
cleanteclib:
	rm -f *~
	@# TECLIB library
	@ if [ -d tecio-$(TECLIB_VERSION) ] ; then \
		cd tecio-$(TECLIB_VERSION) && rm -f tecio.a && rm -f libtecio.a; \
	fi
