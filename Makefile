# Makefile.in generated by automake 1.11.1 from Makefile.am.
# examples/Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
# 2003, 2004, 2005, 2006, 2007, 2008, 2009  Free Software Foundation,
# Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.




pkgdatadir = $(datadir)/tauola-c---interface
pkgincludedir = $(includedir)/tauola-c---interface
pkglibdir = $(libdir)/tauola-c---interface
pkglibexecdir = $(libexecdir)/tauola-c---interface
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = x86_64-unknown-linux-gnu
host_triplet = x86_64-unknown-linux-gnu
example_PROGRAMS = taumain_hepevt_example.exe$(EXEEXT) $(am__EXEEXT_1) \
	$(am__EXEEXT_2) $(am__EXEEXT_3)

### All other examples require HepMC ###
am__append_1 = -R $(HEPMC_DIR)/lib
am__append_2 = -I$(HEPMC_DIR)/include
am__append_3 = -L$(HEPMC_DIR)/lib -lHepMC
am__append_4 = taumain_stand_alone_example.exe

### Taugun & taummk examples (require Pythia8) ###
am__append_5 = -R $(PYTHIA8_DIR)/lib/archive
am__append_6 = -I$(PYTHIA8_DIR)/include
#am__append_7 = -L$(PYTHIA8_DIR)/lib/archive -lpythia8 -llhapdfdummy -lpythia8tohepmc
#am__append_8 = -DPYTHIA8180_OR_LATER
am__append_9 = -L$(PYTHIA8_DIR)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface
am__append_10 = single_tau_gun_example.exe taummk_pythia_example.exe

### Main example (requires Pythia8 and MC-Tester) ###
am__append_11 = -R $(MCTESTER_DIR)/lib  
am__append_12 = -I$(PYTHIA8_DIR)/include -I$(MCTESTER_DIR)/include -I$(ROOTINC)   -I$(USERINC)  
am__append_13 = -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester $(ROOTLIBS) $(USERLIBS)
am__append_14 = taumain_pythia_example.exe
am__append_15 = mypythia_example.exe
subdir = examples
DIST_COMMON = README $(srcdir)/Makefile.am $(srcdir)/Makefile.in
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps = $(top_srcdir)/configure.ac
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
mkinstalldirs = $(install_sh) -d
CONFIG_HEADER = $(top_builddir)/config/config.h
CONFIG_CLEAN_FILES =
CONFIG_CLEAN_VPATH_FILES =
am__EXEEXT_1 =  \
	taumain_stand_alone_example.exe$(EXEEXT)
am__EXEEXT_2 = single_tau_gun_example.exe$(EXEEXT) \
	taummk_pythia_example.exe$(EXEEXT)
am__EXEEXT_3 = taumain_pythia_example.exe$(EXEEXT)
am__EXEEXT_3 = mypythia_example.exe$(EXEEXT)
am__installdirs = "$(DESTDIR)$(exampledir)"
PROGRAMS = $(example_PROGRAMS)
am__single_tau_gun_example_exe_SOURCES_DIST =  \
	single_tau_gun_example.cxx
am_single_tau_gun_example_exe_OBJECTS = single_tau_gun_example.$(OBJEXT)
single_tau_gun_example_exe_OBJECTS =  \
	$(am_single_tau_gun_example_exe_OBJECTS)
single_tau_gun_example_exe_LDADD = $(LDADD)
am__DEPENDENCIES_1 =
am__DEPENDENCIES_2 = $(am__DEPENDENCIES_1)
single_tau_gun_example_exe_DEPENDENCIES =  \
	$(prefix)/lib/libTauolaCxxInterface.so \
	$(prefix)/lib/libTauolaFortran.so $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_1) $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_2)
am_taumain_hepevt_example_exe_OBJECTS =  \
	taumain_hepevt_example.$(OBJEXT)
taumain_hepevt_example_exe_OBJECTS =  \
	$(am_taumain_hepevt_example_exe_OBJECTS)
taumain_hepevt_example_exe_LDADD = $(LDADD)
taumain_hepevt_example_exe_DEPENDENCIES =  \
	$(prefix)/lib/libTauolaCxxInterface.so \
	$(prefix)/lib/libTauolaFortran.so $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_1) $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_2)
am__taumain_pythia_example_exe_SOURCES_DIST =  \
	taumain_pythia_example.cxx
am_taumain_pythia_example_exe_OBJECTS = taumain_pythia_example.$(OBJEXT)
taumain_pythia_example_exe_OBJECTS =  \
	$(am_taumain_pythia_example_exe_OBJECTS)
taumain_pythia_example_exe_LDADD = $(LDADD)
taumain_pythia_example_exe_DEPENDENCIES =  \
	$(prefix)/lib/libTauolaCxxInterface.so \
	$(prefix)/lib/libTauolaFortran.so $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_1) $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_2)
am__mypythia_example_exe_SOURCES_DIST =  \
	mypythia_example.cxx
am_mypythia_example_exe_OBJECTS = mypythia_example.$(OBJEXT)
mypythia_example_exe_OBJECTS =  \
	$(am_mypythia_example_exe_OBJECTS)
mypythia_example_exe_LDADD = $(LDADD)
mypythia_example_exe_DEPENDENCIES =  \
	$(prefix)/lib/libTauolaCxxInterface.so \
	$(prefix)/lib/libTauolaFortran.so $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_1) $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_2)
am__taumain_stand_alone_example_exe_SOURCES_DIST =  \
	taumain_stand_alone_example.cxx
am_taumain_stand_alone_example_exe_OBJECTS =  \
	taumain_stand_alone_example.$(OBJEXT)
taumain_stand_alone_example_exe_OBJECTS =  \
	$(am_taumain_stand_alone_example_exe_OBJECTS)
taumain_stand_alone_example_exe_LDADD = $(LDADD)
taumain_stand_alone_example_exe_DEPENDENCIES =  \
	$(prefix)/lib/libTauolaCxxInterface.so \
	$(prefix)/lib/libTauolaFortran.so $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_1) $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_2)
am__taummk_pythia_example_exe_SOURCES_DIST =  \
	taummk_pythia_example.cxx
am_taummk_pythia_example_exe_OBJECTS = taummk_pythia_example.$(OBJEXT)
taummk_pythia_example_exe_OBJECTS =  \
	$(am_taummk_pythia_example_exe_OBJECTS)
taummk_pythia_example_exe_LDADD = $(LDADD)
taummk_pythia_example_exe_DEPENDENCIES =  \
	$(prefix)/lib/libTauolaCxxInterface.so \
	$(prefix)/lib/libTauolaFortran.so $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_1) $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_2)
DEFAULT_INCLUDES = -I. -I$(top_builddir)/config
depcomp = $(SHELL) $(top_srcdir)/config/depcomp
am__depfiles_maybe = depfiles
am__mv = mv -f
CXXCOMPILE = $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
LTCXXCOMPILE = $(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile $(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
CXXLD = $(CXX)
CXXLINK = $(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) \
	$(LDFLAGS) -o $@
SOURCES = $(single_tau_gun_example_exe_SOURCES) \
	$(taumain_hepevt_example_exe_SOURCES) \
	$(taumain_pythia_example_exe_SOURCES) \
	$(mypythia_example_exe_SOURCES) \
	$(taumain_stand_alone_example_exe_SOURCES) \
	$(taummk_pythia_example_exe_SOURCES)
DIST_SOURCES = $(am__single_tau_gun_example_exe_SOURCES_DIST) \
	$(taumain_hepevt_example_exe_SOURCES) \
	$(am__taumain_pythia_example_exe_SOURCES_DIST) \
	$(am__mypythia_example_exe_SOURCES_DIST) \
	$(am__taumain_stand_alone_example_exe_SOURCES_DIST) \
	$(am__taummk_pythia_example_exe_SOURCES_DIST)
ETAGS = etags
CTAGS = ctags
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
ACLOCAL = ${SHELL} /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/config/missing --run aclocal-1.11
AMTAR = ${SHELL} /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/config/missing --run tar
AR = ar
AUTOCONF = ${SHELL} /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/config/missing --run autoconf
AUTOHEADER = ${SHELL} /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/config/missing --run autoheader
AUTOMAKE = ${SHELL} /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/config/missing --run automake-1.11
AWK = gawk
CC = gcc
CCDEPMODE = depmode=gcc3
CFLAGS = -g -O2
CPPFLAGS = -I/grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/pythia8/176//include -I/grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/lhapdf-5.9.1/workdir//include -I/grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/HepMC-2.06.05/workdir/include  $(am__append_8)
CXX = g++
CXXCPP = g++ -E
CXXDEPMODE = depmode=gcc3
CXXFLAGS = -O2
CYGPATH_W = echo
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
DSYMUTIL = 
DUMPBIN = 
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
F77 = gfortran
FFLAGS = -O2
FGREP = /bin/grep -F
GREP = /bin/grep
HAS_ROOT_CONFIG = yes
HEPMC_DIR = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/HepMC-2.06.05/workdir
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
LD = /usr/bin/ld -m elf_x86_64
LDFLAGS = 
LIBOBJS = 
LIBS = 
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LIPO = 
LN_S = ln -s
LTLIBOBJS = 
MAKEINFO = ${SHELL} /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/config/missing --run makeinfo
MCTESTER_DIR = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/MC-TESTER/
MKDIR_P = /bin/mkdir -p
NM = /usr/bin/nm -B
NMEDIT = 
OBJDUMP = objdump
OBJEXT = o
OTOOL = 
OTOOL64 = 
PACKAGE = tauola-c---interface
PACKAGE_BUGREPORT = tomasz.przedzinski@cern.ch 
PACKAGE_NAME = Tauola C++ Interface
PACKAGE_STRING = Tauola C++ Interface 1.0.6
PACKAGE_TARNAME = tauola-c---interface
PACKAGE_VERSION = 1.0.6
PATH_SEPARATOR = :
PYTHIA8_DIR = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/pythia8/176/
RANLIB = ranlib
ROOTINC = /libcern/root/5.34.18/sl6.3-x86_64/include
ROOTLIBS = -L/libcern/root/5.34.18/sl6.3-x86_64/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
USERINC = UserCodes
USERLIBS = -L./UserCodes -lUserLib
SED = /bin/sed
SET_MAKE = 
SHELL = /bin/sh
STRIP = strip
VERSION = 1.0.6
abs_builddir = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/examples
abs_srcdir = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/examples
abs_top_builddir = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5
abs_top_srcdir = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5
ac_ct_CC = gcc
ac_ct_CXX = g++
ac_ct_DUMPBIN = 
ac_ct_F77 = gfortran
am__include = include
am__leading_dot = .
am__quote = 
am__tar = ${AMTAR} chof - "$$tardir"
am__untar = ${AMTAR} xf -
bindir = ${exec_prefix}/bin
build = x86_64-unknown-linux-gnu
build_alias = 
build_cpu = x86_64
build_os = linux-gnu
build_vendor = unknown
builddir = .
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
dvidir = ${docdir}
exec_prefix = ${prefix}
host = x86_64-unknown-linux-gnu
host_alias = 
host_cpu = x86_64
host_os = linux-gnu
host_vendor = unknown
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = ${SHELL} /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/config/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
lt_ECHO = echo
mandir = ${datarootdir}/man
mkdir_p = /bin/mkdir -p
oldincludedir = /usr/include
pdfdir = ${docdir}
prefix = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/workdir
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
srcdir = .
sysconfdir = ${prefix}/etc
target_alias = 
top_build_prefix = ../
top_builddir = ..
top_srcdir = ..
with_lhapdf = /grid_mnt/home-pbs/vcherepa/taua1/Installation/taudir/tauola++/1.1.5/lhapdf-5.9.1/workdir/
exampledir = $(top_srcdir)/examples
taumain_hepevt_example_exe_SOURCES = taumain_hepevt_example.cxx
INCLUDES = -I$(prefix)/include $(am__append_2) $(am__append_6) \
	$(am__append_12)
AM_LDFLAGS = -R $(prefix)/lib $(am__append_1) $(am__append_5) \
	$(am__append_11)
LDADD = $(FLIBS) $(prefix)/lib/libTauolaCxxInterface.so \
	$(prefix)/lib/libTauolaFortran.so $(am__append_3) \
	$(am__append_7) $(am__append_9) $(am__append_13)
taumain_stand_alone_example_exe_SOURCES = taumain_stand_alone_example.cxx
single_tau_gun_example_exe_SOURCES = single_tau_gun_example.cxx
taummk_pythia_example_exe_SOURCES = taummk_pythia_example.cxx
taumain_pythia_example_exe_SOURCES = taumain_pythia_example.cxx
mypythia_example_exe_SOURCES = mypythia_example.cxx
all: all-am

.SUFFIXES:
.SUFFIXES: .cxx .lo .o .obj
$(srcdir)/Makefile.in:  $(srcdir)/Makefile.am  $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      ( cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh ) \
	        && { if test -f $@; then exit 0; else break; fi; }; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --gnu examples/Makefile'; \
	$(am__cd) $(top_srcdir) && \
	  $(AUTOMAKE) --gnu examples/Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe);; \
	esac;

$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

$(top_srcdir)/configure:  $(am__configure_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(ACLOCAL_M4):  $(am__aclocal_m4_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(am__aclocal_m4_deps):
install-examplePROGRAMS: $(example_PROGRAMS)
	@$(NORMAL_INSTALL)
	test -z "$(exampledir)" || $(MKDIR_P) "$(DESTDIR)$(exampledir)"
	@list='$(example_PROGRAMS)'; test -n "$(exampledir)" || list=; \
	for p in $$list; do echo "$$p $$p"; done | \
	sed 's/$(EXEEXT)$$//' | \
	while read p p1; do if test -f $$p || test -f $$p1; \
	  then echo "$$p"; echo "$$p"; else :; fi; \
	done | \
	sed -e 'p;s,.*/,,;n;h' -e 's|.*|.|' \
	    -e 'p;x;s,.*/,,;s/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/' | \
	sed 'N;N;N;s,\n, ,g' | \
	$(AWK) 'BEGIN { files["."] = ""; dirs["."] = 1 } \
	  { d=$$3; if (dirs[d] != 1) { print "d", d; dirs[d] = 1 } \
	    if ($$2 == $$4) files[d] = files[d] " " $$1; \
	    else { print "f", $$3 "/" $$4, $$1; } } \
	  END { for (d in files) print "f", d, files[d] }' | \
	while read type dir files; do \
	    if test "$$dir" = .; then dir=; else dir=/$$dir; fi; \
	    test -z "$$files" || { \
	    echo " $(INSTALL_PROGRAM_ENV) $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL_PROGRAM) $$files '$(DESTDIR)$(exampledir)$$dir'"; \
	    $(INSTALL_PROGRAM_ENV) $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL_PROGRAM) $$files "$(DESTDIR)$(exampledir)$$dir" || exit $$?; \
	    } \
	; done

uninstall-examplePROGRAMS:
	@$(NORMAL_UNINSTALL)
	@list='$(example_PROGRAMS)'; test -n "$(exampledir)" || list=; \
	files=`for p in $$list; do echo "$$p"; done | \
	  sed -e 'h;s,^.*/,,;s/$(EXEEXT)$$//;$(transform)' \
	      -e 's/$$/$(EXEEXT)/' `; \
	test -n "$$list" || exit 0; \
	echo " ( cd '$(DESTDIR)$(exampledir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(exampledir)" && rm -f $$files

clean-examplePROGRAMS:
	@list='$(example_PROGRAMS)'; test -n "$$list" || exit 0; \
	echo " rm -f" $$list; \
	rm -f $$list || exit $$?; \
	test -n "$(EXEEXT)" || exit 0; \
	list=`for p in $$list; do echo "$$p"; done | sed 's/$(EXEEXT)$$//'`; \
	echo " rm -f" $$list; \
	rm -f $$list
single_tau_gun_example.exe$(EXEEXT): $(single_tau_gun_example_exe_OBJECTS) $(single_tau_gun_example_exe_DEPENDENCIES) 
	@rm -f single_tau_gun_example.exe$(EXEEXT)
	$(CXXLINK) $(single_tau_gun_example_exe_OBJECTS) $(single_tau_gun_example_exe_LDADD) $(LIBS)
taumain_hepevt_example.exe$(EXEEXT): $(taumain_hepevt_example_exe_OBJECTS) $(taumain_hepevt_example_exe_DEPENDENCIES) 
	@rm -f taumain_hepevt_example.exe$(EXEEXT)
	$(CXXLINK) $(taumain_hepevt_example_exe_OBJECTS) $(taumain_hepevt_example_exe_LDADD) $(LIBS)
taumain_pythia_example.exe$(EXEEXT): $(taumain_pythia_example_exe_OBJECTS) $(taumain_pythia_example_exe_DEPENDENCIES) 
	@rm -f taumain_pythia_example.exe$(EXEEXT)
	$(CXXLINK) $(taumain_pythia_example_exe_OBJECTS) $(taumain_pythia_example_exe_LDADD) $(LIBS)
mypythia_example.exe$(EXEEXT): $(mypythia_example_exe_OBJECTS) $(mypythia_example_exe_DEPENDENCIES) 
	@rm -f mypythia_example.exe$(EXEEXT)
	$(CXXLINK) $(mypythia_example_exe_OBJECTS) $(mypythia_example_exe_LDADD) $(LIBS)
taumain_stand_alone_example.exe$(EXEEXT): $(taumain_stand_alone_example_exe_OBJECTS) $(taumain_stand_alone_example_exe_DEPENDENCIES) 
	@rm -f taumain_stand_alone_example.exe$(EXEEXT)
	$(CXXLINK) $(taumain_stand_alone_example_exe_OBJECTS) $(taumain_stand_alone_example_exe_LDADD) $(LIBS)
taummk_pythia_example.exe$(EXEEXT): $(taummk_pythia_example_exe_OBJECTS) $(taummk_pythia_example_exe_DEPENDENCIES) 
	@rm -f taummk_pythia_example.exe$(EXEEXT)
	$(CXXLINK) $(taummk_pythia_example_exe_OBJECTS) $(taummk_pythia_example_exe_LDADD) $(LIBS)

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

include ./$(DEPDIR)/single_tau_gun_example.Po
include ./$(DEPDIR)/taumain_hepevt_example.Po
include ./$(DEPDIR)/taumain_pythia_example.Po
include ./$(DEPDIR)/mypythia_example.Po
include ./$(DEPDIR)/taumain_stand_alone_example.Po
include ./$(DEPDIR)/taummk_pythia_example.Po

.cxx.o:
	$(CXXCOMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
#	source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
#	$(CXXCOMPILE) -c -o $@ $<

.cxx.obj:
	$(CXXCOMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ `$(CYGPATH_W) '$<'`
	$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Po
#	source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
#	$(CXXCOMPILE) -c -o $@ `$(CYGPATH_W) '$<'`

.cxx.lo:
	$(LTCXXCOMPILE) -MT $@ -MD -MP -MF $(DEPDIR)/$*.Tpo -c -o $@ $<
	$(am__mv) $(DEPDIR)/$*.Tpo $(DEPDIR)/$*.Plo
#	source='$<' object='$@' libtool=yes \
#	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
#	$(LTCXXCOMPILE) -c -o $@ $<

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

ID: $(HEADERS) $(SOURCES) $(LISP) $(TAGS_FILES)
	list='$(SOURCES) $(HEADERS) $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	mkid -fID $$unique
tags: TAGS

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	set x; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	shift; \
	if test -z "$(ETAGS_ARGS)$$*$$unique"; then :; else \
	  test -n "$$unique" || unique=$$empty_fix; \
	  if test $$# -gt 0; then \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      "$$@" $$unique; \
	  else \
	    $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	      $$unique; \
	  fi; \
	fi
ctags: CTAGS
CTAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '{ files[$$0] = 1; nonempty = 1; } \
	      END { if (nonempty) { for (i in files) print i; }; }'`; \
	test -z "$(CTAGS_ARGS)$$unique" \
	  || $(CTAGS) $(CTAGSFLAGS) $(AM_CTAGSFLAGS) $(CTAGS_ARGS) \
	     $$unique

GTAGS:
	here=`$(am__cd) $(top_builddir) && pwd` \
	  && $(am__cd) $(top_srcdir) \
	  && gtags -i $(GTAGS_ARGS) "$$here"

distclean-tags:
	-rm -f TAGS ID GTAGS GRTAGS GSYMS GPATH tags

distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
check: check-am
all-am: Makefile $(PROGRAMS)
installdirs:
	for dir in "$(DESTDIR)$(exampledir)"; do \
	  test -z "$$dir" || $(MKDIR_P) "$$dir"; \
	done
install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	  install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	  `test -z '$(STRIP)' || \
	    echo "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'"` install
mostlyclean-generic:

clean-generic:

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)
	-test . = "$(srcdir)" || test -z "$(CONFIG_CLEAN_VPATH_FILES)" || rm -f $(CONFIG_CLEAN_VPATH_FILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-am

clean-am: clean-examplePROGRAMS clean-generic clean-libtool \
	mostlyclean-am

distclean: distclean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
distclean-am: clean-am distclean-compile distclean-generic \
	distclean-tags

dvi: dvi-am

dvi-am:

html: html-am

html-am:

info: info-am

info-am:

install-data-am: install-examplePROGRAMS

install-dvi: install-dvi-am

install-dvi-am:

install-exec-am:

install-html: install-html-am

install-html-am:

install-info: install-info-am

install-info-am:

install-man:

install-pdf: install-pdf-am

install-pdf-am:

install-ps: install-ps-am

install-ps-am:

installcheck-am:

maintainer-clean: maintainer-clean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool

pdf: pdf-am

pdf-am:

ps: ps-am

ps-am:

uninstall-am: uninstall-examplePROGRAMS

.MAKE: install-am install-strip

.PHONY: CTAGS GTAGS all all-am check check-am clean \
	clean-examplePROGRAMS clean-generic clean-libtool ctags \
	distclean distclean-compile distclean-generic \
	distclean-libtool distclean-tags distdir dvi dvi-am html \
	html-am info info-am install install-am install-data \
	install-data-am install-dvi install-dvi-am \
	install-examplePROGRAMS install-exec install-exec-am \
	install-html install-html-am install-info install-info-am \
	install-man install-pdf install-pdf-am install-ps \
	install-ps-am install-strip installcheck installcheck-am \
	installdirs maintainer-clean maintainer-clean-generic \
	mostlyclean mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool pdf pdf-am ps ps-am tags uninstall \
	uninstall-am uninstall-examplePROGRAMS


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
