# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linuxegcs

CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 
DictSuf       = Dict.o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)



ifeq ($(ARCH),linux)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxrh42)
# Linux with gcc 2.7.2.x (RedHat 4.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxdeb)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxdeb2)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxsuse6)
# Linux with gcc 2.7.2.x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxegcs)
# Linux with egcs (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxkcc)
# Linux with the KAI compiler
CXX           = KCC
CXXFLAGS      = -fPIC +K0
LD            = KCC
LDFLAGS       = -O
SOFLAGS       =
endif

ifeq ($(ARCH),linuxppcegcs)
# MkLinux with egcs/glibc
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared -Wl,-soname,
endif

ifeq ($(ARCH),linuxia64gcc)
# Itanium Linux with gcc 2.9x
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxia64sgi)
# Itanium Linux with sgiCC
CXX           = sgiCC
CXXFLAGS      = -O -Wall -fPIC
LD            = gsgiCC
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxalphaegcs)
# Alpha Linux with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxarm)
# ARM Linux with egcs
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

DVCSFLAGS     = -I/afs/in2p3.fr/home/throng/clas/franck/theory/INCLUDE
DVCSLIBS      = -L/afs/in2p3.fr/home/throng/clas/franck/theory/LIB -lmydvcs -lTHEO -lMRST -lGPDQ
#CXXFLAGS     += $(ROOTCFLAGS) $(DVCSFLAGS)
CXXFLAGS     += $(ROOTCFLAGS)
#LIBS          = $(ROOTLIBS) $(SYSLIBS)
LIBS         = -L . -lTGenDVCS -L . -lTGVKelly -L $(ROOTSYS)/lib -lFoam
LIBS	     += $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------


TGVKellyO           =TGVKelly.$(ObjSuf) \
                       TGVKellyDict.$(ObjSuf)

TGVKellyS           =TGVKelly.$(SrcSuf) \
                       TGVKellyDict.$(SrcSuf)

TGVKellySO          =libTGVKelly.$(DllSuf)


TGenDVCSO           =TGenDVCS.$(ObjSuf) \
                       TGenDVCSDict.$(ObjSuf)

TGenDVCSS           =TGenDVCS.$(SrcSuf) \
                       TGenDVCSDict.$(SrcSuf)

TGenDVCSSO          =libTGenDVCS.$(DllSuf)


TGenBaseO           =TGenBase.$(ObjSuf) \
                       TGenBaseDict.$(ObjSuf)

TGenBaseS           =TGenBase.$(SrcSuf) \
                       TGenBaseDict.$(SrcSuf)

TGenBaseSO          =libTGenBase.$(DllSuf)


TGenGeoO           =TGenGeo.$(ObjSuf) \
                       TGenGeoDict.$(ObjSuf)

TGenGeoS           =TGenGeo.$(SrcSuf) \
                       TGenGeoDict.$(SrcSuf)

TGenGeoSO          =libTGenGeo.$(DllSuf)


#InfoO              =Info.$(ObjSuf) \
		        InfoDict.$(ObjSuf)

#InfoS              =Info.$(SrcSuf) \
		        InfoDict.$(SrcSuf)

#InfoSO             =libInfo.$(DllSuf)


TFDISTRO           =TFDISTR.$(ObjSuf) \
                       TFDISTRDict.$(ObjSuf)

TFDISTRS           =TFDISTR.$(SrcSuf) \
                       TFDISTRDict.$(SrcSuf)

TFDISTRSO          =libTFDISTR.$(DllSuf)


MAINO               = foam_demo.$(ObjSuf)
MAINS               = foam_demo.$(SrcSuf)
MAIN                = foam_demo$(ExeSuf)


OBJS          = $(MAINO) $(TFDISTRO) $(TGVKellyO) $(TGenDVCSO)

SLIB          = $(TFDISTRSO) $(TGVKellySO) $(TGenDVCSSO) $(TGenBaseSO) $(TGenGeoSO) 

PROGRAMS      = $(MAIN)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(SLIB)

lib:            $(SLIB)

doc:
	root.exe -b -q make_html.C

%.$(ObjSuf):%.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<

%Dict.$(SrcSuf):%.h
	rootcint -f $@ -c $< 

$(MAIN):  $(MAINO) $(SLIB)
	$(LD) $(CXXFLAGS) $(LIBS)  $^ -o $@ 

lib%.$(DllSuf): %Dict.$(ObjSuf) %.$(ObjSuf)
	$(LD) $(ROOTGLIBS) $(SOFLAGS) $(CXXFLAGS) $^ -o $@

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(SLIB) *Dict.* *.def


























