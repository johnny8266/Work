NAME = topeg

FFLAGS  =  -O3 
CPPFLAGS  =  -O3 -std=c++11

ROOTLIBS = $(shell root-config --libs)
CPPFLAGS += $(shell root-config --cflags)

SRCDIR   = src
OBJDIR   = obj

SRCF  = $(wildcard $(SRCDIR)/*.f)
SRCC  = $(wildcard $(SRCDIR)/*.cc)
OBJ   = $(SRCC:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
OBJ  += $(SRCF:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
rm    = rm -f

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	gfortran -c $(FFLAGS) -o $@ $(SRCDIR)/$*.f
$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	g++ $(CPPFLAGS) -c $(SRCDIR)/$*.cc -o $(OBJDIR)/$*.o

go: mdir $(OBJ)
	g++ $(CPPFLAGS) -o $(NAME) $(OBJ)  -lgfortran -L$(CERN_ROOT)/lib -lmathlib -lpawlib -lpacklib -lkernlib $(ROOTLIBS) -lFoam
#	gfortran  $(FFLAGS) -o $(NAME) $(OBJ)  -L$(CERN_ROOT)/lib -lpawlib -lpacklib -lkernlib \
                                               -L$(CLAS_LIB) -lbosio -lbos -lfpack -lc_bos_io -lrecutl

mdir:
	mkdir -p $(OBJDIR)
	mkdir -p grids

clean:
	rm -f $(NAME) $(OBJ)

