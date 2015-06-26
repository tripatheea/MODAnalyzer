ifeq ($(USER),jthaler)
	-include config_jthaler.mk
	PATH_TO_FASTJET = $(FASTJETLOCATION)/fastjet-config
else 
	PATH_TO_FASTJET = ../fastjet-install/bin/fastjet-config
endif

CXX = g++
CXXFLAGS= -O3 -Wall -Woverloaded-virtual -g -std=c++11

FASTINC = `$(PATH_TO_FASTJET) --cxxflags`
FASTLIB = `$(PATH_TO_FASTJET) --libs --plugins` -lRecursiveTools

ROOTINC = `root-config --cflags --glibs`


OBJDIR=src
EXECDIR=examples
BINDIR=bin
INCDIR=interface
INC= -I$(INCDIR)

_OBJ =calibrated_jet event fractional_jet_multiplicity pfcandidate trigger property condition
OBJ  =$(patsubst %,$(OBJDIR)/%,$(_OBJ:=.o))


_EXEC=skim analyze validate plots
EXEC=$(patsubst %,$(EXECDIR)/%,$(_EXEC:=.o))
BIN=$(patsubst %,$(BINDIR)/%,$(_EXEC))


all: $(BIN)

$(OBJDIR)/%.o : $(OBJDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC) $(ROOTINC)

$(EXECDIR)/%.o : $(EXECDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC) $(ROOTINC)
	
$(BINDIR)/% : $(EXECDIR)/%.o $(OBJ)
	$(CXX) $< $(OBJ) -o $@ $(CXXFLAGS) $(FASTLIB) $(ROOTINC)

.PHONY: clean
.PRECIOUS: $(OBJ) $(EXEC)

clean:
	rm $(OBJ) $(EXEC) $(BIN)


