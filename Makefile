ifeq ($(USER),jthaler)
	-include config_jthaler.mk
	PATH_TO_FASTJET = $(FASTJETLOCATION)/fastjet-config
else 
	PATH_TO_FASTJET = ../fastjet-install/bin/fastjet-config
endif

CXX = g++
CXXFLAGS= -O3 -Wall -Woverloaded-virtual -g -std=c++11

FASTINC = `$(PATH_TO_FASTJET) --cxxflags`
<<<<<<< HEAD
FASTLIB = `$(PATH_TO_FASTJET) --libs --plugins` -lRecursiveTools
=======
FASTLIB = `$(PATH_TO_FASTJET) --libs --plugins`
>>>>>>> fb4fa903be032715b1351db53399a0357c1f7cf4

OBJDIR=src
EXECDIR=examples
BINDIR=bin
INCDIR=interface
INC= -I$(INCDIR)

_OBJ =calibrated_jet event fractional_jet_multiplicity pfcandidate trigger
OBJ  =$(patsubst %,$(OBJDIR)/%,$(_OBJ:=.o))

<<<<<<< HEAD
_EXEC=filter analyze validate corrected_ak5_spectrum corrected_n_tilde analyze_zg invariant_mass
=======
_EXEC=filter analyze validate corrected_ak5_spectrum corrected_n_tilde
>>>>>>> fb4fa903be032715b1351db53399a0357c1f7cf4
EXEC=$(patsubst %,$(EXECDIR)/%,$(_EXEC:=.o))
BIN=$(patsubst %,$(BINDIR)/%,$(_EXEC))


all: $(BIN)

$(OBJDIR)/%.o : $(OBJDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC)

$(EXECDIR)/%.o : $(EXECDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC)
	
$(BINDIR)/% : $(EXECDIR)/%.o $(OBJ)
	$(CXX) $< $(OBJ) -o $@ $(CXXFLAGS) $(FASTLIB)

.PHONY: clean
.PRECIOUS: $(OBJ) $(EXEC)

clean:
	rm $(OBJ) $(EXEC) $(BIN)


