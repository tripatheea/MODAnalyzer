PATH_TO_FASTJET = /usr/local/bin/fastjet-config

CXX = g++
CXXFLAGS= -O3 -Wall -Woverloaded-virtual -g -std=c++11

FASTINC = `$(PATH_TO_FASTJET) --cxxflags`
FASTLIB = `$(PATH_TO_FASTJET) --libs --plugins` -lRecursiveTools

ROOTINC = `root-config --cflags --glibs`


OBJDIR=src
EXECDIR=exec
BINDIR=bin
INCDIR=interface
INC= -I$(INCDIR)

_OBJ = InfoCalibratedJet InfoPFC Event Trigger Property Condition helpers
OBJ  = $(patsubst %,$(OBJDIR)/%,$(_OBJ:=.o))



# _EXEC=skim analyze turn_on convert_to_pristine analyze_data write
_EXEC=analyze_pair_pfc analyze_pfc analyze convert_to_one_jet write analyze_data skim analyze_triggers analyze_lumi triggers analyze_weights
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


