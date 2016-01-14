ifeq ($(USER),jthaler)
	-include config_jthaler.mk
	PATH_TO_FASTJET = $(FASTJETLOCATION)/fastjet-config
else 
	PATH_TO_FASTJET = /usr/local/bin/fastjet-config
endif

CXX = g++
CXXFLAGS= -O3 -Wall -Woverloaded-virtual -g -std=c++11

FASTINC = `$(PATH_TO_FASTJET) --cxxflags`
FASTLIB = `$(PATH_TO_FASTJET) --libs --plugins` -lRecursiveTools

ROOTINC = `root-config --cflags --glibs`
PYTHIA_INC = $(PYTHIA8)/include
PYTHIA_LIB = $(PYTHIA8)/lib

OBJDIR=src
EXECDIR=examples
BINDIR=bin
INCDIR=interface
INC= -I$(INCDIR)

_OBJ = info_calibrated_jet info_pfc event trigger property condition helpers
OBJ  = $(patsubst %,$(OBJDIR)/%,$(_OBJ:=.o))



# _EXEC=skim analyze turn_on convert_to_pristine analyze_data write
_EXEC=analyze convert_to_pristine write analyze_data skim turn_on
EXEC=$(patsubst %,$(EXECDIR)/%,$(_EXEC:=.o))
BIN=$(patsubst %,$(BINDIR)/%,$(_EXEC))


all: $(BIN)

$(OBJDIR)/%.o : $(OBJDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC) -I$(PYTHIA_INC)

$(EXECDIR)/%.o : $(EXECDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC) -I$(PYTHIA_INC)
	
$(BINDIR)/% : $(EXECDIR)/%.o $(OBJ)
	$(CXX) $< $(OBJ) -o $@ $(CXXFLAGS) $(FASTLIB) -lpythia8 -ldl -L$(PYTHIA_LIB)

.PHONY: clean
.PRECIOUS: $(OBJ) $(EXEC)

clean:
	rm $(OBJ) $(EXEC) $(BIN)


