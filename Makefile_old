OBJS = examples/analysis.cc

ifeq ($(USER),jthaler)
	-include config_jthaler.mk
	PATH_TO_FASTJET = $(FASTJETLOCATION)/fastjet-config
else 
	PATH_TO_FASTJET = ../fastjet-install/bin/fastjet-config
endif

CC = g++
DEBUG = -g
STANDARD = -std=c++11
EXTERNAL = `$(PATH_TO_FASTJET) --cxxflags --libs --plugins`
BIN=./bin/

$(BIN)/analysis : $(OBJS)
	$(CC) $(OBJS) $(STANDARD) -o ./bin/analysis $(EXTERNAL)

clean:
	\rm ./bin/analysis
	\rm ./bin/filter