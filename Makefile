OBJS = examples/filter.cc
PATH_TO_FASTJET = ../fastjet-install/bin/fastjet-config
CC = g++
DEBUG = -g
STANDARD = -std=c++11
EXTERNAL = `$(PATH_TO_FASTJET) --cxxflags --libs --plugins`
BIN=./bin/

$(BIN)/filter : $(OBJS)
	$(CC) $(OBJS) $(STANDARD) -o ./bin/filter $(EXTERNAL)

$(BIN)/analysis : $(OBJS)
	$(CC) $(OBJS) $(STANDARD) -o ./bin/analysis $(EXTERNAL)

clean:
	\rm ./bin/analysis
	\rm ./bin/filter