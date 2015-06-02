OBJS = analysis.cc
CC = g++
DEBUG = -g
STANDARD = -std=c++11
EXTERNAL = `../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

analysis : $(OBJS)
	$(CC) $(OBJS) $(STANDARD) -o analysis $(EXTERNAL)

clean:
	\rm *.o *~ p1