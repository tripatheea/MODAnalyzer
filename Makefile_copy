OBJS = cluster.o event.o ntilde.o particle.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 $(DEBUG) -o '../fastjet-install/bin/fastjet-config' $(DEBUG)
LFLAGS = $(DEBUG)

analysis : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o analysis

cluster.o : interface/cluster.h
	$(CC) $(CFLAGS) src/cluster.cc

ntilde.o : interface/ntilde.h
	$(CC) $(CFLAGS) src/ntilde.cc

trigger.o : interface/trigger.h
	$(CC) $(CFLAGS) src/trigger.cc

event.o : interface/trigger.h interface/particle.h
	$(CC) $(CFLAGS) src/event.cc

particle.o : interface/particle.h 
	$(CC) $(CFLAGS) src/particle.cc

clean:
	\rm *.o *~ analysis
