CC      = g++
FLAGS   = -G -Wall
LDFLAGS =
LIBRARY =
INCLUDE =

EXEC = main
SRCS = main.cpp Mat.cpp HouseHolder.cpp functions.cpp
OBJS = main.o Mat.o HouseHolder.o functions.o



all : $(EXEC)

$(EXEC) : $(OBJS)
		$(CC) $(LDFLAGS) -o $@ $^

%.o : %.cpp
		$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean, mrproper

clean :
		rm -rf *.o

mrproper: clean
		rm -rf $(EXEC)
