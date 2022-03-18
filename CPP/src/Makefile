# Makefile!

CC=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-lboost_program_options -lboost_timer
TARGET=main

VERBOSE=TRUE


default: $(TARGET)
all: $(TARGET)

$(TARGET).exe: $(TARGET)


.PHONY: clean doc


run: $(TARGET)
	./$(TARGET) --dt 0.001 --T 10.0 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 50.0


clean:
	rm -rf $(TARGET) *.o

