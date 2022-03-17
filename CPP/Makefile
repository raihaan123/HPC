# Makefile to make main.cpp

CC=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-lboost_program_options -lboost_timer
TARGET=main

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o

.PHONY: clean doc

clean:
	rm -rf $(TARGET) *.o