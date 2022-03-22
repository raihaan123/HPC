# Makefile!

CXX=g++
CXXFLAGS=-g -Wall -O3 -pedantic
LDLIBS=-fopenmp -lboost_program_options -lboost_timer -lblas
TARGET=main
OBJS=main.o ReactionDiffusion.o
HDRS=ReactionDiffusion.h

VERBOSE=TRUE


default: $(TARGET)


%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

main: $(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS)


all: $(TARGET)


.PHONY: clean doc


run: $(TARGET)
	./$(TARGET) --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 50.0
	# ./$(TARGET) --dt 0.001 --T 10 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 13.0


clean:
	rm -rf $(TARGET) *.o

