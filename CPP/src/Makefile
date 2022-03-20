# Makefile!

CXX=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-fopenmp -lboost_program_options -lboost_timer
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
	./$(TARGET) --dt 0.001 --T 10.0 --Nx 11 --Ny 11 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 50.0


clean:
	rm -rf $(TARGET) *.o

