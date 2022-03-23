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
	# (101, 101, 0.75, 0.06, 50.0, 5.0, 0.0)
	# ./$(TARGET) --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 50.0
	# (251, 251, 0.75, 0.06, 13.0, 5.0, 0.0)
	# ./$(TARGET) --dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 13.0
	# (101, 101, 0.5, 0.1, 50.0, 5.0, 0.0)
	# ./$(TARGET) --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.5 --b 0.1 --mu1 5.0 --mu2 0.0 --eps 50.0
	# (151, 81, 0.75, 0.0001, 12.5, 1.0, 0.01)
	./$(TARGET) --dt 0.001 --T 100 --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --mu1 1.0 --mu2 0.01 --eps 12.5

clean:
	rm -rf $(TARGET) *.o