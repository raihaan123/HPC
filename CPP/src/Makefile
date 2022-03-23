# Makefile!

CXX=g++
CXXFLAGS=-Wall -O3 -pedantic
LDLIBS=-fopenmp -lboost_program_options -lboost_timer
TARGET=main
OBJS=main.o ReactionDiffusion.o
HDRS=ReactionDiffusion.h

default: $(TARGET)


# Rule for compiling C++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDLIBS)


# Rule for linking object files
$(TARGET): $(OBJS)			
	$(CXX) -o $@ $^ $(LDLIBS)



# set the OpenMP number of threads variable - OMP_NUM_THREADS - to 16 threads, which empirically is optimal
export OMP_NUM_THREADS := 16

run1: $(TARGET)
	# (101, 101, 0.75, 0.06, 50.0, 5.0, 0.0)
	./$(TARGET) --threads 16 --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 50.0


run2: $(TARGET)
	# (251, 251, 0.75, 0.06, 13.0, 5.0, 0.0)
	./$(TARGET) --threads 16 --dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 13.0


run3: $(TARGET)
	# (101, 101, 0.5, 0.1, 50.0, 5.0, 0.0)
	./$(TARGET) --threads 16 --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.5 --b 0.1 --mu1 5.0 --mu2 0.0 --eps 50.0


run4: $(TARGET)
	# (151, 81, 0.75, 0.0001, 12.5, 1.0, 0.01)
	./$(TARGET) --threads 16 --dt 0.001 --T 100 --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --mu1 1.0 --mu2 0.01 --eps 12.5




# Targets for running the test cases from 1 to 4 in serial
test1: $(TARGET)
	./$(TARGET) --threads 1 --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 50.0


test2: $(TARGET)
	./$(TARGET) --threads 1 --dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 13.0


test3: $(TARGET)
	./$(TARGET) --threads 1 --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.5 --b 0.1 --mu1 5.0 --mu2 0.0 --eps 50.0


test4: $(TARGET)
	./$(TARGET) --threads 1 --dt 0.001 --T 100 --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --mu1 1.0 --mu2 0.01 --eps 12.5



# clean and doc are phony targets
.PHONY: clean doc


clean:
	rm -rf $(TARGET) *.o


doc:
	doxygen Doxyfile