toCompile = votingModelCPI.o votingModel.o vote.o fitCurves.o

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x #-O3

all: vote

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

vote: $(toCompile)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) *.o 

debug:
	gdb --args vote random 1000 100000 100 0.5 2 0.5 0.5

